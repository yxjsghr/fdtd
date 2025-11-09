#include "FDTD2D.h"
#include "OutputManager.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip> 

FDTD2D::FDTD2D() {
    setupGrid();
    setupSource();
    setupMaterialsLorentz();
    setupPML();
    m_out = std::make_unique<OutputManager>("output_2d");
}

FDTD2D::~FDTD2D() = default;

void FDTD2D::setupGrid() {
    // ===== 对齐 MATLAB 参数 =====
    double f0 = 1.5e14;
    double lambda0 = m_c0 / f0;
    int N_lambda = 40;

    m_dx = lambda0 / N_lambda;
    m_dz = m_dx;

    double Lx = 5.0 * lambda0;
    double Lz = 12.0 * lambda0;

    m_Nx = static_cast<int>(std::round(Lx / m_dx));
    m_Nz = static_cast<int>(std::round(Lz / m_dz));

    double nCFL = 2.0; // 稳定性
    m_dt = m_dx / (m_c0 * nCFL);

    m_steps = 1580;

    // 分配内存（严格 Yee 尺寸）
    Dy.assign( (m_Nx+1)*(m_Nz+1), 0.0 );
    Ey = Dy; Ey_n1 = Dy; Ey_n2 = Dy; SEy = Dy; SEy_n1 = Dy; SEy_n2 = Dy;

    Bx.assign( (m_Nx+1)*m_Nz, 0.0 );
    Hx = Bx; Hx_n1 = Bx; Hx_n2 = Bx; SHx = Bx; SHx_n1 = Bx; SHx_n2 = Bx;

    Bz.assign( m_Nx*(m_Nz+1), 0.0 );
    Hz = Bz; Hz_n1 = Bz; Hz_n2 = Bz; SHz = Bz; SHz_n1 = Bz; SHz_n2 = Bz;

    C1.assign( m_Nx*m_Nz, 0.0 );
    C2 = C1; C3 = C1;

    std::cout << "[grid] Nx=" << m_Nx << ", Nz=" << m_Nz
              << ", dx=" << m_dx << ", dz=" << m_dz
              << ", dt=" << m_dt << ", steps=" << m_steps << "\n";
}

void FDTD2D::setupSource() {
    // 线电流源
    int Nx_c = m_Nx / 2;
    int Nz_c = m_Nz / 2;
    int N_lambda = 40;

    i_Jy = Nx_c - N_lambda;
    k_Jy = Nz_c;

    // Jy(t) + 摆线渐升
    double f0 = 1.5e14, omega0 = 2.0*M_PI*f0, T0 = 1.0/f0;
    Jy.resize(m_steps);
    double Jy0 = 1.0;
    for (int n = 0; n < m_steps; ++n) {
        double t = (n+1)*m_dt;
        Jy[n] = Jy0 * std::sin(omega0*t) / (m_dx*m_dz) / (omega0*m_mu0/4.0);
    }
    double ts = 2.0*T0;
    int tsddt = static_cast<int>(std::round(ts/m_dt));
    for (int n = 0; n <= tsddt && n < m_steps; ++n) {
        double s = static_cast<double>(n)/tsddt;
        double ws = s - std::sin(2.0*M_PI*s)/(2.0*M_PI);
        Jy[n] *= ws;
    }
}

void FDTD2D::setupMaterialsLorentz() {
    // 平板透镜区
    int Nx_c = m_Nx / 2;
    int thick_slab = 40; // N_lambda
    int x_left  = Nx_c - thick_slab/2;
    int x_right = Nx_c + thick_slab/2;
    int z_down  = 0, z_up = m_Nz-1;

    // Lorentz
    double f0 = 1.5e14;
    double omega0 = 2.0 * M_PI * f0;
    double wp = omega0;
    double w0 = wp/std::sqrt(2.0);
    double delta0 = 0.0;

    double wpb = wp*m_dt;
    double w0b = w0*m_dt;
    double delta0b = delta0*m_dt;

    for (int ix = x_left; ix <= x_right; ++ix) {
        for (int kz = z_down; kz <= z_up; ++kz) {
            int idc = id_C(ix, kz);
            double denom = (w0b*w0b + 4.0*delta0b + 4.0);
            C1[idc] = (wpb*wpb) / denom;
            C2[idc] = 2.0*(w0b*w0b - 4.0) / denom;
            C3[idc] = (w0b*w0b - 4.0*delta0b + 4.0) / denom;
        }
    }
}

static std::vector<double> rho_neg(int n, double offset) {
    // ((n - offset) : -1 : 0.25)/n
    std::vector<double> r(n);
    for (int j = 0; j < n; ++j) {
        r[j] = ( (n - offset) - j ) / n;
    }
    return r;
}
static std::vector<double> rho_pos(int n, double offset) {
    // (0.25 : 1 : (n - offset)) / n   —— 用 offset=0.75/0.25
    std::vector<double> r(n);
    for (int j = 0; j < n; ++j) {
        r[j] = (offset + j) / n;
    }
    return r;
}

void FDTD2D::setupPML() {
    // 与 MATLAB 一致：四边均 'on'
    const int m_pml = 3;
    const double alpha_max = 0.05;
    const int sigma_factor = 1;
    const double kappa_max = 7.0;
    const double eps_r_pml = 1.0;
    const double mu_r_pml  = 1.0;

    // --- x- ---
    {
        int n = n_pml_xn;
        if (n > 0) {
            auto rho_e = rho_neg(n, 0.75);
            auto rho_m = rho_neg(n, 0.25);
            double sigma_max = sigma_factor * (m_pml + 1) / (150.0*M_PI*std::sqrt(eps_r_pml*mu_r_pml)*m_dx);
            b_e_xn.resize(n); a_e_xn.resize(n); C_Dy_dxn_v.resize(n);
            b_m_xn.resize(n); a_m_xn.resize(n); C_Bz_dxn_v.resize(n);
            for (int j = 0; j < n; ++j) {
                double sig_e = sigma_max * std::pow(rho_e[j], m_pml);
                double sig_m = (m_mu0/m_eps0) * sigma_max * std::pow(rho_m[j], m_pml);
                double kap_e = 1.0 + (kappa_max - 1.0)*std::pow(rho_e[j], m_pml);
                double kap_m = 1.0 + (kappa_max - 1.0)*std::pow(rho_m[j], m_pml);
                double alp_e = alpha_max * (1.0 - rho_e[j]);
                double alp_m = (m_mu0/m_eps0) * alpha_max * (1.0 - rho_m[j]);

                b_e_xn[j] = std::exp( -(sig_e/kap_e + alp_e) * m_dt/m_eps0 );
                b_m_xn[j] = std::exp( -(sig_m/kap_m + alp_m) * m_dt/m_mu0 );

                // 注意：严格按 MATLAB 写法 —— 分母用 dz
                double denom_e = m_dz * kap_e * (sig_e + alp_e*kap_e);
                double denom_m = m_dz * kap_m * (sig_m + alp_m*kap_m);
                a_e_xn[j] = (std::abs(denom_e) < 1e-30) ? 0.0 : (sig_e*(b_e_xn[j]-1.0)/denom_e);
                a_m_xn[j] = (std::abs(denom_m) < 1e-30) ? 0.0 : (sig_m*(b_m_xn[j]-1.0)/denom_m);

                C_Dy_dxn_v[j] = (m_dt/m_dx) / kap_e;
                C_Bz_dxn_v[j] = (m_dt/m_dx) / kap_m;
            }
        }
    }
    // --- x+ ---
    {
        int n = n_pml_xp;
        if (n > 0) {
            auto rho_e = rho_pos(n, 0.25);
            auto rho_m = rho_pos(n, 0.75);
            double sigma_max = sigma_factor * (m_pml + 1) / (150.0*M_PI*std::sqrt(eps_r_pml*mu_r_pml)*m_dx);
            b_e_xp.resize(n); a_e_xp.resize(n); C_Dy_dxp_v.resize(n);
            b_m_xp.resize(n); a_m_xp.resize(n); C_Bz_dxp_v.resize(n);
            for (int j = 0; j < n; ++j) {
                double sig_e = sigma_max * std::pow(rho_e[j], m_pml);
                double sig_m = (m_mu0/m_eps0) * sigma_max * std::pow(rho_m[j], m_pml);
                double kap_e = 1.0 + (kappa_max - 1.0)*std::pow(rho_e[j], m_pml);
                double kap_m = 1.0 + (kappa_max - 1.0)*std::pow(rho_m[j], m_pml);
                double alp_e = alpha_max * (1.0 - rho_e[j]);
                double alp_m = (m_mu0/m_eps0) * alpha_max * (1.0 - rho_m[j]);

                b_e_xp[j] = std::exp( -(sig_e/kap_e + alp_e) * m_dt/m_eps0 );
                b_m_xp[j] = std::exp( -(sig_m/kap_m + alp_m) * m_dt/m_mu0 );

                // 严格按 MATLAB：分母用 dz
                double denom_e = m_dz * kap_e * (sig_e + alp_e*kap_e);
                double denom_m = m_dz * kap_m * (sig_m + alp_m*kap_m);
                a_e_xp[j] = (std::abs(denom_e) < 1e-30) ? 0.0 : (sig_e*(b_e_xp[j]-1.0)/denom_e);
                a_m_xp[j] = (std::abs(denom_m) < 1e-30) ? 0.0 : (sig_m*(b_m_xp[j]-1.0)/denom_m);

                C_Dy_dxp_v[j] = (m_dt/m_dx) / kap_e;
                C_Bz_dxp_v[j] = (m_dt/m_dx) / kap_m;
            }
        }
    }
    // --- z- ---
    {
        int n = n_pml_zn;
        if (n > 0) {
            auto rho_e = rho_neg(n, 0.75);
            auto rho_m = rho_neg(n, 0.25);
            double sigma_max = sigma_factor * (m_pml + 1) / (150.0*M_PI*std::sqrt(eps_r_pml*mu_r_pml)*m_dz);
            b_e_zn.resize(n); a_e_zn.resize(n); C_Dy_dzn_v.resize(n);
            b_m_zn.resize(n); a_m_zn.resize(n); C_Bx_dzn_v.resize(n);
            for (int j = 0; j < n; ++j) {
                double sig_e = sigma_max * std::pow(rho_e[j], m_pml);
                double sig_m = (m_mu0/m_eps0) * sigma_max * std::pow(rho_m[j], m_pml);
                double kap_e = 1.0 + (kappa_max - 1.0)*std::pow(rho_e[j], m_pml);
                double kap_m = 1.0 + (kappa_max - 1.0)*std::pow(rho_m[j], m_pml);
                double alp_e = alpha_max * (1.0 - rho_e[j]);
                double alp_m = (m_mu0/m_eps0) * alpha_max * (1.0 - rho_m[j]);

                b_e_zn[j] = std::exp( -(sig_e/kap_e + alp_e) * m_dt/m_eps0 );
                b_m_zn[j] = std::exp( -(sig_m/kap_m + alp_m) * m_dt/m_mu0 );

                double denom_e = m_dz * kap_e * (sig_e + alp_e*kap_e);
                double denom_m = m_dz * kap_m * (sig_m + alp_m*kap_m);
                a_e_zn[j] = (std::abs(denom_e) < 1e-30) ? 0.0 : (sig_e*(b_e_zn[j]-1.0)/denom_e);
                a_m_zn[j] = (std::abs(denom_m) < 1e-30) ? 0.0 : (sig_m*(b_m_zn[j]-1.0)/denom_m);

                C_Dy_dzn_v[j] = (m_dt/m_dz) / kap_e;
                C_Bx_dzn_v[j] = (m_dt/m_dz) / kap_m;
            }
        }
    }
    // --- z+ ---
    {
        int n = n_pml_zp;
        if (n > 0) {
            auto rho_e = rho_pos(n, 0.25);
            auto rho_m = rho_pos(n, 0.75);
            double sigma_max = sigma_factor * (m_pml + 1) / (150.0*M_PI*std::sqrt(eps_r_pml*mu_r_pml)*m_dz);
            b_e_zp.resize(n); a_e_zp.resize(n); C_Dy_dzp_v.resize(n);
            b_m_zp.resize(n); a_m_zp.resize(n); C_Bx_dzp_v.resize(n);
            for (int j = 0; j < n; ++j) {
                double sig_e = sigma_max * std::pow(rho_e[j], m_pml);
                double sig_m = (m_mu0/m_eps0) * sigma_max * std::pow(rho_m[j], m_pml);
                double kap_e = 1.0 + (kappa_max - 1.0)*std::pow(rho_e[j], m_pml);
                double kap_m = 1.0 + (kappa_max - 1.0)*std::pow(rho_m[j], m_pml);
                double alp_e = alpha_max * (1.0 - rho_e[j]);
                double alp_m = (m_mu0/m_eps0) * alpha_max * (1.0 - rho_m[j]);

                b_e_zp[j] = std::exp( -(sig_e/kap_e + alp_e) * m_dt/m_eps0 );
                b_m_zp[j] = std::exp( -(sig_m/kap_m + alp_m) * m_dt/m_mu0 );

                double denom_e = m_dz * kap_e * (sig_e + alp_e*kap_e);
                double denom_m = m_dz * kap_m * (sig_m + alp_m*kap_m);
                a_e_zp[j] = (std::abs(denom_e) < 1e-30) ? 0.0 : (sig_e*(b_e_zp[j]-1.0)/denom_e);
                a_m_zp[j] = (std::abs(denom_m) < 1e-30) ? 0.0 : (sig_m*(b_m_zp[j]-1.0)/denom_m);

                C_Dy_dzp_v[j] = (m_dt/m_dz) / kap_e;
                C_Bx_dzp_v[j] = (m_dt/m_dz) / kap_m;
            }
        }
    }

    // 分配 ϕ（全 0）
    Phi_ey_zn.assign( std::max(0,m_Nx-1) * n_pml_zn, 0.0 );
    Phi_ey_zp.assign( std::max(0,m_Nx-1) * n_pml_zp, 0.0 );
    Phi_mx_zn.assign( std::max(0,m_Nx-1) * n_pml_zn, 0.0 );
    Phi_mx_zp.assign( std::max(0,m_Nx-1) * n_pml_zp, 0.0 );

    Phi_ey_xn.assign( n_pml_xn * std::max(0,m_Nz-1), 0.0 );
    Phi_ey_xp.assign( n_pml_xp * std::max(0,m_Nz-1), 0.0 );
    Phi_mz_xn.assign( n_pml_xn * std::max(0,m_Nz-1), 0.0 );
    Phi_mz_xp.assign( n_pml_xp * std::max(0,m_Nz-1), 0.0 );
}

void FDTD2D::stepOnce(int n) {
    const double CD_dt = m_dt;
    const double CD_dt_dx = m_dt / m_dx;
    const double CD_dt_dz = m_dt / m_dz;

    const double CB_dt = m_dt;
    const double CB_dt_dx = m_dt / m_dx;
    const double CB_dt_dz = m_dt / m_dz;

    // 快速别名
    const int Nx = m_Nx, Nz = m_Nz;

    // ---------- Dy 更新 ----------
    // z- PML: knD = 2:(n_zn+1) => k = 1..n_zn
    for (int i = 1; i <= Nx-1; ++i) {
        int ix0 = i - 1; // for Phi_? arrays width Nx-1
        for (int l = 0; l < n_pml_zn; ++l) {
            int k = 1 + l;
            if (k-1 < 0 || k >= Nz) continue; // Hx 索引安全
            double curlH = Hx[id_Hx(i, k)] - Hx[id_Hx(i, k-1)];
            int pid = id_x_pml(ix0, l, n_pml_zn);
            Phi_ey_zn[pid] = b_e_zn[l]*Phi_ey_zn[pid] + a_e_zn[l]*curlH;
            Dy[id_Dy(i, k)] += C_Dy_dzn_v[l]*curlH + CD_dt*Phi_ey_zn[pid];
        }
    }
    // z+ PML: kpD = (Nz - n_zp + 1):Nz => k = Nz - n_zp .. Nz-1
    for (int i = 1; i <= Nx-1; ++i) {
        int ix0 = i - 1;
        for (int l = 0; l < n_pml_zp; ++l) {
            int k = Nz - n_pml_zp + l;
            if (k-1 < 0 || k >= Nz) continue;
            double curlH = Hx[id_Hx(i, k)] - Hx[id_Hx(i, k-1)];
            int pid = id_x_pml(ix0, l, n_pml_zp);
            Phi_ey_zp[pid] = b_e_zp[l]*Phi_ey_zp[pid] + a_e_zp[l]*curlH;
            Dy[id_Dy(i, k)] += C_Dy_dzp_v[l]*curlH + CD_dt*Phi_ey_zp[pid];
        }
    }
    // 中间（z）
    for (int i = 1; i <= Nx-1; ++i) {
        int k_begin = std::max(1, n_pml_zn+1);
        int k_end   = std::min(Nz-1, Nz - n_pml_zp - 1);
        if (k_begin <= k_end) {
            for (int k = k_begin; k <= k_end; ++k) {
                double curlH = Hx[id_Hx(i, k)] - Hx[id_Hx(i, k-1)];
                Dy[id_Dy(i, k)] += CD_dt_dz * curlH;
            }
        }
    }
    // x- PML: inD = 2:(n_xn+1) => i = 1..n_xn
    for (int k = 1; k <= Nz-1; ++k) {
        int kz0 = k - 1; // for Phi_? arrays height Nz-1
        for (int l = 0; l < n_pml_xn; ++l) {
            int i = 1 + l;
            if (i-1 < 0 || i >= Nx) continue; // Hz 索引安全（Hz: 0..Nx-1, 0..Nz）
            double curlH = Hz[id_Hz(i, k)] - Hz[id_Hz(i-1, k)];
            int pid = id_z_pml(l, kz0, Nz-1);
            Phi_ey_xn[pid] = b_e_xn[l]*Phi_ey_xn[pid] + a_e_xn[l]*curlH;
            Dy[id_Dy(i, k)] -= C_Dy_dxn_v[l]*curlH + CD_dt*Phi_ey_xn[pid];
        }
    }
    // x+ PML: ipD = (Nx - n_xp + 1):Nx => i = Nx - n_xp .. Nx-1
    for (int k = 1; k <= Nz-1; ++k) {
        int kz0 = k - 1;
        for (int l = 0; l < n_pml_xp; ++l) {
            int i = Nx - n_pml_xp + l;
            if (i-1 < 0 || i >= Nx) continue;
            double curlH = Hz[id_Hz(i, k)] - Hz[id_Hz(i-1, k)];
            int pid = id_z_pml(l, kz0, Nz-1);
            Phi_ey_xp[pid] = b_e_xp[l]*Phi_ey_xp[pid] + a_e_xp[l]*curlH;
            Dy[id_Dy(i, k)] -= C_Dy_dxp_v[l]*curlH + CD_dt*Phi_ey_xp[pid];
        }
    }
    // 中间（x）
    for (int k = 1; k <= Nz-1; ++k) {
        int i_begin = std::max(1, n_pml_xn+1);
        int i_end   = std::min(Nx-1, Nx - n_pml_xp - 1);
        if (i_begin <= i_end) {
            for (int i = i_begin; i <= i_end; ++i) {
                double curlH = Hz[id_Hz(i, k)] - Hz[id_Hz(i-1, k)];
                Dy[id_Dy(i, k)] -= CD_dt_dx * curlH;
            }
        }
    }

    // 源：Dy(i_Jy, k_Jy) -= dt * Jy(nt)
    if (i_Jy >= 0 && i_Jy <= Nx && k_Jy >= 0 && k_Jy <= Nz) {
        Dy[id_Dy(i_Jy, k_Jy)] -= m_dt * Jy[n];
    }

    // ---------- Ey 材料更新（i1=1..Nx, k1=1..Nz） ----------
    for (int i = 0; i <= Nx-1; ++i) {
        for (int k = 0; k <= Nz-1; ++k) {
            int idc = id_C(i, k);
            int ide = id_Ey(i, k);
            double num = (Dy[ide]/m_eps0) - 2.0*C1[idc]*Ey_n1[ide] - C1[idc]*Ey_n2[ide]
                         + C2[idc]*SEy_n1[ide] + C3[idc]*SEy_n2[ide];
            Ey[ide] = num / (1.0 + C1[idc]);
            SEy[ide] = C1[idc]*Ey[ide] + 2.0*C1[idc]*Ey_n1[ide] + C1[idc]*Ey_n2[ide]
                       - C2[idc]*SEy_n1[ide] - C3[idc]*SEy_n2[ide];
        }
    }
    Ey_n2.swap(Ey_n1);
    Ey_n1 = Ey;
    SEy_n2.swap(SEy_n1);
    SEy_n1 = SEy;

    // ---------- Bx 更新（i2=2..Nx => i=1..Nx-1） ----------
    // z-
    for (int i = 1; i <= Nx-1; ++i) {
        int ix0 = i - 1;
        for (int l = 0; l < n_pml_zn; ++l) {
            int k = 0 + l; // knB = 1..n_zn => 0..n-1
            if (k+1 > Nz) continue; // Ey 索引安全（Ey: Nz+1）
            double curlE = Ey[id_Ey(i, k+1)] - Ey[id_Ey(i, k)];
            int pid = id_x_pml(ix0, l, n_pml_zn);
            Phi_mx_zn[pid] = b_m_zn[l]*Phi_mx_zn[pid] + a_m_zn[l]*curlE;
            Bx[id_Bx(i, k)] += C_Bx_dzn_v[l]*curlE + CB_dt*Phi_mx_zn[pid];
        }
    }
    // z+
    for (int i = 1; i <= Nx-1; ++i) {
        int ix0 = i - 1;
        for (int l = 0; l < n_pml_zp; ++l) {
            int k = Nz - n_pml_zp + l;
            if (k+1 > Nz) continue;
            double curlE = Ey[id_Ey(i, k+1)] - Ey[id_Ey(i, k)];
            int pid = id_x_pml(ix0, l, n_pml_zp);
            Phi_mx_zp[pid] = b_m_zp[l]*Phi_mx_zp[pid] + a_m_zp[l]*curlE;
            Bx[id_Bx(i, k)] += C_Bx_dzp_v[l]*curlE + CB_dt*Phi_mx_zp[pid];
        }
    }
    // 中间（z）
    for (int i = 1; i <= Nx-1; ++i) {
        int k_begin = std::max(0, n_pml_zn);
        int k_end   = std::min(Nz-1, Nz - n_pml_zp - 1);
        if (k_begin <= k_end) {
            for (int k = k_begin; k <= k_end; ++k) {
                double curlE = Ey[id_Ey(i, k+1)] - Ey[id_Ey(i, k)];
                Bx[id_Bx(i, k)] += CB_dt_dz * curlE;
            }
        }
    }

    // ---------- Bz 更新 ----------
    // x-
    for (int k = 1; k <= Nz-1; ++k) {
        int kz0 = k - 1;
        for (int l = 0; l < n_pml_xn; ++l) {
            int i = 0 + l; // inB = 1..n_xn => 0..n-1
            if (i+1 > Nx) continue; // Ey 索引安全（Ey: Nx+1）
            double curlE = Ey[id_Ey(i+1, k)] - Ey[id_Ey(i, k)];
            int pid = id_z_pml(l, kz0, Nz-1);
            Phi_mz_xn[pid] = b_m_xn[l]*Phi_mz_xn[pid] + a_m_xn[l]*curlE;
            Bz[id_Bz(i, k)] -= C_Bz_dxn_v[l]*curlE + CB_dt*Phi_mz_xn[pid];
        }
    }
    // x+
    for (int k = 1; k <= Nz-1; ++k) {
        int kz0 = k - 1;
        for (int l = 0; l < n_pml_xp; ++l) {
            int i = Nx - n_pml_xp + l;
            if (i+1 > Nx) continue;
            double curlE = Ey[id_Ey(i+1, k)] - Ey[id_Ey(i, k)];
            int pid = id_z_pml(l, kz0, Nz-1);
            Phi_mz_xp[pid] = b_m_xp[l]*Phi_mz_xp[pid] + a_m_xp[l]*curlE;
            Bz[id_Bz(i, k)] -= C_Bz_dxp_v[l]*curlE + CB_dt*Phi_mz_xp[pid];
        }
    }
    // 中间（x）
    for (int k = 1; k <= Nz-1; ++k) {
        int i_begin = std::max(0, n_pml_xn);
        int i_end   = std::min(Nx-1, Nx - n_pml_xp - 1);
        if (i_begin <= i_end) {
            for (int i = i_begin; i <= i_end; ++i) {
                double curlE = Ey[id_Ey(i+1, k)] - Ey[id_Ey(i, k)];
                Bz[id_Bz(i, k)] -= CB_dt_dx * curlE;
            }
        }
    }

    // ---------- Hx/Hz 材料更新（i1=1..Nx, k1=1..Nz） ----------
    for (int i = 0; i <= Nx-1; ++i) {
        for (int k = 0; k <= Nz-1; ++k) {
            int idc = id_C(i, k);

            // Hx: (Nx+1) x Nz —— 更新 i=0..Nx-1, k=0..Nz-1
            {
                int idh = id_Hx(i, k);
                double num = (Bx[idh]/m_mu0) - 2.0*C1[idc]*Hx_n1[idh] - C1[idc]*Hx_n2[idh]
                             + C2[idc]*SHx_n1[idh] + C3[idc]*SHx_n2[idh];
                Hx[idh] = num / (1.0 + C1[idc]);
                SHx[idh] = C1[idc]*Hx[idh] + 2.0*C1[idc]*Hx_n1[idh] + C1[idc]*Hx_n2[idh]
                           - C2[idc]*SHx_n1[idh] - C3[idc]*SHx_n2[idh];
            }
            // Hz: Nx x (Nz+1) —— 更新 i=0..Nx-1, k=0..Nz-1
            {
                int idh = id_Hz(i, k);
                double num = (Bz[idh]/m_mu0) - 2.0*C1[idc]*Hz_n1[idh] - C1[idc]*Hz_n2[idh]
                             + C2[idc]*SHz_n1[idh] + C3[idc]*SHz_n2[idh];
                Hz[idh] = num / (1.0 + C1[idc]);
                SHz[idh] = C1[idc]*Hz[idh] + 2.0*C1[idc]*Hz_n1[idh] + C1[idc]*Hz_n2[idh]
                           - C2[idc]*SHz_n1[idh] - C3[idc]*SHz_n2[idh];
            }
        }
    }
    Hx_n2.swap(Hx_n1); Hx_n1 = Hx; SHx_n2.swap(SHx_n1); SHx_n1 = SHx;
    Hz_n2.swap(Hz_n1); Hz_n1 = Hz; SHz_n2.swap(SHz_n1); SHz_n1 = SHz;
}

void FDTD2D::dumpFrame(int n) {
    // 为了减小体积，只导出 Ey(i=0..Nx-1, k=0..Nz-1) 的中心区域（与 C1 尺寸一致）
    // 方便可视化：保存为 Nx×Nz 的“表面”场
    std::vector<double> ey_center(m_Nx*m_Nz, 0.0);
    for (int i = 0; i <= m_Nx-1; ++i) {
        for (int k = 0; k <= m_Nz-1; ++k) {
            ey_center[i*m_Nz + k] = Ey[id_Ey(i, k)];
        }
    }
    std::ostringstream oss;
    oss << "ey_t" << std::setw(5) << std::setfill('0') << n << ".csv";
    m_out->saveScalar2D_CSV(oss.str(), ey_center, m_Nx, m_Nz, m_dx, m_dz);

    // 仅首帧保存尺寸
    if (n == 0) m_out->saveDimsTxt2D(m_Nx, m_Nz, m_dx, m_dz);
}

void FDTD2D::run() {
    // 计算半周期 RMS（与 MATLAB 的 Ey_mag 相同思想）
    int N_lambda = 40;
    double nCFL = 2.0;
    int n_T0 = static_cast<int>(std::round(nCFL * N_lambda));
    int n_T0d2 = n_T0 / 2;

    std::vector<double> Ey_mag_sq(m_Nx*m_Nz, 0.0);
    std::vector<double> Ey_mag(m_Nx*m_Nz, 0.0);

    for (int n = 0; n < m_steps; ++n) {
        stepOnce(n);

        // RMS 统计（中心 Nx×Nz）
        if (n % n_T0d2 == 0) std::fill(Ey_mag_sq.begin(), Ey_mag_sq.end(), 0.0);
        for (int i = 0; i <= m_Nx-1; ++i) {
            for (int k = 0; k <= m_Nz-1; ++k) {
                Ey_mag_sq[i*m_Nz + k] += Ey[id_Ey(i,k)]*Ey[id_Ey(i,k)] * (4.0 / n_T0);
            }
        }
        if ((n+1) % n_T0d2 == 0) {
            for (size_t p = 0; p < Ey_mag.size(); ++p) Ey_mag[p] = std::sqrt(Ey_mag_sq[p]);
            std::ostringstream oss;
            oss << "eymag_t" << std::setw(5) << std::setfill('0') << n+1 << ".csv";
            m_out->saveScalar2D_CSV(oss.str(), Ey_mag, m_Nx, m_Nz, m_dx, m_dz);
        }

        if (n % 20 == 0) {
            dumpFrame(n);
            std::cout << "[n=" << n << "] frame dumped.\n";
        }
    }
    dumpFrame(m_steps);
    std::cout << "Done. CSV in ./output_2d\n";
}
