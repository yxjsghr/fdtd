#include "FDTD3D.h"
#include "OutputManager.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>

FDTD3D::FDTD3D() {
    setupGrid();
    setupSource();
    setupMaterials();
    setupPML();
    m_out = std::make_unique<OutputManager>("output_3d");
}

FDTD3D::~FDTD3D() = default;

void FDTD3D::setupGrid() {
    // ===== 对齐 MATLAB 参数 =====
    double lambda0 = 1053e-9; // 1053 nm
    int N_lambda = 20;

    m_dx = m_dy = m_dz = lambda0 / N_lambda;

    double Lx = 10.0 * lambda0;
    double Ly = 10.0 * lambda0;
    double Lz = 10.0 * lambda0;

    m_Nx = static_cast<int>(std::round(Lx / m_dx));
    m_Ny = static_cast<int>(std::round(Ly / m_dy));
    m_Nz = static_cast<int>(std::round(Lz / m_dz));

    double nCFL = 2.0;
    m_dt = m_dx / (m_c0 * nCFL);

    m_steps = 101;

    // 分配内存（严格 Yee 尺寸）
    Dx.assign(m_Nx * (m_Ny + 1) * (m_Nz + 1), 0.0);
    Ex = Dx;

    Dy.assign((m_Nx + 1) * m_Ny * (m_Nz + 1), 0.0);
    Ey = Dy;

    Dz.assign((m_Nx + 1) * (m_Ny + 1) * m_Nz, 0.0);
    Ez = Dz;

    Bx.assign((m_Nx + 1) * m_Ny * m_Nz, 0.0);
    Hx = Bx;

    By.assign(m_Nx * (m_Ny + 1) * m_Nz, 0.0);
    Hy = By;

    Bz.assign(m_Nx * m_Ny * (m_Nz + 1), 0.0);
    Hz = Bz;

    // 材料：定义在 (Nx, Ny, Nz) —— Yee 元胞中心
    m_eps_x.assign(m_Nx * (m_Ny + 1) * (m_Nz + 1), 0.0);
    m_eps_y.assign((m_Nx + 1) * m_Ny * (m_Nz + 1), 0.0);
    m_eps_z.assign((m_Nx + 1) * (m_Ny + 1) * m_Nz, 0.0);

    std::cout << "[grid] Nx=" << m_Nx << ", Ny=" << m_Ny << ", Nz=" << m_Nz
              << ", dx=" << m_dx << ", dt=" << m_dt << ", steps=" << m_steps << "\n";
}

void FDTD3D::setupSource() {
    // 高斯面电流源 Jx at z = k_Jx
    m_k_Jx = 20; // MATLAB: k_Jx = 21 (1-based)

    // 时间波形
    double f0 = m_c0 / 1053e-9;
    double omega0 = 2.0 * M_PI * f0;
    double T0 = 1.0 / f0;
    m_wf_Jx_i.resize(m_steps);
    double E0 = 1.0;
    for (int n = 0; n < m_steps; ++n) {
        double t = (n + 0.5) * m_dt;
        m_wf_Jx_i[n] = (2.0 * E0 / m_Z0) * std::sin(omega0 * t) / m_dz;
    }
    // 3-4-5 渐升
    double ts = 5.0 * T0;
    int tsddt = static_cast<int>(std::round(ts / m_dt));
    for (int n = 0; n <= tsddt && n < m_steps; ++n) {
        double s = static_cast<double>(n) / tsddt;
        double ws = 10.0 * std::pow(s, 3) - 15.0 * std::pow(s, 4) + 6.0 * std::pow(s, 5);
        m_wf_Jx_i[n] *= ws;
    }

    // 空间高斯分布: Jx_i_mag[i][j] for i=0..Nx-1, j=0..Ny
    m_Jx_i_mag.assign(m_Nx * (m_Ny + 1), 0.0);
    double w0 = 1.5 * 1053e-9; // 1.5 lambda0
    int Nx_c = m_Nx / 2;
    int Ny_c = m_Ny / 2;
    for (int i = 0; i < m_Nx; ++i) {
        for (int j = 0; j <= m_Ny; ++j) {
            double x = (i - Nx_c) * m_dx;
            double y = (j - Ny_c) * m_dy;
            m_Jx_i_mag[i * (m_Ny + 1) + j] = E0 * std::exp(-(x * x + y * y) / (w0 * w0));
        }
    }
}

void FDTD3D::setupMaterials() {
    // 背景：熔石英
    double n_fsg = 1.47;
    double eps_r_fsg = n_fsg * n_fsg;
    double eps_fsg = eps_r_fsg * m_eps0;
    std::fill(m_eps_x.begin(), m_eps_x.end(), eps_fsg);
    std::fill(m_eps_y.begin(), m_eps_y.end(), eps_fsg);
    std::fill(m_eps_z.begin(), m_eps_z.end(), eps_fsg);

    // 气泡：中心球形空气
    double n_bub = 1.0;
    double eps_r_bub = n_bub * n_bub;
    double R_bub = 2.0 * 1053e-9;
    double eps_bub = eps_r_bub * m_eps0;

    int Nx_sc = m_Nx / 2;
    int Ny_sc = m_Ny / 2;
    int Nz_sc = m_Nz / 2;

    auto mix_eps = [&](double fraction) {
        return eps_fsg + (eps_bub - eps_fsg) * fraction;
    };

    // eps_x: faces normal to x, located on (i+0.5, j, k)
    for (int i = 0; i < m_Nx; ++i) {
        double x_face = ((i + 0.5) - Nx_sc) * m_dx;
        for (int j = 0; j <= m_Ny; ++j) {
            double y_node = (j - Ny_sc) * m_dy;
            for (int k = 0; k <= m_Nz; ++k) {
                double z_node = (k - Nz_sc) * m_dz;

                int inside_count = 0;
                for (int jj = -1; jj <= 1; ++jj) {
                    double y = y_node + jj * (m_dy * 0.5);
                    for (int kk = -1; kk <= 1; ++kk) {
                        double z = z_node + kk * (m_dz * 0.5);
                        double r = std::sqrt(x_face * x_face + y * y + z * z);
                        if (r <= R_bub) ++inside_count;
                    }
                }

                double fraction = static_cast<double>(inside_count) / 9.0;
                m_eps_x[id_eps_x(i, j, k)] = mix_eps(fraction);
            }
        }
    }

    // eps_y: faces normal to y, located on (i, j+0.5, k)
    for (int i = 0; i <= m_Nx; ++i) {
        double x_node = (i - Nx_sc) * m_dx;
        for (int j = 0; j < m_Ny; ++j) {
            double y_face = ((j + 0.5) - Ny_sc) * m_dy;
            for (int k = 0; k <= m_Nz; ++k) {
                double z_node = (k - Nz_sc) * m_dz;

                int inside_count = 0;
                for (int ii = -1; ii <= 1; ++ii) {
                    double x = x_node + ii * (m_dx * 0.5);
                    for (int kk = -1; kk <= 1; ++kk) {
                        double z = z_node + kk * (m_dz * 0.5);
                        double r = std::sqrt(x * x + y_face * y_face + z * z);
                        if (r <= R_bub) ++inside_count;
                    }
                }

                double fraction = static_cast<double>(inside_count) / 9.0;
                m_eps_y[id_eps_y(i, j, k)] = mix_eps(fraction);
            }
        }
    }

    // eps_z: faces normal to z, located on (i, j, k+0.5)
    for (int i = 0; i <= m_Nx; ++i) {
        double x_node = (i - Nx_sc) * m_dx;
        for (int j = 0; j <= m_Ny; ++j) {
            double y_node = (j - Ny_sc) * m_dy;
            for (int k = 0; k < m_Nz; ++k) {
                double z_face = ((k + 0.5) - Nz_sc) * m_dz;

                int inside_count = 0;
                for (int ii = -1; ii <= 1; ++ii) {
                    double x = x_node + ii * (m_dx * 0.5);
                    for (int jj = -1; jj <= 1; ++jj) {
                        double y = y_node + jj * (m_dy * 0.5);
                        double r = std::sqrt(x * x + y * y + z_face * z_face);
                        if (r <= R_bub) ++inside_count;
                    }
                }

                double fraction = static_cast<double>(inside_count) / 9.0;
                m_eps_z[id_eps_z(i, j, k)] = mix_eps(fraction);
            }
        }
    }
}

// ===== CPML 辅助函数 =====
static std::vector<double> rho_neg(int n, double offset) {
    std::vector<double> r(n);
    for (int j = 0; j < n; ++j) {
        r[j] = ((n - offset) - j) / n;
    }
    return r;
}
static std::vector<double> rho_pos(int n, double offset) {
    std::vector<double> r(n);
    for (int j = 0; j < n; ++j) {
        r[j] = (offset + j) / n;
    }
    return r;
}

void FDTD3D::setupPML() {
    const int m_pml = 3;
    const double alpha_max = 0.05;
    const int sigma_factor = 1;
    const double kappa_max = 7.0;
    const double eps_r_pml = 1.0;
    const double mu_r_pml = 1.0;

    // --- x- ---
    {
        int n = n_pml_xn;
        if (n > 0) {
            auto rho_e = rho_neg(n, 0.75);
            auto rho_m = rho_neg(n, 0.25);
            double sigma_max = sigma_factor * (m_pml + 1) / (150.0 * M_PI * std::sqrt(eps_r_pml * mu_r_pml) * m_dx);
            b_e_xn.resize(n); a_e_xn.resize(n); C_Dy_dxn_v.resize(n); C_Dz_dxn_v.resize(n);
            b_m_xn.resize(n); a_m_xn.resize(n); C_By_dxn_v.resize(n); C_Bz_dxn_v.resize(n);
            for (int j = 0; j < n; ++j) {
                double sig_e = sigma_max * std::pow(rho_e[j], m_pml);
                double sig_m = (m_mu0 / m_eps0) * sigma_max * std::pow(rho_m[j], m_pml);
                double kap_e = 1.0 + (kappa_max - 1.0) * std::pow(rho_e[j], m_pml);
                double kap_m = 1.0 + (kappa_max - 1.0) * std::pow(rho_m[j], m_pml);
                double alp_e = alpha_max * (1.0 - rho_e[j]);
                double alp_m = (m_mu0 / m_eps0) * alpha_max * (1.0 - rho_m[j]);

                b_e_xn[j] = std::exp(-(sig_e / kap_e + alp_e) * m_dt / m_eps0);
                b_m_xn[j] = std::exp(-(sig_m / kap_m + alp_m) * m_dt / m_mu0);

                double denom_e = m_dx * kap_e * (sig_e + alp_e * kap_e);
                double denom_m = m_dx * kap_m * (sig_m + alp_m * kap_m);
                a_e_xn[j] = (std::abs(denom_e) < 1e-30) ? 0.0 : (sig_e * (b_e_xn[j] - 1.0) / denom_e);
                a_m_xn[j] = (std::abs(denom_m) < 1e-30) ? 0.0 : (sig_m * (b_m_xn[j] - 1.0) / denom_m);

                C_Dy_dxn_v[j] = (m_dt / m_dx) / kap_e;
                C_Dz_dxn_v[j] = (m_dt / m_dx) / kap_e;
                C_By_dxn_v[j] = (m_dt / m_dx) / kap_m;
                C_Bz_dxn_v[j] = (m_dt / m_dx) / kap_m;
            }
        }
    }
    // --- x+ ---
    {
        int n = n_pml_xp;
        if (n > 0) {
            auto rho_e = rho_pos(n, 0.25);
            auto rho_m = rho_pos(n, 0.75);
            double sigma_max = sigma_factor * (m_pml + 1) / (150.0 * M_PI * std::sqrt(eps_r_pml * mu_r_pml) * m_dx);
            b_e_xp.resize(n); a_e_xp.resize(n); C_Dy_dxp_v.resize(n); C_Dz_dxp_v.resize(n);
            b_m_xp.resize(n); a_m_xp.resize(n); C_By_dxp_v.resize(n); C_Bz_dxp_v.resize(n);
            for (int j = 0; j < n; ++j) {
                double sig_e = sigma_max * std::pow(rho_e[j], m_pml);
                double sig_m = (m_mu0 / m_eps0) * sigma_max * std::pow(rho_m[j], m_pml);
                double kap_e = 1.0 + (kappa_max - 1.0) * std::pow(rho_e[j], m_pml);
                double kap_m = 1.0 + (kappa_max - 1.0) * std::pow(rho_m[j], m_pml);
                double alp_e = alpha_max * (1.0 - rho_e[j]);
                double alp_m = (m_mu0 / m_eps0) * alpha_max * (1.0 - rho_m[j]);

                b_e_xp[j] = std::exp(-(sig_e / kap_e + alp_e) * m_dt / m_eps0);
                b_m_xp[j] = std::exp(-(sig_m / kap_m + alp_m) * m_dt / m_mu0);

                double denom_e = m_dx * kap_e * (sig_e + alp_e * kap_e);
                double denom_m = m_dx * kap_m * (sig_m + alp_m * kap_m);
                a_e_xp[j] = (std::abs(denom_e) < 1e-30) ? 0.0 : (sig_e * (b_e_xp[j] - 1.0) / denom_e);
                a_m_xp[j] = (std::abs(denom_m) < 1e-30) ? 0.0 : (sig_m * (b_m_xp[j] - 1.0) / denom_m);

                C_Dy_dxp_v[j] = (m_dt / m_dx) / kap_e;
                C_Dz_dxp_v[j] = (m_dt / m_dx) / kap_e;
                C_By_dxp_v[j] = (m_dt / m_dx) / kap_m;
                C_Bz_dxp_v[j] = (m_dt / m_dx) / kap_m;
            }
        }
    }
    // --- y- ---
    {
        int n = n_pml_yn;
        if (n > 0) {
            auto rho_e = rho_neg(n, 0.75);
            auto rho_m = rho_neg(n, 0.25);
            double sigma_max = sigma_factor * (m_pml + 1) / (150.0 * M_PI * std::sqrt(eps_r_pml * mu_r_pml) * m_dy);
            b_e_yn.resize(n); a_e_yn.resize(n); C_Dx_dyn_v.resize(n); C_Dz_dyn_v.resize(n);
            b_m_yn.resize(n); a_m_yn.resize(n); C_Bx_dyn_v.resize(n); C_Bz_dyn_v.resize(n);
            for (int j = 0; j < n; ++j) {
                double sig_e = sigma_max * std::pow(rho_e[j], m_pml);
                double sig_m = (m_mu0 / m_eps0) * sigma_max * std::pow(rho_m[j], m_pml);
                double kap_e = 1.0 + (kappa_max - 1.0) * std::pow(rho_e[j], m_pml);
                double kap_m = 1.0 + (kappa_max - 1.0) * std::pow(rho_m[j], m_pml);
                double alp_e = alpha_max * (1.0 - rho_e[j]);
                double alp_m = (m_mu0 / m_eps0) * alpha_max * (1.0 - rho_m[j]);

                b_e_yn[j] = std::exp(-(sig_e / kap_e + alp_e) * m_dt / m_eps0);
                b_m_yn[j] = std::exp(-(sig_m / kap_m + alp_m) * m_dt / m_mu0);

                double denom_e = m_dy * kap_e * (sig_e + alp_e * kap_e);
                double denom_m = m_dy * kap_m * (sig_m + alp_m * kap_m);
                a_e_yn[j] = (std::abs(denom_e) < 1e-30) ? 0.0 : (sig_e * (b_e_yn[j] - 1.0) / denom_e);
                a_m_yn[j] = (std::abs(denom_m) < 1e-30) ? 0.0 : (sig_m * (b_m_yn[j] - 1.0) / denom_m);

                C_Dx_dyn_v[j] = (m_dt / m_dy) / kap_e;
                C_Dz_dyn_v[j] = (m_dt / m_dy) / kap_e;
                C_Bx_dyn_v[j] = (m_dt / m_dy) / kap_m;
                C_Bz_dyn_v[j] = (m_dt / m_dy) / kap_m;
            }
        }
    }
    // --- y+ ---
    {
        int n = n_pml_yp;
        if (n > 0) {
            auto rho_e = rho_pos(n, 0.25);
            auto rho_m = rho_pos(n, 0.75);
            double sigma_max = sigma_factor * (m_pml + 1) / (150.0 * M_PI * std::sqrt(eps_r_pml * mu_r_pml) * m_dy);
            b_e_yp.resize(n); a_e_yp.resize(n); C_Dx_dyp_v.resize(n); C_Dz_dyp_v.resize(n);
            b_m_yp.resize(n); a_m_yp.resize(n); C_Bx_dyp_v.resize(n); C_Bz_dyp_v.resize(n);
            for (int j = 0; j < n; ++j) {
                double sig_e = sigma_max * std::pow(rho_e[j], m_pml);
                double sig_m = (m_mu0 / m_eps0) * sigma_max * std::pow(rho_m[j], m_pml);
                double kap_e = 1.0 + (kappa_max - 1.0) * std::pow(rho_e[j], m_pml);
                double kap_m = 1.0 + (kappa_max - 1.0) * std::pow(rho_m[j], m_pml);
                double alp_e = alpha_max * (1.0 - rho_e[j]);
                double alp_m = (m_mu0 / m_eps0) * alpha_max * (1.0 - rho_m[j]);

                b_e_yp[j] = std::exp(-(sig_e / kap_e + alp_e) * m_dt / m_eps0);
                b_m_yp[j] = std::exp(-(sig_m / kap_m + alp_m) * m_dt / m_mu0);

                double denom_e = m_dy * kap_e * (sig_e + alp_e * kap_e);
                double denom_m = m_dy * kap_m * (sig_m + alp_m * kap_m);
                a_e_yp[j] = (std::abs(denom_e) < 1e-30) ? 0.0 : (sig_e * (b_e_yp[j] - 1.0) / denom_e);
                a_m_yp[j] = (std::abs(denom_m) < 1e-30) ? 0.0 : (sig_m * (b_m_yp[j] - 1.0) / denom_m);

                C_Dx_dyp_v[j] = (m_dt / m_dy) / kap_e;
                C_Dz_dyp_v[j] = (m_dt / m_dy) / kap_e;
                C_Bx_dyp_v[j] = (m_dt / m_dy) / kap_m;
                C_Bz_dyp_v[j] = (m_dt / m_dy) / kap_m;
            }
        }
    }
    // --- z- ---
    {
        int n = n_pml_zn;
        if (n > 0) {
            auto rho_e = rho_neg(n, 0.75);
            auto rho_m = rho_neg(n, 0.25);
            double sigma_max = sigma_factor * (m_pml + 1) / (150.0 * M_PI * std::sqrt(eps_r_pml * mu_r_pml) * m_dz);
            b_e_zn.resize(n); a_e_zn.resize(n); C_Dx_dzn_v.resize(n); C_Dy_dzn_v.resize(n);
            b_m_zn.resize(n); a_m_zn.resize(n); C_Bx_dzn_v.resize(n); C_By_dzn_v.resize(n);
            for (int j = 0; j < n; ++j) {
                double sig_e = sigma_max * std::pow(rho_e[j], m_pml);
                double sig_m = (m_mu0 / m_eps0) * sigma_max * std::pow(rho_m[j], m_pml);
                double kap_e = 1.0 + (kappa_max - 1.0) * std::pow(rho_e[j], m_pml);
                double kap_m = 1.0 + (kappa_max - 1.0) * std::pow(rho_m[j], m_pml);
                double alp_e = alpha_max * (1.0 - rho_e[j]);
                double alp_m = (m_mu0 / m_eps0) * alpha_max * (1.0 - rho_m[j]);

                b_e_zn[j] = std::exp(-(sig_e / kap_e + alp_e) * m_dt / m_eps0);
                b_m_zn[j] = std::exp(-(sig_m / kap_m + alp_m) * m_dt / m_mu0);

                double denom_e = m_dz * kap_e * (sig_e + alp_e * kap_e);
                double denom_m = m_dz * kap_m * (sig_m + alp_m * kap_m);
                a_e_zn[j] = (std::abs(denom_e) < 1e-30) ? 0.0 : (sig_e * (b_e_zn[j] - 1.0) / denom_e);
                a_m_zn[j] = (std::abs(denom_m) < 1e-30) ? 0.0 : (sig_m * (b_m_zn[j] - 1.0) / denom_m);

                C_Dx_dzn_v[j] = (m_dt / m_dz) / kap_e;
                C_Dy_dzn_v[j] = (m_dt / m_dz) / kap_e;
                C_Bx_dzn_v[j] = (m_dt / m_dz) / kap_m;
                C_By_dzn_v[j] = (m_dt / m_dz) / kap_m;
            }
        }
    }
    // --- z+ ---
    {
        int n = n_pml_zp;
        if (n > 0) {
            auto rho_e = rho_pos(n, 0.25);
            auto rho_m = rho_pos(n, 0.75);
            double sigma_max = sigma_factor * (m_pml + 1) / (150.0 * M_PI * std::sqrt(eps_r_pml * mu_r_pml) * m_dz);
            b_e_zp.resize(n); a_e_zp.resize(n); C_Dx_dzp_v.resize(n); C_Dy_dzp_v.resize(n);
            b_m_zp.resize(n); a_m_zp.resize(n); C_Bx_dzp_v.resize(n); C_By_dzp_v.resize(n);
            for (int j = 0; j < n; ++j) {
                double sig_e = sigma_max * std::pow(rho_e[j], m_pml);
                double sig_m = (m_mu0 / m_eps0) * sigma_max * std::pow(rho_m[j], m_pml);
                double kap_e = 1.0 + (kappa_max - 1.0) * std::pow(rho_e[j], m_pml);
                double kap_m = 1.0 + (kappa_max - 1.0) * std::pow(rho_m[j], m_pml);
                double alp_e = alpha_max * (1.0 - rho_e[j]);
                double alp_m = (m_mu0 / m_eps0) * alpha_max * (1.0 - rho_m[j]);

                b_e_zp[j] = std::exp(-(sig_e / kap_e + alp_e) * m_dt / m_eps0);
                b_m_zp[j] = std::exp(-(sig_m / kap_m + alp_m) * m_dt / m_mu0);

                double denom_e = m_dz * kap_e * (sig_e + alp_e * kap_e);
                double denom_m = m_dz * kap_m * (sig_m + alp_m * kap_m);
                a_e_zp[j] = (std::abs(denom_e) < 1e-30) ? 0.0 : (sig_e * (b_e_zp[j] - 1.0) / denom_e);
                a_m_zp[j] = (std::abs(denom_m) < 1e-30) ? 0.0 : (sig_m * (b_m_zp[j] - 1.0) / denom_m);

                C_Dx_dzp_v[j] = (m_dt / m_dz) / kap_e;
                C_Dy_dzp_v[j] = (m_dt / m_dz) / kap_e;
                C_Bx_dzp_v[j] = (m_dt / m_dz) / kap_m;
                C_By_dzp_v[j] = (m_dt / m_dz) / kap_m;
            }
        }
    }

    // 分配 Φ 场（初始化为 0）
    // E 场 Φ（修正后）
    // Dx: [Nx] x [Ny+1] x [Nz+1]，y方向PML影响x-z平面
    Phi_ex_yn.assign(m_Nx * (m_Nz + 1) * n_pml_yn, 0.0);
    Phi_ex_yp.assign(m_Nx * (m_Nz + 1) * n_pml_yp, 0.0);
    // Dx: z方向PML影响x-y平面
    Phi_ex_zn.assign(m_Nx * (m_Ny + 1) * n_pml_zn, 0.0);
    Phi_ex_zp.assign(m_Nx * (m_Ny + 1) * n_pml_zp, 0.0);

    // Dy: [Nx+1] x [Ny] x [Nz+1]，x方向PML影响y-z平面
    Phi_ey_xn.assign(n_pml_xn * m_Ny * (m_Nz + 1), 0.0);
    Phi_ey_xp.assign(n_pml_xp * m_Ny * (m_Nz + 1), 0.0);
    // Dy: z方向PML影响x-y平面
    Phi_ey_zn.assign((m_Nx + 1) * m_Ny * n_pml_zn, 0.0);
    Phi_ey_zp.assign((m_Nx + 1) * m_Ny * n_pml_zp, 0.0);

    // Dz: [Nx+1] x [Ny+1] x [Nz]，x方向PML影响y-z平面
    Phi_ez_xn.assign(n_pml_xn * (m_Ny + 1) * m_Nz, 0.0);
    Phi_ez_xp.assign(n_pml_xp * (m_Ny + 1) * m_Nz, 0.0);
    // Dz: y方向PML影响x-z平面
    Phi_ez_yn.assign((m_Nx + 1) * n_pml_yn * m_Nz, 0.0);
    Phi_ez_yp.assign((m_Nx + 1) * n_pml_yp * m_Nz, 0.0);

    // B 场 Φ（修正后）
    // Bx: [Nx+1] x [Ny] x [Nz]，y方向PML影响x-z平面
    Phi_mx_yn.assign((m_Nx + 1) * n_pml_yn * m_Nz, 0.0);
    Phi_mx_yp.assign((m_Nx + 1) * n_pml_yp * m_Nz, 0.0);
    // Bx: z方向PML影响x-y平面
    Phi_mx_zn.assign((m_Nx + 1) * m_Ny * n_pml_zn, 0.0);
    Phi_mx_zp.assign((m_Nx + 1) * m_Ny * n_pml_zp, 0.0);

    // By: [Nx] x [Ny+1] x [Nz]，x方向PML影响y-z平面
    Phi_my_xn.assign(n_pml_xn * (m_Ny + 1) * m_Nz, 0.0);
    Phi_my_xp.assign(n_pml_xp * (m_Ny + 1) * m_Nz, 0.0);
    // By: z方向PML影响x-y平面
    Phi_my_zn.assign(m_Nx * (m_Ny + 1) * n_pml_zn, 0.0);
    Phi_my_zp.assign(m_Nx * (m_Ny + 1) * n_pml_zp, 0.0);

    // Bz: [Nx] x [Ny] x [Nz+1]，x方向PML影响y-z平面
    Phi_mz_xn.assign(n_pml_xn * m_Ny * (m_Nz + 1), 0.0);
    Phi_mz_xp.assign(n_pml_xp * m_Ny * (m_Nz + 1), 0.0);
    // Bz: y方向PML影响x-z平面
    Phi_mz_yn.assign(m_Nx * n_pml_yn * (m_Nz + 1), 0.0);
    Phi_mz_yp.assign(m_Nx * n_pml_yp * (m_Nz + 1), 0.0);
}

/* void FDTD3D::applySource(int step) {
    if (step >= static_cast<int>(m_wf_Jx_i.size())) return;
    double Jx_val = m_dt * m_wf_Jx_i[step];
    int i_begin = std::max(0, n_pml_xn);
    int i_end = std::min(m_Nx, m_Nx - n_pml_xp);
    int j_begin = std::max(0, n_pml_yn);
    int j_end = std::min(m_Ny + 1, m_Ny + 1 - n_pml_yp);

    for (int i = i_begin; i < i_end; ++i) {
        for (int j = j_begin; j < j_end; ++j) {
            int id = id_Dx(i, j, m_k_Jx);
            Dx[id] += Jx_val * m_Jx_i_mag[i * (m_Ny + 1) + j];
        }
    }
} */
 void FDTD3D::applySource(int step) {
    if (step >= static_cast<int>(m_wf_Jx_i.size())) return;
    double Jx_val = m_dt * m_wf_Jx_i[step];

    // 施加在整个 x-y 平面 (i=0..Nx-1, j=0..Ny)
    for (int i = 0; i < m_Nx; ++i) {
        for (int j = 0; j <= m_Ny; ++j) {
            // m_k_Jx 是源所在的 z 平面索引
            int id = id_Dx(i, j, m_k_Jx);
            Dx[id] += Jx_val * m_Jx_i_mag[i * (m_Ny + 1) + j];
        }
    }
}

void FDTD3D::updateElectricMaterialRelations() {
    for (int i = 0; i < m_Nx; ++i) {
        for (int j = 0; j <= m_Ny; ++j) {
            for (int k = 0; k <= m_Nz; ++k) {
                int idx = id_Ex(i, j, k);
                Ex[idx] = Dx[idx] / m_eps_x[id_eps_x(i, j, k)];
            }
        }
    }

    for (int i = 0; i <= m_Nx; ++i) {
        for (int j = 0; j < m_Ny; ++j) {
            for (int k = 0; k <= m_Nz; ++k) {
                int idx = id_Ey(i, j, k);
                Ey[idx] = Dy[idx] / m_eps_y[id_eps_y(i, j, k)];
            }
        }
    }

    for (int i = 0; i <= m_Nx; ++i) {
        for (int j = 0; j <= m_Ny; ++j) {
            for (int k = 0; k < m_Nz; ++k) {
                int idx = id_Ez(i, j, k);
                Ez[idx] = Dz[idx] / m_eps_z[id_eps_z(i, j, k)];
            }
        }
    }
}

void FDTD3D::updateMagneticMaterialRelations() {
    for (size_t idx = 0; idx < Hx.size(); ++idx) Hx[idx] = Bx[idx] / m_mu0;
    for (size_t idx = 0; idx < Hy.size(); ++idx) Hy[idx] = By[idx] / m_mu0;
    for (size_t idx = 0; idx < Hz.size(); ++idx) Hz[idx] = Bz[idx] / m_mu0;
}

// ===== 完整的 updateE / updateH 实现 =====
void FDTD3D::updateE(int step) {
    const double CD_dt = m_dt;
    const double CD_dt_dx = m_dt / m_dx;
    const double CD_dt_dy = m_dt / m_dy;
    const double CD_dt_dz = m_dt / m_dz;

    const int Nx = m_Nx, Ny = m_Ny, Nz = m_Nz;

    // ========== Dx 更新：curl H = dHz/dy - dHy/dz ==========
    // z- PML
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int l = 0; l < n_pml_zn; ++l) {
                int k = 1 + l;
                if (k >= Nz) continue;
                double curlH_z = Hy[id_Hy(i, j, k)] - Hy[id_Hy(i, j, k - 1)];
                int pid_z = id_Phi_ex_zn(i, j, l);
                Phi_ex_zn[pid_z] = b_e_zn[l] * Phi_ex_zn[pid_z] + a_e_zn[l] * curlH_z;
                Dx[id_Dx(i, j, k)] -= C_Dx_dzn_v[l] * curlH_z + CD_dt * Phi_ex_zn[pid_z];
            }
        }
    }
    // z+ PML
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int l = 0; l < n_pml_zp; ++l) {
                int k = Nz - n_pml_zp + l;
                if (k < 1 || k >= Nz) continue;
                double curlH_z = Hy[id_Hy(i, j, k)] - Hy[id_Hy(i, j, k - 1)];
                int pid_z = id_Phi_ex_zp(i, j, l);
                Phi_ex_zp[pid_z] = b_e_zp[l] * Phi_ex_zp[pid_z] + a_e_zp[l] * curlH_z;
                Dx[id_Dx(i, j, k)] -= C_Dx_dzp_v[l] * curlH_z + CD_dt * Phi_ex_zp[pid_z];
            }
        }
    }
    // z 中间区域
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            int k_begin = std::max(1, n_pml_zn + 1);
            int k_end = std::min(Nz, Nz - n_pml_zp);
            if (k_begin < k_end) {
                for (int k = k_begin; k < k_end; ++k) {
                    double curlH_z = Hy[id_Hy(i, j, k)] - Hy[id_Hy(i, j, k - 1)];
                    Dx[id_Dx(i, j, k)] -= CD_dt_dz * curlH_z;
                }
            }
        }
    }
    // y- PML
    for (int i = 0; i < Nx; ++i) {
        for (int k = 1; k < Nz; ++k) {
            for (int l = 0; l < n_pml_yn; ++l) {
                int j = 1 + l;
                if (j >= Ny) continue;
                double curlH_y = Hz[id_Hz(i, j, k)] - Hz[id_Hz(i, j - 1, k)];
                int pid_y = id_Phi_ex_yn(i, k - 1, l);
                Phi_ex_yn[pid_y] = b_e_yn[l] * Phi_ex_yn[pid_y] + a_e_yn[l] * curlH_y;
                Dx[id_Dx(i, j, k)] += C_Dx_dyn_v[l] * curlH_y + CD_dt * Phi_ex_yn[pid_y];
            }
        }
    }
    // y+ PML
    for (int i = 0; i < Nx; ++i) {
        for (int k = 1; k < Nz; ++k) {
            for (int l = 0; l < n_pml_yp; ++l) {
                int j = Ny - n_pml_yp + l;
                if (j < 1 || j >= Ny) continue;
                double curlH_y = Hz[id_Hz(i, j, k)] - Hz[id_Hz(i, j - 1, k)];
                int pid_y = id_Phi_ex_yp(i, k - 1, l);
                Phi_ex_yp[pid_y] = b_e_yp[l] * Phi_ex_yp[pid_y] + a_e_yp[l] * curlH_y;
                Dx[id_Dx(i, j, k)] += C_Dx_dyp_v[l] * curlH_y + CD_dt * Phi_ex_yp[pid_y];
            }
        }
    }
    // y 中间区域
    for (int i = 0; i < Nx; ++i) {
        for (int k = 1; k < Nz; ++k) {
            int j_begin = std::max(1, n_pml_yn + 1);
            int j_end = std::min(Ny, Ny - n_pml_yp);
            if (j_begin < j_end) {
                for (int j = j_begin; j < j_end; ++j) {
                    double curlH_y = Hz[id_Hz(i, j, k)] - Hz[id_Hz(i, j - 1, k)];
                    Dx[id_Dx(i, j, k)] += CD_dt_dy * curlH_y;
                }
            }
        }
    }

    

    // ========== Dy 更新：curl H = dHx/dz - dHz/dx ==========
    // z- PML
    for (int i = 1; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int l = 0; l < n_pml_zn; ++l) {
                int k = 1 + l;
                if (k >= Nz) continue;
                double curlH_z = Hx[id_Hx(i, j, k)] - Hx[id_Hx(i, j, k - 1)];
                int pid_z = id_Phi_ey_zn(i - 1, j, l);
                Phi_ey_zn[pid_z] = b_e_zn[l] * Phi_ey_zn[pid_z] + a_e_zn[l] * curlH_z;
                Dy[id_Dy(i, j, k)] += C_Dy_dzn_v[l] * curlH_z + CD_dt * Phi_ey_zn[pid_z];
            }
        }
    }
    // z+ PML
    for (int i = 1; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int l = 0; l < n_pml_zp; ++l) {
                int k = Nz - n_pml_zp + l;
                if (k < 1 || k >= Nz) continue;
                double curlH_z = Hx[id_Hx(i, j, k)] - Hx[id_Hx(i, j, k - 1)];
                int pid_z = id_Phi_ey_zp(i - 1, j, l);
                Phi_ey_zp[pid_z] = b_e_zp[l] * Phi_ey_zp[pid_z] + a_e_zp[l] * curlH_z;
                Dy[id_Dy(i, j, k)] += C_Dy_dzp_v[l] * curlH_z + CD_dt * Phi_ey_zp[pid_z];
            }
        }
    }
    // z 中间区域
    for (int i = 1; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            int k_begin = std::max(1, n_pml_zn + 1);
            int k_end = std::min(Nz, Nz - n_pml_zp);
            if (k_begin < k_end) {
                for (int k = k_begin; k < k_end; ++k) {
                    double curlH_z = Hx[id_Hx(i, j, k)] - Hx[id_Hx(i, j, k - 1)];
                    Dy[id_Dy(i, j, k)] += CD_dt_dz * curlH_z;
                }
            }
        }
    }
    // x- PML
    for (int j = 0; j < Ny; ++j) {
        for (int k = 1; k < Nz; ++k) {
            for (int l = 0; l < n_pml_xn; ++l) {
                int i = 1 + l;
                if (i >= Nx) continue;
                double curlH_x = Hz[id_Hz(i, j, k)] - Hz[id_Hz(i - 1, j, k)];
                int pid_x = id_Phi_ey_xn(l, j, k - 1);
                Phi_ey_xn[pid_x] = b_e_xn[l] * Phi_ey_xn[pid_x] + a_e_xn[l] * curlH_x;
                Dy[id_Dy(i, j, k)] -= C_Dy_dxn_v[l] * curlH_x + CD_dt * Phi_ey_xn[pid_x];
            }
        }
    }
    // x+ PML
    for (int j = 0; j < Ny; ++j) {
        for (int k = 1; k < Nz; ++k) {
            for (int l = 0; l < n_pml_xp; ++l) {
                int i = Nx - n_pml_xp + l;
                if (i < 1 || i >= Nx) continue;
                double curlH_x = Hz[id_Hz(i, j, k)] - Hz[id_Hz(i - 1, j, k)];
                int pid_x = id_Phi_ey_xp(l, j, k - 1);
                Phi_ey_xp[pid_x] = b_e_xp[l] * Phi_ey_xp[pid_x] + a_e_xp[l] * curlH_x;
                Dy[id_Dy(i, j, k)] -= C_Dy_dxp_v[l] * curlH_x + CD_dt * Phi_ey_xp[pid_x];
            }
        }
    }
    // x 中间区域
    for (int j = 0; j < Ny; ++j) {
        for (int k = 1; k < Nz; ++k) {
            int i_begin = std::max(1, n_pml_xn + 1);
            int i_end = std::min(Nx, Nx - n_pml_xp);
            if (i_begin < i_end) {
                for (int i = i_begin; i < i_end; ++i) {
                    double curlH_x = Hz[id_Hz(i, j, k)] - Hz[id_Hz(i - 1, j, k)];
                    Dy[id_Dy(i, j, k)] -= CD_dt_dx * curlH_x;
                }
            }
        }
    }

    // ========== Dz 更新：curl H = dHy/dx - dHx/dy ==========
    // x- PML
    for (int j = 1; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            for (int l = 0; l < n_pml_xn; ++l) {
                int i = 1 + l;
                if (i >= Nx) continue;
                double curlH_x = Hy[id_Hy(i, j, k)] - Hy[id_Hy(i - 1, j, k)];
                int pid_x = id_Phi_ez_xn(l, j - 1, k);
                Phi_ez_xn[pid_x] = b_e_xn[l] * Phi_ez_xn[pid_x] + a_e_xn[l] * curlH_x;
                Dz[id_Dz(i, j, k)] += C_Dz_dxn_v[l] * curlH_x + CD_dt * Phi_ez_xn[pid_x];
            }
        }
    }
    // x+ PML
    for (int j = 1; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            for (int l = 0; l < n_pml_xp; ++l) {
                int i = Nx - n_pml_xp + l;
                if (i < 1 || i >= Nx) continue;
                double curlH_x = Hy[id_Hy(i, j, k)] - Hy[id_Hy(i - 1, j, k)];
                int pid_x = id_Phi_ez_xp(l, j - 1, k);
                Phi_ez_xp[pid_x] = b_e_xp[l] * Phi_ez_xp[pid_x] + a_e_xp[l] * curlH_x;
                Dz[id_Dz(i, j, k)] += C_Dz_dxp_v[l] * curlH_x + CD_dt * Phi_ez_xp[pid_x];
            }
        }
    }
    // x 中间区域
    for (int j = 1; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            int i_begin = std::max(1, n_pml_xn + 1);
            int i_end = std::min(Nx, Nx - n_pml_xp);
            if (i_begin < i_end) {
                for (int i = i_begin; i < i_end; ++i) {
                    double curlH_x = Hy[id_Hy(i, j, k)] - Hy[id_Hy(i - 1, j, k)];
                    Dz[id_Dz(i, j, k)] += CD_dt_dx * curlH_x;
                }
            }
        }
    }
    // y- PML
    for (int i = 1; i < Nx; ++i) {
        for (int k = 0; k < Nz; ++k) {
            for (int l = 0; l < n_pml_yn; ++l) {
                int j = 1 + l;
                if (j >= Ny) continue;
                double curlH_y = Hx[id_Hx(i, j, k)] - Hx[id_Hx(i, j - 1, k)];
                int pid_y = id_Phi_ez_yn(i - 1, l, k);
                Phi_ez_yn[pid_y] = b_e_yn[l] * Phi_ez_yn[pid_y] + a_e_yn[l] * curlH_y;
                Dz[id_Dz(i, j, k)] -= C_Dz_dyn_v[l] * curlH_y + CD_dt * Phi_ez_yn[pid_y];
            }
        }
    }
    // y+ PML
    for (int i = 1; i < Nx; ++i) {
        for (int k = 0; k < Nz; ++k) {
            for (int l = 0; l < n_pml_yp; ++l) {
                int j = Ny - n_pml_yp + l;
                if (j < 1 || j >= Ny) continue;
                double curlH_y = Hx[id_Hx(i, j, k)] - Hx[id_Hx(i, j - 1, k)];
                int pid_y = id_Phi_ez_yp(i - 1, l, k);
                Phi_ez_yp[pid_y] = b_e_yp[l] * Phi_ez_yp[pid_y] + a_e_yp[l] * curlH_y;
                Dz[id_Dz(i, j, k)] -= C_Dz_dyp_v[l] * curlH_y + CD_dt * Phi_ez_yp[pid_y];
            }
        }
    }
    // y 中间区域
    for (int i = 1; i < Nx; ++i) {
        for (int k = 0; k < Nz; ++k) {
            int j_begin = std::max(1, n_pml_yn + 1);
            int j_end = std::min(Ny, Ny - n_pml_yp);
            if (j_begin < j_end) {
                for (int j = j_begin; j < j_end; ++j) {
                    double curlH_y = Hx[id_Hx(i, j, k)] - Hx[id_Hx(i, j - 1, k)];
                    Dz[id_Dz(i, j, k)] -= CD_dt_dy * curlH_y;
                }
            }
        }
    }

    applySource(step);
    updateElectricMaterialRelations();
}

void FDTD3D::updateH() {
    const double CB_dt = m_dt;
    const double CB_dt_dx = m_dt / m_dx;
    const double CB_dt_dy = m_dt / m_dy;
    const double CB_dt_dz = m_dt / m_dz;

    const int Nx = m_Nx, Ny = m_Ny, Nz = m_Nz;

    // ========== Bx 更新：curl E = dEy/dz - dEz/dy ==========
    // z- PML
    for (int i = 1; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int l = 0; l < n_pml_zn; ++l) {
                int k = 0 + l;
                if (k + 1 > Nz) continue;
                double curlE_z = Ey[id_Ey(i, j, k + 1)] - Ey[id_Ey(i, j, k)];
                int pid_z = id_Phi_mx_zn(i - 1, j, l);
                Phi_mx_zn[pid_z] = b_m_zn[l] * Phi_mx_zn[pid_z] + a_m_zn[l] * curlE_z;
                Bx[id_Bx(i, j, k)] += C_Bx_dzn_v[l] * curlE_z + CB_dt * Phi_mx_zn[pid_z];
            }
        }
    }
    // z+ PML
    for (int i = 1; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int l = 0; l < n_pml_zp; ++l) {
                int k = Nz - n_pml_zp + l;
                if (k + 1 > Nz) continue;
                double curlE_z = Ey[id_Ey(i, j, k + 1)] - Ey[id_Ey(i, j, k)];
                int pid_z = id_Phi_mx_zp(i - 1, j, l);
                Phi_mx_zp[pid_z] = b_m_zp[l] * Phi_mx_zp[pid_z] + a_m_zp[l] * curlE_z;
                Bx[id_Bx(i, j, k)] += C_Bx_dzp_v[l] * curlE_z + CB_dt * Phi_mx_zp[pid_z];
            }
        }
    }
    // z 中间区域
    for (int i = 1; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            int k_begin = std::max(0, n_pml_zn);
            int k_end = std::min(Nz - 1, Nz - n_pml_zp);
            for (int k = k_begin; k < k_end; ++k) {
                double curlE_z = Ey[id_Ey(i, j, k + 1)] - Ey[id_Ey(i, j, k)];
                Bx[id_Bx(i, j, k)] += CB_dt_dz * curlE_z;
            }
        }
    }
    // y- PML
    for (int i = 1; i < Nx; ++i) {
        for (int k = 0; k < Nz; ++k) {
            for (int l = 0; l < n_pml_yn; ++l) {
                int j = 0 + l;
                if (j + 1 > Ny) continue;
                double curlE_y = Ez[id_Ez(i, j + 1, k)] - Ez[id_Ez(i, j, k)];
                int pid_y = id_Phi_mx_yn(i - 1, l, k);
                Phi_mx_yn[pid_y] = b_m_yn[l] * Phi_mx_yn[pid_y] + a_m_yn[l] * curlE_y;
                Bx[id_Bx(i, j, k)] -= C_Bx_dyn_v[l] * curlE_y + CB_dt * Phi_mx_yn[pid_y];
            }
        }
    }
    // y+ PML
    for (int i = 1; i < Nx; ++i) {
        for (int k = 0; k < Nz; ++k) {
            for (int l = 0; l < n_pml_yp; ++l) {
                int j = Ny - n_pml_yp + l;
                if (j + 1 > Ny) continue;
                double curlE_y = Ez[id_Ez(i, j + 1, k)] - Ez[id_Ez(i, j, k)];
                int pid_y = id_Phi_mx_yp(i - 1, l, k);
                Phi_mx_yp[pid_y] = b_m_yp[l] * Phi_mx_yp[pid_y] + a_m_yp[l] * curlE_y;
                Bx[id_Bx(i, j, k)] -= C_Bx_dyp_v[l] * curlE_y + CB_dt * Phi_mx_yp[pid_y];
            }
        }
    }
    // y 中间区域
    for (int i = 1; i < Nx; ++i) {
        for (int k = 0; k < Nz; ++k) {
            int j_begin = std::max(0, n_pml_yn);
            int j_end = std::min(Ny - 1, Ny - n_pml_yp);
            for (int j = j_begin; j < j_end; ++j) {
                double curlE_y = Ez[id_Ez(i, j + 1, k)] - Ez[id_Ez(i, j, k)];
                Bx[id_Bx(i, j, k)] -= CB_dt_dy * curlE_y;
            }
        }
    }

    // ========== By 更新：curl E = dEz/dx - dEx/dz ==========
    // x- PML
    for (int j = 1; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            for (int l = 0; l < n_pml_xn; ++l) {
                int i = 0 + l;
                if (i + 1 > Nx) continue;
                double curlE_x = Ez[id_Ez(i + 1, j, k)] - Ez[id_Ez(i, j, k)];
                int pid_x = id_Phi_my_xn(l, j - 1, k);
                Phi_my_xn[pid_x] = b_m_xn[l] * Phi_my_xn[pid_x] + a_m_xn[l] * curlE_x;
                By[id_By(i, j, k)] += C_By_dxn_v[l] * curlE_x + CB_dt * Phi_my_xn[pid_x];
            }
        }
    }
    // x+ PML
    for (int j = 1; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            for (int l = 0; l < n_pml_xp; ++l) {
                int i = Nx - n_pml_xp + l;
                if (i + 1 > Nx) continue;
                double curlE_x = Ez[id_Ez(i + 1, j, k)] - Ez[id_Ez(i, j, k)];
                int pid_x = id_Phi_my_xp(l, j - 1, k);
                Phi_my_xp[pid_x] = b_m_xp[l] * Phi_my_xp[pid_x] + a_m_xp[l] * curlE_x;
                By[id_By(i, j, k)] += C_By_dxp_v[l] * curlE_x + CB_dt * Phi_my_xp[pid_x];
            }
        }
    }
    // x 中间区域
    for (int j = 1; j < Ny; ++j) {
        for (int k = 0; k < Nz; ++k) {
            int i_begin = std::max(0, n_pml_xn);
            int i_end = std::min(Nx - 1, Nx - n_pml_xp);
            for (int i = i_begin; i < i_end; ++i) {
                double curlE_x = Ez[id_Ez(i + 1, j, k)] - Ez[id_Ez(i, j, k)];
                By[id_By(i, j, k)] += CB_dt_dx * curlE_x;
            }
        }
    }
    // z- PML
    for (int i = 0; i < Nx; ++i) {
        for (int j = 1; j < Ny; ++j) {
            for (int l = 0; l < n_pml_zn; ++l) {
                int k = 0 + l;
                if (k + 1 > Nz) continue;
                double curlE_z = Ex[id_Ex(i, j, k + 1)] - Ex[id_Ex(i, j, k)];
                int pid_z = id_Phi_my_zn(i, j - 1, l);
                Phi_my_zn[pid_z] = b_m_zn[l] * Phi_my_zn[pid_z] + a_m_zn[l] * curlE_z;
                By[id_By(i, j, k)] -= C_By_dzn_v[l] * curlE_z + CB_dt * Phi_my_zn[pid_z];
            }
        }
    }
    // z+ PML
    for (int i = 0; i < Nx; ++i) {
        for (int j = 1; j < Ny; ++j) {
            for (int l = 0; l < n_pml_zp; ++l) {
                int k = Nz - n_pml_zp + l;
                if (k + 1 > Nz) continue;
                double curlE_z = Ex[id_Ex(i, j, k + 1)] - Ex[id_Ex(i, j, k)];
                int pid_z = id_Phi_my_zp(i, j - 1, l);
                Phi_my_zp[pid_z] = b_m_zp[l] * Phi_my_zp[pid_z] + a_m_zp[l] * curlE_z;
                By[id_By(i, j, k)] -= C_By_dzp_v[l] * curlE_z + CB_dt * Phi_my_zp[pid_z];
            }
        }
    }
    // z 中间区域
    for (int i = 0; i < Nx; ++i) {
        for (int j = 1; j < Ny; ++j) {
            int k_begin = std::max(0, n_pml_zn);
            int k_end = std::min(Nz - 1, Nz - n_pml_zp);
            for (int k = k_begin; k < k_end; ++k) {
                double curlE_z = Ex[id_Ex(i, j, k + 1)] - Ex[id_Ex(i, j, k)];
                By[id_By(i, j, k)] -= CB_dt_dz * curlE_z;
            }
        }
    }

    // ========== Bz 更新：curl E = dEx/dy - dEy/dx ==========
    // y- PML
    for (int i = 0; i < Nx; ++i) {
        for (int k = 1; k < Nz; ++k) {
            for (int l = 0; l < n_pml_yn; ++l) {
                int j = 0 + l;
                if (j + 1 > Ny) continue;
                double curlE_y = Ex[id_Ex(i, j + 1, k)] - Ex[id_Ex(i, j, k)];
                int pid_y = id_Phi_mz_yn(i, l, k - 1);
                Phi_mz_yn[pid_y] = b_m_yn[l] * Phi_mz_yn[pid_y] + a_m_yn[l] * curlE_y;
                Bz[id_Bz(i, j, k)] += C_Bz_dyn_v[l] * curlE_y + CB_dt * Phi_mz_yn[pid_y];
            }
        }
    }
    // y+ PML
    for (int i = 0; i < Nx; ++i) {
        for (int k = 1; k < Nz; ++k) {
            for (int l = 0; l < n_pml_yp; ++l) {
                int j = Ny - n_pml_yp + l;
                if (j + 1 > Ny) continue;
                double curlE_y = Ex[id_Ex(i, j + 1, k)] - Ex[id_Ex(i, j, k)];
                int pid_y = id_Phi_mz_yp(i, l, k - 1);
                Phi_mz_yp[pid_y] = b_m_yp[l] * Phi_mz_yp[pid_y] + a_m_yp[l] * curlE_y;
                Bz[id_Bz(i, j, k)] += C_Bz_dyp_v[l] * curlE_y + CB_dt * Phi_mz_yp[pid_y];
            }
        }
    }
    // y 中间区域
    for (int i = 0; i < Nx; ++i) {
        for (int k = 1; k < Nz; ++k) {
            int j_begin = std::max(0, n_pml_yn);
            int j_end = std::min(Ny - 1, Ny - n_pml_yp);
            for (int j = j_begin; j < j_end; ++j) {
                double curlE_y = Ex[id_Ex(i, j + 1, k)] - Ex[id_Ex(i, j, k)];
                Bz[id_Bz(i, j, k)] += CB_dt_dy * curlE_y;
            }
        }
    }
    // x- PML
    for (int j = 0; j < Ny; ++j) {
        for (int k = 1; k < Nz; ++k) {
            for (int l = 0; l < n_pml_xn; ++l) {
                int i = 0 + l;
                if (i + 1 > Nx) continue;
                double curlE_x = Ey[id_Ey(i + 1, j, k)] - Ey[id_Ey(i, j, k)];
                int pid_x = id_Phi_mz_xn(l, j, k - 1);
                Phi_mz_xn[pid_x] = b_m_xn[l] * Phi_mz_xn[pid_x] + a_m_xn[l] * curlE_x;
                Bz[id_Bz(i, j, k)] -= C_Bz_dxn_v[l] * curlE_x + CB_dt * Phi_mz_xn[pid_x];
            }
        }
    }
    // x+ PML
    for (int j = 0; j < Ny; ++j) {
        for (int k = 1; k < Nz; ++k) {
            for (int l = 0; l < n_pml_xp; ++l) {
                int i = Nx - n_pml_xp + l;
                if (i + 1 > Nx) continue;
                double curlE_x = Ey[id_Ey(i + 1, j, k)] - Ey[id_Ey(i, j, k)];
                int pid_x = id_Phi_mz_xp(l, j, k - 1);
                Phi_mz_xp[pid_x] = b_m_xp[l] * Phi_mz_xp[pid_x] + a_m_xp[l] * curlE_x;
                Bz[id_Bz(i, j, k)] -= C_Bz_dxp_v[l] * curlE_x + CB_dt * Phi_mz_xp[pid_x];
            }
        }
    }
    // x 中间区域
    for (int j = 0; j < Ny; ++j) {
        for (int k = 1; k < Nz; ++k) {
            int i_begin = std::max(0, n_pml_xn);
            int i_end = std::min(Nx - 1, Nx - n_pml_xp);
            for (int i = i_begin; i < i_end; ++i) {
                double curlE_x = Ey[id_Ey(i + 1, j, k)] - Ey[id_Ey(i, j, k)];
                Bz[id_Bz(i, j, k)] -= CB_dt_dx * curlE_x;
            }
        }
    }

    updateMagneticMaterialRelations();
}

void FDTD3D::dumpFrame(int step) {
    std::vector<double> E_mag(m_Nx * m_Ny * m_Nz, 0.0);
    const double inv4 = 0.25;

    auto center_e_sq = [&](int i, int j, int k) {
        auto avg_sq = [&](double a, double b, double c, double d) {
            return inv4 * (a * a + b * b + c * c + d * d);
        };

        double ex_sq = avg_sq(Ex[id_Ex(i, j, k)],     Ex[id_Ex(i, j + 1, k)],
                              Ex[id_Ex(i, j, k + 1)], Ex[id_Ex(i, j + 1, k + 1)]);
        double ey_sq = avg_sq(Ey[id_Ey(i, j, k)],     Ey[id_Ey(i + 1, j, k)],
                              Ey[id_Ey(i, j, k + 1)], Ey[id_Ey(i + 1, j, k + 1)]);
        double ez_sq = avg_sq(Ez[id_Ez(i, j, k)],     Ez[id_Ez(i + 1, j, k)],
                              Ez[id_Ez(i, j + 1, k)], Ez[id_Ez(i + 1, j + 1, k)]);
        return ex_sq + ey_sq + ez_sq;
    };

    for (int i = 0; i < m_Nx; ++i) {
        for (int j = 0; j < m_Ny; ++j) {
            for (int k = 0; k < m_Nz; ++k) {
                E_mag[i * m_Ny * m_Nz + j * m_Nz + k] = std::sqrt(center_e_sq(i, j, k));
            }
        }
    }

    m_out->saveScalar3D_VTI("E_mag", E_mag, m_Nx, m_Ny, m_Nz, m_dx, m_dy, m_dz, step);
}

void FDTD3D::run() {
    int N_lambda = 20;
    int nCFL = 2;
    int n_T0 = static_cast<int>(std::round(nCFL * N_lambda));
    int n_T0d2 = n_T0 / 2;

    std::vector<double> E_mag_sq(m_Nx * m_Ny * m_Nz, 0.0);
    const double inv16 = 1.0 / 16.0;

    auto center_e_sq = [&](int i, int j, int k) {
        double ex_sum = Ex[id_Ex(i, j, k)] + Ex[id_Ex(i, j + 1, k)]
                      + Ex[id_Ex(i, j, k + 1)] + Ex[id_Ex(i, j + 1, k + 1)];
        double ey_sum = Ey[id_Ey(i, j, k)] + Ey[id_Ey(i + 1, j, k)]
                      + Ey[id_Ey(i, j, k + 1)] + Ey[id_Ey(i + 1, j, k + 1)];
        double ez_sum = Ez[id_Ez(i, j, k)] + Ez[id_Ez(i + 1, j, k)]
                      + Ez[id_Ez(i, j + 1, k)] + Ez[id_Ez(i + 1, j + 1, k)];
        return inv16 * (ex_sum * ex_sum + ey_sum * ey_sum + ez_sum * ez_sum);
    };

    for (int n = 0; n < m_steps; ++n) {
        updateE(n);
        updateH();

        // RMS 统计
        if (n % n_T0d2 == 0) std::fill(E_mag_sq.begin(), E_mag_sq.end(), 0.0);
        for (int i = 0; i < m_Nx; ++i) {
            for (int j = 0; j < m_Ny; ++j) {
                for (int k = 0; k < m_Nz; ++k) {
                    double e_sq = center_e_sq(i, j, k);
                    E_mag_sq[i * m_Ny * m_Nz + j * m_Nz + k] += e_sq * (4.0 / n_T0);
                }
            }
        }
        if ((n + 1) % n_T0d2 == 0) {
            std::vector<double> E_mag(m_Nx * m_Ny * m_Nz);
            for (size_t p = 0; p < E_mag.size(); ++p) E_mag[p] = std::sqrt(E_mag_sq[p]);
            m_out->saveScalar3D_VTI("E_mag_rms", E_mag, m_Nx, m_Ny, m_Nz, m_dx, m_dy, m_dz, n + 1);
        }

        if (n % 20 == 0) {
            dumpFrame(n);
            std::cout << "[n=" << n << "] frame dumped.\n";
        }
    }
    std::cout << "Done. VTI in ./output_3d\n";
}