#pragma once

#include <memory>
#include <string>
#include <vector>
#include <cmath>

class OutputManager;

class FDTD3D {
public:
    FDTD3D();
    ~FDTD3D(); // 只声明，在 cpp 中定义（避免不完全类型）

    void run();

private:
    // ========== 基本物理常量 ==========
    double m_c0   = 299792458.0;
    double m_mu0  = 4.0 * M_PI * 1e-7;
    double m_eps0 = 1.0 / (m_c0 * m_c0 * m_mu0);
    double m_Z0   = std::sqrt(m_mu0 / m_eps0);

    // ========== 网格参数 ==========
    int m_Nx = 0, m_Ny = 0, m_Nz = 0;          // Yee 网格单元数（电场 D/E 定义域）
    double m_dx = 0.0, m_dy = 0.0, m_dz = 0.0;
    int m_steps = 0;
    double m_dt = 0.0;

    // ========== 材料系数（与各向 Yee 面/边对齐） ==========
    // eps_x: [Nx] x [Ny+1] x [Nz+1]
    std::vector<double> m_eps_x;
    // eps_y: [Nx+1] x [Ny] x [Nz+1]
    std::vector<double> m_eps_y;
    // eps_z: [Nx+1] x [Ny+1] x [Nz]
    std::vector<double> m_eps_z;

    // ========== 场量（扁平化：x 快变，y 中变，z 慢变） ==========
    // Dx, Ex: [Nx] x [Ny+1] x [Nz+1]
    std::vector<double> Dx, Ex;
    // Dy, Ey: [Nx+1] x [Ny] x [Nz+1]
    std::vector<double> Dy, Ey;
    // Dz, Ez: [Nx+1] x [Ny+1] x [Nz]
    std::vector<double> Dz, Ez;

    // Bx, Hx: [Nx+1] x [Ny] x [Nz]
    std::vector<double> Bx, Hx;
    // By, Hy: [Nx] x [Ny+1] x [Nz]
    std::vector<double> By, Hy;
    // Bz, Hz: [Nx] x [Ny] x [Nz+1]
    std::vector<double> Bz, Hz;

    // ========== CPML 参数 ==========
    int n_pml_xn = 20, n_pml_xp = 20;
    int n_pml_yn = 20, n_pml_yp = 20;
    int n_pml_zn = 20, n_pml_zp = 20;

    // --- 电场 CPML 系数（1D 剖面）---
    // y 方向 PML（影响 Dx, Dz）
    std::vector<double> b_e_yn, a_e_yn, C_Dx_dyn_v, C_Dz_dyn_v;
    std::vector<double> b_e_yp, a_e_yp, C_Dx_dyp_v, C_Dz_dyp_v;
    // z 方向 PML（影响 Dx, Dy）
    std::vector<double> b_e_zn, a_e_zn, C_Dx_dzn_v, C_Dy_dzn_v;
    std::vector<double> b_e_zp, a_e_zp, C_Dx_dzp_v, C_Dy_dzp_v;
    // x 方向 PML（影响 Dy, Dz）
    std::vector<double> b_e_xn, a_e_xn, C_Dy_dxn_v, C_Dz_dxn_v;
    std::vector<double> b_e_xp, a_e_xp, C_Dy_dxp_v, C_Dz_dxp_v;

    // --- 磁场 CPML 系数（1D 剖面）---
    // y 方向 PML（影响 Bx, Bz）
    std::vector<double> b_m_yn, a_m_yn, C_Bx_dyn_v, C_Bz_dyn_v;
    std::vector<double> b_m_yp, a_m_yp, C_Bx_dyp_v, C_Bz_dyp_v;
    // z 方向 PML（影响 Bx, By）
    std::vector<double> b_m_zn, a_m_zn, C_Bx_dzn_v, C_By_dzn_v;
    std::vector<double> b_m_zp, a_m_zp, C_Bx_dzp_v, C_By_dzp_v;
    // x 方向 PML（影响 By, Bz）
    std::vector<double> b_m_xn, a_m_xn, C_By_dxn_v, C_Bz_dxn_v;
    std::vector<double> b_m_xp, a_m_xp, C_By_dxp_v, C_Bz_dxp_v;

    // --- CPML 辅助场 Φ（扁平存储）---
    // Dx 的 Φ（来自 y 和 z 方向）
    std::vector<double> Phi_ex_yn, Phi_ex_yp; // [Nx] x [n_pml_y*] x [Nz-1] → 扁平为 Nx*(Nz-1)*n_pml
    std::vector<double> Phi_ex_zn, Phi_ex_zp; // [Nx] x [Ny-1] x [n_pml_z*]

    // Dy 的 Φ（来自 z 和 x 方向）
    std::vector<double> Phi_ey_zn, Phi_ey_zp; // [Nx-1] x [Ny] x [n_pml_z*]
    std::vector<double> Phi_ey_xn, Phi_ey_xp; // [n_pml_x*] x [Ny] x [Nz-1]

    // Dz 的 Φ（来自 x 和 y 方向）
    std::vector<double> Phi_ez_xn, Phi_ez_xp; // [n_pml_x*] x [Ny-1] x [Nz]
    std::vector<double> Phi_ez_yn, Phi_ez_yp; // [Nx-1] x [n_pml_y*] x [Nz]

    // Bx 的 Φ（来自 y 和 z 方向）
    std::vector<double> Phi_mx_yn, Phi_mx_yp; // [Nx-1] x [n_pml_y*] x [Nz]
    std::vector<double> Phi_mx_zn, Phi_mx_zp; // [Nx-1] x [Ny] x [n_pml_z*]

    // By 的 Φ（来自 x 和 z 方向）
    std::vector<double> Phi_my_xn, Phi_my_xp; // [n_pml_x*] x [Ny-1] x [Nz]
    std::vector<double> Phi_my_zn, Phi_my_zp; // [Nx] x [Ny-1] x [n_pml_z*]

    // Bz 的 Φ（来自 x 和 y 方向）
    std::vector<double> Phi_mz_xn, Phi_mz_xp; // [n_pml_x*] x [Ny] x [Nz-1]
    std::vector<double> Phi_mz_yn, Phi_mz_yp; // [Nx] x [n_pml_y*] x [Nz-1]

    // ========== 源 ==========
    int m_k_Jx = 0; // z index of Jx source plane (k_Jx in MATLAB)
    std::vector<double> m_Jx_i_mag; // [Nx] x [Ny+1] → 扁平
    std::vector<double> m_wf_Jx_i;  // time waveform, size = m_steps

    // ========== 输出管理 ==========
    std::unique_ptr<OutputManager> m_out;

private:
    // ========== 扁平索引函数（Yee 网格）==========
    inline int id_Dx(int ix, int jy, int kz) const { return ix*(m_Ny+1)*(m_Nz+1) + jy*(m_Nz+1) + kz; } // ix:0..Nx-1, jy:0..Ny, kz:0..Nz
    inline int id_Dy(int ix, int jy, int kz) const { return ix*m_Ny*(m_Nz+1) + jy*(m_Nz+1) + kz; }     // ix:0..Nx, jy:0..Ny-1, kz:0..Nz
    inline int id_Dz(int ix, int jy, int kz) const { return ix*(m_Ny+1)*m_Nz + jy*m_Nz + kz; }         // ix:0..Nx, jy:0..Ny, kz:0..Nz-1

    inline int id_Ex(int ix, int jy, int kz) const { return id_Dx(ix, jy, kz); }
    inline int id_Ey(int ix, int jy, int kz) const { return id_Dy(ix, jy, kz); }
    inline int id_Ez(int ix, int jy, int kz) const { return id_Dz(ix, jy, kz); }

    inline int id_Bx(int ix, int jy, int kz) const { return ix*m_Ny*m_Nz + jy*m_Nz + kz; }             // ix:0..Nx, jy:0..Ny-1, kz:0..Nz-1
    inline int id_By(int ix, int jy, int kz) const { return ix*(m_Ny+1)*m_Nz + jy*m_Nz + kz; }         // ix:0..Nx-1, jy:0..Ny, kz:0..Nz-1
    inline int id_Bz(int ix, int jy, int kz) const { return ix*m_Ny*(m_Nz+1) + jy*(m_Nz+1) + kz; }     // ix:0..Nx-1, jy:0..Ny-1, kz:0..Nz

    inline int id_Hx(int ix, int jy, int kz) const { return id_Bx(ix, jy, kz); }
    inline int id_Hy(int ix, int jy, int kz) const { return id_By(ix, jy, kz); }
    inline int id_Hz(int ix, int jy, int kz) const { return id_Bz(ix, jy, kz); }

    inline int id_eps_x(int ix, int jy, int kz) const { return id_Dx(ix, jy, kz); }                  // Nx * (Ny+1) * (Nz+1)
    inline int id_eps_y(int ix, int jy, int kz) const { return id_Dy(ix, jy, kz); }                  // (Nx+1) * Ny * (Nz+1)
    inline int id_eps_z(int ix, int jy, int kz) const { return id_Dz(ix, jy, kz); }                  // (Nx+1) * (Ny+1) * Nz

    // ========== CPML Φ 场索引（修正后，匹配新的内存分配）==========
    // Dx 的 Φ 索引
    // Phi_ex_yn: [Nx] x [Nz+1] x [n_pml_yn]
    inline int id_Phi_ex_yn(int ix, int kz, int l) const { return (ix*(m_Nz+1) + kz)*n_pml_yn + l; }
    inline int id_Phi_ex_yp(int ix, int kz, int l) const { return (ix*(m_Nz+1) + kz)*n_pml_yp + l; }
    // Phi_ex_zn: [Nx] x [Ny+1] x [n_pml_zn]
    inline int id_Phi_ex_zn(int ix, int jy, int l) const { return (ix*(m_Ny+1) + jy)*n_pml_zn + l; }
    inline int id_Phi_ex_zp(int ix, int jy, int l) const { return (ix*(m_Ny+1) + jy)*n_pml_zp + l; }

    // Dy 的 Φ 索引
    // Phi_ey_zn: [Nx+1] x [Ny] x [n_pml_zn]
    inline int id_Phi_ey_zn(int ix, int jy, int l) const { return (ix*m_Ny + jy)*n_pml_zn + l; }
    inline int id_Phi_ey_zp(int ix, int jy, int l) const { return (ix*m_Ny + jy)*n_pml_zp + l; }
    // Phi_ey_xn: [n_pml_xn] x [Ny] x [Nz+1]
    inline int id_Phi_ey_xn(int l, int jy, int kz) const { return (l*m_Ny*(m_Nz+1) + jy*(m_Nz+1) + kz); }
    inline int id_Phi_ey_xp(int l, int jy, int kz) const { return (l*m_Ny*(m_Nz+1) + jy*(m_Nz+1) + kz); }

    // Dz 的 Φ 索引
    // Phi_ez_xn: [n_pml_xn] x [Ny+1] x [Nz]
    inline int id_Phi_ez_xn(int l, int jy, int kz) const { return (l*(m_Ny+1)*m_Nz + jy*m_Nz + kz); }
    inline int id_Phi_ez_xp(int l, int jy, int kz) const { return (l*(m_Ny+1)*m_Nz + jy*m_Nz + kz); }
    // Phi_ez_yn: [Nx+1] x [n_pml_yn] x [Nz]
    inline int id_Phi_ez_yn(int ix, int l, int kz) const { return (ix*n_pml_yn*m_Nz + l*m_Nz + kz); }
    inline int id_Phi_ez_yp(int ix, int l, int kz) const { return (ix*n_pml_yp*m_Nz + l*m_Nz + kz); }

    // Bx 的 Φ 索引
    // Phi_mx_yn: [Nx+1] x [n_pml_yn] x [Nz]
    inline int id_Phi_mx_yn(int ix, int l, int kz) const { return (ix*n_pml_yn*m_Nz + l*m_Nz + kz); }
    inline int id_Phi_mx_yp(int ix, int l, int kz) const { return (ix*n_pml_yp*m_Nz + l*m_Nz + kz); }
    // Phi_mx_zn: [Nx+1] x [Ny] x [n_pml_zn]
    inline int id_Phi_mx_zn(int ix, int jy, int l) const { return (ix*m_Ny*n_pml_zn + jy*n_pml_zn + l); }
    inline int id_Phi_mx_zp(int ix, int jy, int l) const { return (ix*m_Ny*n_pml_zp + jy*n_pml_zp + l); }

    // By 的 Φ 索引
    // Phi_my_xn: [n_pml_xn] x [Ny+1] x [Nz]
    inline int id_Phi_my_xn(int l, int jy, int kz) const { return (l*(m_Ny+1)*m_Nz + jy*m_Nz + kz); }
    inline int id_Phi_my_xp(int l, int jy, int kz) const { return (l*(m_Ny+1)*m_Nz + jy*m_Nz + kz); }
    // Phi_my_zn: [Nx] x [Ny+1] x [n_pml_zn]
    inline int id_Phi_my_zn(int ix, int jy, int l) const { return (ix*(m_Ny+1)*n_pml_zn + jy*n_pml_zn + l); }
    inline int id_Phi_my_zp(int ix, int jy, int l) const { return (ix*(m_Ny+1)*n_pml_zp + jy*n_pml_zp + l); }

    // Bz 的 Φ 索引
    // Phi_mz_xn: [n_pml_xn] x [Ny] x [Nz+1]
    inline int id_Phi_mz_xn(int l, int jy, int kz) const { return (l*m_Ny*(m_Nz+1) + jy*(m_Nz+1) + kz); }
    inline int id_Phi_mz_xp(int l, int jy, int kz) const { return (l*m_Ny*(m_Nz+1) + jy*(m_Nz+1) + kz); }
    // Phi_mz_yn: [Nx] x [n_pml_yn] x [Nz+1]
    inline int id_Phi_mz_yn(int ix, int l, int kz) const { return (ix*n_pml_yn*(m_Nz+1) + l*(m_Nz+1) + kz); }
    inline int id_Phi_mz_yp(int ix, int l, int kz) const { return (ix*n_pml_yp*(m_Nz+1) + l*(m_Nz+1) + kz); }

    // ========== 初始化 ==========
    void setupGrid();
    void setupSource();
    void setupMaterials();
    void setupPML();

    // ========== 时间推进 ==========
    void updateE(int step);
    void updateH();
    void applySource(int step);
    void updateElectricMaterialRelations();
    void updateMagneticMaterialRelations();

    // ========== 输出 ==========
    void dumpFrame(int step);
};