#pragma once
#include <memory>
#include <string>
#include <vector>
#include <cmath>

class OutputManager;

class FDTD2D {
public:
    FDTD2D();
    ~FDTD2D(); // 只声明，在 cpp 中定义（避免不完全类型）

    void run();

private:
    // ========== 常量 ==========
    double m_c0  = 299792458.0;
    double m_mu0 = 4.0 * M_PI * 1e-7;
    double m_eps0 = 1.0 / (m_c0*m_c0*m_mu0);

    // ========== 网格 ==========
    int m_Nx = 0, m_Nz = 0;      // 单元数（与 MATLAB 对齐）
    double m_dx = 0.0, m_dz = 0.0;
    int m_steps = 0;
    double m_dt = 0.0;

    // ========== Lorentz 系数（定义在 i1=1..Nx, k1=1..Nz 上） ==========
    std::vector<double> C1, C2, C3; // 尺寸 Nx * Nz

    // ========== 场 ==========
    // 扁平化索引（x 快变，z 慢变）
    // Dy, Ey: (Nx+1) x (Nz+1)
    std::vector<double> Dy, Ey, Ey_n1, Ey_n2, SEy, SEy_n1, SEy_n2;
    // Bx, Hx: (Nx+1) x Nz
    std::vector<double> Bx, Hx, Hx_n1, Hx_n2, SHx, SHx_n1, SHx_n2;
    // Bz, Hz: Nx x (Nz+1)
    std::vector<double> Bz, Hz, Hz_n1, Hz_n2, SHz, SHz_n1, SHz_n2;

    // ========== CPML 参数 ==========
    int n_pml_xn = 20, n_pml_xp = 20, n_pml_zn = 20, n_pml_zp = 20;
    // 1D 侧向剖面系数（按 MATLAB 写法）
    std::vector<double> b_e_zn, a_e_zn, C_Dy_dzn_v;
    std::vector<double> b_e_zp, a_e_zp, C_Dy_dzp_v;
    std::vector<double> b_e_xn, a_e_xn, C_Dy_dxn_v;
    std::vector<double> b_e_xp, a_e_xp, C_Dy_dxp_v;

    std::vector<double> b_m_zn, a_m_zn, C_Bx_dzn_v;
    std::vector<double> b_m_zp, a_m_zp, C_Bx_dzp_v;
    std::vector<double> b_m_xn, a_m_xn, C_Bz_dxn_v;
    std::vector<double> b_m_xp, a_m_xp, C_Bz_dxp_v;

    // ϕ 辅助量：形状与 MATLAB 一致（扁平存储）
    //  z向 PML —— 宽度：Nx-1（i2），高度：n_pml_z*
    std::vector<double> Phi_ey_zn, Phi_ey_zp; // (Nx-1)*n_pml_z*
    std::vector<double> Phi_mx_zn, Phi_mx_zp; // (Nx-1)*n_pml_z*
    //  x向 PML —— 宽度：n_pml_x*，高度：Nz-1（k2）
    std::vector<double> Phi_ey_xn, Phi_ey_xp; // (n_pml_x*)*(Nz-1)
    std::vector<double> Phi_mz_xn, Phi_mz_xp; // (n_pml_x*)*(Nz-1)

    // ========== 源 ==========
    int i_Jy = 0, k_Jy = 0;
    std::vector<double> Jy;

    // ========== 输出 ==========
    std::unique_ptr<OutputManager> m_out;

private:
    // 扁平索引
    inline int id_Dy(int ix, int kz) const { return ix*(m_Nz+1) + kz; }     // 0..Nx, 0..Nz
    inline int id_Ey(int ix, int kz) const { return ix*(m_Nz+1) + kz; }
    inline int id_Bx(int ix, int kz) const { return ix*m_Nz + kz; }         // 0..Nx, 0..Nz-1
    inline int id_Hx(int ix, int kz) const { return ix*m_Nz + kz; }
    inline int id_Bz(int ix, int kz) const { return ix*(m_Nz+1) + kz; }     // 0..Nx-1, 0..Nz
    inline int id_Hz(int ix, int kz) const { return ix*(m_Nz+1) + kz; }
    inline int id_C (int ix, int kz) const { return ix*m_Nz + kz; }         // 0..Nx-1, 0..Nz-1

    // CPML ϕ 索引
    inline int id_x_pml(int i0, int l, int nlayer) const { return i0*nlayer + l; }       // (Nx-1) x nlayer
    inline int id_z_pml(int l, int k0, int Nz_1)  const { return l*Nz_1 + k0; }          // nlayer x (Nz-1)

    // 初始化
    void setupGrid();
    void setupSource();
    void setupMaterialsLorentz();
    void setupPML();

    // 单步推进
    void stepOnce(int n);

    // 输出
    void dumpFrame(int n);
};
