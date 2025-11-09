
#include <string>
#include <vector>

class OutputManager {
public:
    explicit OutputManager(const std::string& directory);

    // ========== 3D 输出：VTK ImageData (.vti) ==========
    /**
     * @brief 保存 3D 标量场为 VTK .vti 文件（ParaView 可直接打开）
     * 
     * 场量 field 的逻辑维度为 [Nx][Ny][Nz]，按 x-fast, y-medium, z-slow 扁平化存储
     * 即：index = (k * Ny + j) * Nx + i  对应 (x=i, y=j, z=k)
     * 
     * @param fname 输出文件名（不含路径，自动加 .vti）
     * @param field 扁平化的 3D 标量场数据（size = Nx * Ny * Nz）
     * @param Nx, Ny, Nz 网格点数（物理单元数）
     * @param dx, dy, dz 网格步长
     * @param step 时间步（用于命名或元数据）
     */
    void saveScalar3D_VTI(const std::string& fname,
                          const std::vector<double>& field,
                          int Nx, int Ny, int Nz,
                          double dx, double dy, double dz,
                          int step = -1);

    // ========== 兼容性接口（避免链接错误，即使不使用）==========
    // 2D: 保存标量场（每行：x,z,val）
    void saveScalar2D_CSV(const std::string& fname,
                          const std::vector<double>& field,
                          int Nx, int Nz, double dx, double dz);

    // 2D: 保存尺寸
    void saveDimsTxt2D(int Nx, int Nz, double dx, double dz);

    // 1D 兼容（即使不使用，也保留以避免链接缺符号）
    void saveEzCSV(int step, const std::vector<double>& ez, double ds);
    void saveExMagCSV(int step, const std::vector<double>& exmag, double ds);

private:
    std::string m_base_directory;
};