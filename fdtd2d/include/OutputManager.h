#pragma once
#include <string>
#include <vector>

class OutputManager {
public:
    explicit OutputManager(const std::string& directory);

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
