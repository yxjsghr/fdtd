#include "OutputManager.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

OutputManager::OutputManager(const std::string& directory) {
    namespace fs = std::filesystem;
    fs::path abs = fs::absolute(directory);
    m_base_directory = abs.string();
    try {
        if (fs::exists(m_base_directory)) {
            fs::remove_all(m_base_directory);
        }
        fs::create_directory(m_base_directory);
    } catch (const fs::filesystem_error& e) {
        std::cerr << "[OutputManager] Filesystem error: " << e.what() << std::endl;
    }
}

void OutputManager::saveScalar3D_VTI(const std::string& fname,
                                     const std::vector<double>& field,
                                     int Nx, int Ny, int Nz,
                                     double dx, double dy, double dz,
                                     int step) {
    // 构建文件名
    std::ostringstream oss;
    if (step >= 0) {
        oss << m_base_directory << "/" << fname << "_t"
            << std::setw(5) << std::setfill('0') << step << ".vti";
    } else {
        oss << m_base_directory << "/" << fname << ".vti";
    }
    std::string full_path = oss.str();

    std::ofstream ofs(full_path);
    if (!ofs.is_open()) {
        std::cerr << "[OutputManager] Failed to open file: " << full_path << std::endl;
        return;
    }

    // VTK ImageData 要求：WholeExtent 是点索引范围 [x0, x1, y0, y1, z0, z1]
    // 对于 Nx × Ny × Nz 个点，extent 为 [0, Nx-1, 0, Ny-1, 0, Nz-1]
    int x0 = 0, x1 = Nx - 1;
    int y0 = 0, y1 = Ny - 1;
    int z0 = 0, z1 = Nz - 1;

    ofs << "<?xml version=\"1.0\"?>\n";
    ofs << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    ofs << "  <ImageData WholeExtent=\"" 
        << x0 << " " << x1 << " "
        << y0 << " " << y1 << " "
        << z0 << " " << z1 << "\" "
        << "Origin=\"0 0 0\" "
        << "Spacing=\"" << dx << " " << dy << " " << dz << "\">\n";
    ofs << "    <Piece Extent=\"" 
        << x0 << " " << x1 << " "
        << y0 << " " << y1 << " "
        << z0 << " " << z1 << "\">\n";
    ofs << "      <PointData Scalars=\"" << fname << "\">\n";
    ofs << "        <DataArray type=\"Float64\" Name=\"" << fname 
        << "\" format=\"ascii\" NumberOfComponents=\"1\">\n";

    // 写入数据：VTK 要求按 z-slow, y-medium, x-fast 顺序（即 Fortran 列优先）
    // 但我们假设输入 field 是 C 风格：index = ((k * Ny) + j) * Nx + i （x 快变）
    // 这与 VTK 的 PointData 存储顺序一致（VTK 也按 x-fast, y, z 存储点数据）
    ofs << std::setprecision(16);
    for (size_t idx = 0; idx < field.size(); ++idx) {
        if (idx > 0 && idx % 4 == 0) ofs << "\n"; // 每行 4 个数，美观
        ofs << field[idx] << " ";
    }
    if (!field.empty()) ofs << "\n";

    ofs << "        </DataArray>\n";
    ofs << "      </PointData>\n";
    ofs << "    </Piece>\n";
    ofs << "  </ImageData>\n";
    ofs << "</VTKFile>\n";

    ofs.close();
    std::cout << "[OutputManager] Saved VTI: " << full_path << std::endl;
}

// ========== 2D 兼容接口（保留）==========
void OutputManager::saveScalar2D_CSV(const std::string& fname,
                                     const std::vector<double>& field,
                                     int Nx, int Nz, double dx, double dz) {
    std::ofstream ofs(m_base_directory + "/" + fname);
    if (!ofs.is_open()) return;
    ofs << "# x,z,val\n";
    ofs << std::setprecision(16);
    for (int ix = 0; ix < Nx; ++ix) {
        for (int kz = 0; kz < Nz; ++kz) {
            int id = ix * Nz + kz;
            if (static_cast<size_t>(id) < field.size()) {
                ofs << (ix * dx) << "," << (kz * dz) << "," << field[id] << "\n";
            }
        }
    }
}

void OutputManager::saveDimsTxt2D(int Nx, int Nz, double dx, double dz) {
    std::ofstream ofs(m_base_directory + "/dims.txt");
    if (!ofs.is_open()) return;
    ofs << "Nx " << Nx << "\n";
    ofs << "Nz " << Nz << "\n";
    ofs << std::setprecision(16);
    ofs << "dx " << dx << "\n";
    ofs << "dz " << dz << "\n";
}

// ========== 1D 兼容接口（保留）==========
void OutputManager::saveEzCSV(int step, const std::vector<double>& ez, double ds) {
    std::ostringstream oss;
    oss << m_base_directory << "/ez_t" << std::setw(5) << std::setfill('0') << step << ".csv";
    std::ofstream ofs(oss.str());
    if (!ofs.is_open()) return;
    ofs << "# z,Ez\n";
    ofs << std::setprecision(16);
    for (size_t i = 0; i < ez.size(); ++i) {
        ofs << (i * ds) << "," << ez[i] << "\n";
    }
}

void OutputManager::saveExMagCSV(int step, const std::vector<double>& exmag, double ds) {
    std::ostringstream oss;
    oss << m_base_directory << "/exmag_t" << std::setw(5) << std::setfill('0') << step << ".csv";
    std::ofstream ofs(oss.str());
    if (!ofs.is_open()) return;
    ofs << "# z,Ex_mag\n";
    ofs << std::setprecision(16);
    for (size_t i = 0; i < exmag.size(); ++i) {
        ofs << (i * ds) << "," << exmag[i] << "\n";
    }
}
