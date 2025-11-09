#include "OutputManager.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

OutputManager::OutputManager(const std::string& directory) {
    namespace fs = std::filesystem;
    fs::path abs = fs::absolute(directory);
    m_base_directory = abs.string();
    try {
        if (fs::exists(m_base_directory)) fs::remove_all(m_base_directory);
        fs::create_directory(m_base_directory);
    } catch (const fs::filesystem_error& e) {
        std::cerr << "[OutputManager] Filesystem error: " << e.what() << std::endl;
    }
}

void OutputManager::saveScalar2D_CSV(const std::string& fname,
                                     const std::vector<double>& field,
                                     int Nx, int Nz, double dx, double dz) {
    std::ofstream ofs(m_base_directory + "/" + fname);
    if (!ofs.is_open()) return;
    ofs << "# x,z,val\n";
    ofs << std::setprecision(16);
    for (int ix = 0; ix < Nx; ++ix) {
        for (int kz = 0; kz < Nz; ++kz) {
            int id = ix*Nz + kz;  // 调用方保证 field 的布局是 [Nx][Nz]
            ofs << (ix*dx) << "," << (kz*dz) << "," << field[id] << "\n";
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

/************ 1D 兼容：不使用也保留 ************/
void OutputManager::saveEzCSV(int step, const std::vector<double>& ez, double ds) {
    std::ostringstream oss;
    oss << m_base_directory << "/ez_t" << std::setw(5) << std::setfill('0') << step << ".csv";
    std::ofstream ofs(oss.str());
    if (!ofs.is_open()) return;
    ofs << "# z,Ez\n";
    ofs << std::setprecision(16);
    for (size_t i = 0; i < ez.size(); ++i) {
        ofs << (i*ds) << "," << ez[i] << "\n";
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
        ofs << (i*ds) << "," << exmag[i] << "\n";
    }
}
