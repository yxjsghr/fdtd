#!/bin/bash

# FDTD3D 编译和运行脚本

set -e  # 遇到错误立即退出

# 获取脚本所在目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=========================================="
echo "FDTD3D 编译和运行脚本"
echo "=========================================="

# 创建构建目录
BUILD_DIR="build"
if [ -d "$BUILD_DIR" ]; then
    echo "删除旧的构建目录: $BUILD_DIR"
    rm -rf "$BUILD_DIR"
fi

echo "创建构建目录: $BUILD_DIR"
mkdir -p "$BUILD_DIR"

cd "$BUILD_DIR"

# 运行 CMake 配置
echo ""
echo "配置 CMake..."
cmake .. -DCMAKE_BUILD_TYPE=Release

# 编译
echo ""
echo "开始编译..."
make -j$(nproc)

# 检查可执行文件是否存在
if [ ! -f "fdtd3d" ]; then
    echo "错误: 编译失败，未找到可执行文件 fdtd3d"
    exit 1
fi

echo ""
echo "=========================================="
echo "编译成功！"
echo "=========================================="
echo ""

# 运行程序
echo "开始运行 FDTD3D 仿真..."
echo ""
./fdtd3d

echo ""
echo "=========================================="
echo "运行完成！"
echo "=========================================="
echo "输出文件位于: $BUILD_DIR/output_3d/"

