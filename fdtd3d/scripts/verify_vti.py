#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FDTD3D VTI文件验证脚本
- 读取VTI文件并检查数据合理性
- 显示统计信息和可视化
"""
import argparse
import os
import re
import sys
import numpy as np
try:
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    HAS_VTK = True
except ImportError:
    HAS_VTK = False
    print("[警告] 未安装vtk，将使用XML解析方式")

try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("[警告] 未安装matplotlib，将跳过可视化")


TIME_STEP_PATTERN = re.compile(r'_t(\d+)\.vti$')
MAX_SERIES_STEPS = 1500


def extract_time_step(path):
    """提取文件名中的时间步编号，若不存在则返回None"""
    match = TIME_STEP_PATTERN.search(os.path.basename(path))
    return int(match.group(1)) if match else None


def sort_vti_files(files):
    """按照时间步编号对VTI文件排序，无法解析的排在最后"""
    def sort_key(path):
        time_step = extract_time_step(path)
        return (0, time_step) if time_step is not None else (1, os.path.basename(path))

    return sorted(files, key=sort_key)

def read_vti_xml(vti_path):
    """直接从XML格式的VTI文件读取数据（不依赖vtk库）"""
    import xml.etree.ElementTree as ET
    
    tree = ET.parse(vti_path)
    root = tree.getroot()
    
    # 获取维度信息
    image_data = root.find('ImageData')
    if image_data is None:
        raise ValueError("未找到ImageData元素")
    
    whole_extent = [int(x) for x in image_data.get('WholeExtent', '').split()]
    origin = [float(x) for x in image_data.get('Origin', '0 0 0').split()]
    spacing = [float(x) for x in image_data.get('Spacing', '1 1 1').split()]
    
    # 计算维度
    dims = [whole_extent[1] - whole_extent[0] + 1,
            whole_extent[3] - whole_extent[2] + 1,
            whole_extent[5] - whole_extent[4] + 1]
    
    # 查找数据数组
    data_array = root.find('.//DataArray')
    if data_array is None:
        raise ValueError("未找到DataArray元素")
    
    # 读取数据
    data_str = data_array.text.strip()
    data = np.fromstring(data_str, sep=' ')
    
    # 重塑为3D数组
    data = data.reshape(dims[2], dims[1], dims[0])  # VTK使用(z,y,x)顺序
    data = np.transpose(data, (2, 1, 0))  # 转换为(x,y,z)顺序
    
    return data, dims, origin, spacing

def read_vti_vtk(vti_path):
    """使用VTK库读取VTI文件"""
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(vti_path)
    reader.Update()
    
    image_data = reader.GetOutput()
    dims = image_data.GetDimensions()
    origin = image_data.GetOrigin()
    spacing = image_data.GetSpacing()
    
    # 获取数据
    point_data = image_data.GetPointData()
    if point_data.GetNumberOfArrays() == 0:
        raise ValueError("VTI文件中没有数据数组")
    
    array = point_data.GetArray(0)
    data = vtk_to_numpy(array)
    
    # 重塑为3D数组 (VTK使用Fortran顺序)
    data = data.reshape(dims[2], dims[1], dims[0])
    data = np.transpose(data, (2, 1, 0))  # 转换为(x,y,z)顺序
    
    return data, dims, origin, spacing

def read_vti(vti_path):
    """读取VTI文件（自动选择方法）"""
    if HAS_VTK:
        try:
            return read_vti_vtk(vti_path)
        except Exception as e:
            print(f"[警告] VTK读取失败: {e}，尝试XML解析")
            return read_vti_xml(vti_path)
    else:
        return read_vti_xml(vti_path)

def analyze_data(data, name="数据"):
    """分析数据并打印统计信息"""
    print(f"\n{'='*60}")
    print(f"{name} 统计信息:")
    print(f"{'='*60}")
    print(f"  形状 (Nx, Ny, Nz): {data.shape}")
    print(f"  总点数: {data.size:,}")
    print(f"  数据类型: {data.dtype}")
    print(f"\n  数值范围:")
    print(f"    最小值: {np.min(data):.6e}")
    print(f"    最大值: {np.max(data):.6e}")
    print(f"    平均值: {np.mean(data):.6e}")
    print(f"    标准差: {np.std(data):.6e}")
    print(f"\n  分位数:")
    print(f"    1%:   {np.percentile(data, 1):.6e}")
    print(f"    25%:  {np.percentile(data, 25):.6e}")
    print(f"    50%:  {np.percentile(data, 50):.6e}")
    print(f"    75%:  {np.percentile(data, 75):.6e}")
    print(f"    99%:  {np.percentile(data, 99):.6e}")
    
    # 检查异常值
    finite_mask = np.isfinite(data)
    if not np.all(finite_mask):
        n_inf = np.sum(~finite_mask)
        print(f"\n  [警告] 发现 {n_inf} 个非有限值 (NaN/Inf)")
    
    # 检查能量守恒（RMS应该随时间增长或稳定）
    total_energy = np.sum(data**2)
    print(f"\n  总能量 (sum(E²)): {total_energy:.6e}")
    
    return {
        'shape': data.shape,
        'min': np.min(data),
        'max': np.max(data),
        'mean': np.mean(data),
        'std': np.std(data),
        'total_energy': total_energy
    }

def plot_slices(data, dims, output_dir=None, prefix="slice"):
    """绘制2D切片"""
    if not HAS_MATPLOTLIB:
        return
    
    Nx, Ny, Nz = dims
    
    # 选择中心切片
    mid_x, mid_y, mid_z = Nx//2, Ny//2, Nz//2
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # XY平面 (z=mid_z)
    im1 = axes[0].imshow(data[:, :, mid_z], origin='lower', cmap='jet', aspect='auto')
    axes[0].set_title(f'XY平面 (z={mid_z}/{Nz-1})')
    axes[0].set_xlabel('Y')
    axes[0].set_ylabel('X')
    plt.colorbar(im1, ax=axes[0])
    
    # XZ平面 (y=mid_y)
    im2 = axes[1].imshow(data[:, mid_y, :], origin='lower', cmap='jet', aspect='auto')
    axes[1].set_title(f'XZ平面 (y={mid_y}/{Ny-1})')
    axes[1].set_xlabel('Z')
    axes[1].set_ylabel('X')
    plt.colorbar(im2, ax=axes[1])
    
    # YZ平面 (x=mid_x)
    im3 = axes[2].imshow(data[mid_x, :, :], origin='lower', cmap='jet', aspect='auto')
    axes[2].set_title(f'YZ平面 (x={mid_x}/{Nx-1})')
    axes[2].set_xlabel('Z')
    axes[2].set_ylabel('Y')
    plt.colorbar(im3, ax=axes[2])
    
    plt.tight_layout()
    
    if output_dir:
        output_path = os.path.join(output_dir, f"{prefix}_slices.png")
        plt.savefig(output_path, dpi=150)
        print(f"\n[*] 已保存切片图: {output_path}")
    else:
        plt.show()
    
    plt.close()

def verify_consistency(files, output_dir):
    """验证多个时间步的一致性"""
    print(f"\n{'='*60}")
    print("验证时间序列一致性...")
    print(f"{'='*60}")
    
    stats_list = []
    times = []
    
    sorted_files = sort_vti_files(files)
    max_steps = min(MAX_SERIES_STEPS, len(sorted_files))
    
    for vti_file in sorted_files[:max_steps]:
        # 提取时间步
        match = re.search(r'_t(\d+)\.vti', vti_file)
        if match:
            time_step = int(match.group(1))
            times.append(time_step)
        else:
            continue
        
        try:
            data, dims, origin, spacing = read_vti(vti_file)
            stats = analyze_data(data, f"时间步 {time_step}")
            stats['time'] = time_step
            stats_list.append(stats)
        except Exception as e:
            print(f"[错误] 读取 {os.path.basename(vti_file)} 失败: {e}")
            continue
    
    # 绘制能量演化
    if HAS_MATPLOTLIB and stats_list:
        energies = [s['total_energy'] for s in stats_list]
        time_steps = [s['time'] for s in stats_list]
        
        plt.figure(figsize=(10, 6))
        plt.plot(time_steps, energies, 'b-o', linewidth=2, markersize=6)
        plt.xlabel('时间步')
        plt.ylabel('总能量 (sum(E²))')
        plt.title('能量演化')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        if output_dir:
            energy_plot = os.path.join(output_dir, "energy_evolution.png")
            plt.savefig(energy_plot, dpi=150)
            print(f"\n[*] 已保存能量演化图: {energy_plot}")
        else:
            plt.show()
        
        plt.close()
        
        # 检查能量是否合理增长
        if len(energies) > 1:
            energy_growth = energies[-1] / energies[0] if energies[0] > 0 else 0
            print(f"\n  能量增长比 (最后/最初): {energy_growth:.4f}")
            if energy_growth > 100:
                print(f"  [警告] 能量增长过快，可能有问题")
            elif energy_growth < 0.1:
                print(f"  [警告] 能量衰减过快，可能有问题")

def main():
    parser = argparse.ArgumentParser(description="FDTD3D VTI文件验证工具")
    parser.add_argument("--file", type=str, help="单个VTI文件路径")
    parser.add_argument("--dir", type=str, help="包含VTI文件的目录")
    parser.add_argument("--pattern", type=str, default="E_mag.*\\.vti", 
                       help="文件匹配模式（正则表达式）")
    parser.add_argument("--plot", action="store_true", help="生成可视化图像")
    parser.add_argument("--verify-series", action="store_true", 
                       help="验证时间序列一致性")
    
    args = parser.parse_args()
    
    if args.file:
        # 验证单个文件
        if not os.path.isfile(args.file):
            print(f"[错误] 文件不存在: {args.file}")
            return 1
        
        print(f"读取文件: {args.file}")
        try:
            data, dims, origin, spacing = read_vti(args.file)
            analyze_data(data, os.path.basename(args.file))
            
            if args.plot:
                output_dir = os.path.dirname(args.file)
                plot_slices(data, dims, output_dir, 
                          prefix=os.path.splitext(os.path.basename(args.file))[0])
        except Exception as e:
            print(f"[错误] 读取失败: {e}")
            import traceback
            traceback.print_exc()
            return 1
    
    elif args.dir:
        # 验证目录中的文件
        if not os.path.isdir(args.dir):
            print(f"[错误] 目录不存在: {args.dir}")
            return 1
        
        # 查找匹配的文件
        pattern = re.compile(args.pattern)
        files = [os.path.join(args.dir, f) for f in os.listdir(args.dir) 
                if pattern.match(f) and f.endswith('.vti')]
        files = sort_vti_files(files)
        
        if not files:
            print(f"[错误] 在 {args.dir} 中未找到匹配的VTI文件")
            return 1
        
        print(f"找到 {len(files)} 个VTI文件")
        
        # 验证第一个文件
        print(f"\n验证第一个文件: {os.path.basename(files[0])}")
        try:
            data, dims, origin, spacing = read_vti(files[0])
            analyze_data(data, os.path.basename(files[0]))
            
            print(f"\n网格信息:")
            print(f"  维度: {dims}")
            print(f"  原点: {origin}")
            print(f"  间距: {spacing}")
            
            if args.plot:
                plot_slices(data, dims, args.dir, 
                          prefix=os.path.splitext(os.path.basename(files[0]))[0])
            
            # 验证时间序列
            if args.verify_series:
                verify_consistency(files, args.dir)
        except Exception as e:
            print(f"[错误] 处理失败: {e}")
            import traceback
            traceback.print_exc()
            return 1
    else:
        parser.print_help()
        return 1
    
    print(f"\n{'='*60}")
    print("验证完成！")
    print(f"{'='*60}")
    return 0

if __name__ == "__main__":
    sys.exit(main())

