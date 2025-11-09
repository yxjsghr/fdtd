#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""FDTD3D 可视化脚本（MATLAB 样式 2x2 切片图）

- 支持从 VTI 文件读取单个时间步的 Ex 与 E_mag（振幅）场
- 参考 MATLAB `S6/S7` 文件的布局生成四个子图
- 允许通过命令行参数灵活设置 PML 厚度、色标、视角、刻度等

示例:
    python visualize_slices.py \
        --ex-file build/output_3d/Ex_t00540.vti \
        --emag-file build/output_3d/E_mag_rms_t00540.vti \
        --output build/output_3d/slice_layout_t00540.png
"""

import argparse
import glob
import os
import re
import shutil
import sys
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

try:
    import matplotlib.pyplot as plt
    from matplotlib import animation as mpl_animation
    from matplotlib.colors import Normalize
    HAS_MPL = True
except ImportError:  # pragma: no cover - 可视化环境缺失
    HAS_MPL = False
    plt = None
    mpl_animation = None

try:
    import vtk  # type: ignore
    from vtk.util.numpy_support import vtk_to_numpy  # type: ignore
    HAS_VTK = True
except ImportError:
    HAS_VTK = False

import xml.etree.ElementTree as ET


DEFAULT_EXPORT_DIR = os.path.join('build', 'visualizations')


# ---------------------------------------------------------------------------
# VTI 读取工具（与 scripts/verify_vti.py 中的逻辑保持一致，但独立拷贝，避免 import 副作用）
# ---------------------------------------------------------------------------

def _read_vti_xml(path: str) -> Tuple[np.ndarray, Tuple[int, int, int], Tuple[float, float, float]]:
    tree = ET.parse(path)
    root = tree.getroot()
    image_data = root.find('ImageData')
    if image_data is None:
        raise ValueError(f"{os.path.basename(path)} 缺少 ImageData 节点")

    whole_extent = [int(x) for x in image_data.get('WholeExtent', '').split()]
    spacing = tuple(float(x) for x in image_data.get('Spacing', '1 1 1').split())

    dims = (
        whole_extent[1] - whole_extent[0] + 1,
        whole_extent[3] - whole_extent[2] + 1,
        whole_extent[5] - whole_extent[4] + 1,
    )

    data_array = root.find('.//DataArray')
    if data_array is None:
        raise ValueError(f"{os.path.basename(path)} 缺少 DataArray 数据")

    data_str = data_array.text.strip()
    data = np.fromstring(data_str, sep=' ')
    data = data.reshape(dims[2], dims[1], dims[0])  # vtk:z,y,x
    data = np.transpose(data, (2, 1, 0))  # -> (x,y,z)
    return data, dims, spacing


def _read_vti_vtk(path: str) -> Tuple[np.ndarray, Tuple[int, int, int], Tuple[float, float, float]]:
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(path)
    reader.Update()

    image_data = reader.GetOutput()
    dims = image_data.GetDimensions()
    spacing = image_data.GetSpacing()

    point_data = image_data.GetPointData()
    if point_data.GetNumberOfArrays() == 0:
        raise ValueError(f"{os.path.basename(path)} 不包含数据数组")

    array = point_data.GetArray(0)
    data = vtk_to_numpy(array)
    data = data.reshape(dims[2], dims[1], dims[0])
    data = np.transpose(data, (2, 1, 0))
    return data, dims, spacing


def read_vti(path: str) -> Tuple[np.ndarray, Tuple[int, int, int], Tuple[float, float, float]]:
    if HAS_VTK:
        try:
            return _read_vti_vtk(path)
        except Exception as exc:  # pragma: no cover - 运行时备援
            print(f"[警告] VTK 解析 {os.path.basename(path)} 失败 ({exc})，回退到 XML 读取")
    return _read_vti_xml(path)


# ---------------------------------------------------------------------------
# 可视化核心逻辑
# ---------------------------------------------------------------------------

def _slice_bounds(length: int, n_neg: int, n_pos: int) -> slice:
    start = max(0, n_neg)
    stop = max(start + 1, length - max(0, n_pos))
    return slice(start, stop)


def _center_index(bounds: slice) -> int:
    return bounds.start + (bounds.stop - bounds.start) // 2


def _extract_plane_xz(data: np.ndarray, y_idx: int, x_bounds: slice, z_bounds: slice) -> np.ndarray:
    plane = data[x_bounds, y_idx, z_bounds]
    return np.array(plane)


def _extract_plane_yz(data: np.ndarray, x_idx: int, y_bounds: slice, z_bounds: slice) -> np.ndarray:
    plane = data[x_idx, y_bounds, z_bounds]
    return np.array(plane)


def _extent_from_bounds(bounds: slice) -> Tuple[float, float]:
    return 0.0, float(bounds.stop - bounds.start)


def _apply_lambda_ticks(ax, length: int, axis: str, n_lambda: int, center_zero: bool,
                        spacing: float = 1.0, offset: float = 0.0) -> None:
    if n_lambda <= 0:
        return
    ticks = np.arange(0, length + 1, n_lambda)
    if not len(ticks):
        return
    tick_positions = (ticks + offset) * spacing
    if center_zero:
        center = (offset + length / 2.0) * spacing
        labels = [(pos - center) / (n_lambda * spacing) for pos in tick_positions]
    else:
        labels = [pos / (n_lambda * spacing) for pos in tick_positions]
    label_str = [f"{val:g}" for val in labels]
    if axis == 'x':
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(label_str)
    else:
        ax.set_yticks(tick_positions)
        ax.set_yticklabels(label_str)


TIME_STEP_PATTERN = re.compile(r'_t(\d+)\.vti$')


def _extract_time_step(path: str) -> Optional[int]:
    match = TIME_STEP_PATTERN.search(os.path.basename(path))
    return int(match.group(1)) if match else None


def _collect_series_map(pattern: Optional[str]) -> Dict[int, str]:
    if not pattern:
        return {}
    mapping: Dict[int, str] = {}
    for path in sorted(glob.glob(pattern)):
        step = _extract_time_step(path)
        if step is not None:
            mapping[step] = path
    return mapping


def _prepare_slices(ex_data: Optional[np.ndarray], emag_data: Optional[np.ndarray],
                    range_x: slice, range_y: slice, range_z: slice) -> Dict[str, Optional[np.ndarray]]:
    slices: Dict[str, Optional[np.ndarray]] = {
        'ex_xz': None,
        'emag_xz': None,
        'ex_yz': None,
        'emag_yz': None,
    }

    mid_y = _center_index(range_y)
    mid_x = _center_index(range_x)

    if ex_data is not None:
        slices['ex_xz'] = _extract_plane_xz(ex_data, mid_y, range_x, range_z)
        slices['ex_yz'] = _extract_plane_yz(ex_data, mid_x, range_y, range_z)

    if emag_data is not None:
        slices['emag_xz'] = _extract_plane_xz(emag_data, mid_y, range_x, range_z)
        slices['emag_yz'] = _extract_plane_yz(emag_data, mid_x, range_y, range_z)
    elif ex_data is not None:
        mag_data = np.abs(ex_data)
        slices['emag_xz'] = _extract_plane_xz(mag_data, mid_y, range_x, range_z)
        slices['emag_yz'] = _extract_plane_yz(mag_data, mid_x, range_y, range_z)

    return slices


def _select_animation_steps(steps: List[int], interval: int, max_frames: int) -> List[int]:
    if not steps:
        return []
    selected: List[int] = []
    last = None
    for step in steps:
        if last is None or interval <= 1 or step - last >= interval:
            selected.append(step)
            last = step
        if len(selected) >= max_frames:
            break
    if not selected:
        # fallback: at least keep first element if interval too strict
        return steps[:min(len(steps), max_frames)]
    return selected


def _default_animation_path(args: argparse.Namespace, fmt: str) -> str:
    base_dir_candidates = [
        os.path.dirname(args.animate_output or ''),
        os.path.dirname(args.output or ''),
        os.path.dirname(args.ex_file or ''),
        os.path.dirname(args.emag_file or ''),
        DEFAULT_EXPORT_DIR,
    ]
    base_dir = next((d for d in base_dir_candidates if d), DEFAULT_EXPORT_DIR)
    if not os.path.isdir(base_dir):
        os.makedirs(base_dir, exist_ok=True)
    base_name = 'fdtd3d_slices'
    suffix = '.mp4' if fmt == 'mp4' else '.gif'
    return os.path.join(base_dir, base_name + suffix)


def _default_static_path(args: argparse.Namespace) -> str:
    os.makedirs(DEFAULT_EXPORT_DIR, exist_ok=True)
    step = None
    for path in (args.ex_file, args.emag_file):
        if path:
            step = _extract_time_step(path)
            if step is not None:
                break
    suffix = f'_t{step:05d}' if step is not None else ''
    return os.path.join(DEFAULT_EXPORT_DIR, f'slice_layout{suffix}.png')


def render_animation(ex_series: Dict[int, str], emag_series: Dict[int, str],
                     bounds: Tuple[slice, slice, slice], args: argparse.Namespace) -> None:
    if not HAS_MPL:
        raise RuntimeError('matplotlib 未安装，无法生成动画。')

    all_steps = sorted(set(ex_series.keys()) | set(emag_series.keys()))
    if not all_steps:
        print('[WARN] --animate requires --ex-series or --emag-series; nothing matched.')
        return

    steps = _select_animation_steps(all_steps, max(args.animate_interval, 1), args.animate_max_frames)
    print(f'[*] animation: {len(steps)} frames selected (interval={args.animate_interval}, max={args.animate_max_frames})')
    first_step = steps[0]

    spacing_cache: Optional[Tuple[float, float, float]] = None

    def load_step(step: int) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        print(f'[debug] loading step {step} ...', flush=True)
        ex_data = None
        emag_data = None
        candidate_spacing: Optional[Tuple[float, float, float]] = None
        if step in ex_series:
            ex_data, _, ex_spacing = read_vti(ex_series[step])
            candidate_spacing = ex_spacing
        if step in emag_series:
            emag_data, _, emag_spacing = read_vti(emag_series[step])
            if candidate_spacing is None:
                candidate_spacing = emag_spacing
            elif candidate_spacing != emag_spacing:
                print(f'[WARN] spacing mismatch between Ex 和 E_mag at step {step}: {candidate_spacing} vs {emag_spacing}')
        if candidate_spacing is not None:
            nonlocal spacing_cache
            if spacing_cache is None:
                spacing_cache = candidate_spacing
            elif spacing_cache != candidate_spacing:
                print(f'[WARN] spacing mismatch across steps: cached {spacing_cache} vs {candidate_spacing}')
        print(f'[debug] loaded step {step}: ex={ex_data is not None}, emag={emag_data is not None}', flush=True)
        return ex_data, emag_data

    ex_first, emag_first = load_step(first_step)
    if spacing_cache is None:
        spacing_cache = (1.0, 1.0, 1.0)
        print('[WARN] 无法获取 VTI spacing，使用默认 (1,1,1)。')
    fig, images = build_figure(ex_first, emag_first, bounds, spacing_cache, args)
    print('[debug] initial figure built', flush=True)

    def update(frame_idx: int):
        step = steps[frame_idx]
        print(f'[debug] updating frame index {frame_idx} (step {step})', flush=True)
        ex_frame, emag_frame = load_step(step)
        slices = _prepare_slices(ex_frame, emag_frame, *bounds)
        artists = []
        for key, im in images.items():
            if im is None or slices[key] is None:
                continue
            im.set_data(slices[key])
            artists.append(im)
        fig.suptitle(f't = {step}')
        print(f'[debug] frame {frame_idx} ready', flush=True)
        return artists

    ani = mpl_animation.FuncAnimation(fig, update, frames=len(steps), interval=1000 // max(args.animate_fps, 1), blit=False)

    out_fmt = args.animate_format
    if out_fmt == 'mp4' and shutil.which('ffmpeg') is None:
        print('[WARN] ffmpeg not found; switching animation format to GIF')
        out_fmt = 'gif'
    out_path = args.animate_output
    if not out_path:
        out_path = _default_animation_path(args, out_fmt)
    else:
        out_dir = os.path.dirname(out_path)
        if out_dir and not os.path.isdir(out_dir):
            os.makedirs(out_dir, exist_ok=True)

    try:
        print(f'[debug] saving animation to {out_path} (format={out_fmt})', flush=True)
        if out_fmt == 'mp4':
            writer = mpl_animation.FFMpegWriter(fps=args.animate_fps)
            ani.save(out_path, writer=writer, dpi=args.dpi)
        else:
            writer = mpl_animation.PillowWriter(fps=args.animate_fps)
            ani.save(out_path, writer=writer)
        print(f'[*] animation saved: {out_path}')
    except KeyboardInterrupt:
        print('[WARN] animation generation interrupted by user; partial output may be incomplete.')
    except Exception as exc:
        if out_fmt == 'mp4':
            fallback = os.path.splitext(out_path)[0] + '.gif'
            print(f'[!] failed to save MP4 ({exc}); trying GIF: {fallback}')
            try:
                writer = mpl_animation.PillowWriter(fps=args.animate_fps)
                ani.save(fallback, writer=writer)
                print(f'[*] animation saved: {fallback}')
            except KeyboardInterrupt:
                print('[WARN] animation generation interrupted during GIF fallback; no animation saved.')
            except Exception as gif_exc:
                print(f'[ERROR] failed to save GIF fallback: {gif_exc}')
        else:
            print(f'[ERROR] failed to save animation: {exc}')
    finally:
        plt.close(fig)

def _format_axes(ax, xlabel: str, ylabel: str, title: str, grid: bool,
                 ticks_config: Tuple[int, int, bool]) -> None:
    n_lambda_axis, n_lambda_other, center_zero = ticks_config
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(grid)
    _apply_lambda_ticks(ax, int(ax.get_xlim()[1]), 'x', n_lambda_axis, center_zero)
    _apply_lambda_ticks(ax, int(ax.get_ylim()[1]), 'y', n_lambda_other, center_zero)
    ax.tick_params(axis='both', which='major', labelsize=11)


def build_figure(ex_data: Optional[np.ndarray], emag_data: Optional[np.ndarray],
                 bounds: Tuple[slice, slice, slice], spacing: Tuple[float, float, float],
                 args: argparse.Namespace) -> Tuple[plt.Figure, Dict[str, Optional[Any]]]:
    if not HAS_MPL:
        raise RuntimeError("matplotlib 未安装，无法生成图形。请先 pip install matplotlib。")

    range_x, range_y, range_z = bounds
    len_x = range_x.stop - range_x.start
    len_y = range_y.stop - range_y.start
    len_z = range_z.stop - range_z.start

    dx, dy, dz = spacing

    fig, axes = plt.subplots(2, 2, figsize=args.figsize, dpi=args.dpi)
    axes = axes.reshape(2, 2)

    cmap = args.cmap
    grid = args.grid

    # 共享色标设置
    ref_scale = None
    if args.e0 is not None:
        ref_scale = abs(args.e0)
    elif ex_data is not None:
        ref_scale = float(np.max(np.abs(ex_data)))
    elif emag_data is not None:
        ref_scale = float(np.max(emag_data))
    else:
        ref_scale = 1.0
    if ref_scale <= 0.0:
        ref_scale = 1.0
    ex_limit = args.max_scale * ref_scale
    mag_limit = args.max_scale * ref_scale

    def add_imshow(ax, data, title, symmetric, xlabel, ylabel,
                    extent: Tuple[float, float, float, float],
                    x_len: int, y_len: int,
                    x_spacing: float, y_spacing: float,
                    x_offset: float, y_offset: float,
                    center_zero_ticks: bool):
        if data is None:
            ax.axis('off')
            ax.text(0.5, 0.5, '数据缺失', ha='center', va='center', fontsize=12)
            return None
        norm = Normalize(vmin=-ex_limit, vmax=ex_limit) if symmetric else Normalize(vmin=0.0, vmax=mag_limit)
        im = ax.imshow(data, origin='lower', aspect='auto', cmap=cmap, norm=norm, extent=extent)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.formatter.set_powerlimits((0, 0))
        cb.ax.tick_params(labelsize=10)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontweight='bold')
        if args.lambda_ticks:
            _apply_lambda_ticks(ax, x_len, 'x', args.lambda_ticks, center_zero_ticks, x_spacing, x_offset)
            _apply_lambda_ticks(ax, y_len, 'y', args.lambda_ticks, center_zero_ticks, y_spacing, y_offset)
        if grid:
            ax.grid(True, linestyle='--', alpha=0.3)
        return im

    def axis_label(axis_key: str) -> str:
        if args.lambda_ticks:
            return f'{axis_key}/λ0'
        return f'索引 ({axis_key})'

    slices = _prepare_slices(ex_data, emag_data, range_x, range_y, range_z)
    images: Dict[str, Optional[Any]] = {}

    axis_mode = getattr(args, 'axis_mode', 'physical')
    if axis_mode == 'index':
        x_extent = (0.0, float(len_x))
        y_extent = (0.0, float(len_y))
        z_extent = (-0.5 * float(len_z), 0.5 * float(len_z))
        x_spacing = y_spacing = z_spacing = 1.0
        x_offset = 0.0
        y_offset = 0.0
        z_offset = -0.5 * float(len_z)
    else:
        x_extent = (0.0, len_x * dx)
        y_extent = (0.0, len_y * dy)
        z_extent = (-0.5 * len_z * dz, 0.5 * len_z * dz)
        x_spacing, y_spacing, z_spacing = dx, dy, dz
        x_offset = 0.0
        y_offset = 0.0
        z_offset = -0.5 * float(len_z)

    def plot_plane(plane: Optional[np.ndarray]) -> Optional[np.ndarray]:
        if plane is None:
            return None
        return plane.T

    images['ex_xz'] = add_imshow(
        axes[0, 0],
        plot_plane(slices['ex_xz']),
        '(a) Ex(x, 0, z, t)',
        symmetric=True,
        xlabel=axis_label('x'),
        ylabel=axis_label('z'),
        extent=(x_extent[0], x_extent[1], z_extent[0], z_extent[1]),
        x_len=len_x,
        y_len=len_z,
        x_spacing=x_spacing,
        y_spacing=z_spacing,
        x_offset=x_offset,
        y_offset=z_offset,
        center_zero_ticks=True,
    )
    images['emag_xz'] = add_imshow(
        axes[0, 1],
        plot_plane(slices['emag_xz']),
        '(b) |E|(x, 0, z)',
        symmetric=False,
        xlabel=axis_label('x'),
        ylabel=axis_label('z'),
        extent=(x_extent[0], x_extent[1], z_extent[0], z_extent[1]),
        x_len=len_x,
        y_len=len_z,
        x_spacing=x_spacing,
        y_spacing=z_spacing,
        x_offset=x_offset,
        y_offset=z_offset,
        center_zero_ticks=True,
    )
    images['ex_yz'] = add_imshow(
        axes[1, 0],
        plot_plane(slices['ex_yz']),
        '(c) Ex(0, y, z, t)',
        symmetric=True,
        xlabel=axis_label('y'),
        ylabel=axis_label('z'),
        extent=(y_extent[0], y_extent[1], z_extent[0], z_extent[1]),
        x_len=len_y,
        y_len=len_z,
        x_spacing=y_spacing,
        y_spacing=z_spacing,
        x_offset=y_offset,
        y_offset=z_offset,
        center_zero_ticks=True,
    )
    images['emag_yz'] = add_imshow(
        axes[1, 1],
        plot_plane(slices['emag_yz']),
        '(d) |E|(0, y, z)',
        symmetric=False,
        xlabel=axis_label('y'),
        ylabel=axis_label('z'),
        extent=(y_extent[0], y_extent[1], z_extent[0], z_extent[1]),
        x_len=len_y,
        y_len=len_z,
        x_spacing=y_spacing,
        y_spacing=z_spacing,
        x_offset=y_offset,
        y_offset=z_offset,
        center_zero_ticks=True,
    )

    fig.tight_layout()
    return fig, images


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: Optional[list] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='FDTD3D 切片可视化 (MATLAB 风格)')
    parser.add_argument('--ex-file', type=str, default=None, help='Ex VTI 文件路径 (瞬时值)')
    parser.add_argument('--emag-file', type=str, default=None, help='E 振幅 / RMS VTI 文件路径')
    parser.add_argument('--output', type=str, default=None, help='输出 PNG 路径 (默认显示)')
    parser.add_argument('--figsize', type=float, nargs=2, default=(8.0, 6.5), help='图形尺寸 (英寸)')
    parser.add_argument('--dpi', type=int, default=150, help='图像 DPI')
    parser.add_argument('--cmap', type=str, default='jet', help='色标方案 (默认 jet)')
    parser.add_argument('--grid', action='store_true', help='显示网格线')

    parser.add_argument('--max-scale', type=float, default=1.5, help='色标最大值系数 (默认 1.5×E0)')
    parser.add_argument('--e0', type=float, default=None, help='手动指定参考场幅 E0 (不指定则自动统计)')

    parser.add_argument('--n-pml-xn', type=int, default=20, help='负向 x PML 单元数')
    parser.add_argument('--n-pml-xp', type=int, default=20, help='正向 x PML 单元数')
    parser.add_argument('--n-pml-yn', type=int, default=20, help='负向 y PML 单元数')
    parser.add_argument('--n-pml-yp', type=int, default=20, help='正向 y PML 单元数')
    parser.add_argument('--n-pml-zn', type=int, default=20, help='负向 z PML 单元数')
    parser.add_argument('--n-pml-zp', type=int, default=20, help='正向 z PML 单元数')

    parser.add_argument('--lambda-ticks', type=int, default=0,
                        help='若 >0，则按 λ 刻度显示坐标轴 (值为每个 λ 所含网格数)')

    parser.add_argument('--animate', action='store_true', help='生成多时间步动画（需要 --*_series）')
    parser.add_argument('--animate-format', choices=['mp4', 'gif'], default='mp4', help='动画输出格式 (默认 mp4，若失败自动尝试 gif)')
    parser.add_argument('--animate-interval', type=int, default=20, help='动画帧间步数 (默认 20)')
    parser.add_argument('--animate-fps', type=int, default=15, help='动画帧率 (默认 15 fps)')
    parser.add_argument('--animate-output', type=str, default=None, help='动画输出文件路径 (默认与输出目录同名)')
    parser.add_argument('--ex-series', type=str, default=None, help='Ex VTI 序列 glob 模式 (例如 output/Ex_t*.vti)')
    parser.add_argument('--emag-series', type=str, default=None, help='E_mag VTI 序列 glob 模式')
    parser.add_argument('--animate-max-frames', type=int, default=200, help='动画最多渲染帧数，防止过大 (默认 200)')
    parser.add_argument('--axis-mode', choices=['physical', 'index'], default='physical', help='坐标轴刻度模式：physical 表示真实空间长度，index 表示原始网格索引')
    parser.add_argument('--show', action='store_true', help='生成后在窗口中展示 (默认仅保存/退出)')
    return parser.parse_args(argv)


def main(argv: Optional[list] = None) -> int:
    args = parse_args(argv)

    os.makedirs(DEFAULT_EXPORT_DIR, exist_ok=True)

    ex_series = _collect_series_map(args.ex_series)
    emag_series = _collect_series_map(args.emag_series)

    if args.ex_file is None and args.emag_file is None and not args.animate:
        print('[ERROR] require at least one of --ex-file or --emag-file.')
        return 1

    def load_optional(path: Optional[str]) -> Tuple[Optional[np.ndarray], Optional[Tuple[int, int, int]], Optional[Tuple[float, float, float]]]:
        if path is None:
            return None, None, None
        if not os.path.isfile(path):
            raise FileNotFoundError(f'未找到文件: {path}')
        data, dims, spacing = read_vti(path)
        return data, dims, spacing

    try:
        ex_data, ex_dims, ex_spacing = load_optional(args.ex_file)
        emag_data, emag_dims, emag_spacing = load_optional(args.emag_file)
    except Exception as exc:
        print(f'[ERROR] failed to load data: {exc}')
        return 1

    dims = None
    spacing: Optional[Tuple[float, float, float]] = None
    if ex_dims is not None:
        dims = ex_dims
    if ex_spacing is not None:
        spacing = ex_spacing
    if emag_dims is not None:
        if dims is None:
            dims = emag_dims
        elif dims != emag_dims:
            print('[WARN] Ex and E_mag dimensions differ; using Ex dimensions.')
    if emag_spacing is not None and spacing is None:
        spacing = emag_spacing
    if dims is None:
        for series in (ex_series, emag_series):
            if series:
                first_step = sorted(series.keys())[0]
                try:
                    _, series_dims, series_spacing = read_vti(series[first_step])
                except Exception as exc:
                    print(f'[ERROR] failed to read animation series: {exc}')
                    return 1
                dims = series_dims
                if spacing is None:
                    spacing = series_spacing
                break
    if dims is None:
        print('[ERROR] unable to determine grid size.')
        return 1
    if spacing is None:
        spacing = (1.0, 1.0, 1.0)

    Nx, Ny, Nz = dims
    range_x = _slice_bounds(Nx, args.n_pml_xn, args.n_pml_xp)
    range_y = _slice_bounds(Ny, args.n_pml_yn, args.n_pml_yp)
    range_z = _slice_bounds(Nz, args.n_pml_zn, args.n_pml_zp)

    fig: Optional[plt.Figure] = None
    if ex_data is not None or emag_data is not None:
        fig, _ = build_figure(ex_data, emag_data, (range_x, range_y, range_z), spacing, args)
        out_path = args.output or _default_static_path(args)
        out_dir = os.path.dirname(out_path)
        if out_dir and not os.path.isdir(out_dir):
            os.makedirs(out_dir, exist_ok=True)
        fig.savefig(out_path, dpi=args.dpi)
        print(f'[*] saved image: {out_path}')
        if args.show:
            fig.show()
            try:
                plt.show()
            except Exception:
                pass
        plt.close(fig)
    elif args.show:
        print('[WARN] --show requires --ex-file or --emag-file; ignored.')

    if args.animate:
        render_animation(ex_series, emag_series, (range_x, range_y, range_z), args)

    return 0


if __name__ == '__main__':
    sys.exit(main())
