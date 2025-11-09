#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# FDTD3D 可视化批处理脚本
# 直接运行本脚本即可调用 visualize_slices.py 生成动画或静态图。
# 如需调整参数，请修改下方变量；每个变量旁的中文注释说明了调整方法。
# -----------------------------------------------------------------------------

set -euo pipefail

# 取得脚本所在目录（fdtd3d 根目录），确保无论从哪个路径调用都能找到 Python 脚本。
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

# ======================== 基本路径设置 ========================
# EX_SERIES: Ex 序列 VTI 文件的通配符路径；若只生成静态图，可改为单个文件路径 (配合 EX_FILE 使用)。
EX_SERIES='build/output_3d/Ex_t*.vti'
# EMAG_SERIES: E_mag 或 E_mag_rms 序列 VTI 文件；若无此数据，可置空字符串。
EMAG_SERIES='build/output_3d/E_mag_rms_t*.vti'
# EX_FILE / EMAG_FILE: 单帧静态图使用的文件；动画模式下通常留空即可。
EX_FILE=''      # 例如 'build/output_3d/Ex_t00540.vti'
EMAG_FILE=''    # 例如 'build/output_3d/E_mag_rms_t00540.vti'

# ======================== 动画输出控制 ========================
# ENABLE_ANIMATION: 是否生成动画；设为 'true' 则启用动画流程，'false' 只生成静态图。
ENABLE_ANIMATION='true'
# ANIMATE_FORMAT: 动画格式，可选 'mp4' 或 'gif'。
ANIMATE_FORMAT='mp4'
# ANIMATE_INTERVAL: 帧间步数，每隔多少个时间步取一帧。
ANIMATE_INTERVAL=20
# ANIMATE_MAX_FRAMES: 最多渲染多少帧，防止动画过大。
ANIMATE_MAX_FRAMES=200
# ANIMATE_FPS: 动画帧率 (frames per second)。
ANIMATE_FPS=10
# ANIMATE_OUTPUT: 动画输出路径；留空则自动保存到 build/visualizations/ 下。
ANIMATE_OUTPUT=''

# ======================== 静态图输出控制 ========================
# OUTPUT_IMAGE: 静态 PNG 输出路径；留空则保存在 build/visualizations/ 目录。
OUTPUT_IMAGE=''
# SHOW_FIG: 是否在生成后弹出窗口显示 (True/False)；服务器环境建议保持 False。
SHOW_FIG='False'

# ======================== 绘图外观设置 ========================
# AXIS_MODE: 坐标轴模式。'physical' 表示显示真实空间长度，'index' 表示网格索引。
AXIS_MODE='physical'
# LAMBDA_TICKS: 若 >0，则按给定网格数作为 1λ 刻度；设置为 0 则使用默认坐标刻度。
LAMBDA_TICKS=20
# FIGSIZE: 图像尺寸 (单位英寸)；例 "8 6.5"。
FIGSIZE='8 6.5'
# DPI: 输出分辨率；数值越大图像越清晰但文件更大。
DPI=150
# CMAP: 颜色映射表，例如 'jet'、'viridis' 等。
CMAP='jet'
# GRID: 是否显示网格线；填 'true' 显示，'false' 隐藏。
GRID='false'
# MAX_SCALE: 色标最大值系数，默认 1.5 表示以 (1.5×参考幅值) 为最大值。
MAX_SCALE=1.5
# E0: 手动指定参考场幅值；留空表示自动统计最大值。
E0=''

# ======================== PML (吸收边界) 裁剪 ========================
# 下列参数控制可视化时剔除前后 PML 单元数，可根据仿真配置调整。
N_PML_XN=20
N_PML_XP=20
N_PML_YN=20
N_PML_YP=20
N_PML_ZN=20
N_PML_ZP=20

# -----------------------------------------------------------------------------
# 组装命令
# -----------------------------------------------------------------------------
PY_SCRIPT="${SCRIPT_DIR}/scripts/visualize_slices.py"
ARGS=("${PY_SCRIPT}")

# 静态图文件
if [[ -n "${EX_FILE}" ]]; then
  ARGS+=("--ex-file" "${EX_FILE}")
fi
if [[ -n "${EMAG_FILE}" ]]; then
  ARGS+=("--emag-file" "${EMAG_FILE}")
fi
if [[ -n "${OUTPUT_IMAGE}" ]]; then
  ARGS+=("--output" "${OUTPUT_IMAGE}")
fi
if [[ "${SHOW_FIG}" == 'True' ]]; then
  ARGS+=("--show")
fi

# 动画相关
if [[ "${ENABLE_ANIMATION}" == 'true' ]]; then
  ARGS+=("--animate")
  ARGS+=("--animate-format" "${ANIMATE_FORMAT}")
  ARGS+=("--animate-interval" "${ANIMATE_INTERVAL}")
  ARGS+=("--animate-fps" "${ANIMATE_FPS}")
  ARGS+=("--animate-max-frames" "${ANIMATE_MAX_FRAMES}")
  if [[ -n "${ANIMATE_OUTPUT}" ]]; then
    ARGS+=("--animate-output" "${ANIMATE_OUTPUT}")
  fi
fi

if [[ -n "${EX_SERIES}" ]]; then
  ARGS+=("--ex-series" "${EX_SERIES}")
fi
if [[ -n "${EMAG_SERIES}" ]]; then
  ARGS+=("--emag-series" "${EMAG_SERIES}")
fi

# 绘图外观
ARGS+=("--axis-mode" "${AXIS_MODE}")
ARGS+=("--figsize" ${FIGSIZE})
ARGS+=("--dpi" "${DPI}")
ARGS+=("--cmap" "${CMAP}")
if [[ "${GRID}" == 'true' ]]; then
  ARGS+=("--grid")
fi
ARGS+=("--max-scale" "${MAX_SCALE}")
if [[ -n "${E0}" ]]; then
  ARGS+=("--e0" "${E0}")
fi
if [[ ${LAMBDA_TICKS} -gt 0 ]]; then
  ARGS+=("--lambda-ticks" "${LAMBDA_TICKS}")
fi

# PML 裁剪
ARGS+=("--n-pml-xn" "${N_PML_XN}")
ARGS+=("--n-pml-xp" "${N_PML_XP}")
ARGS+=("--n-pml-yn" "${N_PML_YN}")
ARGS+=("--n-pml-yp" "${N_PML_YP}")
ARGS+=("--n-pml-zn" "${N_PML_ZN}")
ARGS+=("--n-pml-zp" "${N_PML_ZP}")

# -----------------------------------------------------------------------------
# 执行命令
# -----------------------------------------------------------------------------
echo "[info] 调用命令:" "python3 ${ARGS[*]}"
python3 "${ARGS[@]}"
