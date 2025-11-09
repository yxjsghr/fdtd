#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FDTD 2D 可视化（漂亮稳健版）
- 读取 OutputManager 的 CSV: "# x,z,val"，id = ix*Nz + kz
- 横轴 x、纵轴 z；默认 transpose=True（arr=(Nx,Nz)->(Nz,Nx)）
- Ey: 对称色标 [-A,+A]；Eymag: 默认 [0,A]（更不易“全蓝”），可 --sym-mag
- 色标优先：--jy0*--max-scale；否则全局分位数 (--clip-ptile)；若仍太小，逐帧自动
- 输出：动画（mp4/gif）+ 最后一帧 PNG
"""
import argparse, os, re, sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def read_dims(outdir):
    p = os.path.join(outdir, "dims.txt")
    if not os.path.isfile(p): raise FileNotFoundError(f"未找到 {p}")
    d={}
    for ln in open(p):
        s=ln.split()
        if len(s)==2:
            k,v=s
            d[k]=int(v) if k in ("Nx","Nz") else float(v)
    for k in ("Nx","Nz","dx","dz"):
        if k not in d: raise ValueError(f"dims.txt 缺少 {k}")
    return d["Nx"], d["Nz"], d["dx"], d["dz"]

def find_files(outdir, prefix):
    pat=re.compile(rf"^{re.escape(prefix)}_t(\d+)\.csv$")
    lst=[]
    for fn in os.listdir(outdir):
        m=pat.match(fn)
        if m: lst.append((int(m.group(1)), os.path.join(outdir,fn)))
    lst.sort(key=lambda x:x[0])
    return [p for _,p in lst]

def load_frame(path, Nx, Nz, transpose=True):
    data=np.loadtxt(path, delimiter=",", comments="#")
    if data.ndim==1: data=data[None,:]
    if data.shape[1]<3: raise ValueError(f"{os.path.basename(path)} 列不足3")
    vals=data[:,2]
    if vals.size!=Nx*Nz: raise ValueError(f"{os.path.basename(path)} 点数 {vals.size} != {Nx*Nz}")
    arr=vals.reshape(Nx, Nz)     # (Nx,Nz): x 行、z 列
    if transpose: arr=arr.T      # -> (Nz,Nx): 行= z，列= x
    return arr

def robust_percentile(files, Nx, Nz, transpose, ptile=99.0, step=10):
    vmax=0.0
    for i in range(0, len(files), max(1,int(step))):
        a=load_frame(files[i], Nx, Nz, transpose)
        aa=np.abs(a[np.isfinite(a)])
        if aa.size:
            vmax=max(vmax, float(np.percentile(aa, ptile)))
    return vmax

def choose_clim(field, files, Nx, Nz, transpose, args):
    # 1) 明确给了 jy0
    if args.jy0 is not None:
        A=float(args.max_scale)*float(args.jy0)
        if field=="ey":   return (-A, +A)
        return ((-A, +A) if args.sym_mag else (0.0, A))

    # 2) 手工 vmin/vmax
    if args.vmin is not None or args.vmax is not None:
        vmin = args.vmin if args.vmin is not None else (0.0 if (field=="eymag" and not args.sym_mag) else None)
        if vmin is None:  # ey 或 eymag 对称
            # 估一个全局上限
            A=robust_percentile(files, Nx, Nz, args.transpose, args.clip_ptile, args.scan_step)
            vmin=-A
        vmax = args.vmax if args.vmax is not None else robust_percentile(files, Nx, Nz, args.transpose, args.clip_ptile, args.scan_step)
        return (vmin, vmax)

    # 3) 自动：全局分位数
    A=robust_percentile(files, Nx, Nz, args.transpose, args.clip_ptile, args.scan_step)
    if A<=1e-15:  # 极小，交给逐帧
        return None
    if field=="ey":   return (-A,+A)
    return ((-A,+A) if args.sym_mag else (0.0, A))

def set_lambda_ticks(ax, Nx, Nz, dx, dz, args):
    # 以索引单位标 tick，但标签按 λ0
    if args.lambda0 is not None:
        nlam_x=max(1,int(round(args.lambda0/dx)))
        nlam_z=max(1,int(round(args.lambda0/dz)))
    else:
        nlam_x=nlam_z=int(args.n_lambda)

    xs=[0]+[k for k in range(nlam_x, Nx+1, nlam_x)]
    zs=[0]+[k for k in range(nlam_z, Nz+1, nlam_z)]
    ax.set_xticks(xs); ax.set_yticks(zs)
    ax.set_xticklabels([str(k//nlam_x) if k>0 else "0" for k in xs])
    ax.set_yticklabels([str(k//nlam_z) if k>0 else "0" for k in zs])
    ax.set_xlabel(r"$x/\lambda_0$"); ax.set_ylabel(r"$z/\lambda_0$")

def render_field(outdir, field, args):
    Nx,Nz,dx,dz = read_dims(outdir)
    files = find_files(outdir, "ey" if field=="ey" else "eymag")
    if not files:
        print(f"[!] {outdir} 下没有 {field}_t*.csv"); return

    A0=load_frame(files[0], Nx, Nz, args.transpose)
    extent=[0, Nx*dx, 0, Nz*dz] if args.extent_real else [0, Nx, 0, Nz]

    fig, ax = plt.subplots(figsize=(9,6), dpi=args.dpi)
    fig.patch.set_facecolor("white")
    im = ax.imshow(A0, origin="lower", extent=extent, cmap=args.cmap, interpolation="nearest", aspect="equal")
    if not args.grid: ax.grid(False)

    set_lambda_ticks(ax, Nx, Nz, dx, dz, args)  # λ 刻度（更像你的 MATLAB）
    title = "Ey(x,z,t)" if field=="ey" else "E_y magnitude"
    ax.set_title(title)

    # 颜色条 & 色标
    clim = choose_clim(field, files, Nx, Nz, args.transpose, args)
    if clim is not None:
        im.set_clim(*clim)
        print(f"[clim-{field}] vmin={clim[0]:.3g}, vmax={clim[1]:.3g}  (fixed)")
    else:
        print(f"[clim-{field}] 全局幅度过小，启用逐帧自适应（避免空白）")

    if args.colorbar:
        cb=fig.colorbar(im, ax=ax)
        cb.set_label(field.upper())

    plt.tight_layout()

    def update(i):
        A=load_frame(files[i], Nx, Nz, args.transpose)
        im.set_data(A)
        if clim is None:  # 逐帧自适应
            amax=float(np.nanmax(np.abs(A))) if A.size else 1.0
            if field=="ey":
                im.set_clim(-amax, +amax)
            else:
                im.set_clim(((-amax,+amax) if args.sym_mag else (0.0, amax)))
        ax.set_title(f"{title} (frame {i})")
        return (im,)

    ani = animation.FuncAnimation(fig, update, frames=len(files),
                                  interval=args.interval_ms, blit=False, repeat=False)

    if args.save:
        outmp4=os.path.join(outdir, f"{field}.mp4")
        try:
            ani.save(outmp4, writer="ffmpeg", fps=args.fps, dpi=args.dpi)
            print(f"[*] 已保存: {outmp4}")
        except Exception as e:
            outgif=os.path.join(outdir, f"{field}.gif")
            print(f"[!] MP4 失败: {e}\n[*] 回退 GIF: {outgif}")
            ani.save(outgif, writer="pillow", fps=args.fps)

    if args.snapshot_last:
        last=len(files)-1
        update(last)
        png=os.path.join(outdir, f"{field}_frame{last:05d}.png")
        fig.savefig(png, dpi=args.dpi)
        print(f"[*] 已保存快照: {png}")

    if not args.save and not args.snapshot_last:
        try: plt.show()
        except Exception: pass
    plt.close(fig)

def main():
    ap=argparse.ArgumentParser(description="FDTD 2D 可视化（漂亮稳健版）")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--field", default="both", choices=["ey","eymag","both"])
    ap.add_argument("--save", action="store_true")
    ap.add_argument("--snapshot-last", action="store_true")
    ap.add_argument("--cmap", default="jet")
    ap.add_argument("--grid", action="store_true", help="显示网格（默认关闭）")
    ap.add_argument("--colorbar", action="store_true", default=True)
    ap.add_argument("--no-colorbar", dest="colorbar", action="store_false")
    ap.add_argument("--transpose", action="store_true", default=True, help="默认转置到 (z,x)。若仍感觉颠倒，可加 --no-transpose")
    ap.add_argument("--no-transpose", dest="transpose", action="store_false")

    # 色标
    ap.add_argument("--jy0", type=float, default=None)
    ap.add_argument("--max-scale", type=float, default=1.25)
    ap.add_argument("--clip-ptile", type=float, default=99.0)
    ap.add_argument("--scan-step", type=int, default=10)
    ap.add_argument("--vmin", type=float, default=None)
    ap.add_argument("--vmax", type=float, default=None)
    ap.add_argument("--sym-mag", action="store_true", default=False, help="Eymag 用对称色标（默认关闭，使用[0,max]）")

    # 动画参数
    ap.add_argument("--interval-ms", type=int, default=30)
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--dpi", type=int, default=150)

    # 轴刻度
    ap.add_argument("--n-lambda", type=int, default=40)
    ap.add_argument("--lambda0", type=float, default=None)
    ap.add_argument("--extent-real", action="store_true", help="extent 用米（否则用索引）")

    args=ap.parse_args()
    if "MPLBACKEND" not in os.environ: matplotlib.use("Agg")
    fields=["ey","eymag"] if args.field=="both" else [args.field]
    for f in fields:
        print(f"[*] 渲染 {f} …")
        render_field(args.outdir, f, args)
    print("[✓] 完成。")

if __name__=="__main__":
    main()
