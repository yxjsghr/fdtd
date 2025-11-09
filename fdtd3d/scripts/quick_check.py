#!/usr/bin/env python3
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))
from verify_vti import read_vti
import numpy as np

# 检查前几个时间步
for t in [0, 20, 40, 60, 80, 100]:
    fname = f"build/output_3d/E_mag_t{t:05d}.vti"
    if os.path.exists(fname):
        data, dims, origin, spacing = read_vti(fname)
        finite = np.isfinite(data)
        n_finite = np.sum(finite)
        max_val = np.max(data[finite]) if n_finite > 0 else 0
        total_energy = np.sum(data[finite]**2) if n_finite > 0 else 0
        print(f"t={t:5d}: max={max_val:.6e}, energy={total_energy:.6e}, finite={n_finite}/{data.size}")

