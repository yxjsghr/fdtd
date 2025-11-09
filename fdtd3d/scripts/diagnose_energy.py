#!/usr/bin/env python3
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))
from verify_vti import read_vti
import numpy as np

# 检查能量增长模式
print("检查能量增长模式...")
print("="*60)

for t in [0, 20, 40, 60, 80, 100, 120, 140, 160]:
    fname = f"build/output_3d/E_mag_t{t:05d}.vti"
    if not os.path.exists(fname):
        continue
    try:
        data, dims, origin, spacing = read_vti(fname)
        finite_mask = np.isfinite(data)
        if np.all(finite_mask):
            max_val = np.max(data)
            total_energy = np.sum(data**2)
            # 检查能量分布
            energy_density = data**2
            center_energy = np.sum(energy_density[90:110, 90:110, 90:110])
            edge_energy = (np.sum(energy_density[:20, :, :]) + 
                          np.sum(energy_density[-20:, :, :]) +
                          np.sum(energy_density[:, :20, :]) +
                          np.sum(energy_density[:, -20:, :]) +
                          np.sum(energy_density[:, :, :20]) +
                          np.sum(energy_density[:, :, -20:]))
            print(f"t={t:5d}: max={max_val:.6e}, total={total_energy:.6e}, "
                  f"center={center_energy:.6e}, edge={edge_energy:.6e}, "
                  f"ratio={edge_energy/total_energy if total_energy>0 else 0:.4f}")
        else:
            n_inf = np.sum(~finite_mask)
            print(f"t={t:5d}: [ERROR] {n_inf} non-finite values")
    except Exception as e:
        print(f"t={t:5d}: [ERROR] {e}")

