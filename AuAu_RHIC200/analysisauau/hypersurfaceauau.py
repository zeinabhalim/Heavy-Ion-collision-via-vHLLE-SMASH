#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ============================================================
# 1. Read OSCAR2013 particle file
# ============================================================
def read_oscar(filename):
    t, px, py, pz, E, pid = [], [], [], [], [], []
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 12:
                continue
            t.append(float(parts[0]))
            px.append(float(parts[6]))
            py.append(float(parts[7]))
            pz.append(float(parts[8]))
            E.append(float(parts[5]))
            pid.append(int(parts[9]))
    return np.array(t), np.array(px), np.array(py), np.array(pz), np.array(E), np.array(pid)

# ============================================================
# 2. Load particle data
# ============================================================
oscar_file = "../sampler.out/cent0_5/particle_lists_0.oscar"
t_part, px, py, pz, E, pid = read_oscar(oscar_file)
print(f"Loaded {len(t_part)} particles")

# ============================================================
# 3. Sample for plotting
# ============================================================
Nmax_config = 12000
Nmax_mom = 8000

idx_config = np.random.choice(len(px), Nmax_config, replace=False) if len(px) > Nmax_config else slice(None)
idx_mom = np.random.choice(len(px), Nmax_mom, replace=False) if len(px) > Nmax_mom else np.arange(len(px))

# ============================================================
# 4. Create horizontal row layout figure
# ============================================================
fig = plt.figure(figsize=(36, 10))  # slightly wider figure

# Axes positions: [left, bottom, width, height]
ax_width = 0.28
ax_height = 0.85
ax_spacing = 0.03

ax_hypersurface = fig.add_axes([0.03, 0.07, ax_width, ax_height], projection='3d')
ax_config = fig.add_axes([0.03 + ax_width + ax_spacing, 0.07, ax_width, ax_height], projection='3d')
ax_momentum = fig.add_axes([0.03 + 2*(ax_width + ax_spacing), 0.07, ax_width, ax_height], projection='3d')

# ------------------------------
# (a) Centered Freeze-out hypersurface
# ------------------------------
eta = np.linspace(-3, 3, 60)
tau = np.linspace(1, 12, 60)
ETA, TAU = np.meshgrid(eta, tau)

T = TAU * np.cosh(ETA)
T_centered = T - np.mean(T)
Z = TAU * np.sinh(ETA)
Z_centered = Z - np.mean(Z)

surf = ax_hypersurface.plot_surface(Z_centered, TAU, T_centered, cmap="viridis",
                                    alpha=0.85, edgecolor='k', linewidth=0.2)
ax_hypersurface.set_xlabel("z [fm]")
ax_hypersurface.set_ylabel("τ [fm]")
ax_hypersurface.set_zlabel("t [fm]")
ax_hypersurface.set_title("(a) Freeze-out Hypersurface", fontweight="bold")
ax_hypersurface.view_init(elev=30, azim=-60)

cbar_h = fig.colorbar(surf, ax=ax_hypersurface, fraction=0.05, pad=0.04)
cbar_h.set_label(r"Coordinate $\langle t \rangle$ [fm]", rotation=270, labelpad=18)

# ------------------------------
# (b) Freeze-out configuration space
# ------------------------------
x_f = px[idx_config]
y_f = py[idx_config]
z_f = pz[idx_config]
t_f = t_part[idx_config]

sc_config = ax_config.scatter(
    x_f, y_f, z_f,
    c=t_f,
    cmap="plasma",
    s=20,
    alpha=0.85
)
ax_config.set_xlabel(r"$x_f$ [fm]")
ax_config.set_ylabel(r"$y_f$ [fm]")
ax_config.set_zlabel(r"$z_f$ [fm]")
ax_config.set_title("(b) Freeze-out Configuration Space", fontweight="bold")
ax_config.view_init(elev=25, azim=-45)

# ------------------------------
# (c) 3D Momentum-space colored by freeze-out time
# ------------------------------
px_sel = px[idx_mom] - np.mean(px[idx_mom])
py_sel = py[idx_mom] - np.mean(py[idx_mom])
pz_sel = pz[idx_mom] - np.mean(pz[idx_mom])
t_sel = t_part[idx_mom]

sc_mom = ax_momentum.scatter(
    px_sel, py_sel, pz_sel,
    c=t_sel,
    cmap="plasma",
    s=20,
    alpha=0.85
)
ax_momentum.set_xlabel(r"$p_x$ [GeV]")
ax_momentum.set_ylabel(r"$p_y$ [GeV]")
ax_momentum.set_zlabel(r"$p_z$ [GeV]")
ax_momentum.set_title("(c) Momentum Coordinates Distribution", fontweight="bold")
ax_momentum.view_init(elev=25, azim=-45)

fig.suptitle(
    "Hybrid Model: vHLLE Hydrodynamics + SMASH Hadronic Transport\n"
    r"Au+Au at $\sqrt{s_{NN}}=200$ GeV, Centrality 0–5%",
    fontsize=22,
    fontweight="bold",
    alpha=0.9
)


cbar_m = fig.colorbar(sc_mom, ax=ax_momentum, fraction=0.05, pad=0.15)
cbar_m.set_label(r"Freeze-out time $t_f$ [fm]", rotation=270, labelpad=18)

# ------------------------------
# Finalize
# ------------------------------
#fig.suptitle(
#    "Hybrid Model: vHLLE Hydrodynamics + SMASH Hadronic Transport",
#    fontsize=28, fontweight="bold"
#)
plt.show()


