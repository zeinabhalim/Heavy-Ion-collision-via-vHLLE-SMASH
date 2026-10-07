#!/usr/bin/env python3
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ============================================================
# 0. Publication-quality font & output settings
# ============================================================
# Serif font family matching journal body text (Times-like).
# Fallback chain because "Times New Roman" is often not installed on Linux.
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Nimbus Roman", "Liberation Serif",
                   "DejaVu Serif", "STIXGeneral"],
    "mathtext.fontset": "stix",       # math glyphs matching a Times-like body font
    "font.size": 16,
    "axes.titlesize": 18,
    "axes.labelsize": 16,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "legend.fontsize": 12,
    "figure.titlesize": 20,
    # Embed fonts as Type 42 (TrueType outlines) instead of matplotlib's
    # default Type 3 (bitmap), which many journals reject or misprint.
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

# ============================================================
# 1. Read OSCAR2013 particle file
#    Reads t, x, y, z (spatial freeze-out coordinates), px, py, pz, E, pid.
# ============================================================
def read_oscar(filename):
    t, x, y, z, px, py, pz, E, pid = [], [], [], [], [], [], [], [], []
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 12:
                continue
            t.append(float(parts[0]))
            x.append(float(parts[1]))
            y.append(float(parts[2]))
            z.append(float(parts[3]))
            px.append(float(parts[6]))
            py.append(float(parts[7]))
            pz.append(float(parts[8]))
            E.append(float(parts[5]))
            pid.append(int(parts[9]))
    return (np.array(t), np.array(x), np.array(y), np.array(z),
            np.array(px), np.array(py), np.array(pz),
            np.array(E), np.array(pid))

# ============================================================
# 2. Load particle data
# ============================================================
oscar_file = "../sampler.out/cent0_5/particle_lists_0.oscar"
t_part, x_part, y_part, z_part, px, py, pz, E, pid = read_oscar(oscar_file)
print(f"Loaded {len(t_part)} particles")

# ============================================================
# 2b. Reconstruct the REAL Milne coordinates (tau, eta_s) for
#     every particle from its actual (t, z) freeze-out point.
# ============================================================

valid = t_part > np.abs(z_part)
n_dropped = np.sum(~valid)
if n_dropped:
    print(f"Dropping {n_dropped} particles with |z| >= t "
          "(outside the future light cone; likely numerical edge cases)")

t_v = t_part[valid]
z_v = z_part[valid]

tau_part = np.sqrt(t_v**2 - z_v**2)
eta_part = 0.5 * np.log((t_v + z_v) / (t_v - z_v))

# ============================================================
# 3. Sample for plotting
# ============================================================
Nmax_hyper = 15000
idx_hyper = (np.random.choice(len(t_v), Nmax_hyper, replace=False)
             if len(t_v) > Nmax_hyper else np.arange(len(t_v)))

Nmax_config = 15000
Nmax_mom = 8000

idx_config = (np.random.choice(len(px), Nmax_config, replace=False)
              if len(px) > Nmax_config else np.arange(len(px)))
idx_mom = (np.random.choice(len(px), Nmax_mom, replace=False)
           if len(px) > Nmax_mom else np.arange(len(px)))
# ============================================================
# 4. Create horizontal row layout figure
# ============================================================

fig = plt.figure(figsize=(18, 6))

ax_width = 0.28
ax_height = 0.85
ax_spacing = 0.03

ax_hypersurface = fig.add_axes([0.03, 0.07, ax_width, ax_height], projection='3d')
ax_config       = fig.add_axes([0.03 + ax_width + ax_spacing, 0.07, ax_width, ax_height], projection='3d')
ax_momentum     = fig.add_axes([0.03 + 2*(ax_width + ax_spacing), 0.07, ax_width, ax_height], projection='3d')



# (a) REAL freeze-out hypersurface, reconstructed from the actual
#     (t, z) of every sampled particle -- no synthetic grid.
# ------------------------------
z_sel = z_v[idx_hyper]
tau_sel = tau_part[idx_hyper]
t_sel_hyper = t_v[idx_hyper]

surf = ax_hypersurface.scatter(
    z_sel, tau_sel, t_sel_hyper,
    c=t_sel_hyper,
    cmap="plasma", s=20, alpha=0.85
)
ax_hypersurface.set_xlabel("z [fm]")
ax_hypersurface.set_ylabel(r"$\tau$ [fm/$c$]")
ax_hypersurface.set_zlabel(r"$t$ [fm/$c$]")

ax_hypersurface.set_title("(a) Freeze-out Hypersurface ", fontweight="bold")
ax_hypersurface.view_init(elev=25, azim=-45)

#cbar_h = fig.colorbar(surf, ax=ax_hypersurface, fraction=0.05, pad=0.04)
#cbar_h.set_label(r"Freeze-out time $t_f$ [fm]", rotation=270, labelpad=18)

# ------------------------------
# (b) Freeze-out configuration space (true positions)
# ------------------------------
x_f = x_part[idx_config]
y_f = y_part[idx_config]
z_f = z_part[idx_config]
t_f = t_part[idx_config]

sc_config = ax_config.scatter(x_f, y_f, z_f, c=t_f, cmap="plasma", s=20, alpha=0.85)
#cbar_config = fig.colorbar(sc_config, ax=ax_config, fraction=0.05, pad=0.15)
#cbar_config.set_label(r"Freeze-out time $t_f$ [fm]", rotation=270, labelpad=18)
ax_config.set_xlabel(r"$x_f$ [fm]")
ax_config.set_ylabel(r"$y_f$ [fm]")
ax_config.set_zlabel(r"$z_f$ [fm]")
ax_config.set_title("(b) Freeze-out Configuration Space", fontweight="bold")
ax_config.view_init(elev=25, azim=-45)

# ------------------------------
# (c) 3D Momentum-space colored by freeze-out time
# ------------------------------

px_sel = px[idx_config]
py_sel = py[idx_config]
pz_sel = pz[idx_config]
t_sel = t_part[idx_config]


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

cbar_m = fig.colorbar(sc_mom, ax=ax_momentum, fraction=0.05, pad=0.15)
cbar_m.set_label(r"Freeze-out time $t_f$ [fm/$c$]", rotation=270, labelpad=20)
cbar_m.set_label(r"Lab-frame freeze-out time $t_f$ [fm/$c$]", rotation=270, labelpad=20)


# ============================================================
# 5. Save for the paper (vector PDF + high-res PNG preview)
# ============================================================
out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(out_dir, exist_ok=True)

fig.savefig(os.path.join(out_dir, "vhlle_AuAu2000.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(out_dir, "vhlle_AuAu2000.png"), dpi=300, bbox_inches="tight")
print(f"Saved figures to: {out_dir}")

plt.show()
