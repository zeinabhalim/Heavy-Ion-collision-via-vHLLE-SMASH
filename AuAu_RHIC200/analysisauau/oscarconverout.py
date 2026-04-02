#!/usr/bin/env python3

import csv
import numpy as np
import matplotlib.pyplot as plt

# ==================== CONFIGURATION ====================
input_file  = "../sampler.out/cent0_5/particle_lists_0.oscar"
output_file = "particles_processed.csv"

species = {
    "eta":       221,
    "eta_prime": 331,
    "Lambda":    3122,
    "K0S":       310
}
# =======================================================

# Initialize containers
t_species = {key: [] for key in species}
pid_counts = {}
particle_data = []          # will hold [event_id, pid, px, py, pz, E, x_f, y_f, z_f, t_f]
# For global momentum distributions (all particles)
all_px = []
all_py = []
all_pz = []

current_event = 0
n_particles = 0

with open(input_file, 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        
        if line.startswith('#'):
            if 'event' in line:
                parts = line.split()
                if 'event' in parts:
                    current_event = int(parts[parts.index('event') + 1])
            continue
        
        parts = line.split()
        if len(parts) < 12:
            continue
        
        # OSCAR2013 columns (adjust if your format differs)
        t_f = float(parts[10])
        pid = int(parts[9])
        px = float(parts[6])
        py = float(parts[7])
        pz = float(parts[8])
        E  = float(parts[5])
        x  = float(parts[1])
        y  = float(parts[2])
        z  = float(parts[3])
        
        n_particles += 1
        pid_counts[pid] = pid_counts.get(pid, 0) + 1
        
        # Store for species‑specific histograms
        for name, code in species.items():
            if pid == code:
                t_species[name].append(t_f)
        
        # Store for global momentum distributions
        all_px.append(px)
        all_py.append(py)
        all_pz.append(pz)
        
        # Store full record for CSV
        particle_data.append([
            current_event, pid,
            px, py, pz, E,
            x, y, z, t_f
        ])
        
        if n_particles % 100000 == 0:
            print(f"  processed {n_particles} particles...")

print(f"Total particles processed: {n_particles}")

# ==================== Identify most common particle ====================
if pid_counts:
    most_common_pid = max(pid_counts, key=pid_counts.get)
    print(f"Most abundant PID: {most_common_pid} (count: {pid_counts[most_common_pid]})")

# ==================== Write CSV ====================
with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow([
        "event_id", "pid",
        "px", "py", "pz", "E",
        "x_f", "y_f", "z_f", "t_f"
    ])
    writer.writerows(particle_data)

print(f"CSV written to {output_file}")

# ==================== Plot 1: Freeze‑out time distributions ====================
plt.figure(figsize=(10,6))
colors = {"eta": "red", "eta_prime": "blue", "Lambda": "green", "K0S": "orange"}
for name, times in t_species.items():
    if times:
        plt.hist(times, bins=100, alpha=0.5, label=f"{name} (PID {species[name]})", color=colors.get(name, "gray"))
plt.xlabel("Freeze-out time t_f [fm/c]")
plt.ylabel("Number of particles")
plt.legend()
plt.grid(alpha=0.3)
plt.show()

# ==================== NEW: Plot 2 – px, py, pz distributions (all particles) ====================
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Common bins for comparison (optional)
bins = np.linspace(-10, 10, 100)   # adjust range based on your data

# px
axes[0].hist(all_px, bins=bins, color='blue', alpha=0.7, edgecolor='black')
axes[0].set_xlabel('p_x [GeV/c]')
axes[0].set_ylabel('Number of particles')
axes[0].set_title('Transverse momentum x-component')
axes[0].grid(alpha=0.3)

# py
axes[1].hist(all_py, bins=bins, color='green', alpha=0.7, edgecolor='black')
axes[1].set_xlabel('p_y [GeV/c]')
axes[1].set_ylabel('Number of particles')
axes[1].set_title('Transverse momentum y-component')
axes[1].grid(alpha=0.3)

# pz
axes[2].hist(all_pz, bins=bins, color='red', alpha=0.7, edgecolor='black')
axes[2].set_xlabel('p_z [GeV/c]')
axes[2].set_ylabel('Number of particles')
axes[2].set_title('Longitudinal momentum')
axes[2].grid(alpha=0.3)

plt.tight_layout()
plt.show()


print("\nMomentum statistics (all particles):")
print(f"p_x: mean = {np.mean(all_px):.4f}, RMS = {np.std(all_px):.4f}")
print(f"p_y: mean = {np.mean(all_py):.4f}, RMS = {np.std(all_py):.4f}")
print(f"p_z: mean = {np.mean(all_pz):.4f}, RMS = {np.std(all_pz):.4f}")
