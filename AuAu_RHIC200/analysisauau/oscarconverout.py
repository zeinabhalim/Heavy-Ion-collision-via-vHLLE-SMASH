#!/usr/bin/env python3

import csv
import os

input_file  = "../sampler.out/cent0_5/particle_lists_0.oscar"
output_file = "final_stateauau.csv"

print(f"Reading: {input_file}")
print(f"Writing: {output_file}")

rows = []
current_event = 0
n_particles = 0

with open(input_file, "r") as f:
    for line in f:
        line = line.strip()

        # Skip comments
        if not line or line.startswith("#"):
            # Extract event number if present
            if "event" in line:
                parts = line.split()
                if "event" in parts:
                    current_event = int(parts[parts.index("event") + 1])
            continue

        parts = line.split()
        if len(parts) < 12:
            continue

        # OSCAR2013 format (your file!)
        t  = float(parts[0])
        x  = float(parts[1])
        y  = float(parts[2])
        z  = float(parts[3])
        m  = float(parts[4])
        E  = float(parts[5])
        px = float(parts[6])
        py = float(parts[7])
        pz = float(parts[8])
        pid = int(parts[9])

        rows.append([
            current_event,
            pid,
            px, py, pz, E,
            x, y, z, t
        ])

        n_particles += 1
        if n_particles % 100000 == 0:
            print(f"  processed {n_particles} particles...")

print(f"Total particles: {n_particles}")

# Write CSV
with open(output_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "event_id", "pid",
        "px", "py", "pz", "E",
        "x_f", "y_f", "z_f", "t_f"
    ])
    writer.writerows(rows)

print("Done.")


