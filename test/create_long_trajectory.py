#!/usr/bin/env python
"""
Create a longer trajectory by combining multiple trajectory files.
"""
from ase.io import read, write
from pathlib import Path
import glob

# Read multiple trajectories and combine them
traj_dir = "data/ocp/traj/example"
traj_files = sorted(glob.glob(f"{traj_dir}/101_*.traj"))[:10]  # First 10 trajectories

all_frames = []
for traj_file in traj_files:
    frames = read(traj_file, index=':')
    if isinstance(frames, list):
        all_frames.extend(frames)
    else:
        all_frames.append(frames)

print(f"Combined {len(traj_files)} trajectories into {len(all_frames)} frames")

# Save combined trajectory
output_file = "test/combined_trajectory.traj"
write(output_file, all_frames)
print(f"Saved to: {output_file}")
