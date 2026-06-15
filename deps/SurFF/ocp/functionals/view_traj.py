from ase.io import read, Trajectory, BundleTrajectory

# Replace 'your_file.traj' with the path to your .traj file
trajectory = Trajectory(r'D:\Pyprojects\ocp\ocp\traj\test\102_11.traj')

# 'trajectory' is now a list of ASE Atoms objects, one for each frame in the trajectory.
# You can iterate over it or access individual frames.

# For example, to access the first frame:
first_frame = trajectory[0]
print(first_frame)

# To access the last frame:
last_frame = trajectory[-1]
print(last_frame)
