import lmdb

lmdb_pth = r"D:\OneDrive\OneDrive - National University of Singapore\Pyprojects\Intermetallic_Crystal_Nanoparticle_Model\ocp\data\example\surface_relaxation.lmdb"

env = lmdb.open(lmdb_pth, map_async=True, subdir=False, readonly=False)
print(env.stat())