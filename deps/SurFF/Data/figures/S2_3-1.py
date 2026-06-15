# visualize POSCAR of crystal slab and oriented unit cell
from ase.io import read

slab_pth = r"Data/figures/45_4_slab"
ouc_pth = r"Data/figures/45_4_bulk"

slab = read(slab_pth, format='vasp')
ouc = read(ouc_pth, format='vasp')

# plot
from ase.visualize import view
view(ouc)
view(slab)



