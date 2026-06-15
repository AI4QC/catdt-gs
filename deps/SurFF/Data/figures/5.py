import ase
from ase.io import read
from pymatgen.core import Structure
from pymatgen.core.surface import SlabGenerator

CIF_path = r"Data/figures/CuPd.cif"

# read CIF file in pymatgen
structure = Structure.from_file(CIF_path)

# create a slab generator for the 001 surface
slab_gen = SlabGenerator(initial_structure=structure,
                         miller_index=(1, 1, 0),
                         min_slab_size=10.0,  # thickness of the slab in Ångströms
                         min_vacuum_size=10.0,  # thickness of the vacuum layer in Ångströms
                         center_slab=True)  # center the slab in the vacuum layer

# generate the slab
slab = slab_gen.get_slab()

# make it 4*4*1
slab.make_supercell([4, 4, 1])

# Convert to ASE Atoms object for visualization
from pymatgen.io.ase import AseAtomsAdaptor

ase_atoms = AseAtomsAdaptor().get_atoms(slab)

# Visualize the slab
from ase.visualize import view

view(ase_atoms)
