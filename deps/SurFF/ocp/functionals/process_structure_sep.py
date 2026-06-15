import os

import numpy as np
import pandas as pd
from pymatgen.core import Structure, Lattice
from pymatgen.core.surface import generate_all_slabs
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def fix_middle_atom(structure: Structure, fix_distance: float = 3.0) -> Poscar:
    selective_dynamics = []

    A = structure.lattice.matrix[0]
    B = structure.lattice.matrix[1]
    C = structure.lattice.matrix[2]

    N = np.cross(A, B)
    N = N / np.linalg.norm(N)

    # Get the z-coordinate of all sites
    # z_coordinates = [site.z for site in structure.sites]

    distance = [np.abs(np.dot(N, site)) for site in structure.cart_coords]

    # Find the maximum and minimum z-coordinate that corresponds to atomic positions
    distance_max = max(distance)
    distance_min = min(distance)

    # print(z_max, z_min)

    for d in distance:
        # Calculate the distance of the atom from the top and bottom of the slab
        distance_from_top = distance_max - d
        distance_from_bottom = d - distance_min

        # If the atom is more than "distance" away from both the top and bottom, fix the atom
        if distance_from_top > fix_distance and distance_from_bottom > fix_distance:
            selective_dynamics.append([False, False, False])
        else:
            selective_dynamics.append([True, True, True])

    return Poscar(structure, selective_dynamics=selective_dynamics)


def generate_unique_slabs(bulk):
    """
    Generates unique slabs for a given crystal structure.

    Args:
        bulk (Structure): Pymatgen Structure object.

    Returns:
        list: List of unique slabs.
    """

    min_slab_size = 10.0

    try:
        slabs = generate_all_slabs(bulk, max_index=2, min_slab_size=min_slab_size,
                                   min_vacuum_size=15.0, primitive=False, center_slab=True)
    except:
        print(f'\nFailed to generate slabs for {bulk.formula} at thickness {min_slab_size}')
        return []

    return slabs


def process_structure_sep(crystal_id, slab_dir, bulk_dir, source_dir):
    slab_data = []

    POSCAR = Poscar.from_file(os.path.join(source_dir, f"{crystal_id}"))
    structure = POSCAR.structure
    # get conventional cell
    structure = SpacegroupAnalyzer(structure, symprec=0.01).get_conventional_standard_structure()
    fomula = structure.formula

    unique_slabs = generate_unique_slabs(structure)

    for j, slab in enumerate(unique_slabs):
        miller_index = slab.miller_index
        index = ''.join(map(str, miller_index))
        shift = slab.shift
        num_atom = len(slab.sites)
        # slab_thickness = slab.cart_coords[:, 2].max() - slab.cart_coords[:, 2].min()
        # print(f"slab_thickness {slab_thickness}")

        slab_filename = os.path.join(slab_dir, f"{crystal_id}_{j}")
        # slab.to(filename=slab_filename, fmt='poscar')
        with open(slab_filename, 'w') as f:
            p = fix_middle_atom(slab)
            f.write(str(p))

        slab_info = [f"{crystal_id}_{j}", fomula, index, shift, num_atom]
        slab_data.append(slab_info)

        bulk_filename = os.path.join(bulk_dir, f"{crystal_id}_{j}")

        ouc = slab.oriented_unit_cell

        ouc_lattice = Lattice.from_parameters(
            ouc.lattice.a,
            ouc.lattice.b,
            ouc.lattice.c,
            ouc.lattice.alpha,
            ouc.lattice.beta,
            ouc.lattice.gamma, )

        ouc = Structure(ouc_lattice, ouc.species, ouc.frac_coords, coords_are_cartesian=False, )

        ouc.to(filename=bulk_filename, fmt='poscar')

    return slab_data
