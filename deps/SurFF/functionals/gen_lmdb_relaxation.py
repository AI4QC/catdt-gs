import os
import pdb
import pickle

import lmdb
import numpy as np
import pandas as pd
import torch
from pymatgen.io.vasp import Poscar
from torch_geometric.data import Data
from tqdm import tqdm


def get_POSCAR_data(pth_enegy_csv_lst, large_dif_id=None):
    energy_csv_lst = []
    for pth_energy_csv in pth_enegy_csv_lst:
        energy_csv_lst.append(pd.read_csv(pth_energy_csv))

    energy_csv = pd.concat(energy_csv_lst, axis=0, ignore_index=True)
    # energy_csv = energy_csv.drop(energy_csv[energy_csv["converge"] == False].index)

    # Drop rows where 'chem' column contains "Fe" or "Cr"
    # energy_csv = energy_csv[~energy_csv['formula'].str.contains('Cr|Mn|Fe|Ni|Co|Cu')]

    # energy_csv = energy_csv[~energy_csv['slab_id'].isin(large_dif_id)]
    # energy_csv.to_csv(r'Data/Workspace/test_all_clean.csv', index=False)

    # energy_csv = energy_csv.sample(frac=1, random_state=11).reset_index(drop=True)
    # energy_csv = energy_csv[int(len(energy_csv)*0.6):]

    # slab_index = pd.read_csv(pth_slab_index)

    slab_id = energy_csv["slab_id"].values
    slab_cif_pth = [os.path.join(POSCAR_pth, f"{n}") for n in slab_id]

    energy_csv["POSCAR_pth"] = slab_cif_pth

    return energy_csv

def time_it(func):
    def wrapper(*args, **kwargs):
        import time
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        print(f"Time elapsed: {end - start}")
        return result

    return wrapper

# @time_it
def process_item(i: int, dataset):
    """
    Process a dataset item given its index and return a PyTorch Geometric Data object
    or None if there was an error.

    Args:
        i (int): Index of the item in the dataset.
        dataset (Dataset): The dataset containing the item.

    Returns:
        Union[Data, None]: A PyTorch Geometric Data object containing the item's
                           information or None if there was an error.
    """

    try:
        # Retrieve features, neighbor features, neighbor feature indices, and target
        data = dataset[i]
        POSCAR_pth = data["POSCAR_pth"]

        sid = data["slab_id"]
        frame_number = 0

        crystal = Poscar.from_file(POSCAR_pth).structure
        natoms = len(crystal)

        pos = crystal.cart_coords
        pos = torch.Tensor(pos.tolist())

        try:
            tags = crystal.site_properties['selective_dynamics']
            tags = list(map(lambda x: 1 if x[0] else 0, tags))
            tags = torch.LongTensor(tags)
        except:
            tags = torch.LongTensor([1]*natoms)

        fixed = (tags == 0).float()

        lattice = crystal.lattice.matrix
        lattice = torch.Tensor([lattice.tolist()])

        # Calculate atom features
        atom_fea = np.vstack([crystal[i].specie.number for i in range(len(crystal))])
        atom_fea = torch.Tensor(atom_fea).squeeze()

        data = Data(pos=pos, cell=lattice, atomic_numbers=atom_fea, natoms=natoms, tags=tags,
                    fixed=fixed, sid=sid, fid=int(frame_number))

    except Exception as e:
        print(f"Error processing item {i}: {e}")
        return None

    # Create a PyTorch Geometric Data object from the retrieved information
    return data


def gen_lmdb(dataset, DB_path: str, map_size: int = 10 ** 9):
    """
    Generate an LMDB database from CIF or LMDB data using multiple workers.

    Args:
        data (str or array): Path to the input data (CIF or LMDB).
        DB_path (str): Path to store the generated LMDB database.
        data_type (str, optional): Type of input data, either "CIF" or "LMDB".
                                   Defaults to None.
        map_size (int, optional): Maximum size of the LMDB database. Defaults to 10 ** 9.
    """

    qty = len(dataset)
    db_syn_freq = 100

    # Open LMDB environment
    db = lmdb.open(DB_path, map_size=map_size, subdir=False, meminit=False, map_async=True)

    db_idx = 0
    txn = db.begin(write=True)
    for index in tqdm(range(db_idx, qty), smoothing=0.1):

        try:
            # Process the item at the given index in the dataset
            data_value = process_item(index, dataset)

            if data_value is None:
                continue

            # Store the processed item in the LMDB environment
            txn.put(key=f"{db_idx}".encode("ascii"), value=pickle.dumps(data_value, protocol=-1))
            db_idx += 1

            # Commit and synchronize the LMDB environment periodically
            if db_idx % db_syn_freq == 0:
                txn.commit()
                db.sync()
                txn = db.begin(write=True)

        except Exception as e:
            print(f"Error processing item index={index}, db_idx={db_idx}: {e}")
            break

    # Commit and synchronize the LMDB environment one last time
    txn.commit()
    db.sync()

    # Store the length of the dataset in the LMDB environment
    txn = db.begin(write=True)
    L = db_idx
    txn.put(key='length'.encode("ascii"), value=pickle.dumps(L, protocol=-1))
    txn.commit()
    db.sync()
    db.close()

    print(f"save to {DB_path}")


if __name__=="__main__":
    pth_energy_lst = [
        # r"D:\Pyprojects\cat_particle_surface\Data\workspace\element_2_4.csv",
        # r"D:\Pyprojects\cat_particle_surface\Data\Data_index\slab_data_3.csv"
        # r"D:\Pyprojects\cat_particle_surface\Data\workspace\element_3_4.csv"
        # r"D:\Pyprojects\cat_particle_surface\Data\workspace\element_3_7.csv"
        # r"D:\Temp_Data\surface_model\surface\slab_element_2\slab_data_2_exp_processed.csv"
        r"D:\Temp_Data\surface_model\surface\all\slab_data_all.csv"
    ]

    DB_pth = r'traj/is2re/relax_all.lmdb'
    POSCAR_pth = r"D:\Temp_Data\surface_model\surface\all\slab"

    POSCAR_data = get_POSCAR_data(pth_energy_lst)
    gen_lmdb(POSCAR_data, DB_pth, map_size=1 * 10 ** 9)
