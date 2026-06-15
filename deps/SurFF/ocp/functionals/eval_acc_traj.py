import os
from ase.io import Trajectory
import numpy as np
import pandas as pd
from tqdm import tqdm


def collect_prediction(dir:str):

    # get all .traj file in dir
    traj_file_list = []
    for file in os.listdir(dir):
        if file.endswith(".traj"):
            traj_file_list.append(os.path.join(dir, file))

    slab_id_list = []
    crystal_id_list = []
    surface_energy_list = []

    # read .traj file
    for traj_file in tqdm(traj_file_list):

        slab_id = os.path.basename(traj_file).split(".")[0]
        crystal_id = slab_id.split("_")[0]

        traj = Trajectory(traj_file)
        traj_relax = traj[-1]
        energy = traj_relax.get_potential_energy()
        cell = traj_relax.get_cell()
        # get area by cross product of cell vectors
        area = np.linalg.norm(np.cross(cell[0,:], cell[1,:]))

        surface_energy = energy / (2 * area)

        slab_id_list.append(slab_id)
        crystal_id_list.append(crystal_id)
        surface_energy_list.append(surface_energy)

    return pd.DataFrame({"slab_id": slab_id_list, "surface_energy": surface_energy_list, "crystal_id": crystal_id_list})


def collect_true(dir: str):
    # read all .csv file in the dir
    csv_file_list = []
    for file in os.listdir(dir):
        if file.endswith(".csv"):
            csv_file_list.append(os.path.join(dir, file))

    # read all .csv file and concat
    df_list = []
    for csv_file in csv_file_list:
        df_list.append(pd.read_csv(csv_file))
    df = pd.concat(df_list, axis=0, ignore_index=True)

    return df

def eval_acc(pred, true):
    """

    :param pred: data_frame with columns: slab_id, surface_energy
    :param true: data_frame with columns: slab_id, surface_energy
    :return: mae of surface energy prediction
    """
    df = pd.merge(pred, true, on="slab_id", how="inner")
    mae = np.mean(np.abs(df["surface_energy_x"] - df["surface_energy_y"]))
    print(f"mae: {mae}")

    return mae


