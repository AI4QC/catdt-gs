import os
import pdb
import pickle

import pandas as pd
from scipy.stats import gaussian_kde

from functionals.eval_acc_traj import collect_prediction, collect_true, eval_acc
from functionals.wulff_results import main_wulff, Analyzer


def run_wulff(predict_traj_dir: str, info_dir: str = None,
              crystal_dir: str = r"D:\Temp_Data\surface_model\CIF_MP\bulk_opt_exp",
              save_dir: str = r"results/wulff_results", ):
    # create dir if not exist
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    print(f"predict_traj_dir: {predict_traj_dir}\n"
          f"info_dir: {info_dir}\n"
          f"crystal_dir: {crystal_dir}\n"
          f"save_dir: {save_dir}\n")

    print("Start collect prediction traj data...")

    prediction = collect_prediction(predict_traj_dir)
    info = collect_true(info_dir)

    # if surface_energy is not in info, use as true
    provided_true = True if "surface_energy" in info.columns else False

    if provided_true:
        # combine df
        df = pd.merge(prediction, info, on="slab_id", how="inner", suffixes=("_pred", "_true"))
        eval_acc(prediction, info)
    else:
        df = pd.merge(prediction, info, on="slab_id", how="inner")
        # change column name surface_energy to surface_energy_pred
        df = df.rename(columns={"surface_energy": "surface_energy_pred"})

    results = main_wulff(df, crystal_dir=crystal_dir)

    save_pth = os.path.join(save_dir, "wulff_results.pkl")

    with open(save_pth, 'wb') as f:
        pickle.dump(results, f, pickle.HIGHEST_PROTOCOL)

    print(f"Save results to {save_pth}")

    # also save some wulff_result as csv
    crystal_ids = list(results.keys())

    crystal_id_list = []
    slab_id_list = []
    miller_list = []
    shift_list = []
    surface_energy_true_list = []
    surface_energy_pred_list = []
    area_true_list = []
    area_pred_list = []

    for crystal_id in crystal_ids:
        data = results[crystal_id]

        slab_id_list.extend(data["slab_id"])
        crystal_id_list.extend([crystal_id] * len(data["slab_id"]))
        miller_list.extend(data["miller_index"])
        shift_list.extend(data["shift"])
        surface_energy_pred_list.extend(data["surface_energy_pred"])
        area_pred_list.extend(data["pred"]["color_area"] / sum(data["pred"]["color_area"]))

        if provided_true:
            surface_energy_true_list.extend(data["surface_energy_true"])
            area_true_list.extend(data["true"]["color_area"] / sum(data["true"]["color_area"]))
        else:
            # fill with nan
            surface_energy_true_list.extend([None])
            area_true_list.extend([None])

    if provided_true:
        df = pd.DataFrame({"crystal_id": crystal_id_list,
                           "slab_id": slab_id_list,
                           "miller": miller_list,
                           "shift": shift_list,
                           "surface_energy_true": surface_energy_true_list,
                           "surface_energy_pred": surface_energy_pred_list,
                           "area_true": area_true_list,
                           "area_pred": area_pred_list})
    else:
        df = pd.DataFrame({"crystal_id": crystal_id_list,
                           "slab_id": slab_id_list,
                           "miller": miller_list,
                           "shift": shift_list,
                           "surface_energy_pred": surface_energy_pred_list,
                           "area_pred": area_pred_list})

    # save csv, keep .5f
    df.to_csv(os.path.join(save_dir, "wulff_results.csv"), index=False,
              float_format='%.5f')

    print(f"Save wulff results to {os.path.join(save_dir, 'wulff_results.csv')}")

    pickle_pth = save_pth

    analyser = Analyzer(pickle_pth)

    if provided_true:
        metric = analyser.get_metric()
        metric.to_csv(os.path.join(save_dir, "metric.csv"), index=False)

    print(f"Start save wulff shape to {os.path.join(save_dir, 'wulff_shape')}")

    analyser.save_wulff_shape(crystal_dir=crystal_dir,
                              wulff_save_dir=os.path.join(save_dir, "wulff_shape"))

    print("Complete!")


if __name__ == "__main__":
    predict_dir = r"ocp/traj/001/element3_exp"
    info_dir = r"ocp/traj/test_true"
    # info_dir = r"D:\Pyprojects\ocp\functionals\crystal_surface_id_info"

    run_wulff(predict_dir, info_dir=info_dir,
              crystal_dir=r"D:\Temp_Data\surface_model\CIF_MP\bulk_opt_exp_3",
              save_dir=r"results/element3_exp")
