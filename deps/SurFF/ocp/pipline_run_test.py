import os
from ase.io import read
import pandas as pd
from ase.visualize import view
from functionals.gen_bulk_slab import process_structure_sep
from functionals.gen_lmdb_relaxation import gen_lmdb
from main import Runner
import yaml
import json
import argparse
import subprocess
from ase.io import read
from ase.visualize import view
from results.run_wulff import run_wulff
import pandas as pd


def run_data_preprocess_test(folder: str, slab_dir: str, bulk_dir: str, lmdb_path: str):
    "pipline running test"

    files = os.listdir(folder)
    print(files, "\n")
    # plot the structure
    structure = read(os.path.join(folder, files[0]), format='vasp')
    # view(structure)
    # print(structure)

    check_dir(slab_dir)
    check_dir(bulk_dir)
    info = []
    for poscar in files:
        slab_info = process_structure_sep(crystal_id=int(poscar),  # a unique id for each crystal
                                          slab_dir=slab_dir,
                                          bulk_dir=bulk_dir,
                                          source_dir=folder)  # the folder where the crystal structures are stored
        info.extend(slab_info)

    # check the generated surface structures
    info = pd.DataFrame(info, columns=['slab_id', 'formula', 'miller_index', 'shift', 'num_atom'])
    print(info, "\n")

    # # check the generated surface structures
    # info = pd.DataFrame(info, columns=['slab_id', 'formula', 'miller_index', 'shift', 'num_atom'])
    # print(info, "\n")

    # Let us check the generated surface structures
    slab_files = os.listdir(slab_dir)
    # print(slab_files,"\n")

    # Let first visualize the one structure using `ase` package
    # structure = read(os.path.join(slab_dir, slab_files[0]), format='vasp')
    # print(structure, "n/")
    # view(structure)

    # surface relaxation
    # We first compile all the surface structures into a single lmdb file

    # prepare a list of slabs info
    dataset = [{'slab_id': files, "POSCAR_pth": os.path.join(slab_dir, files)} for files in slab_files]

    gen_lmdb(dataset=dataset,
             DB_path=lmdb_path,   # the path to save the lmdb file
             map_size=int(1e8))


    print("--------------------pre-process over! ------------------")
    return info


def check_dir(dirs):
    if not os.path.exists(dirs):
        os.makedirs(dirs)
        print(f"!!! mkdir {dirs} !!! \n")


# def save_dict_to_yaml(dict_value: dict, save_path: str):
#     """dict保存为yaml"""
#     with open(save_path, 'w') as file:
#         file.write(yaml.dump(dict_value, allow_unicode=True))


# def read_yaml_to_dict(yaml_path: str, ):
#     with open(yaml_path) as file:
#         dict_value = yaml.load(file.read(), Loader=yaml.FullLoader)
#         return dict_value

# def json_to_dict(args_path:str):
#     with open(args_path, 'r') as f:  
#         args_dict = json.load(f)

#     return args_dict

# def network_infer(args_dict:str, config_dir:str):
#
#     # args_dict = read_yaml_to_dict(args_dict)
#     config = read_yaml_to_dict(config_dir)
#     # args_dict = json_to_dict(args_path=args_dict)
#     # args = argparse.Namespace(**args_dict)
#     args: argparse.Namespace
#
#     Runner()(config)

def run_infer_bash(bash_path: str = "relax_example.sh"):
    print("-----------------Start Infer------------- \n")
    print("Please wait for the inference to finish...")
    try:
        result = subprocess.run(['bash', bash_path], capture_output=True, text=True)
        print("-----------------Finish Infer------------- \n")
        print(result.stdout)
        return result
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running the script: {e}")
        print(f"Error output: {e.stderr}")
        return e


def data_post_process(predicted_data_path: str, preprocess_info, folder: str,
                      post_results_path: str = "results/example",
                      isshow_img: bool = False):
    # predict_dir = r"traj/example"

    # results_saved_path = r"results/example"
    run_wulff(predicted_data_path, info=preprocess_info,
              crystal_dir=folder,
              save_dir=post_results_path)

    results_csv = os.path.join(post_results_path, "wulff_results.csv")
    # results_csv = r"results/example/wulff_results.csv"
    results = pd.read_csv(results_csv)
    print(results)
    img_dir = r"results/example/wulff_shape"

    if isshow_img:
        # let view the predicted shape and surface exposure

        images = os.listdir(img_dir)
        import matplotlib.pyplot as plt
        import time
        for image in images:
            img = plt.imread(os.path.join(img_dir, image))
            plt.imshow(img)
            plt.show()

        time.sleep(10)
        plt.close("all")
    print("------------------- post-process over! ----------------------")

    return results, img_dir


if __name__ == "__main__":
    folder = r'../Data/example/crystal_structures'  # Data/example/crystal_structures
    # Need to sync in the confing
    slab_dir = r"../Data/example/slabs/test_slabs"
    # the folder to save the surface slab structures "Data/example/slabs"
    bulk_dir = r"../Data/example/bulks/test_bulks"
    # the folder to save the corresponding OUC bulk structures (not used in this tutorial)
    # "Data/example/bulks" 
    lmdb_path = r"data/example/surface_relaxation.lmdb"
    info = run_data_preprocess_test(folder=folder, slab_dir=slab_dir, bulk_dir=bulk_dir, lmdb_path=lmdb_path)

    # set this folder which save the inference result, need to sync with relax_example and config__.yml
    traj_folder = r"traj/example"  # ocp/traj/example

    # parser = argparse.ArgumentParser(description='test argparse')
    # # 向创建的解析器对象中添加参数
    # parser.add_argument('--mode', dest='mode', type=str, default='run-relaxations', help='mode')
    # parser.add_argument('-l', '--learning_rate', dest='lr', type=float, default=0.001, help='the learning rate in training')
    # # 参数解析
    # args = parser.parse_args()

    # args_dir = "ocp/configs/config_test/args_web_test.yaml"
    # config_dir = "ocp/configs/config_test/config_web_test.yaml"
    # network_infer(args_dict=args_dir, config_dir=config_dir)

    # network inference
    infer_bash_path = "relax_example.sh"
    run_infer_bash(infer_bash_path)

    #infered data post-process
    trajs = os.listdir(traj_folder)
    print(trajs)

    traj = read(os.path.join(traj_folder,
                             trajs[0]),
                format='traj', index=':')
    # view(traj)
    print(traj[0])

    # post process on the network prediciton results
    predict_dir = r"traj/example"
    results_saved_path = r"results/example"
    results, _ = data_post_process(predict_dir, info, folder, results_saved_path)

    # run_wulff(predict_dir, info=info,
    #             crystal_dir=folder,
    #             save_dir=results_saved_path)

    # results_csv = r"results/example/wulff_results.csv"
    # results = pd.read_csv(results_csv)
    # print(results)

    # # let view the predicted shape and surface exposure
    # image_dir = r"results/example/wulff_shape"
    # images = os.listdir(image_dir)

    print("--------------- Successfuly Running! -------------------")
