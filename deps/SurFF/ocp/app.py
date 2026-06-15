from flask import Flask, render_template, request, jsonify
import os
from pipline_run_test import run_data_preprocess_test, run_infer_bash, data_post_process
import shutil

app = Flask(__name__)

example_data = {
    "Test001": {'crystal_id': ["-"], 'slab_id': ["-"],
                'miller': ["-"], 'shift': ["-"],
                'surface_energy_pred': ["-"], 'area_pred': ["-"]},
}

example_images = {
    "Test001": 'Picture2.png',
}


# # example_data.clear()
# # example_images.clear()
@app.route('/')
def index():
    keys = list(example_data.keys())
    if keys:
        initial_data = example_data[keys[0]]
        initial_image = example_images[keys[0]]
    else:
        initial_data = example_data[keys[0]]
        initial_image = example_images[keys[0]]
    return render_template('index.html', keys=keys, initial_data=initial_data, initial_image=initial_image)


@app.route('/data/<crystal_id>')
def get_data(crystal_id):
    crystal_id_ = int(crystal_id)
    if crystal_id_ not in example_data.keys():
        return jsonify(error="Invalid crystal ID"), 404
    return jsonify(data=example_data[crystal_id_], image=example_images[crystal_id_])


@app.route('/upload-directory', methods=['POST'])
def upload_directory():
    directory_path = request.form.get('directory')
    if not directory_path or not os.path.isdir(directory_path):
        return jsonify(error="Invalid directory path"), 400

    # 处理目录中的文件
    new_data, new_images = run_pipeline(directory_path)
    global example_data, example_images
    example_data.clear()
    example_data.update(new_data)
    example_images.clear()
    example_images.update(new_images)

    return jsonify(success=True)


def mycopy(srcpath, dstpath):
    if not os.path.exists(srcpath):
        print("srcpath not exist!")
    if not os.path.exists(dstpath):
        print("dstpath not exist!")
        os.mkdir(dstpath)
    for root, dirs, files in os.walk(srcpath, True):
        for eachfile in files:
            shutil.copy(os.path.join(root, eachfile), dstpath)


def run_pipeline(input_path: str = r'../Data/example/crystal_structures'):
    # update_num = os.listdir(parent_directory)

    # perform data pre-process 
    # folder = r'../Data/example/crystal_structures' # Data/example/crystal_structures
    folder = input_path
    # Need to sync in the confing
    slab_dir = r"../Data/example/slabs"
    # the folder to save the surface slab structures "Data/example/slabs"
    bulk_dir = r"../Data/example/bulks"
    # the folder to save the corresponding OUC bulk structures (not used in this tutorial)

    # create if not exist and empty the folders before running the pipeline
    if not os.path.exists(slab_dir):
        os.makedirs(slab_dir)
    else:
        shutil.rmtree(slab_dir)
        os.makedirs(slab_dir)

    if not os.path.exists(bulk_dir):
        os.makedirs(bulk_dir)
    else:
        shutil.rmtree(bulk_dir)
        os.makedirs(bulk_dir)

    lmdb_path = r"./data/example/surface_relaxation.lmdb"
    # delete the lmdb file if exist
    if os.path.exists(lmdb_path):
        os.remove(lmdb_path)
    
    lmdb_folder = os.path.dirname(lmdb_path)
    if not os.path.exists(lmdb_folder):
        os.makedirs(lmdb_folder)

    pre_process_info = run_data_preprocess_test(folder=folder, slab_dir=slab_dir, bulk_dir=bulk_dir,
                                                lmdb_path=lmdb_path)

    # set this folder which save the inference result, need to sync with relax_example and config__.yml
    traj_folder = r"traj/example"  # ocp/traj/example
    # clear the folder if exist
    if os.path.exists(traj_folder):
        shutil.rmtree(traj_folder)
    os.makedirs(traj_folder)

    results_saved_path = r"results/example"
    if not os.path.exists(results_saved_path):
        os.makedirs(results_saved_path)
    # perfom network inference
    run_infer_bash()

    # perform post-process
    post_process_resutls, img_list = data_post_process(traj_folder, pre_process_info, folder, results_saved_path)
    output_dict = post_process_resutls.groupby('crystal_id').apply(
        lambda x: x.to_dict(orient='list')).to_dict()
    output_id = list(output_dict.keys())
    # have to put the image in the staic content, otherwise it would not show up
    mycopy('results/example/wulff_shape', 'static/example/wullf_shape')
    img_path_list = os.listdir("static/example/wullf_shape")
    img_path_list.sort()
    img_path_dict = {key: value for key, value in zip(output_id, img_path_list)}

    return output_dict, img_path_dict


# def run_pipeline(directory_path:str):
#     # 这里替换为实际的 pipeline 处理逻辑
#     # 遍历目录中的文件并处理
#     # standard funciton output format
#     new_data = {
#         103: {'crystal_id': [103, 103, 103], 'slab_id': ['101_1', '101_2', '101_3'], 'miller': ['(1, 0, 0)', '(1, 1, 0)', '(1, 1, 1)'], 'shift': [0.1, 0.2, 0.3], 'surface_energy_pred': [0.12, 0.13, 0.14], 'area_pred': [0.01, 0.02, 0.03]},
#         104: {'crystal_id': [104, 104, 104], 'slab_id': ['102_1', '102_2', '102_3'], 'miller': ['(2, 0, 0)', '(2, 2, 0)', '(2, 2, 2)'], 'shift': [0.15, 0.25, 0.35], 'surface_energy_pred': [0.15, 0.16, 0.17], 'area_pred': [0.04, 0.05, 0.06]},
#         105: {'crystal_id': [105, 105, 105], 'slab_id': ['102_1', '102_2', '102_3'], 'miller': ['(2, 0, 0)', '(2, 2, 0)', '(2, 2, 2)'], 'shift': [0.15, 0.25, 0.35], 'surface_energy_pred': [0.15, 0.16, 0.17], 'area_pred': [0.04, 0.05, 0.06]}    
#     }
#     new_images = {
#         103: '103_pred.png',
#         104: '104_pred.png',
#         105: '105_pred.png'
#     }
#     return new_data, new_images

if __name__ == '__main__':
    app.run(debug=True)
