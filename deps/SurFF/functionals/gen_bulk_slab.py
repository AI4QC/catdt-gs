import os

import pandas as pd

import time

from concurrent.futures import ProcessPoolExecutor

from functionals.process_structure_sep import process_structure_sep

# data_pth = r"D:\Pyprojects\cat_particle\Data\Crystal_CIF\crystal_structure.pkl"
bulk_dir = r"Data/example/bulks"
slab_dir = r"Data/example/slabs"

# create dir if not exist
if not os.path.exists(bulk_dir):
    os.makedirs(bulk_dir)

if not os.path.exists(slab_dir):
    os.makedirs(slab_dir)


class ProgressBar:
    def __init__(self, total):
        self.total = total
        self.start_time = time.time()
        self.current = 0

    def update(self, current=None):
        if current is not None:
            self.current = current
        else:
            self.current += 1

        progress = self.current / self.total
        bar_length = 60
        block = int(round(bar_length * progress))

        elapsed_time = time.time() - self.start_time
        estimated_total_time = elapsed_time / (self.current + 1e-5) * self.total
        remaining_time = estimated_total_time - elapsed_time

        text = "\rProgress: [{0}] {1}/{2} -- Elapsed Time: {3}s -- Estimated Total Time: {4}s -- Remaining Time: {5}s".format(
            "#" * block + "-" * (bar_length - block),
            self.current, self.total,
            int(elapsed_time),
            int(estimated_total_time),
            int(remaining_time))

        print(text, end="")

def process_data(data, slab_dir, bulk_dir):
    slab_data = []

    data_size = len(data)
    progress_bar = ProgressBar(data_size)

    with ProcessPoolExecutor(max_workers=8) as executor:
        futures = [executor.submit(process_structure_sep, data, i, slab_dir, bulk_dir)
                   for i in range(data_size)]

        for future in futures:
            progress_bar.update()
            slab_data.extend(future.result())  # This will raise ValueError("Something went wrong!")

    return slab_data


def save_slab_data(slab_data, output_file):
    """
    Saves the slab data to a CSV file.

    Args:
        slab_data (list): List of slab data.
        output_file (str): Path to the output CSV file.
    """

    df = pd.DataFrame(slab_data,
                      columns=['slab_id', 'mp_id', "formula", 'crystal_id', 'miller_index', 'shift', "num_atom", "test",
                               'energy'])
    df.to_csv(output_file, index=False)


def main():
    """
    Main function that loads the data, processes it,
    and saves the resulting slab data to a CSV file.
    """
    data_csv = pd.read_csv(r"D:\Pyprojects\cat_particle_surface\Data\Data_index\crystal_id_all.csv")

    # drop unconverged row
    data_csv = data_csv[data_csv['convergence'] == True]

    slab_data = process_data(data_csv, slab_dir, bulk_dir)
    print("\nStart saving data info")

    save_slab_data(slab_data, r'"D:\Pyprojects\cat_particle_surface\Data\Data_index\slab_data_all.csv"')


if __name__ == "__main__":
    main()
