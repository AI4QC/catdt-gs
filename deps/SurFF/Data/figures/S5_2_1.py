# get infomation based on mp-id and plot the figure
import pandas as pd
from matplotlib import pyplot as plt
from mp_api.client import MPRester


def fetch_material_data(mp_ids):
    with MPRester("jBmMd3ojMrsTVBZsDRtddygpko6T20Eq") as mpr:
        # do stuff with mpr...
        docs = mpr.summary.search(
            material_ids=mp_ids,
            fields=['symmetry', 'elements', 'energy_above_hull'],
        )

    return docs

if __name__ == "__main__":

    all_crystal = r"Data/Crystal/crystal_bulk_all.csv"
    crystal_id_df = pd.read_csv(all_crystal)
    mp_id_list = crystal_id_df['mp_id'].tolist()

    docs = fetch_material_data(mp_id_list)

    save_dir = r"Data/figures/fig_save"

    # change global matplotlib settings
    # change dpi to 600, tight layout, and font size to 12
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['savefig.bbox'] = 'tight'
    plt.rcParams['font.size'] = 12


    # figure S5_2.1 hostogram of the number of elements in the materials (bar plot)
    num_elements = [len(d.elements) for d in docs]
    categories = list(set(num_elements))
    counts = [num_elements.count(c) for c in categories]
    # adjust figure width based on the number of categories
    plt.figure(figsize=(len(categories)*0.5 + 2, 5))
    plt.bar(categories, counts, color='gray')
    plt.xlabel("Number of elements")
    plt.ylabel("Number of materials")
    plt.title("Number of elements in crystals")
    # show the number of materials in each bar
    for i in range(len(categories)):
        plt.text(categories[i], counts[i], counts[i], ha='center', va='bottom')

    plt.savefig(f"{save_dir}/S5_2-1.png")

    # figure S2.1-2, hostogram of the crystal system in the materials (bar plot)
    crystal_system = [d.symmetry.crystal_system.value for d in docs]
    categories = list(set(crystal_system))
    counts = [crystal_system.count(c) for c in categories]

    # adjust figure width based on the number of categories
    plt.figure(figsize=(len(categories)*0.5 + 2, 5))
    plt.bar(categories, counts, color='gray')
    plt.xlabel("Crystal system")
    plt.ylabel("Number of materials")
    plt.title("Crystal system in crystals")
    # add actual number of materials in each bar
    for i in range(len(categories)):
        plt.text(categories[i], counts[i], counts[i], ha='center', va='bottom')

    # change x label to vertical layout
    plt.xticks(rotation=30)
    plt.savefig(f"{save_dir}/S5_2-2.png")

    # figure S2.1-3, hostogram of the space group in the materials (bar plot)
    point_group = [d.symmetry.point_group for d in docs]
    categories = list(set(point_group))
    counts = [point_group.count(c) for c in categories]
    # adjust figure width based on the number of categories
    plt.figure(figsize=(len(categories)*0.5 + 2, 5))
    plt.bar(categories, counts, color='gray')
    plt.xlabel("Point group")
    plt.ylabel("Number of materials")
    plt.title("Point group in crystals")
    # add actual number of materials in each bar
    for i in range(len(categories)):
        plt.text(categories[i], counts[i], counts[i], ha='center', va='bottom')
    plt.savefig(f"{save_dir}/S5_2-3.png")

    # figure S2.1-4, hostogram of elements occurance in the materials (bar plot)
    font_size = 28
    elements = [element.value for d in docs for element in d.elements]
    categories = list(set(elements))
    counts = [elements.count(c) for c in categories]
    # adjust figure width based on the number of categories
    plt.figure(figsize=(len(categories)*0.5 + 2, 5))
    # change font size to 20

    plt.bar(categories, counts, color='gray')
    plt.xlabel("Element", fontdict={'fontsize': font_size})
    plt.ylabel("Number of materials", fontdict={'fontsize': font_size})
    plt.title("Element occurance in crystals", fontdict={'fontsize': font_size})
    # add actual number of materials in each bar
    for i in range(len(categories)):
        plt.text(categories[i], counts[i], counts[i], ha='center', va='bottom')

    plt.savefig(f"{save_dir}/S5_2-4.png")

    # figure S2.1-5, hostogram of the energy above hull in the materials, just histogram
    plt.close()
    energy_above_hull = [d.energy_above_hull for d in docs]

    plt.hist(energy_above_hull, bins=20, color='gray')
    plt.xlabel("Energy above hull (eV/atom)")
    plt.ylabel("Number of materials")
    plt.title("Energy above hull in crystals")
    plt.savefig(f"{save_dir}/S5_2-5.png")




