import itertools
import math
import os
import pdb
import re


import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt, animation
from pymatgen.analysis.wulff import WulffShape
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from tqdm import tqdm
# from wulffpack import SingleCrystal


class NewWulffShape(WulffShape):
    def get_plot(
            self,
            color_set="PuBu",
            grid_off=True,
            axis_off=True,
            show_area=False,
            alpha=1,
            off_color="red",
            direction=None,
            bar_pos=(0.75, 0.15, 0.05, 0.65),
            bar_on=False,
            units_in_JPERM2=True,
            legend_on=True,
            aspect_ratio=(10, 8),
            custom_colors=None,
    ):
        """
        Get the Wulff shape plot.

        Args:
            color_set: default is 'PuBu'
            grid_off (bool): default is True
            axis_off (bool): default is True
            show_area (bool): default is False
            alpha (float): chosen from 0 to 1 (float), default is 1
            off_color: Default color for facets not present on the Wulff shape.
            direction: default is (1, 1, 1)
            bar_pos: default is [0.75, 0.15, 0.05, 0.65]
            bar_on (bool): default is False
            legend_on (bool): default is True
            aspect_ratio: default is (8, 8)
            custom_colors ({(h,k,l}: [r,g,b,alpha]}): Customize color of each
                facet with a dictionary. The key is the corresponding Miller
                index and value is the color. Undefined facets will use default
                color site. Note: If you decide to set your own colors, it
                probably won't make any sense to have the color bar on.
            units_in_JPERM2 (bool): Units of surface energy, defaults to
                Joules per square meter (True)

        Return:
            (matplotlib.pyplot)
        """
        from mpl_toolkits.mplot3d import art3d

        colors = self._get_colors(color_set, alpha, off_color, custom_colors=custom_colors or {})
        color_list, color_proxy, color_proxy_on_wulff, miller_on_wulff, e_surf_on_wulff = colors

        if not direction:
            # If direction is not specified, use the miller indices of
            # maximum area.
            direction = max(self.area_fraction_dict.items(), key=lambda x: x[1])[0]

        fig = plt.figure()
        fig.set_size_inches(aspect_ratio[0], aspect_ratio[1])
        azim, elev = 0, 0.15

        wulff_pt_list = self.wulff_pt_list

        ax_3d = fig.add_subplot(projection="3d")
        ax_3d.view_init(azim=azim, elev=elev)
        fig.add_axes(ax_3d)

        for plane in self.facets:
            # check whether [pts] is empty
            if len(plane.points) < 1:
                # empty, plane is not on_wulff.
                continue
            # assign the color for on_wulff facets according to its
            # index and the color_list for on_wulff
            plane_color = color_list[plane.index]
            pt = self.get_line_in_facet(plane)
            # plot from the sorted pts from [simpx]
            tri = art3d.Poly3DCollection([pt])
            tri.set_color(plane_color)
            tri.set_edgecolor("#808080")
            ax_3d.add_collection3d(tri)

        # set ranges of x, y, z
        # find the largest distance between on_wulff pts and the origin,
        # to ensure complete and consistent display for all directions
        r_range = max(np.linalg.norm(x) for x in wulff_pt_list)
        ax_3d.set_xlim([-r_range * 1.1, r_range * 1.1])
        ax_3d.set_ylim([-r_range * 1.1, r_range * 1.1])
        ax_3d.set_zlim([-r_range * 1.1, r_range * 1.1])
        # add legend
        if legend_on:
            if show_area:
                legend_list = []

                miller_index = list(self.area_fraction_dict.keys())
                miller_area = list(self.area_fraction_dict.values())
                exposed = ["H" if x > 0.1 else "M" if x > 0.01 else "L" for x in miller_area]

                for i in range(len(miller_index)):
                    legend_list.append(
                        f"{miller_index[i]} {miller_area[i]:.3f} ({exposed[i]})")

                ax_3d.legend(
                    color_proxy,
                    legend_list,
                    loc="center left",
                    bbox_to_anchor=(-0.2, 0.5),
                    fancybox=True,
                    shadow=False,
                )
            else:
                ax_3d.legend(
                    color_proxy_on_wulff,
                    miller_on_wulff,
                    loc="upper center",
                    bbox_to_anchor=(0.5, 1),
                    ncol=3,
                    fancybox=True,
                    shadow=False,
                )
        ax_3d.set(xlabel="x", ylabel="y", zlabel="z")

        # Add color bar
        if bar_on:
            cmap = plt.get_cmap(color_set)
            cmap.set_over("0.25")
            cmap.set_under("0.75")
            bounds = [round(e, 2) for e in e_surf_on_wulff]
            bounds.append(1.2 * bounds[-1])
            norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
            # display surface energies
            ax1 = fig.add_axes(bar_pos)
            cbar = mpl.colorbar.ColorbarBase(
                ax1,
                cmap=cmap,
                norm=norm,
                boundaries=[0, *bounds, 10],
                extend="both",
                ticks=bounds[:-1],
                spacing="proportional",
                orientation="vertical",
            )
            units = "$J/m^2$" if units_in_JPERM2 else r"$eV/\AA^2$"
            cbar.set_label(f"Surface Energies ({units})", fontsize=25)

        if grid_off:
            ax_3d.grid("off")
        if axis_off:
            ax_3d.axis("off")
        return plt


def main_wulff(data, crystal_dir=None):
    results_dict = {}
    crystal_id_list = get_crystal_id_list(data)

    print(f"total {len(crystal_id_list)} crystals, {len(data)} slabs, start Wulff analysis")

    for crystal_id in tqdm(crystal_id_list):
        wulff_info = get_wulff_info(data, crystal_id, crystal_dir)
        results_dict[crystal_id] = wulff_info

    return results_dict


def get_wulff_info(data, crystal_id, crystal_dir):
    crystal_pth = os.path.join(crystal_dir, f"{crystal_id}")
    crystal_struc = Poscar.from_file(crystal_pth).structure

    df = data[data["crystal_id"] == crystal_id]

    # if surface_energy_true provided
    if "surface_energy_true" in df.columns:

        miller_index = df["miller_index"].values.tolist()
        miller_index = [str_to_tuple(i) for i in miller_index]

        surface_energy_true = df["surface_energy_true"].values.tolist()
        surface_energy_pred = df["surface_energy_pred"].values.tolist()
        wulffshape_true = NewWulffShape(crystal_struc.lattice, miller_index, surface_energy_true)
        wulffshape_pred = NewWulffShape(crystal_struc.lattice, miller_index, surface_energy_pred)

        new_df = df.copy()
        # create two new columns true_area and pred_area to new_df
        new_df["true_area"] = np.zeros((len(new_df), 1))
        new_df["pred_area"] = np.zeros((len(new_df), 1))
        new_df["miller_index"] = miller_index

        # load area from wulffshape_true
        for i in range(len(wulffshape_true.color_area)):
            surface_energy = wulffshape_true.e_surf_list[i]
            miller = wulffshape_true.miller_list[i]
            # # tuple (1,1,1) to int 111
            # miller = "".join(str(x) for x in miller)
            # find row indx based on miller and surface_energy, give some threshold to float comparison

            idx = new_df[(new_df["miller_index"] == miller) & (new_df["surface_energy_true"] - surface_energy < 1e-5)].index

            new_df.loc[idx, "true_area"] = wulffshape_true.color_area[i]

        # load area from wulffshape_pred
        for i in range(len(wulffshape_pred.color_area)):
            surface_energy = wulffshape_pred.e_surf_list[i]
            miller = wulffshape_pred.miller_list[i]
            # tuple (1,1,1) to int 111
            # miller = "".join(str(x) for x in miller)
            # find row indx based on miller and surface_energy, give some threshold to float comparison

            idx = new_df[(new_df["miller_index"] == miller) & (new_df["surface_energy_pred"] - surface_energy < 1e-5)].index

            new_df.loc[idx, "pred_area"] = wulffshape_pred.color_area[i]

        return {"true": wulffshape_true.__dict__, "pred": wulffshape_pred.__dict__,
                "miller_index": new_df["miller_index"].values.tolist(),
                "surface_energy_true": new_df["surface_energy_true"].values.tolist(),
                "surface_energy_pred": new_df["surface_energy_pred"].values.tolist(),
                "true_area": new_df["true_area"].values.tolist(),
                "pred_area": new_df["pred_area"].values.tolist(),
                "slab_id": new_df["slab_id"].values.tolist(),
                "shift": new_df["shift"].values.tolist()
                }

    # if surface_energy_true not provided
    else:
        miller_index = df["miller_index"].values.tolist()
        miller_index = [str_to_tuple(i) for i in miller_index]
        surface_energy_pred = df["surface_energy_pred"].values.tolist()
        wulffshape_pred = NewWulffShape(crystal_struc.lattice, miller_index, surface_energy_pred)

        new_df = df.copy()
        # create new column pred_area to new_df
        new_df["pred_area"] = np.zeros((len(new_df), 1))
        new_df["miller_index"] = miller_index

        # load area from wulffshape_pred
        for i in range(len(wulffshape_pred.color_area)):
            surface_energy = wulffshape_pred.e_surf_list[i]
            miller = wulffshape_pred.miller_list[i]
            # tuple (1,1,1) to int 111
            # miller = "".join(str(x) for x in miller)
            # find row indx based on miller and surface_energy, give some threshold to float comparison

            idx = new_df[(new_df["miller_index"] == miller) & (new_df["surface_energy_pred"] - surface_energy < 1e-5)].index

            new_df.loc[idx, "pred_area"] = wulffshape_pred.color_area[i]


        return {"pred": wulffshape_pred.__dict__,
                "miller_index": new_df["miller_index"].values.tolist(),
                "surface_energy_pred": new_df["surface_energy_pred"].values.tolist(),
                "pred_area": new_df["pred_area"].values.tolist(),
                "slab_id": new_df["slab_id"].values.tolist(),
                "shift": new_df["shift"].values.tolist()
                }


def str_to_tuple(s):
    matches = re.findall(r'-?\d', s)
    numbers = [int(num) for num in matches]
    return tuple(numbers)


def get_crystal_id_list(data):
    # get all unique crystal id in the dataset

    all_id = data["crystal_id"].values
    all_id = list(set(all_id))
    return all_id


class Analyzer:
    def __init__(self, pickle_pth):
        self.data: dict = self.load_pickle(pickle_pth)
        self.crystal_id_list = self.data.keys()
        self.table_data = self.get_table_data()

    def load_pickle(self, pickle_pth):
        import pickle
        with open(pickle_pth, "rb") as f:
            data = pickle.load(f)
        return data

    def get_table_data(self):
        table_data = []

        # if surface_energy_true provided
        if "true" in self.data[list(self.data.keys())[0]].keys():

            for crystal_id in self.crystal_id_list:
                results = self.data[crystal_id]
                pred = results["pred"]
                true = results["true"]

                df1 = pd.DataFrame({'miller_index': list(true["miller_list"]),
                                    "true": true["color_area"]/sum(true["color_area"])})
                df2 = pd.DataFrame({'miller_index': list(pred["miller_list"]),
                                    "pred": pred["color_area"]/sum(pred["color_area"])})

                # remove duplicated miller_index in df, keep the higher surface area
                miller_index, surface_energy_true = remove_duplicated_miller(df1["miller_index"], df1["true"])
                df1 = pd.DataFrame({'miller_index': miller_index,
                                    "true": surface_energy_true})

                miller_index, surface_energy_pred = remove_duplicated_miller(df2["miller_index"], df2["pred"])
                df2 = pd.DataFrame({'miller_index': miller_index,
                                    "pred": surface_energy_pred})

                # Merge the DataFrames using outer join to include NaN values
                results = df1.merge(df2, on='miller_index', how='outer')
                results["crystal_id"] = [crystal_id] * len(results)

                table_data.append(results)

        # if surface_energy_true not provided
        else:
            for crystal_id in self.crystal_id_list:
                results = self.data[crystal_id]
                pred = results["pred"]

                df1 = pd.DataFrame({'miller_index': list(pred["miller_list"]),
                                    "pred": pred["color_area"]})

                results = df1
                results["crystal_id"] = [crystal_id] * len(results)

                table_data.append(results)

        table_data = pd.concat(table_data)

        return table_data

    def get_metric(self, expose_threshold=0.1):
        # get overall_acc, precision and recall

        table_data = self.table_data

        # l "low", m "medium", h "high"
        # l<1e-3, 1e-3<=m<expose_threshold, h>=expose_threshold
        table_data["type_true"] = "l"
        table_data.loc[table_data["true"] > 1e-3, "type_true"] = "m"
        table_data.loc[table_data["true"] > expose_threshold, "type_true"] = "h"

        table_data["type_pred"] = "l"
        table_data.loc[table_data["pred"] > 1e-3, "type_pred"] = "m"
        table_data.loc[table_data["pred"] > expose_threshold, "type_pred"] = "h"

        l_acc = table_data[table_data["type_true"] == "l"]["type_pred"].value_counts()["l"] / \
                table_data["type_true"].value_counts()["l"]
        m_acc = table_data[table_data["type_true"] == "m"]["type_pred"].value_counts()["m"] / \
                table_data["type_true"].value_counts()["m"]
        h_acc = table_data[table_data["type_true"] == "h"]["type_pred"].value_counts()["h"] / \
                table_data["type_true"].value_counts()["h"]

        lm_rate = table_data[table_data["type_true"] == "l"]["type_pred"].value_counts()["m"] / \
                  table_data["type_true"].value_counts()["l"]
        hm_rate = table_data[table_data["type_true"] == "h"]["type_pred"].value_counts()["m"] / \
                  table_data["type_true"].value_counts()["h"]

        mae = np.mean(np.abs(table_data["true"] - table_data["pred"]))

        overall_acc = table_data[table_data["type_true"] == table_data["type_pred"]]["type_pred"].count() / len(
            table_data)


        # get high N rate

        correct_5 = 0
        all_5 = 0

        correct_3 = 0
        all_3 = 0

        for crystal_id in self.crystal_id_list:

            data = self.table_data[self.table_data["crystal_id"] == crystal_id]
            # get miller index with highest 5 true surface area
            miller_index_high_5_true = data.sort_values(by="true", ascending=False)["miller_index"].values[:5].tolist()
            miller_index_high_5_pred = data.sort_values(by="pred", ascending=False)["miller_index"].values[:5].tolist()

            # pred in true
            for miller in miller_index_high_5_pred:
                if miller in miller_index_high_5_true:
                    correct_5 += 1

            all_5 += 5

            # get miller index with highest 3 true surface area
            miller_index_high_3_true = data.sort_values(by="true", ascending=False)["miller_index"].values[:3].tolist()
            miller_index_high_3_pred = data.sort_values(by="pred", ascending=False)["miller_index"].values[:3].tolist()

            # pred in true
            for miller in miller_index_high_3_pred:
                if miller in miller_index_high_3_true:
                    correct_3 += 1

            all_3 += 3


        acc = pd.DataFrame({"overall_acc": [overall_acc],
                            "l_acc": [l_acc],
                            "lm_rate": [lm_rate + l_acc],
                            "m_acc": [m_acc],
                            "h_acc": [h_acc],
                            "hm_rate": [hm_rate + h_acc],
                            "mae": [mae],
                            "high_5_rate": [correct_5 / all_5],
                            "high_3_rate": [correct_3 / all_3]})
        print(acc)

        return acc

    def save_wulff_shape(self, crystal_dir, wulff_save_dir):

        if not os.path.exists(wulff_save_dir):
            os.makedirs(wulff_save_dir)

        color_dict = self.get_color_dict()

        for crystal_id in tqdm(self.crystal_id_list):
            crystal_pth = os.path.join(crystal_dir, f"{crystal_id}")
            crystal_struc = Poscar.from_file(crystal_pth).structure
            crystal_struc = SpacegroupAnalyzer(crystal_struc).get_conventional_standard_structure()

            results = self.data[crystal_id]

            miller_index = results["miller_index"]
            surface_energy_pred = results["surface_energy_pred"]

            miller_index, surface_energy_pred = remove_duplicated_miller(miller_index, surface_energy_pred)

            wulffshape_pred = NewWulffShape(crystal_struc.lattice, miller_index, surface_energy_pred)

            p = wulffshape_pred.get_plot(custom_colors=color_dict, show_area=True, off_color="black")
            p.savefig(os.path.join(wulff_save_dir, f"{crystal_id}_pred.png"))
            plt.close()

            if "surface_energy_true" in results.keys():
                surface_energy_true = results["surface_energy_true"]

                miller_index, surface_energy_true = remove_duplicated_miller(miller_index, surface_energy_true)

                wulffshape_true = NewWulffShape(crystal_struc.lattice, miller_index, surface_energy_true)
                p = wulffshape_true.get_plot(custom_colors=color_dict, show_area=True, off_color="black")
                p.savefig(os.path.join(wulff_save_dir, f"{crystal_id}_true.png"))
                plt.close()

        print(f"save wulff shape to {wulff_save_dir}")

        pass

    def save_wulff_shape_gif(self, crystal_dir, wulff_save_dir, frame_num=36):
        import imageio
        for crystal_id in tqdm(self.crystal_id_list):
            crystal_pth = os.path.join(crystal_dir, f"{crystal_id}")
            crystal_struc = Poscar.from_file(crystal_pth).structure
            crystal_struc = SpacegroupAnalyzer(crystal_struc).get_conventional_standard_structure()

            results = self.data[crystal_id]

            miller_index = results["miller_index"]

            surface_energy_pred = results["surface_energy_pred"]
            miller_index, surface_energy_pred = remove_duplicated_miller(miller_index, surface_energy_pred)

            if "surface_energy_true" in results.keys():
                surface_energy_true = results["surface_energy_true"]
                miller_index, surface_energy_true = remove_duplicated_miller(miller_index, surface_energy_true)

            if not os.path.exists(wulff_save_dir):
                os.makedirs(wulff_save_dir)

            color_dict = self.get_color_dict()

            pred_png_list = []
            true_png_list = []

            for i in range(frame_num):
                wulffshape_pred = NewWulffShape(crystal_struc.lattice, miller_index, surface_energy_pred)
                p = wulffshape_pred.get_plot(custom_colors=color_dict, legend_on=False,
                                             azim=(i + 1) * (360.0 / frame_num), elev=30)
                pred_png = os.path.join(wulff_save_dir, f"{crystal_id}_{i}_pred.png")
                pred_png_list.append(pred_png)
                p.savefig(pred_png)
                plt.close()

                if "surface_energy_true" in results.keys():

                    wulffshape_true = NewWulffShape(crystal_struc.lattice, miller_index, surface_energy_true)
                    p = wulffshape_true.get_plot(custom_colors=color_dict, legend_on=False,
                                                 azim=(i + 1) * (360.0 / frame_num), elev=30)
                    true_png = os.path.join(wulff_save_dir, f"{crystal_id}_{i}_true.png")
                    true_png_list.append(true_png)
                    p.savefig(true_png)
                    plt.close()

            imageio.mimsave(os.path.join(wulff_save_dir, f"{crystal_id}_pred.gif"),
                            [imageio.imread(file) for file in pred_png_list], loop=0, fps=int(frame_num / 4))
            if "surface_energy_true" in results.keys():
                imageio.mimsave(os.path.join(wulff_save_dir, f"{crystal_id}_true.gif"),
                                [imageio.imread(file) for file in true_png_list], loop=0, fps=int(frame_num / 4))

            if "surface_energy_true" in results.keys():
                for file in true_png_list:
                    os.remove(file)

            for file in pred_png_list:
                os.remove(file)

        print(f"save wulff shape to {wulff_save_dir}")

    def get_color_dict(self, cmap="YlGnBu"):

        def has_common_factor(t):
            g = math.gcd(t[0], t[1])
            g = math.gcd(g, t[2])
            return g != 1

        numbers = (-2, -1, 0, 1, 2)
        combinations = list(itertools.product(numbers, repeat=3))
        combinations = [t for t in combinations if not has_common_factor(t)]
        sorted_combinations = sorted(combinations, key=lambda x: sum(map(abs, x)))
        cmap = plt.cm.get_cmap(cmap, len(sorted_combinations))
        colors = [cmap(i)[:3] for i in range(len(sorted_combinations))]

        color_dict = {sorted_combinations[i]: colors[i] for i in range(len(sorted_combinations))}
        return color_dict


def remove_duplicated_miller(miller_index:[tuple], surface_energy: [float]):
    new_miller_index = []
    new_surface_energy = []

    for i in range(len(miller_index)):
        miller = miller_index[i]
        energy = surface_energy[i]

        if miller not in new_miller_index:
            new_miller_index.append(miller)
            new_surface_energy.append(energy)
        else:
            # keep the lower surface energy
            idx = new_miller_index.index(miller)
            if energy > new_surface_energy[idx]:
                new_surface_energy[idx] = energy


    return new_miller_index, new_surface_energy