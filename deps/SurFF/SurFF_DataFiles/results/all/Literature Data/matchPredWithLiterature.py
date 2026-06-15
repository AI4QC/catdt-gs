import pandas

pred_csv_pth = "results/all/wulff_results.csv"
literature_excel_pth = "results/all/all.csv"

pred_df = pandas.read_csv(pred_csv_pth)
literature_df = pandas.read_csv(literature_excel_pth)

expose_dict = {}

# for each row in literature_df, all rows in pred_df with the same crystal_id/mp_id
for i in range(len(literature_df)):
    mpid = literature_df.loc[i, "mp_id"]
    # find all surface from pred_df with the same crystal_id
    surface_row = pred_df[pred_df["crystal_id"] == mpid]
    # check if surface_row is empty
    if surface_row.empty:
        continue

    positive_surf = surface_row
    # convert miller index, shift, surface energy, and area_pred to a dict
    miller = positive_surf["miller"].tolist()
    shift = positive_surf["shift"].tolist()
    surface_energy_pred = positive_surf["surface_energy_pred"].tolist()
    area_pred = positive_surf["area_pred"].tolist()

    # combine miller, shift, and area_pred into a dict
    data_list =  [[miller[i], shift[i], surface_energy_pred[i], area_pred[i]] for i in range(len(miller))]
    # write data_list as string
    data_string = f"Miller_index, shift, surface_energy_pred, area_pred\n"
    for data in data_list:
        data_string += f"{data[0]}, {data[1]}, {data[2]}, {data[3]}\n"

    expose_dict[str(mpid)] = data_string

# write expose_dict to literature_df in a new column, "predicted"
literature_df["predicted"] = literature_df["mp_id"].map(expose_dict, na_action="ignore")

# save literature_df to csv
literature_df.to_csv("results/all/comparison.csv", index=False)

