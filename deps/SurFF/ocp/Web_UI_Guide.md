# Web UI Guide
Please first install the required packages and download the core data files before running the web UI.

## Crystal Files Preparation
Please prepare the crystal files as POSCAR files in a folder.
The POSCAR files should be named using unique integer numbers.

![img0](app_img/img.png) 

![img1](app_img/img_1.png) 


## Start the Web UI
Run the following command in the terminal:
```bash
cd ocp
python app.py
```
Then simply click the browser link provided in the terminal.

## Usage

Upload the POSCAR files by input the folder path and click the "Submit" button.
One example folder is provided in the [../Data/example/crystal_structures](../Data/example/crystal_structures) folder.

To run the example, you should download the SurFF_CoreDataFiles.zip in `SurFF/`, and run following command:

```bash
unzip SurFF_CoreDataFiles.zip
```

![img2](app_img/img_2.png)

Wait for the prediction results. The results will be displayed in the table.