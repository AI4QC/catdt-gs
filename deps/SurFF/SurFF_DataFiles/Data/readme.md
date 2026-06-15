Please Download the Data by running the following command:
```bash
python download_data.py
```

## Data
The data is stored in the following format:
```txt
Data                        
├── Crystal
│   ├── Opt_crystal
│   │   ├── All
│   │   ├── ID
│   │   └── OOD
│   │       ├── element2
│   │       └── element3
│   └── Raw_crystal
│       ├── All
│       ├── ID
│       └── OOD
│           ├── element2
│           └── element3
├── dataset_generation
│   ├── LMDB
│   └── Model
│       └── save_dir
│           ├── element_1
│           ├── element_2_0
│           ├── element_2_1
│           ├── element_2_2
│           ├── element_2_3
│           ├── element_2_4
│           ├── element_3_0
│           ├── element_3_1
│           ├── element_3_2
│           ├── element_3_3
│           ├── element_3_4
│           ├── element_3_5
│           ├── element_3_6
│           └── element_3_7
├── figures
│   └── fig_save
└── Surface
    ├── Surface_Energy
    │   ├── data
    │   │   ├── ID_testset
    │   │   ├── OOD_testset
    │   │   └── trainset
    │   └── original_data
    │       ├── slab_element_1
    │       │   ├── bulk_1
    │       │   └── slab_1
    │       ├── slab_element_2
    │       │   ├── bulk_2_exp
    │       │   └── slab_2_exp
    │       └── slab_element_3
    │           ├── bulk_3_exp
    │           └── slab_3_exp
    └── Traj
```