# ph_2

This module provides a workflow for processing phenotype microscopy files, performing cell segmentation, and extracting phenotypic features. It is designed to handle phenotype acquisitions and integrate with the previous preprocessing step.

## Contents

1. `ph_2_smk_test.ipynb`: Python script for ensuring correct PH image loading and processing at the tile level.
2. `ph_2.sh`: Bash script to execute the Snakemake workflow.
3. `ph_2.smk`: Main Snakemake file that includes rules for phenotype image processing, cell segmentation, and feature extraction.
4. `ph_2_eval.py`: Python script for evaluating phenotype results and generating quality control plots.

## Usage

1. Test input patterns and processing:
   - Run `ph_2_smk_test.ipynb` to ensure PH images are loaded and processed correctly.
   - Modify the input patterns, channels, and other parameters as needed for your specific setup.

2. Run phenotype processing workflow:
   - Adjust the parameters in `ph_2.smk` based on your specific setup and requirements.
   - Execute `ph_2.sh` to run the Snakemake workflow.

3. Evaluate results:
   - Run `ph_2_eval.py`, which generates quality control plots and statistics for the phenotype processing results.

## Key Features

- Aligns PH images across channels
- Performs illumination correction
- Segments cells and nuclei and identifies cytoplasmic regions
- Extracts phenotypic features using a CellProfiler emulator
- Generates quality control plots and statistics

## Dependencies

- Python libraries: numpy, pandas, matplotlib, seaborn, tifffile
- ops package modules:
  - ops.firesnake: For specific phenotype processing functions used in the Snakemake workflow
      - ops.features: For feature extraction
      - ops.morphology_features: For morphological feature calculations
      - ops.cp_emulator: For emulating CellProfiler functionality
      - ops.cellpose: For cell segmentation using Cellpose
      - ops.process: For image processing operations
  - ops.qc: For quality control functions

## Notes

Given different formats of saving data in either a well based (single or multichannel) or tile based (multichannel) format, we provide example analysis for the two respective cases. These workflows differ largely only their alignment and image correction step.

- Ensure all paths and parameters are correctly set in the Snakemake and Python files before running the workflow.
- The evaluation step is crucial to confirm the quality and correctness of the phenotype processing results.
- Before running the full workflow, it is recommended to test the pipeline on a single well or cycle using `ph_2_smk_test.ipynb`.

## File structure for running

```
plate/
├── ph_2.smk
├── input_ph_tif/
├── ph_2/
│   ├── ph_2_smk_test.ipynb
│   ├── ph_2.sh
│   ├── ph_2_eval.py
│   ├── ph_2_eval.sh
│   ├── csv/*
│   ├── tif/*
│   ├── hdf/*
│   └── qc/*
└── process_ph/
    ├── images/*
    └── tables/*
```

Directories marked with an asterisk contain files generated by this workflow.

Note: 
- Ensure that your input TIFF files from the preprocessing step are in the `input_ph_tif/` and `illumination_correction/` directories before running the workflow.
- The test workflow will generate example output files for one tile in the `ph_2/csv/`, `ph_2/tif/` directories.
- The Snakemake run will generate output files in `process_ph/images/` and `process_ph/tables/` directories.
- The evaluation step will generate output files in `ph_2/hdf/` and `ph_2/qc/` directories.
- Maintaining this file structure is crucial for the correct execution of the workflow. Any changes to the structure may require corresponding adjustments in the scripts.
