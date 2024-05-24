# SlideCNA Analysis

This repository contains vignettes to guide you in reproducing SlideCNA results from our manuscript and to get you started on implementing SlideCNA on your own Slide-seq-like spatial transcriptomics datasets. We demonstrate SlideCNA usage on two Slide-seq datasets on three modes: spatial, non-spatial (for non-spatial transcriptomics single cell-like data), and spatial with bead splitting via [TACCO](https://www.nature.com/articles/s41587-023-01657-3).

For more information about the SlideCNA workflow, please refer to documentation here: https://github.com/dkzhang777/SlideCNA .

### Data
Data for the SlideCNA vignettes may be downloaded via [Zenodo](https://doi.org/10.5281/zenodo.10658097).

### Conda environments
Create a clean environment with conda using the SlideCNA_env.yml file to run the SlideCNA implementation vignettes:
```
conda env create -f "https://github.com/dkzhang777/SlideCNA_Analysis/envs/SlideCNA_env.yml"
```

and the bead_split_env.yml file to run the TACCO bead splitting vignette:
```
conda env create -f "https://github.com/dkzhang777/SlideCNA_Analysis/bead_split_env.yml"
```

(For older versions of conda one needs to download the environment.yml and use the local file for installation.)

### Example Notebooks
1. SlideCNA implementation on Slide-seq metastatic breast cancer (sample HTAPP-895-SMP-7359)
2. SlideCNA implementation on Single-nucleus RNA-seq metastatic breast cancer (sample HTAPP-895-SMP-7359)
3. Bead splitting of Slide-seq metastatic breast cancer with [TACCO](https://www.nature.com/articles/s41587-023-01657-3) (sample HTAPP-895-SMP-7359)
4. SlideCNA implementation on bead split Slide-seq metastatic breast cancer (sample HTAPP-895-SMP-7359)