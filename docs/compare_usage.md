
# Example: Comparing Global Datasets with agoge_viz.py

Assuming you have two directories with simulation outputs (each containing HDF5 files or a global VDS):
- Experiment 1: /path/to/experiment1
- Experiment 2: /path/to/experiment2

Use the following command in your terminal to compare the field "rho":

```bash
python3 /Users/smc/Library/CloudStorage/Dropbox/Teaching/CMSE822/cmse822-codex-private/agoge/viz/agoge_viz.py --compare /path/to/experiment1 /path/to/experiment2 --field rho
```

The script will compute and display:
- Maximum and mean absolute error
- Maximum and mean relative error
