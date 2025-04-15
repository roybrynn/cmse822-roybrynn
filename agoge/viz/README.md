# Agoge Visualization Module

This module provides an object‐oriented framework to load, analyze, and visualize Agoge simulation output. It is designed to work with both individual simulation snapshots and time‐series data. All main functionalities are wrapped in classes such as `AgogeOutput` (for single snapshots), `AgogeTimeSeries` (for analyzing sequential outputs), and `AgogeViz` (a higher‑level class for handling multiple outputs and custom comparisons).

> **Note:** The module relies on a virtual data set (VDS) created from multiple MPI rank files. If the VDS file does not exist, it will be built automatically using the included `vds_builder.py`.

## Prerequisites

Make sure you have Python 3 installed along with the following packages:
- h5py
- numpy
- matplotlib

Install them via pip if needed:

```bash
pip install h5py numpy matplotlib
```

## Module Structure

- **agoge_viz.py**: Contains class definitions:
  - `AgogeOutput`: Loads a simulation output (from a directory containing rank files or a VDS file), lists available fields, loads individual field data, and offers plotting methods (e.g. slice plots, line plots, and field comparisons).
  - `AgogeTimeSeries`: Scans a base directory for sequential outputs, loads outputs by index, and provides methods for animations and time evolution plots.
  - `AgogeViz`: (Skeleton) A wrapper class intended to aggregate multiple outputs and offer higher-level custom comparisons or visualizations.

## Command Line Examples

You can execute the `agoge_viz.py` script directly to visualize or compare Agoge simulation outputs.

### Visualize a Field Slice

To view a slice of a field (for example, `'rho'`) along a chosen axis:
```bash
python agoge_viz.py /path/to/output_directory --field rho --axis z --index 10
```
- `/path/to/output_directory`: Directory containing HDF5 files (or the VDS file).
- `--field`: Name of the field to visualize.
- `--axis`: The axis along which to slice ('x', 'y', or 'z').
- `--index`: Specific slice index (if omitted, the middle slice is used).

### Compare Two Outputs

To compare the same field between two separate output directories:
```bash
python agoge_viz.py --compare /path/to/output1 /path/to/output2 --field rho --axis z --index 10
```
This command compares the `'rho'` field between the outputs in the two directories and displays the difference (with optionally, logarithmic error scale if implemented).

## Usage in a Jupyter Notebook

You can also use this module interactively. Here are some examples:

### 1. Working with a Single Output

```python
import matplotlib.pyplot as plt
from agoge_viz import AgogeOutput

# Initialize AgogeOutput from a directory
output = AgogeOutput('/path/to/output_directory')

# List available fields
fields = output.list_fields()
print("Available fields:", fields)

# Plot a slice of the 'rho' field along the z-axis
fig, ax = output.plot_slice('rho', axis='z', index=10)
plt.show()
```

### 2. Comparing Two Outputs

```python
from agoge_viz import AgogeOutput
import matplotlib.pyplot as plt

# Load two different outputs
output1 = AgogeOutput('/path/to/output1')
output2 = AgogeOutput('/path/to/output2')

# Compare the 'rho' field along the z-axis at the desired slice index
fig, axes = output1.compare_with(output2, field_name='rho', axis='z', index=10)
plt.show()
```

### 3. Analyzing a Time Series

```python
from agoge_viz import AgogeTimeSeries
import matplotlib.pyplot as plt

# Initialize a time series from a base directory containing multiple outputs
time_series = AgogeTimeSeries('/path/to/base_directory')

# Plot the time evolution of the mean value of the 'rho' field
fig, ax = time_series.plot_time_evolution('rho', stat='mean')
plt.show()

# For an animated visualization of a field over time (e.g., along z-axis)
anim = time_series.animate_field('rho', axis='z', index=10)
# To display the animation in a Jupyter Notebook (requires HTML5 video support)
from IPython.display import HTML
HTML(anim.to_html5_video())
```

### 4. Using the High-Level AgogeViz Class

*(Note: This class provides a scaffold for more advanced visualizations and comparisons.)*

```python
from agoge_viz import AgogeViz, AgogeOutput

# Create a list of outputs (e.g., load two outputs for comparison)
outputs = [AgogeOutput('/path/to/output1'), AgogeOutput('/path/to/output2')]

# Initialize AgogeViz (optionally, a custom compare module can be provided)
viz = AgogeViz(outputs)

# (Placeholder) Use methods like plot_slice or compare_outputs when fully implemented.
# viz.plot_slice('rho', slice_axis='z', slice_index=10)
# viz.compare_outputs(0, 1, 'rho')
```

## Contributing and Extensions

The refactored code follows an object‑oriented design to improve modularity and maintainability. Additional improvements and features (such as more advanced animation, statistical analysis, and custom comparison routines) are welcome. Feel free to extend the `AgogeViz` class to suit your analysis needs.

## License

This project is licensed under the MIT License.

For further questions or contributions, please visit the project repository.
