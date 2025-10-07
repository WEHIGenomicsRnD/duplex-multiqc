# WEHI MultiQC Plugin
This is a custom MultiQC plugin developed at WEHI for parsing `.tsv` input files with `Sample`, `Metric`, and `Value` columns. 
It adds per-metric plots to your MultiQC reports.

---

## Features

- Parses `.tsv` files with 3-column structure (`Sample`, `Metric`, `Value`)
- Dynamically generates per-metric plots (e.g., bar plots)
- Works with MultiQC v1.14

---

## Installation and Usage

First, make sure you have [Python](https://www.python.org/) and [MultiQC](https://multiqc.info/) installed.

```bash
pip install multiqc==1.14

Then install the plugin locally:

git clone https://github.com/YOUR-ORG/wehi_multiqc_plugin.git
cd wehi_multiqc_plugin
pip install .

This registers the plugin with MultiQC automatically.

The plugin currently expects tab-delimited .tsv files with the following structure:
Sample    Metric      Value
Note: Column headers must be exactly: Sample, Metric, and Value.

Usage:

multiqc path/to/input/

- To enable verbose/debug mode:
multiqc -v path/to/input/
This helps identify issues with file detection or parsing.
You can also --force to overwrite any existing reports.

Output:

The report will be generated as multiqc_report.html in the current or specified directory.
Additional data will be saved in the multiqc_data/ folder.

Notes:

The plugin automatically creates per-metric plots.
Metric-specific plotting styles can be customized in test_module.py.
