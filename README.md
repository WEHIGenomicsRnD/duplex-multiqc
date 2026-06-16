# duplex-multiqc

A [MultiQC](https://multiqc.info/) plugin for parsing and visualising duplex
sequencing QC metrics. It reads per-sample, per-metric tabular data (CSV or
TSV) and generates interactive bar plots in your MultiQC report. The pipeline
works on output from the
[calculate-duplex-metrics](https://github.com/WEHIGenomicsRnD/calculate-duplex-metrics)
tool that has been integrated into the WEHI Genomics R&D
[duplex-seq-pipeline](https://github.com/WEHIGenomicsRnD/duplex-seq-pipeline).

---

## Features

- Parses `.csv` and `.tsv` files with a `sample`, `metric`, `value` column structure
- Column headers are case-insensitive (`Sample` and `sample` both work)
- Dynamically generates per-metric bar plots
- Groups related metrics into combined plots:
  - **On-target rate** — `on_target_rate_raw`, `on_target_rate_duplex`
  - **GC Metrics** — `gc_single`, `gc_both`
  - **Within-family stats** — `family_max`, `family_median`, `family_mean`
  - **Family size comparison** — `families_gt1`, `single_families`,
    `paired_families`, `paired_and_gt1`, `total_families`
- Adds all metrics to the MultiQC general statistics table

---

## Requirements

- Python ≥ 3.8
- MultiQC == 1.14

---

## Installation

### With uv (recommended)

[uv](https://docs.astral.sh/uv/) provides fast, reproducible environment management.

```bash
# Clone the repository
git clone https://github.com/WEHIGenomicsRnD/duplex-multiqc.git
cd duplex-multiqc

# Create a virtual environment and install the plugin (and its dependencies)
uv sync

# Run MultiQC using the managed environment
uv run multiqc path/to/input/
```

To include development dependencies (e.g. for running tests):

```bash
uv sync --group dev
```

### With pip

```bash
git clone https://github.com/WEHIGenomicsRnD/duplex-multiqc.git
cd duplex-multiqc
pip install .
```

---

## Input Format

The plugin matches files named `*.csv` or `*.tsv` whose first line contains the
headers `sample,metric,value` (CSV) or `sample\tmetric\tvalue` (TSV). Column
names are **case-insensitive**.

### CSV example

```
sample,metric,value
SampleA,efficiency,0.055
SampleA,drop_out_rate,0.192
SampleB,efficiency,0.043
```

### TSV example

```
Sample	Metric	Value
SampleA	Efficiency	0.055
SampleA	drop_out_rate	0.192
```

### Supported metrics

| Metric | Description |
|---|---|
| `efficiency` | Duplex efficiency |
| `drop_out_rate` | Read drop-out rate |
| `on_target_rate_raw` | Fraction of reads on target before duplex filtering (optional) |
| `on_target_rate_duplex` | Fraction of reads on target after duplex filtering (optional) |
| `gc_single` | GC content (single-strand families) |
| `gc_both` | GC content (duplex families) |
| `gc_deviation` | GC deviation between strands |
| `total_families` | Total number of families |
| `family_mean` | Mean family size |
| `family_median` | Median family size |
| `family_max` | Maximum family size |
| `families_gt1` | Families with more than one read |
| `single_families` | Single-read families |
| `paired_families` | Paired families |
| `paired_and_gt1` | Paired families with more than one read |
| `frac_singletons` | Fraction of singleton reads |

Any additional metrics present in the input files will be plotted as individual bar charts.

---

## Usage

```bash
# Basic usage — MultiQC will discover matching files automatically
multiqc path/to/input/

# Force overwrite an existing report
multiqc --force path/to/input/

# Verbose/debug output (useful for diagnosing file detection issues)
multiqc -v path/to/input/
```

The report is written to `multiqc_report.html` in the current directory.
Additional data files are saved in `multiqc_data/`.

---

## Development

### Running tests

```bash
# With uv
uv run pytest tests/

# With pip-installed environment
pytest tests/
```

### Project structure

```
duplex_multiqc/
  custom_code.py              # MultiQC hook — registers search patterns
  modules/
    duplex_seq/
      duplex_seq.py           # Main module: parsing and plot generation
      __init__.py
tests/
  test.csv                    # Sample test fixture (CSV, lowercase headers)
  sampledata.tsv              # Sample test fixture (TSV, uppercase headers)
  test_plugin.py              # pytest test suite
pyproject.toml
```

---

## Contributing

1. Fork the repository and create a feature branch.
2. Install the dev environment: `uv sync --group dev`
3. Make your changes and add tests in `tests/test_plugin.py`.
4. Run `uv run pytest tests/` and ensure all tests pass.
5. Open a pull request against `main`.
