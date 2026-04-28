from multiqc.modules.base_module import BaseMultiqcModule
import logging
import csv
import io
from multiqc.plots import bargraph
from collections import OrderedDict
from matplotlib import colormaps, colors as mcolors

log = logging.getLogger(__name__)

# Metrics that are grouped into combined plots rather than individual plots
_GC_METRICS = {"gc_single", "gc_both", "gc_deviation"}
_FAMILY_METRICS = {"family_max", "family_median", "family_mean", "single_families", "total_families"}
_FAMILY_SIZE_METRICS = {"families_gt1", "single_families", "paired_families", "paired_and_gt1"}
_GROUPED_METRICS = _GC_METRICS | _FAMILY_METRICS | _FAMILY_SIZE_METRICS

_TABLE_COLORS = ["PuBu", "BuPu", "BuGn", "Oranges", "RdPu"]
_METRIC_LABELS = {
    "efficiency": "Efficiency",
    "drop_out_rate": "Drop Out Rate",
    "frac_singletons": "Frac Singletons",
    "gc_single": "GC Single",
    "gc_both": "GC Both",
    "gc_deviation": "GC Deviation",
    "total_families": "Total Families",
    "family_mean": "Family Mean",
    "family_median": "Family Median",
    "family_max": "Family Max",
    "families_gt1": "Families Gt1",
    "single_families": "Single Families",
    "paired_families": "Paired Families",
    "paired_and_gt1": "Paired And Gt1",
}


def _metric_label(metric_name):
    return _METRIC_LABELS.get(metric_name.lower(), metric_name.replace("_", " ").title())


def _build_cats(metric_names):
    """
    Build MultiQC category metadata using a ColorBrewer-style palette.
    """
    metric_names = list(metric_names)
    if not metric_names:
        return OrderedDict()
    if len(metric_names) == 1:
        metric_name = metric_names[0]
        return OrderedDict(
            {
                metric_name: {
                    "name": _metric_label(metric_name),
                    "color": "#4d4d4d",
                }
            }
        )

    cmap_name = "Set2" if len(metric_names) <= 8 else "tab20"
    cmap = colormaps[cmap_name].resampled(len(metric_names))
    cats = OrderedDict()
    for idx, metric_name in enumerate(metric_names):
        cats[metric_name] = {
            "name": _metric_label(metric_name),
            "color": mcolors.to_hex(cmap(idx)),
        }
    return cats


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super().__init__(
            name="Duplex Sequencing Metrics",
            anchor="duplex_metrics",
            info="MultiQC module for plotting duplex sequencing QC metrics",
        )

        # sample_data[sample][metric] = value
        self.sample_data = {}
        # metrics_dict[metric][sample] = value
        metrics_dict = {}

        matched_files = self.find_log_files("duplex_seq")

        if not matched_files:
            log.info("[duplex_seq] No matched files found.")
            return

        for f in matched_files:
            log.info(f"[duplex_seq] Processing file: {f['fn']}")
            delimiter = "\t" if f["fn"].endswith(".tsv") else ","
            try:
                reader = csv.DictReader(io.StringIO(f["f"]), delimiter=delimiter)
                for row in reader:
                    # Normalise headers to lowercase so both 'Sample' and 'sample' work
                    row = {k.lower(): v for k, v in row.items()}
                    sample = row.get("sample")
                    metric = row.get("metric")
                    value = row.get("value")

                    if not sample or not metric or value is None:
                        log.warning(f"[duplex_seq] Missing required field in row: {row}")
                        continue

                    try:
                        value = float(value)
                    except ValueError:
                        log.warning(
                            f"[duplex_seq] Non-numeric value for sample '{sample}' "
                            f"metric '{metric}': {value}"
                        )
                        continue

                    if sample not in self.sample_data:
                        self.sample_data[sample] = {}
                    self.sample_data[sample][metric] = value

                    if metric not in metrics_dict:
                        metrics_dict[metric] = {}
                    metrics_dict[metric][sample] = value

            except Exception as e:
                log.error(f"[duplex_seq] Failed to read {f['fn']}: {e}")
                continue

        if not self.sample_data:
            log.warning("[duplex_seq] No data rows found to display.")
            return

        # Build general stats table headers once, after all files are parsed
        table_headers = {
            metric: {
                "title": _metric_label(metric),
                "description": f"{_metric_label(metric)} column",
                "format": "{:.5f}",
                "scale": _TABLE_COLORS[i % len(_TABLE_COLORS)],
            }
            for i, metric in enumerate(metrics_dict)
        }
        self.general_stats_addcols(self.sample_data, headers=table_headers)

        # --- Individual metric sections ---
        for metric, sample_values in metrics_dict.items():
            if metric not in _GROUPED_METRICS:
                self.add_section(
                    name=_metric_label(metric),
                    anchor=f"duplex_seq_{metric.lower()}_plot",
                    description=f"Plot for metric: {_metric_label(metric)}",
                    plot=self.plot_bargraph(sample_values, metric),
                )

        # --- Grouped GC section ---
        gc_metrics = {k: v for k, v in metrics_dict.items() if k in _GC_METRICS}
        if gc_metrics:
            self.add_section(
                name="GC Metrics",
                anchor="duplex_seq_gc_metrics",
                description="Grouped GC content metrics per sample",
                plot=self.plot_grouped_bargraph(gc_metrics, "GC"),
            )

        # --- Grouped family section ---
        family_metrics = {k: v for k, v in metrics_dict.items() if k in _FAMILY_METRICS}
        if family_metrics:
            self.add_section(
                name="Family Metrics",
                anchor="duplex_seq_family_metrics",
                description="Family size metrics per sample",
                plot=self.plot_family_metrics(family_metrics, "family"),
            )

        # --- Grouped family size section ---
        family_size_metrics = {k: v for k, v in metrics_dict.items() if k in _FAMILY_SIZE_METRICS}
        if family_size_metrics:
            self.add_section(
                name="Family Size Metrics",
                anchor="duplex_seq_family_size_metrics",
                description="Grouped family size metrics per sample",
                plot=self.plot_grouped_bargraph(family_size_metrics, "size"),
            )

    def plot_bargraph(self, data_dict, metric):
        data = [{sample: {metric: value} for sample, value in data_dict.items()}]
        cats = _build_cats([metric])
        pconfig = {
            "id": f"duplex_seq_barplot_{metric.lower()}",
            "title": f"{_metric_label(metric)} Bar Plot",
            "xlab": "Sample",
            "ylab": _metric_label(metric),
            "xmin": 0,
            "tt_decimals": 6,
            "cpswitch": False,
        }
        return bargraph.plot(data, cats, pconfig=pconfig)

    def plot_grouped_bargraph(self, metric_dict, metric_title):
        """
        metric_dict format: {'gc_single': {'S1': 10}, 'gc_both': {'S1': 20}}
        """
        plot_data = {}
        for metric_name, samples in metric_dict.items():
            for sample_name, val in samples.items():
                if sample_name not in plot_data:
                    plot_data[sample_name] = {}
                plot_data[sample_name][metric_name] = val

        cats = _build_cats(metric_dict.keys())

        pconfig = {
            "id": f"duplex_seq_{metric_title.lower()}_grouped_plot",
            "title": f"{metric_title} Comparison",
            "ylab": "Value",
            "tt_decimals": 6,
            "stacking": None,
            "cpswitch": False,
        }
        return bargraph.plot(plot_data, cats, pconfig=pconfig)

    def plot_family_metrics(self, metric_dict, metric_title):
        """Creates a grouped bar chart for family size metrics."""
        plot_data = {}
        for metric_name, samples in metric_dict.items():
            for sample_name, val in samples.items():
                if sample_name not in plot_data:
                    plot_data[sample_name] = {}
                plot_data[sample_name][metric_name] = val

        priority_order = ["single_families", "family_max", "family_mean", "family_median"]
        ordered_metrics = [m for m in priority_order if m in metric_dict]
        ordered_metrics.extend([m for m in metric_dict if m not in ordered_metrics])
        cats = _build_cats(ordered_metrics)

        pconfig = {
            "id": f"duplex_seq_{metric_title.lower()}_stacked_bar",
            "title": f"{metric_title} Metrics Stacked",
            "xlab": "Samples",
            "ylab": "Family Size",
            "stacking": None,
            "tt_decimals": 6,
            "cpswitch": True,
            "logswitch": True,
            "logswitch_active": True,
        }
        return bargraph.plot(plot_data, cats, pconfig=pconfig)
