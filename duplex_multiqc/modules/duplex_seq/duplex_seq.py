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
_ON_TARGET_RATE_METRICS = {"on_target_rate_raw", "on_target_rate_duplex"}
_ON_TARGET_COVERAGE_METRICS = {"on_target_coverage_raw", "on_target_coverage_duplex"}
_ON_TARGET_RATIO_METRICS = {"on_target_duplex_ratio"}
_SINGLETON_FRAC_METRICS = {"frac_singletons"}
_WITHIN_FAMILY_STATS = {"family_max", "family_median", "family_mean"}
_FAMILY_SIZE_METRICS = {
    "families_gt1",
    "single_families",
    "paired_families",
    "paired_and_gt1",
    "total_families",
}
_GROUPED_METRICS = _GC_METRICS | _ON_TARGET_RATE_METRICS | _ON_TARGET_COVERAGE_METRICS | _WITHIN_FAMILY_STATS | _FAMILY_SIZE_METRICS

# Preferred display order for individual (non-grouped) metrics
_INDIVIDUAL_METRIC_ORDER = ["efficiency", "drop_out_rate"]

_TABLE_COLORS = ["PuBu", "BuPu", "BuGn", "Oranges", "RdPu"]
_METRIC_LABELS = {
    "efficiency": "Efficiency",
    "drop_out_rate": "Drop Out Rate",
    "on_target_rate_raw": "On Target Raw",
    "on_target_rate_duplex": "On Target Duplex",
    "on_target_coverage_raw": "On Target Coverage Raw",
    "on_target_coverage_duplex": "On Target Coverage Duplex",
    "on_target_duplex_ratio": "On Target Duplex Ratio",
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
        # Render known metrics in preferred order first, then any remaining unknowns.
        individual_metrics = {
            k: v for k, v in metrics_dict.items()
            if k not in _GROUPED_METRICS and k not in _SINGLETON_FRAC_METRICS and k not in _ON_TARGET_RATIO_METRICS
        }
        ordered = [m for m in _INDIVIDUAL_METRIC_ORDER if m in individual_metrics]
        ordered += [m for m in individual_metrics if m not in _INDIVIDUAL_METRIC_ORDER]
        for metric in ordered:
            self.add_section(
                name=_metric_label(metric),
                anchor=f"duplex_seq_{metric.lower()}_plot",
                description=f"Plot for metric: {_metric_label(metric)}",
                plot=self.plot_single_bargraph(individual_metrics[metric], metric),
            )

        # --- Grouped on-target rate section ---
        on_target_rate_metrics = {k: v for k, v in metrics_dict.items() if k in _ON_TARGET_RATE_METRICS}
        if on_target_rate_metrics:
            self.add_section(
                name="On-target Rate",
                anchor="duplex_seq_on_target_rate",
                description="On-target rate per sample (raw and duplex)",
                plot=self.plot_grouped_bargraph(on_target_rate_metrics, "On-target Rate"),
            )

        # --- Grouped on-target coverage section ---
        on_target_coverage_metrics = {k: v for k, v in metrics_dict.items() if k in _ON_TARGET_COVERAGE_METRICS}
        if on_target_coverage_metrics:
            self.add_section(
                name="On-target Coverage",
                anchor="duplex_seq_on_target_coverage",
                description="On-target coverage per sample (raw and duplex)",
                plot=self.plot_grouped_bargraph(on_target_coverage_metrics, "On-target Coverage"),
            )

        # --- On-target duplex ratio metric section ---
        on_target_ratio_metrics = {k: v for k, v in metrics_dict.items() if k in _ON_TARGET_RATIO_METRICS}
        for metric, sample_values in on_target_ratio_metrics.items():
            self.add_section(
                name=_metric_label(metric),
                anchor=f"duplex_seq_{metric.lower()}_plot",
                description=f"Plot for metric: {_metric_label(metric)}",
                plot=self.plot_single_bargraph(sample_values, metric),
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

        # --- Within-family stats section ---
        within_family_stats = {k: v for k, v in metrics_dict.items() if k in _WITHIN_FAMILY_STATS}
        if within_family_stats:
            self.add_section(
                name="Within-family stats",
                anchor="duplex_seq_within_family_stats",
                description="Within-family statistics per sample",
                plot=self.plot_family_metrics(within_family_stats, "Within-family stats"),
            )

        # --- Grouped family size section ---
        family_size_metrics = {k: v for k, v in metrics_dict.items() if k in _FAMILY_SIZE_METRICS}
        if family_size_metrics:
            self.add_section(
                name="Family size comparison",
                anchor="duplex_seq_family_size_comparison",
                description="Comparison of family size metrics per sample",
                plot=self.plot_grouped_bargraph(family_size_metrics, "Family size"),
            )

        # --- Singleton fraction metric section ---
        singleton_frac_metrics = {k: v for k, v in metrics_dict.items() if k in _SINGLETON_FRAC_METRICS}
        for metric, sample_values in singleton_frac_metrics.items():
            self.add_section(
                name=_metric_label(metric),
                anchor=f"duplex_seq_{metric.lower()}_plot",
                description=f"Plot for metric: {_metric_label(metric)}",
                plot=self.plot_single_bargraph(sample_values, metric),
            )

    def plot_single_bargraph(self, data_dict, metric):
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
        metric_dict format: {'group1': {'S1': 10}, 'group2': {'S1': 20}}
        """
        plot_data = {}
        for metric_name, samples in metric_dict.items():
            for sample_name, val in samples.items():
                if sample_name not in plot_data:
                    plot_data[sample_name] = {}
                plot_data[sample_name][metric_name] = val

        cats = _build_cats(metric_dict.keys())
        plot_slug = metric_title.lower().replace(" ", "_").replace("-", "_")

        pconfig = {
            "id": f"duplex_seq_{plot_slug}_grouped_plot",
            "title": f"{metric_title}",
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

        priority_order = ["family_max", "family_median", "family_mean"]
        ordered_metrics = [m for m in priority_order if m in metric_dict]
        ordered_metrics.extend([m for m in metric_dict if m not in ordered_metrics])
        cats = _build_cats(ordered_metrics)

        plot_slug = metric_title.lower().replace(" ", "_").replace("-", "_")
        pconfig = {
            "id": f"duplex_seq_{plot_slug}_stacked_bar",
            "title": metric_title,
            "xlab": "Samples",
            "ylab": "Value",
            "stacking": None,
            "tt_decimals": 6,
            "cpswitch": False,
        }
        return bargraph.plot(plot_data, cats, pconfig=pconfig)
