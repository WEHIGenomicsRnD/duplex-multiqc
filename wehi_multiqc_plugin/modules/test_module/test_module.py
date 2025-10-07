from multiqc.modules.base_module import BaseMultiqcModule
import os
import logging
import csv

from multiqc.utils import config
from multiqc.plots import bargraph

log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super().__init__(
            name="Efficiency E.Coli K12 Analyses",
            anchor="test_module",
            info="..."
        )

        self.data = {}
        table_headers = {}
        plot_data = {}

        # Metric → Plot Data
        # Format: { metric1: { sample1: value1, sample2: value2 }, ... }
        metrics_dict = {}

        matched_files = self.find_log_files('test_module')
        if not matched_files:
            log.info("[TestModule] No matched files found.")
            return

        for f in matched_files:
            filepath = os.path.join(f['root'], f['fn'])
            log.info(f"[TestModule] Processing file: {filepath}")

            try:
                with open(filepath, 'r') as handle:
                    reader = csv.DictReader(handle, delimiter='\t')
                    for row in reader:
                        sample = row.get('Sample')
                        metric = row.get('Metric')
                        value = row.get('Value')

                        if not sample or not metric or value is None:
                            log.warning(f"Missing required field in row: {row}")
                            continue

                        try:
                            value = float(value)
                        except ValueError:
                            log.warning(f"Non-numeric value for sample '{sample}' metric '{metric}': {value}")
                            continue

                        # Add to general stats table
                        if sample not in self.data:
                            self.data[sample] = {}
                        self.data[sample][metric] = value

                        # Prepare plot data by metric
                        if metric not in metrics_dict:
                            metrics_dict[metric] = {}
                        metrics_dict[metric][sample] = value

                    # Store table headers only once
                    if reader.fieldnames:
                        for col in reader.fieldnames:
                            table_headers[col] = {
                                "title": col,
                                "description": f"{col} column"
                            }

            except Exception as e:
                log.error(f"[TestModule] Failed to read {filepath}: {e}")
                continue

        if not self.data:
            log.warning("[TestModule] No data rows found to display.")
            return

        # Add general stats table
        self.general_stats_addcols(self.data, table_headers)

        # Define metric-to-plot function mapping
        metric_plot_mapping = {
            "Efficiency": self.plot_bargraph,
            # Add more mappings here
        }

        # plot section for each metric
        for metric, sample_values in metrics_dict.items():
            plot_func = metric_plot_mapping.get(metric, self.plot_bargraph)  # default = bar
            plot_config = plot_func(sample_values, metric)

            self.add_section(
                name=f"{metric} Plot",
                anchor=f"test_module_{metric.lower()}_plot",
                description=f"Plot for metric: {metric}",
                plot=plot_config
            )
     # Uses MultiQC's built-in bargraph plotting utility
    def plot_bargraph(self, data_dict, metric):
    # Convert {sample: value} dict to a list with one dict as required by MultiQC
      data = [
        {sample: {"count": value} for sample, value in data_dict.items()}
        ]
      pconfig = {
        'id': f'test_module_barplot_{metric.lower()}',
        'title': f'{metric} Bar Plot',
        'xlab': 'Sample',
        'ylab': metric,
        'xmin': 0,
        }
      return bargraph.plot(data, pconfig=pconfig)


        
        

