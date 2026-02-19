from multiqc.modules.base_module import BaseMultiqcModule
import os
import logging
import csv
from multiqc.plots import bargraph
from collections import OrderedDict

log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super().__init__(
            name="Duplex Sequencing Metrics",
            anchor="duplex_metrics",
            info="MultiQC module for plotting duplex sequencing QC metrics"
        )

        self.data = {}
        table_headers = {}

        # Metric → Plot Data
        # Format: { metric1: { sample1: value1, sample2: value2 }, ... }
        metrics_dict = {}

        matched_csv_files = self.find_log_files('test_module')

        if not matched_csv_files:
            log.info("[TestModule] No matched files found.")
            return

        for f in matched_csv_files:
            
            filepath = os.path.join(f['root'], f['fn'])
            log.warning(f"[TestModule] Processing file: {f['fn']}")
            log.info(f"[TestModule] Processing file: {filepath}")

            try:
                with open(filepath, 'r') as handle:
                    reader = csv.DictReader(handle, delimiter=',')
                    for row in reader:
                        sample = row.get('sample')
                        metric = row.get('metric')
                        value = row.get('value')

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
                    colors = ["PuBu", "BuPu", "BuGn", "Oranges", "RdPu"]
                    if metrics_dict:
                        for i, metric in enumerate(metrics_dict):
                            table_headers[metric] = {"title" : metric, 
                                                  "description" : f"{metric} column",
                                                  "format": "{:.5f}",
                                                  "scale":  colors[i % len(colors)]
                                                  }

            except Exception as e:
                log.error(f"[TestModule] Failed to read {filepath}: {e}")
                continue

        if not self.data:
            log.warning("[TestModule] No data rows found to display.")
            return
        
        # Add general stats table
        self.general_stats_addcols(self.data, headers=table_headers)
        
        grouped_keys = []

        gc_metrics = {}
        for k, v in metrics_dict.items():
            if k in ["gc_single", "gc_both"]:
               grouped_keys.append(k)
               gc_metrics[k]= v
        

        family_metrics = {}

        for k, v in metrics_dict.items():
            if k in ["family_max", "family_median","family_mean","single_families"]:
                family_metrics[k] = v
                grouped_keys.append(k)
        

        family_size_metrics = {}
        for k, v in metrics_dict.items():
            if k in ["families_gt1", "single_families","paired_families","paired_and_gt1"]:
                family_size_metrics[k]=v
                grouped_keys.append(k)
        


        # Define metric-to-plot function mapping
        metric_plot_mapping = {
            "Efficiency": self.plot_bargraph,
            "drop_out_rate" : self.plot_bargraph,
            "frac_singletons": self.plot_bargraph,
            "gc_single": self.plot_bargraph,
            "gc_both": self.plot_bargraph,
            "gc_deviation": self.plot_bargraph,
            "total_families": self.plot_bargraph,
            "family_mean" : self.plot_bargraph,
            "family_median" : self.plot_bargraph,
            "family_max" : self.plot_bargraph,
            "families_gt1" : self.plot_bargraph,
            "single_families" : self.plot_bargraph,
            "paired_families": self.plot_bargraph,
            "paired_and_gt1": self.plot_bargraph,
            
            # Add more mappings here
        }


        # plot section for each metric
        for metric, sample_values in metrics_dict.items():
            if metric not in grouped_keys:
                plot_func = metric_plot_mapping.get(metric, self.plot_bargraph)  # default = bar
                plot_config = plot_func(sample_values, metric)

                self.add_section(
                    name=f"{metric}",
                    anchor=f"test_module_{metric.lower()}_plot",
                    description=f"Plot for metric: {metric}",
                    plot=plot_config
                )

        # plot gc section 
        if gc_metrics:
            self.add_section(
                name="GC Metrics",
                anchor="test_module_gc_metrics",
                description="Grouped GC content metrics per sample",
                plot=self.plot_grouped_bargraph(gc_metrics, "GC" )
            )

        if family_metrics:
            self.add_section(
                name="family Metrics",
                anchor="test_module_family_metrics",
                description=" family metrics per sample",
                plot=self.plot_family_metrics(family_metrics, "family" )
            )

        if family_size_metrics:
            self.add_section(
                name="family size Metrics",
                anchor="test_module_gc_metrics",
                description="Grouped size content metrics per sample",
                plot=self.plot_grouped_bargraph(family_size_metrics, "size" )
            )


    def plot_bargraph(self, data_dict, metric):
      
   
      data = [
        {sample: { metric: value} for sample, value in data_dict.items()}
        ]
      
      cats = OrderedDict()
    
      cats[metric] = {
            "name": metric,
            "color": "#4d4d4d" 
        }
      pconfig = {
        'id': f'test_module_barplot_{metric.lower()}',
        'title': f'{metric} Bar Plot',
        'xlab': 'Sample',
        'ylab': metric,
        'xmin': 0,
        "tt_decimals": 6,
        "cpswitch": False,
        }
      
      return bargraph.plot(data,cats, pconfig=pconfig)
    

    def plot_grouped_bargraph(self, metric_dict, metric_title):
        """
        metric_dict format: {'gc_single': {'S1': 10}, 'gc_both': {'S1': 20}}
        """
        plot_data = {}
        
        # 1. Re-structure the data
        for metric_name, samples in metric_dict.items():
            for sample_name, val in samples.items():
                if sample_name not in plot_data:
                    plot_data[sample_name] = {}
                
                plot_data[sample_name][metric_name] = val

        # 2. Configure the plot
        pconfig = {
            'id': f'test_module_{metric_title.lower()}_grouped_plot',
            'title': f'{metric_title} Comparison',
            'ylab': 'Value',
            "tt_decimals": 6,
            'stacking': None,   
            'cpswitch': False,   \
        }

        # Pass the re-structured dict directly to the bargraph utility
        return bargraph.plot(plot_data, pconfig=pconfig)

    
    def plot_family_metrics(self, metric_dict, metric_title):
            """
            Creates a group bar graph starting with single_families at the bottom.
            """
            plot_data = {}
            
            # 1. Re-structure the data for bargraph.plot
            for metric_name, samples in metric_dict.items():
                for sample_name, val in samples.items():
                    if sample_name not in plot_data:
                        plot_data[sample_name] = {}
                    plot_data[sample_name][metric_name] = val

            # 2. Define category colors with explicit order
            cats = OrderedDict()
            
            # Define the priority order: single_families first
            priority_order = ["single_families","family_max","family_mean", "family_median" ]
            
            metric_colors = {
                "single_families": "#f7a35c", # Orange
                "family_mean": "#7cb5ec",      # Blue
                "family_median": "#434348",    # Dark grey
                "family_max": "#90ed7d"        # Green
            }
            
            # Add priority metrics first to set the stack order
            for metric_name in priority_order:
                if metric_name in metric_dict:
                    cats[metric_name] = {
                        "name": metric_name,
                        "color": metric_colors.get(metric_name, "#4d4d4d")
                    }
            
            # Add any remaining metrics that weren't in the priority list
            for metric_name in metric_dict.keys():
                if metric_name not in cats:
                    cats[metric_name] = {
                        "name": metric_name,
                        "color": "#4d4d4d"
                    }

            # 3. Configure the plot
            pconfig = {
                "id": f"test_module_{metric_title.lower()}_stacked_bar",
                "title": f"{metric_title} Metrics Stacked",
                "xlab": "Samples",
                "ylab": "Family Size",
                "stacking": None,
                "tt_decimals": 6,
                "cpswitch": True,
                "logswitch": True,
                "logswitch_active": True, # Loads with Log10 scale active
            }

            return bargraph.plot(plot_data, cats, pconfig=pconfig)