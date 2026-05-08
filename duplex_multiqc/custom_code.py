#!/usr/bin/env python
"""
This file defines hooks and utility functions for the MultiQC plugin.
"""

import logging
import importlib.metadata
from multiqc.utils import config

log = logging.getLogger("multiqc")


def execute():
    """
    Code to execute after the config files have been parsed.
    """
    log.info("Running execute() from duplex_multiqc plugin")

    version = importlib.metadata.version("duplex_multiqc")
    log.info("Running MultiQC Plugin v{}".format(version))

    search_patterns = {
        "duplex_seq": [
            {
                "fn": "*.csv",
                "contents": "sample,metric,value",
                "num_lines": 1,
            },
            {
                "fn": "*.tsv",
                "contents": "sample\tmetric\tvalue",
                "num_lines": 1,
            },
        ]
    }

    for pattern_name, pattern in search_patterns.items():
        if pattern_name not in config.sp:
            config.update_dict(config.sp, {pattern_name: pattern})
            log.debug("Added {} to the search patterns".format(pattern_name))
        else:
            log.debug(
                "Not adding {} to the search patterns as it is already set".format(
                    pattern_name
                )
            )

