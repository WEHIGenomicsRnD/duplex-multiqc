#!/usr/bin/env python
"""
This file defines hooks and utility functions for the MultiQC plugin.
"""

from __future__ import print_function
import logging
from multiqc.utils import config, report
import importlib_metadata

# Initialising main MultiQC logger
log = logging.getLogger('multiqc')

def execute():
    """
    Code to execute after the config files have been parsed.
    """
    log.info("Running execute() from custom MultiQC plugin")

    # Plugin's version number defined in pyproject.toml:
    version = importlib_metadata.version("duplex_multiqc")
    log.info("Running MultiQC Plugin v{}".format(version))
    search_patterns = {
            'test_module': {
                'fn': '*.csv',
                'contents' : "sample,metric,value",
                'num_lines': 1
            }
        }
    # Add a custom search pattern for 'test_module' if it doesn't already exist
    # Add to the search patterns used by modules
    for pattern_name, pattern in search_patterns.items():
        if pattern_name not in config.sp:
            config.update_dict( config.sp, { pattern_name: pattern } )
            log.debug("Added {} to the search patterns".format(pattern_name))
        else:
            log.debug("Not adding {} to the search patterns as it is already set".format(pattern_name))
    
    # Some additional filename cleaning
    config.fn_clean_exts.extend([".my_tool_extension", ".removeMetoo"])
    
