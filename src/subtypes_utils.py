#!/usr/bin/env python3
"""
Utility functions for the subtypes project
"""

def gene_to_filename(gene):
    """
    Convert a formal gene name to be suitable for file names
    """
    return gene.lower().replace(' ', '')