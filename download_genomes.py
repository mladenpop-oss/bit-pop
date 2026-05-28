#!/usr/bin/env python3
"""Download 70 bacterial genomes from NCBI for benchmark."""

import os
import time
import argparse
from Bio import Entrez, SeqIO

Entrez.email = "mladen.popovic@gmail.com"

GENOMES = [
    # E. coli (5)
    ("NC_000913.3", "EC_K12"),
    ("NC_002695.2", "EC_O157"),
    ("NC_004431.1", "EC_CFT073"),
    ("NC_007950.1", "EC_042"),
    ("NC_011108.1", "EC_IHD180"),
    # S. aureus (4)
    ("NC_002741.2", "SA_N315"),
    ("NC_003745.2", "SA MU5"),
    ("NC_007780.1", "SA_E
