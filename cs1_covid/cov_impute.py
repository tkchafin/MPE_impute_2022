#!/usr/bin/env python

# Standard library imports
import argparse
import sys

import numpy as np
import pandas as pd
import scipy.stats as stats

import impute

from impute.read_input.read_input import GenotypeData
from impute.read_input.simgenodata import SimGenotypeData
from impute.simple_imputers.simple_imputers import ImputePhylo
from impute.simple_imputers.simple_imputers import ImputeNMF
from impute.simple_imputers.simple_imputers import ImputeAlleleFreq

sys.setrecursionlimit(10**7)

def main():
	data = GenotypeData(
		filename="alignment_0.phylip",
		filetype="phylip",
		siterates_iqtree="alignment_0.phylip.rate",
		qmatrix_iqtree="alignment_0.phylip.iqtree",
		guidetree="alignment_0.phylip.treefile"
	)
	sim=SimGenotypeData(data, prop_missing=0.2, strategy="random")
	imputed=ImputePhylo(genotype_data=sim, prefix="alignment_0.phylip")
	print(sim.accuracy(imputed))
	print(sim.accuracy_by_site(imputed))
	print(sim.missingness_by_site())
	print(sim.missingness_by_site(mask=True))

if __name__ == "__main__":
    main()
