# Mactop

![](https://github.com/ydduanran/Mactop/blob/main/Mactop_Overview.jpg)

## Overview
Mactop is designed for identifying chromatin domains, domain communities, and chromunities in Hi-C data and higher-order interaction data.

Mactop, a Markov clustering-based tool to accurately identify TADs and provide biologically important classifications of TADs and their boundaries. Mactop distinguishes stable and dynamic boundaries based on biological significance. More importantly, leveraging spatial interactions among TADs, Mactop uncovers TAD communities characterized by chromatin accessibility and enriched histone modifications. Mactop unveils the ‘chromunity’ within TADs in high-order interaction data, showing that interactions within TADs are diverse rather than uniform. In short, Mactop is a versatile, accurate, robust tool for deciphering chromatin domain, domain community, and chromunity for 3D genome maps.


## Getting started
See [Documentation and Tutorials](https://stagate.readthedocs.io/en/latest/index.html).

## Software dependencies
scanpy

tensorflow==1.15.0

## Installation
cd STAGATE-main

python setup.py build

python setup.py install

## Citation
Dong, Kangning, and Shihua Zhang. “Deciphering spatial domains from spatially resolved transcriptomics with adaptive graph attention auto-encoder.” bioRxiv (2021). doi: https://doi.org/10.1101/2021.08.21.457240
