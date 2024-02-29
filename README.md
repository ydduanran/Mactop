# Mactop
Mactop is designed for identifying chromatin domains, domain communities, and chromunities in Hi-C data and higher-order interaction data.

## Overview
Mactop, a Markov clustering-based tool to accurately identify TADs and provide biologically important classifications of TADs and their boundaries. Mactop distinguishes stable and dynamic boundaries based on biological significance. More importantly, leveraging spatial interactions among TADs, Mactop uncovers TAD communities characterized by chromatin accessibility and enriched histone modifications. Mactop unveils the ‘chromunity’ within TADs in high-order interaction data, showing that interactions within TADs are diverse rather than uniform. In short, Mactop is a versatile, accurate, robust tool for deciphering chromatin domain, domain community, and chromunity for 3D genome maps.
![](https://github.com/ydduanran/Mactop/blob/main/Mactop_Overview.jpg)


## Getting started
See [Documentation and Tutorials](https://stagate.readthedocs.io/en/latest/index.html).


## Getting start

### Installation
It's recommended to create a conda environment:

```shell
conda create -n mactop python=3.6
conda activate mactop
```

### dependencies
```shell
cooler >= 0.8.11
numpy  >= 1.17.3 
pandas >= 0.24.1 
scipy  >= 1.1.0 
iced   >= 0.5.10
networkx >= 2.5.1
markov_clustering >= 0.0.6
pickle >= 0.7.5
```
