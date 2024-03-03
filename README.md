# Mactop
Mactop is designed for identifying chromatin domains, domain communities, and chromunities in Hi-C data and higher-order interaction data.

![](https://github.com/ydduanran/Mactop/blob/main/Mactop_Overview.jpg)

## Overview
Mactop, a Markov clustering-based tool to accurately identify TADs and provide biologically important classifications of TADs and their boundaries. Mactop distinguishes stable and dynamic boundaries based on biological significance. More importantly, leveraging spatial interactions among TADs, Mactop uncovers TAD communities characterized by chromatin accessibility and enriched histone modifications. Mactop unveils the ‘chromunity’ within TADs in high-order interaction data, showing that interactions within TADs are diverse rather than uniform. In short, Mactop is a versatile, accurate, robust tool for deciphering chromatin domain, domain community, and chromunity for 3D genome maps.

## Getting start

### Installation
It's recommended to create a conda environment:

```shell
conda create -n mactop python=3.6
conda activate mactop
```

### Installation dependency
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
### Installation Mactop

Install from PyPI:

```shell
pip install mactop3D
```

install from source

```shell
git clone https://github.com/ydduanran/Mactop.git
cd mactop
python setup.py build
python setup.py install
```
## Documentation and Tutorials
You can download the mcool and mutiway data used in our tutorial at [NCBI(GSE149117)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149117).
### For Command Line Interface(CLI) user
The parameter file for mactop needs to be prepared according to the file in [the example](./Tutorial/CLI/paramaters.txt).
```shell
cd mactop
python mactop_CLI.py [path/to/you/paramaters.txt]
```

### For Interactive user
See [Mactop_Tutorials.ipynb](./Tutorial/Mactop_Tutorials.ipynb).


### Support

If you are having issues, please let us know. We have a mailing list located at:

* duanran@mail.ynu.edu.cn
