## SnapATAC
Single Nuclesus Analysis Package for ATAC-seq. 

## Introduction
SnapATAC is fast, accurate and truely unbiased method for analyzing single cell ATAC-seq datasets. Compared to previous methods, SnapATAC 1) overcomes the bias introduced by reliance on accessibility peaks defined by aggragte signals; 2) as a results, it reveals novel cis-elements only acitive in rare populations; 3) adjusts for differing sequencing depth effect between cells; 4) scales up to millions of cells using a novel dimensionality reduction method. 

SnapATAC is composed of two components: [Snaptools](https://github.com/r3fang/SnapTools) and [SnapATAC](https://github.com/r3fang/SnapATAC). 

* Snaptools is a module for working with snap files in Python. 
* SnapATAC is a R package for the downstream analysis. 

Together, SnapATAC represents an end-to-end solution of single cell ATAC-seq analysis.

## Requirements 
* Python (2.7)
* R (>= 3.4.0)

## Install SnapTools
Install snaptools from PyPI

```bash
$ pip install snaptools==1.2.3 --user
```

Alteratively, install snaptools from source code

```bash
$ git clone https://github.com/r3fang/snaptools.git
$ cd snaptools
$ python setup.py install --user
$ ./bin/snaptools

usage: snaptools [-h]  ...

Program: snaptools (A module for working with snap files in Python)
Version: 1.2.3
Contact: Rongxin Fang
E-mail:  r4fang@gmail.com

optional arguments:
  -h, --help        show this help message and exit

functions:

    index-genome    Index reference genome.
    align-paired-end
                    Align paired-end reads.
    align-single-end
                    Align single-end reads.
    snap-pre        Create a snap file from bam or bed file.
    snap-add-bmat   Add cell x bin count matrix to snap file.
    snap-add-pmat   Add cell x peak count matrix to snap file.
    snap-add-gmat   Add cell x gene count matrix to snap file.
    dump-fragment   Dump fragments of selected barcodes from a snap file.
    dump-barcode    Dump barcodes from a snap file.
    call-peak       Call peak using selected barcodes.
    louvain         Louvain communities finding.
```

See how to install snaptools on [MAC OS](https://github.com/r3fang/SnapATAC/wiki/SnapTools-Installation). 

## Install SnapATAC

```
$ R
> install.packages("devtools")
> library(devtools)
> install_github("r3fang/SnapATAC");
> library(SnapATAC);
```

## Galleries 
[<img src="./images/Fang_2019.png" width="250" height="290" />](./examples/Fang_2019/Fang_2019.md)
[<img src="./images/Cusanovich_2018.png" width="250" height="290" />](./examples/Cusanovich_2018/Cusanovich_2018.md)
[<img src="./images/Lake_2018.png" width="250" height="290" />](./examples/Lake_2018/Lake_2018.md)
[<img src="./images/Schep_2017.png" width="250" height="290" />](./examples/Schep_2017/Schep_2017.md)
[<img src="./images/Habib_2017.png" width="250" height="290" />](./examples/Habib_2017/Habib_2017.md)
[<img src="./images/Simulated_2019.png" width="250" height="290" />](./examples/Simulated_2019/Simulated_2019.md)
