## SnapATAC
Single Nuclesus Analysis Package for ATAC-seq. 

## Introduction
SnapATAC is fast, accurate and unbiased method for analyzing single cell ATAC-seq datasets. Compared to previous methods, SnapATAC 1) overcomes limitation of reliance on open chromatin peaks defined by aggregate/bulk signal; 2) reveals novel cis-elements active in rare populations; 3) adjusts for differing sequencing depth between cells; 4) scales up to millions of cells.

SnapATAC has two components: [Snaptools](https://github.com/r3fang/SnapTools) and [SnapATAC](https://github.com/r3fang/SnapATAC). 

* Snaptools is a module for working with [snap](https://github.com/r3fang/SnapATAC/wiki/What-is-a-snap-file%3F) file in Python. 
* SnapATAC is a R package for the downstream analysis. 

## Requirements  
* Python (2.7)
* R (>= 3.4.0)

## Install SnapTools
Install snaptools from PyPI

```bash
$ pip install snaptools==1.2.5 --user
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

## Galleries & Tutorials
[<img src="./images/Fang_2019.png" width="250" height="290" />](./examples/Fang_2019/Fang_2019.md)
[<img src="./images/Cusanovich_2018.png" width="250" height="290" />](./examples/Cusanovich_2018/Cusanovich_2018.md)
[<img src="./images/Lake_2018.png" width="250" height="290" />](./examples/Lake_2018/Lake_2018.md)
[<img src="./images/Schep_2017.png" width="250" height="290" />](./examples/Schep_2017/Schep_2017.md)
[<img src="./images/Habib_2017.png" width="250" height="290" />](./examples/Habib_2017/Habib_2017.md)
[<img src="./images/10X_2018.png" width="250" height="290" />](./examples/10X_2018/10X_2018.md)

