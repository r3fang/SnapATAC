## SnapATAC (internal testing)
Single Nuclesus Analysis Package for ATAC-seq. 

## Introduction
SnapATAC is fast, accurate and unbiased method for analyzing single cell ATAC-seq datasets. Compared to previous methods, SnapATAC 1) overcomes limitation of reliance on open chromatin peaks defined by aggregate/bulk signal; 2) reveals novel cis-elements active in rare populations; 3) adjusts for differing sequencing depth between cells; 4) scales up to millions of cells.

SnapATAC has two components: [Snaptools](https://github.com/r3fang/SnapTools) and [SnapATAC](https://github.com/r3fang/SnapATAC). 

* Snaptools is a module for working with [snap](https://github.com/r3fang/SnapATAC/wiki/FAQs) file in Python. 
* SnapATAC is a R package for the downstream analysis. 

## Requirements  
* Python (2.7)
* R (>= 3.4.0)

## Install SnapTools
Install snaptools from PyPI

```bash
$ pip install snaptools==1.2.5 --user
```

See how to install snaptools on [FAQs](https://github.com/r3fang/SnapATAC/wiki/FAQs). 

## Install SnapATAC

```
$ R
> install.packages("devtools")
> library(devtools)
> install_github("r3fang/SnapATAC");
> library(SnapATAC);
```

## Get Started (see more details in MOs 2k in tutorials)

```R
> library(SnapATAC);
> data(mos);
> mos
> mos = calJaccard(mos);
> mos = runPCA(mos, pc.num=20);
> mos = runCluster(mos, k=15, resolution=0.5);
> mos = runViz(mos, dims=2, method="Rtsne");
```

## Galleries & Tutorials (click on the image for details)
[<img src="./images/MOS_2k.png" width="275" height="315" />](./examples/MOS_2k/MOS_2k.md)
[<img src="./images/Fang_2019.png" width="275" height="315" />](./examples/Fang_2019/Fang_2019.md)
[<img src="./images/10X_2018.png" width="275" height="315" />](./examples/10X_2018/10X_2018.md)
[<img src="./images/Cusanovich_2018.png" width="275" height="315" />](./examples/Cusanovich_2018/Cusanovich_2018.md)
[<img src="./images/Lake_2018.png" width="275" height="315" />](./examples/Lake_2018/Lake_2018.md)
[<img src="./images/Schep_2017.png" width="275" height="315" />](./examples/Schep_2017/Schep_2017.md)
[<img src="./images/Habib_2017.png" width="275" height="315" />](./examples/Habib_2017/Habib_2017.md)

