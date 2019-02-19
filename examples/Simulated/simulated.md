## Analysis of simulated scATAC-seq

**Step 1. Download Snap files of simulated scATAC-seq datasets**. 

```
$ wget http://renlab.sdsc.edu/r3fang/share/Fang_2019/simulated_scATAC/500.snap
$ wget http://renlab.sdsc.edu/r3fang/share/Fang_2019/simulated_scATAC/1000.snap
$ wget http://renlab.sdsc.edu/r3fang/share/Fang_2019/simulated_scATAC/2500.snap
$ wget http://renlab.sdsc.edu/r3fang/share/Fang_2019/simulated_scATAC/5000.snap
$ wget http://renlab.sdsc.edu/r3fang/share/Fang_2019/simulated_scATAC/10000.snap
$ wget http://renlab.sdsc.edu/r3fang/share/Fang_2019/simulated_scATAC/10000.offpeak.snap
```

**Step 2. Create a snap object**

```R
$ R
> library(SnapATAC);
> x.sp = createSnap("10000.offpeak.snap");
> x.sp = addBmat(x.sp, "10000.offpeak.snap", binSize=5000);
> x.sp

number of barcodes: 5500
number of bins: 545114
number of peaks: 0
number of genes: 0
==========================
meta data            (metData) :  TRUE
cellxbin matrix      (bmat)    :  TRUE
cellxpeak matrix     (pmat)    :  FALSE
cellxgene matrix     (gmat)    :  FALSE
jaccard matrix       (jmat)    :  FALSE
normalization        (nmat)    :  FALSE
PCA:                 (smat)    :  FALSE
cluster:             (cluster) :  FALSE
t-sne:               (tsne)    :  FALSE
umap:                (umap)    :  FALSE

```

**Step 3. Randomly sample 2,000 cells as described in the paper**

```R
> set.seed(1)
> idx = sort(sample(seq(nrow(x.sp)), 2000));
> x.sp = x.sp[idx,]
```

**Step 3. Convert cell-by-bin count matrix to binary matrix**

```R
> x.sp = makeBinary(x.sp);
```

**Step 4. Jaccard matrix & normalization**

```R
> x.sp = calJaccard(x.sp, mat = "bmat");
```

**Step 5. PCA**

```R
> x.sp = runPCA(
	x.sp,
	pc.num=30,
	input.mat="nmat",
	method="svd",
	weight.by.sd = TRUE,
	center=TRUE,
	scale=FALSE,
	seed.use=10
	);
```

**Step 6. Cluster**

```R
> x.sp = runCluster(
	x.sp,
	pca_dims=1:10,
	k=30,
	resolution=1.0,
	method="jaccard_louvain",
	path_to_louvain="../../../github/snapR/bin/findCommunityLouvain"
	);
```

**Step 7. visulization**

```
# visulization using tsne and umap
> x.sp = runViz(
	x.sp, 
	pca_dims=1:10, 
	dims=2, 
	method="umap",
	);
> x.sp = runViz(
	x.sp, 
	pca_dims=1:10, 
	dims=2, 
	method="Rtsne",
	init_dims=x.sp@tsne
	);

# visulization
> plotViz(x.sp, method="tsne", pch=19, cex=0.5);
> plotViz(x.sp, method="umap", pch=19, cex=0.5);

```

<img src="./Viz_tsne.PNG" width="300" height="300" /> <img src="./Viz_umap.PNG" width="300" height="300" />


**Step 8. Comparision with cell label**

```R
# compare with cell label
> Y = factor(do.call(c, lapply(strsplit(x.sp@barcode, "[.]"), function(x) x[1])));
> compare(Y, x.sp@cluster, method="nmi");
[1] 0.998709
```


