# StemID and RaceID2 algorithms

RaceID2 is an advanced version of RaceID, an algorithm for the identification of rare and abundant cell types from single cell transcriptome data. The method is based on transcript counts obtained with unique molecular identifies.

StemID is an algorithm for the derivation of cell lineage trees based on RaceID2 results and predicts multipotent cell identites.

RaceID2 and StemID are written in the R computing language.

## Methods
* **initialize**. Creates a SCseq object. <br />
 As input we need data frame of transcript counts, columns are cells, rows are genes. Run as:
  + sc <- SCseq(inputdata)

* **filterdata**.
Filters data. <br /> 
Input parameters and default values are: 
  1. _mintotal=1000_ (discards cells with less than mintotal reads)
  2. _minexpr=5, minnumber=1_ (discards genes with less than minexpr transcripts in at least minnumber cells)
  3. _maxexpr=Inf_ (discards genes with more than maxexpr transcripts in at least one cell)
  4. _downsample=FALSE_ (logical; when TRUE data is downsampled to mintotal transcripts per cell, otherwise it is median normalized)
  5. _dsn=1_ (number of downsamplings; output is an average over dsn downsamplings)
  6. _rseed=17000_ (seed used for downsampling)
  7. _dsversion="JCB"_ (downsampling function version) <br />
  
  Input parameters are stored in slot sc@filterparameters. 
  The method first median normalizes or downsamples (dependeing of _downsample_) transcripts across cells with more than _mintotal_ transcripts and stores the result in slot sc@ndata.
  Then removes genes according to _minexpr_, _minnumber_ and _maxexpr_ and stores resulting data.frame into sc@fdata. 

  + sc <- filterdata(sc, mintotal=1000, minexpr=5, minnumber=1, maxexpr=Inf, downsample=FALSE, dsn=1, rseed=17000, dsversion = 'JCB')
  + sc <- filterdata(sc) -- runs function with default values.

* **clustexp**. Clusters data using kmedoids. <br/>
  Input parameters and default values are: 
  1. _clustnr=20_ (Number of clusters. Must be greater than 1.)
  2. _bootnr=50_ (Maximum number of clusters for the computation of the gap statistics or the derivation of the cluster number by saturation criterion.)
  3. _metric="pearson"_ (Metric to compute distance between cells. Options are:  "spearman","pearson","kendall","euclidean","maximum","manhattan","canberra","binary","minkowski". Check function dist.gen for more information. Distances are stored in sc@distances.)
  4. _do.gap=TRUE_ (If set to TRUE, the number of clusters is determined using gap statistics. Default is TRUE.)
  5. _sat=FALSE_ (incorporated in RaceID2, computes the number of clusters using saturation criterion.)
  6. _SE.method="Tibs2001SEmax"_ ()
  7. _SE.factor=.25_ ()
  8. _B.gap=50_ (Number of bootstrap runs for the gap statistics.)
  9. _cln=0_ (Number of clusters for clustering. In case it is 0, will be determined by either gap statistics of saturation criterion.)
  10. _rseed=17000_ (Seed for random number generator used in case of gap statistics and for posterior clustering.)
  11. _FUNcluster="kmeans"_ (incorporated in RaceID2, this can be kmeans, hclust or kmedoids. ) <br />
  
 Input parameters are stored in slot sc@clusterpar. Default is taken when no specified. <br/>
 Data in sc@fdata in clustered using clustfun functions. Results are stored in sc@cluster, sc@distances and sc@fcol. Go to _Functions_ section to learn more. <br/>
  Run as:

  + sc <- clustexp(sc, clustnr=20, bootnr=50, metric="pearson", do.gap=FALSE, sat=TRUE, SE.method="Tibs2001SEmax", SE.factor=0.25, B.gap=50, cln=0, rseed=17000, FUNcluster="kmedoids")
  + sc <- clustexp(sc) -- runs function with default values

## Functions
* **downsample**. Downsamples inputdata. <br />
Transcript data is converted to integer data and random sampling is done _dsn_ times and averaged. A peudocount equal to 0.1 is added to the resulting data.frame. 
There are two versions (DG and JCB, written by Dominic Gr"un and Jean-Charles Boisset respectively). By default the functions uses JCB version. To choose another one use _dsversion_ in method _filterdata_.

* **clustfun**. Clusters sc@fdata. <br/>
Version 2, from RaceID2. Computes distance between cells (using dist.gen function) using specified metric. Determines cluster number if required using gap statistics or saturation criterion. Then clusters data (using clusGapExt function) using the specified method -kmedoids, kmeans or hclust-. <br />

* **dist.gen**. Distance between cells. <br/>
Computes and returns the distance matrix computed by using the specified distance (mmetric) measure to compute the distances between the cells. In case of metric "spearman", "pearson", or "kendall", the function takes 1 - correlation as a distance, and takes the direct measurement of the distance for metric "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".

* **clusGapExt**. Clustering.

* **clusterboot kmeansCBI**. kmeans clustering. fpc package.
* **clusterboot with pamkCBI**. kmedoids clustering. fpc package.
* **clusterboot with hclusterCBI**. hclust clustering. fpc package.


The following files are provided:

StemID/RaceID2 class definition: RaceID2_StemID_class.R 
StemID/RaceID2 sample code: RaceID2_StemID_sample.R
StemID/RaceID2 reference manual: Reference_manual_RaceID2_StemID.pdf
StemID/RaceID2 sample data: transcript_counts_intestine_5days_YFP.xls

