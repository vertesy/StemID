# StemID and RaceID2 algorithms

RaceID2 is an advanced version of RaceID, an algorithm for the identification of rare and abundant cell types from single cell transcriptome data. The method is based on transcript counts obtained with unique molecular identifies.

StemID is an algorithm for the derivation of cell lineage trees based on RaceID2 results and predicts multipotent cell identites.

RaceID2 and StemID are written in the R computing language.

## RaceID2 methods
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

  + sc <- filterdata(sc, mintotal=1000, minexpr=5, minnumber=1, maxexpr=Inf, downsample=FALSE, dsn=1, rseed=17000)
  + sc <- filterdata(sc) -- runs function with default values.

* **clustexp**. Clusters data using kmedoids. <br/>
  Input parameters and default values are: 
  1. _clustnr=20_ ()
  2. _bootnr=50_ ()
  3. _metric="pearson"_ ()
  4. _do.gap=TRUE_ ()
  5. _sat=FALSE_ (incorporated in RaceID2, )
  6. _SE.method="Tibs2001SEmax"_ ()
  7. _SE.factor=.25_ ()
  8. _B.gap=50_ ()
  9. _cln=0_ ()
  10. _rseed=17000_ ()
  11. _FUNcluster="kmeans"_ (incorporated in RaceID2, )
  12. _version = 2_ (version of RaceID) <br />
  
  Clusters sc@fdata <br/>
  Inut parameters are stored in slot sc@clusterpar. <br />
  Run as:

  + sc <- clustexp(sc, clustnr=30, bootnr=50, metric="pearson", do.gap=FALSE, sat=TRUE, SE.method="Tibs2001SEmax", SE.factor=0.25, B.gap=50, cln=0, rseed=17000, FUNcluster="kmedoids", version = 2)
  + sc <- clustexp(sc) -- runs function with default values


## RaceID2 functions
* **downsample**. Downsamples inputdata. <br />
Transcript data is converted to integer data and random sampling is done _dsn_ times and averaged. A peudocount equal to 0.1 is added to the resulting data.frame. 
There are two versions (DG and JCB, written by Dominic Gr"un and Jean-Charles Boisset respectively). By default the functions uses JCB version. To choose another one use _dsversion_ in method _filterdata_.

The following files are provided:

StemID/RaceID2 class definition: RaceID2_StemID_class.R 
StemID/RaceID2 sample code: RaceID2_StemID_sample.R
StemID/RaceID2 reference manual: Reference_manual_RaceID2_StemID.pdf
StemID/RaceID2 sample data: transcript_counts_intestine_5days_YFP.xls

