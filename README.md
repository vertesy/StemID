# StemID and RaceID2 algorithms

RaceID2 is an advanced version of RaceID, an algorithm for the identification of rare and abundant cell types from single cell transcriptome data. The method is based on transcript counts obtained with unique molecular identifies.

StemID is an algorithm for the derivation of cell lineage trees based on RaceID2 results and predicts multipotent cell identites.

RaceID2 and StemID are written in the R computing language.

## RaceID2 methods
* **initialize**.
To create a SCseq object we need an input data frame of transcript counts, columns are cells, rows are genes.
  + sc <- SCseq(inputdata)

* **filterdata**.
Filters data. Input parameters and default values are: mintotal=3000 (discards cells with less than mintotal reads), minexpr=5, minnumber=1 (discards genes with less than minexpr transcripts in at least minnumber cells), maxexpr=Inf (discards genes with more than maxexpr transcripts in at least one cell), downsample=TRUE (when downsample is set to TRUE data is downsampled to mindownsample transcripts per cell, otherwise it is only median normalized), dsn=1 (number of downsamplings, ), rseed=17000 (used in case of downsample). <br />
Input parameters are stored in slot sc@filterparameters. 
The script first normalizes transcripts across cells with more than mintotal transcripts and stores the result in slot sc@ndata.
Then removes genes according to given filters and stores resulting data to sc@fdata. 

  + sc <- filterdata(sc, mintotal=minreadspercell, minexpr=5, minnumber=1, maxexpr=maxexprpergene, downsample=dodownsample, dsn=1, rseed=17000)

* **clustexp**. 


## RaceID2 functions

The following files are provided:

StemID/RaceID2 class definition: RaceID2_StemID_class.R 
StemID/RaceID2 sample code: RaceID2_StemID_sample.R
StemID/RaceID2 reference manual: Reference_manual_RaceID2_StemID.pdf
StemID/RaceID2 sample data: transcript_counts_intestine_5days_YFP.xls

