# STARMAPs
STARMAPs is a statistical tool for a comparison between two microbiota datasets. It provides inference for alike and opposite matchings of the microbial changes associated with treatments/factors in two different datasets. 

# Implementation
STARMAPs is currently implemented as an R source file ("STARMAPs_v1.R"). We are currently working to develop it into an R package.

# Usage
In general, use the following for two datasets with taxon-by-sample count matrices data1 and data2 with meta data files meta1 and meta2 with matching sample ids. Names of the columns in the meta files that specify the groups to compare need to be provided to fcol1 and fcol2.
~~~
source("STARMAPs_v1.R")
res <- starmaps(data1, meta1, fcol1, data2, meta2, fcol2)
plot.starmaps(res, "res.pdf")
tbl <- summary.starmaps(res)
~~~

Two additional examples are provided in the file "Usage_Examples.R".

Note that, by default, STARMAPs expects taxon names that uses ";" to seperate each taxonomic level beginning with Phylum (see get.taxa() function in "STARMAPs_v1.R").

# Benchmarking
Codes and results of benchmarking evaluation of STARMAPs, as described in Jiang et al (citation below), is provided in the "Benchmarking_Evaluations.zip" file.

# Contact: 
Peng Jiang (peng.jiang@northwestern.edu)

# Citation: 
Jiang P, Green SJ, Chlipala GE, Turek FW, Vitaterna MH (submitted) Reproducible changes in the gut microbiome reveal a shift in microbial and host metabolism during spaceflight.

# Required packages
- MCMCpack
- vegan
- compositions
- ggplot2
- gridExtra
- parallel

Benchmarking also used
- phyloseq
- plotROC
