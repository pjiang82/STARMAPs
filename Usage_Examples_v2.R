#### STARMAPs Usage examples using	human lung microbiome data from Segal et al. PMIDs: 27572644 and 27855721

dir.create("STARMAPs_example")
setwd("STARMAPs_example")

source("../STARMAPs_v2.R")


#### Example 1: Each dataset involves a factor of two levels (when more than two levels, need to specify the pair of comparison)

#get data
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2242nnn/GSM2242724/suppl/GSM2242724_otu_tableRA.Sarc.Cont.biom.txt.gz", 
              destfile="GSE84608_otu.gz")
counts1 <- read.delim(gzfile("GSE84608_otu.gz"), skip = 1, comment.char = "", check.names = F, row.names = 1)
counts1 <- aggregate(. ~ taxonomy, data = counts1, FUN = sum)
row.names(counts1) <- counts1$taxonomy
counts1 <- counts1[,-1]
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84608/suppl/GSE84608_complete_metadata.txt.gz", 
              destfile="GSE84608_meta.gz")
meta1 <- read.delim(gzfile("GSE84608_meta.gz"), skip = 9, row.names = 1, nrows = 58)
names(meta1)[4] <- "disease"

meta.ra <- meta1[c(1:19, which(meta1$disease == "RA")), ]
counts.ra <- counts1[, row.names(meta.ra)]
meta.sc <- meta1[c(20:28, which(meta1$disease == "Sarcoid")), ]
counts.sc <- counts1[, row.names(meta.sc)]

#STARMAPs using 1 CPU thread
res1 <- starmaps(data1 = counts.ra, meta1 = meta.ra, fcol1 = "disease", sbst1 = c("Control", "RA"),
                 data2 = counts.sc, meta2 = meta.sc, fcol2 = "disease", 
                 mc = 1)
plot.starmaps(res1, "Example1.pdf")
tbl1 <- summary.starmaps(res1)
write.csv(tbl1, "Example1.csv")


#### Example 2: variable(s) of interest is numeric, e.g., microbiota in response to a grandient of doses  

#get antohter set of data
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1919nnn/GSM1919500/suppl/GSM1919500_otu_table_MultiC.Helahty.s.ID.biom.txt.gz", 
              destfile="GSE74395_otu.gz")
counts2 <- read.delim(gzfile("GSE74395_otu.gz"), skip = 1, comment.char = "", check.names = F, row.names = 1)
counts2 <- aggregate(. ~ taxonomy, data = counts2, FUN = sum)
row.names(counts2) <- counts2$taxonomy
counts2 <- counts2[,-1]
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74395/suppl/GSE74395_Pneumotype.sep.Map.A1.txt.gz", 
              destfile="GSE74395_meta.gz")
meta2 <- read.delim(gzfile("GSE74395_meta.gz"), skip = 7, row.names = 1, nrows = 100)
meta2 <- meta2[meta2$source.name == "BAL", ]
counts2 <- counts2[, row.names(meta2)]

#STARMAPs using 16 CPU threads
res2 <- starmaps(data1 = counts1, meta1 = meta1, fcol1 = "disease", sbst1 = c("Control", "RA"),
                 counts2, meta2, "Fractalkine", 
                 mc = 16)
plot.starmaps(res2, "Example2.pdf")
tbl2 <- summary.starmaps(res2)
write.csv(tbl2, "Example2.csv")
