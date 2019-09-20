##########################################
#### STARMAPs v2.0
#### 
#### Performs a comparison between two microbiota datasets
#### Linking the microbial changes associated with different treatments/factors in two datasets based on Princpal Component Analysis
#### Provides inference for alike and opposite matchings
#### Permutation-based tests
#### Includes three main functions: starmaps(), plot.starmaps(), and summary.starmaps()
####
#### Author: Peng Jiang (peng.jiang@northwestern.edu)
#### Reference: Jiang P, Green SJ, Chlipala GE, Turek FW, Vitaterna MH. Reproducible changes in the gut microbiome reveal a shift in microbial and host metabolism during spaceflight. Microbiome. 2019.Jiang P, Green SJ, Chlipala GE, Turek FW, Vitaterna MH. Reproducible changes in the gut microbiome reveal a shift in microbial and host metabolism during spaceflight. Microbiome. 2019 Aug 9;7(1):113. doi: 10.1186/s40168-019-0724-4
####
#### Last edited: 2019-03-28 
####
##########################################

library(MCMCpack)
library(vegan)
library(compositions)
library(parallel)
library(ggplot2)
library(gridExtra)

#### starmaps:  The main function of STARMAPs; Testing at each of the matched taxonomic level whether the group-segregations in two microbome datasets are similar/opposite .
####
#### Arguments:
#### data1/2:    count matrices (or data.frame) with samples in columns and features in rows.
####             Feature names (rownames) needs to be in the same format with each taxonomic level (beginning from Kingdom) seperated by ";".
####             Unasigned counts should be removed.
#### meta1/2:    A data.frame containing experimental conditions in columns, samples in rows matching samples in data1/data2.
####             Experimental conditions can be either factor/character or integer/numeric.
#### fcol1/2:    A character giving the name of the column in meta1/meta2 that specifies the experimental factor of interest.
#### mc:         Number of CPUs/threads allowed for multi-core processing, default = 1.
#### nperm:      Number of permutation in PERMANOVA (both sample-wise and feture wise).
#### nboot:      Number of bootstrap to test the similarity in the dirction of changes.
#### sbst1/2:    A character vector with length of 2 to specifiy subsets of groups to be compared, when cond has more than 2 levels.
####             Can be NULL if the cond is numeric or a factor (or character) of 2 levels.
####             Note: For data1, providing the full set of cond1 and specifying sbst1 is prefered than providing subsetted data, as the former capture the full picture in PCA (i.e., PC loading would be different)
####             For data2, providing the full set of cond2 and specifying sbst2 is not different from providing subsetted data in the test, other than the figure.
#### PCs:        A numeric vector speciflying the selected PC space for PERMANOVA; default: PC1 & PC2 
#### all.levels: TRUE/FALSE; whether to perform tests at higher taxonomic levels; 
####             if FALSE, only results at the input taxonomic level will be computed; default = TRUE   
####
#### Returns a list containing: 
#### metadata:         A list containing experimental varibles, variable names, selected PCs, etc.
#### Species - Phylum: Lists of results at each taxonomic level, depending lowest levels are provided in the input
####                   Each contains the following:
####                     overlap:     The proportion of features in data1 can be matched to data2
####                     counts1/2:   Taxa matched counts table
####                     ilr1/2:      ilr-transformed values
####                     ilr1.PCs:    Values of data1 PCs (all PCs without selection)
####                     ilr2.Prj:    Projected values of data2 on the same PC space
####                     pamanova1/2: aov tables from PERMANOVA in the selected PC space
####                                  parmanova2 gives both sample and feature permutation P values (sp.Pr(>F) and fp.Pr(>F), repectively).
####                     direction:   Comparasion of the directions of changes in each dataset.
####                                    Comparason1/2: Groups being compared
####                                    Distance1/2:   The estimated distance between the centers of groups or per 1 unit of numeric change in each dataset.
####                                    cos.theta:     cosine of the angel between the directions of changes in two datasets.
####                                                   cos.theta = 1, exact same direction; = -1 exact opposite direction; 0, perpendicular
####                                    P:             Bootstrap-determined P value of cos.theta != 0             
starmaps <- function(data1, meta1, fcol1, data2, meta2, fcol2, mc = 1, nperm = 1000, nboot = nperm,
                  PCs = NULL, sbst1 = NULL, sbst2 = NULL, all.levels = TRUE){
  
  cond1 <- meta1[, fcol1]
  cond2 <- meta2[, fcol2]
  
  na1 <- is.na(cond1)
  data1 <- as.matrix(data1[, !na1])
  cond1 <- cond1[!na1]
  
  na2 <- is.na(cond2)
  data2 <- as.matrix(data2[, !na2])
  cond2 <- cond2[!na2]
  
  if(!is.numeric(cond1)){
    cond1 <- droplevels(factor(cond1))
    if(is.null(sbst1)){
      if(length(levels(cond1)) > 2)
        stop("sbst1 cannot be NULL when cond1 involves more than 2 levels.") 
      else 
        sbst1 = levels(cond1)
    }
  }
  
  if(!is.numeric(cond2)){
    cond2 <- droplevels(factor(cond2))
    if(is.null(sbst2)){
      if(length(levels(cond2)) > 2)
        stop("sbst2 cannot be NULL when cond2 involves more than 2 levels.") 
      else 
        sbst2 = levels(cond2)
    }
  }
  
  #data1 <- data1[rowSums(data1) != 0, ]
  #data2 <- data2[rowSums(data2) != 0, ]
  
  output <- list()
  output[["metadata"]] <- list()
  output$metadata[["cond1"]] <- cond1
  output$metadata[["cond2"]] <- cond2
  output$metadata[["cond1.name"]] <- fcol1  
  output$metadata[["cond2.name"]] <- fcol2
  output$metadata[["sbst1"]] <- sbst1
  output$metadata[["sbst2"]] <- sbst2
  
  ### some cleaning to ensure best matching between the two datasets 
  taxa1 <- get.taxa(data1)
  taxa2 <- get.taxa(data2)
  
  lowLv <- min(ncol(taxa1), ncol(taxa2))
  
  taxa1 <- taxa1[, 1:lowLv]
  taxa2 <- taxa2[, 1:lowLv]
  
  data1 <- taxa.agg(data1, taxa1)
  data2 <- taxa.agg(data2, taxa2)

  taxa1 <- get.taxa(data1)
  taxa2 <- get.taxa(data2)
  
  ### taxa matching
  cat("Matching taxa...", "\n")
  taxa0 <- match.taxa(taxa1, taxa2)
  
  nr1 <- nrow(data1)
  nr2 <- nrow(data2)
  
  data1 <- rbind(data1, matrix(0, nr2, ncol(data1)))
  data2 <- rbind(matrix(0, nr1, ncol(data2)), data2)

  ### test at each level
  cat("Testing at the level of", "\n")
  lowLv <- ncol(taxa0)
  highLv <- ifelse(all.levels, 2, lowLv)
  
  for (i in lowLv:highLv){
    #i = lowLv
    LVi <- colnames(taxa0)[i]
    cat(LVi, "\n")
    output[[LVi]] <- list()
    
    f1 <- unique(apply(taxa0[1:nr1, 1:i], 1, paste, collapse = ";"))
    f2 <- unique(apply(taxa0[-(1:nr1), 1:i], 1, paste, collapse = ";"))
    fc <- intersect(f1, f2)
    output[[LVi]][["overlap"]] <- length(fc)/length(f1)
    
    data1.i <- taxa.agg(data1, taxa = taxa0, lv = i)
    data2.i <- taxa.agg(data2, taxa = taxa0, lv = i)
    
    output[[LVi]][["counts1"]] <- data1.i
    output[[LVi]][["counts2"]] <- data2.i
    
    ## ilr
    ilr1 <- apply(data1.i, 2, PCilr)
    ilr2 <- apply(data2.i, 2, PCilr)
    
    output[[LVi]][["ilr1"]] <- ilr1
    output[[LVi]][["ilr2"]] <- ilr2
    
    ## PCA via singular value decomposition
    # centering to the center of ilr1
    ctr<- rowMeans(ilr1) 
    ilr1 <- t(ilr1 - ctr)
    ilr2 <- t(ilr2 - ctr)
    
    # svd
    svd1 <- svd(ilr1)
    # PCs from svd
    ilr1.PCs <- ilr1 %*% svd1$v
    # apply the same loadings to irl2
    ilr2.Prj <- ilr2 %*% svd1$v
    
    colnames(ilr1.PCs) <- paste0("PC", 1:ncol(ilr1.PCs))
    colnames(ilr2.Prj) <- paste0("PC", 1:ncol(ilr2.Prj))
    
    output[[LVi]][["ilr1.PCs"]] <- ilr1.PCs
    output[[LVi]][["ilr2.Prj"]] <- ilr2.Prj

    ## PERMANOVA
    # if PC space is selected, use selected, otherwise use all
    if(is.null(PCs)) PCspace <- 1:ncol(ilr1.PCs) else PCspace <- PCs
    output$metadata[["PCs"]] <- PCspace
    
    if(is.numeric(cond1)){
      ilr1.PCs_ <- ilr1.PCs[, PCspace]
      cond1_ <- cond1
    }else{
      sel1 <- which(cond1 %in% sbst1)
      ilr1.PCs_ <- ilr1.PCs[sel1, PCspace]
      cond1_ <- factor(cond1[sel1], levels=sbst1)
    }
    
    permanova1 <- vegan::adonis(ilr1.PCs_ ~ cond1_, permutations = nperm, method = "euclidean", 
                                contr.unordered = "contr.treatment")
    
    aov.tbl1 <- permanova1$aov.tab
    row.names(aov.tbl1)[1] <- ifelse(is.null(sbst1), fcol1, paste(sbst1, collapse = " vs. "))
    
    if(is.numeric(cond2)){
      ilr2.Prj_ <- ilr2.Prj[, PCspace]
      cond2_ <- cond2
      sel2 <- 1:length(cond2)
    }else{
      sel2 <- which(cond2 %in% sbst2)
      ilr2.Prj_ <- ilr2.Prj[sel2, PCspace]
      cond2_ <- factor(cond2[sel2], levels=sbst2)
    }
    
    permanova2 <- vegan::adonis(ilr2.Prj_ ~ cond2_, permutations = nperm, method = "euclidean", 
                                contr.unordered = "contr.treatment")
    aov.tbl2 <- permanova2$aov.tab
    row.names(aov.tbl2)[1] <- ifelse(is.null(sbst2), fcol2, paste(sbst2, collapse = " vs. "))
    
    # "feature" permutation of ilr-transformed data
    if(mc > 1){
      cl <- parallel::makeCluster(mc, type="PSOCK")
      fp.F <- parallel::parSapply(cl, 1:nperm, fperm.adonis2, ilr = ilr2, V = svd1$v, sel = sel2, PCs = PCspace, cond = cond2_)
      parallel::stopCluster(cl);rm(cl)
    }else{
      fp.F <- sapply(1:nperm, fperm.adonis2, ilr = ilr2, V = svd1$v, sel = sel2, PCs = PCspace, cond = cond2_)
    }
    
    fp.P <- sum(fp.F > aov.tbl2[1, 4])/(1 + nperm)
    
    aov.tbl2[,"fp.Pr(>F)"] <- c(fp.P, NA, NA)
    names(aov.tbl2)[6] <- "sp.Pr(>F)"
    row.names(aov.tbl2)[1] <- ifelse(is.null(sbst2), fcol2, paste(sbst2, collapse = " vs. "))
    
    output[[LVi]][["permanova1"]] <- aov.tbl1
    output[[LVi]][["permanova2"]] <- aov.tbl2

    # vectors of the group effect in each of the datasets
    vec1 <- permanova1$coefficients[2, ]
    vec2 <- permanova2$coefficients[2, ]
    # angle between vectors
    cos.theta <- sum(vec1 * vec2) / (sqrt(sum(vec1 * vec1)) * sqrt(sum(vec2 * vec2)))
    
    # bootstrap
    if(mc > 1){
      cl <- parallel::makeCluster(mc, type="PSOCK")
      bt.cos <- parallel::parSapply(cl, 1:nboot, boot.angle, prj = ilr2.Prj_, cond = cond2_, vec1 = vec1)
      parallel::stopCluster(cl);rm(cl)
    }else{
      bt.cos <- sapply(1:nboot, boot.angle, prj = ilr2.Prj_, cond = cond2_, vec1 = vec1)
    }
    
    bt.cos <- bt.cos - cos.theta
    bt.P <- mean(abs(bt.cos) > abs(cos.theta))
    
    direction <- data.frame(Comparison1 = ifelse(is.null(sbst1), paste("Increase in", fcol1), paste(sbst1, collapse = " vs. ")), 
                            Distance1 = sqrt(sum(vec1 * vec1)),
                            Comparison2 = ifelse(is.null(sbst2), paste("Increase in", fcol2), paste(sbst2, collapse = " vs. ")),
                            Distance2 = sqrt(sum(vec2 * vec2)),
                            cos.theta = cos.theta,
                            P = bt.P)
    
    output[[LVi]][["direction"]] <- direction
    
  }# loop through taxa levels
  cat("Done!", "\n")
  return(output)
}


#### plot.starmaps: A function to plot results of starmaps.
####
#### Arguments:
#### starmap:    The output object from the starmaps function.
#### pc.plane:   The PC plane to plot; default: PC1 & PC2; currently the graph is 2D only.
#### pdf:        Name and path of the figure output file (pdf).
#### col1/col2:  A vector of colors mapped to experimental groups (pass to scale_colour_manual()).
####             If the experimental condition/factor is numeric, col1/col2 can be used to specify n-color gradient.
#### pch1/pch2:  Plotting character/symbol to use (same as pch in points{graphics}).
####             Note that col2 is mapped to colour() and col1 is mapped to fill(), thus pch1 needs to be from 21:25.
#### cex1/cex2:  Size of the plotting symbol.
#### alpha:      Set the fading/transparency of data1 points on the overlap plot.
####
#### Produces a graph at each taxa level found in the [starmap] object; each graph contains 5 panels.
#### panel 1: (top-right): PC values of data1
#### panel 2: (top-left):  projected values of data2
#### panel 3: (mid-right): PERMANOVA table for data1 on the plotted PC plane
#### panel 4: (mid-left):  PERMANOVA table for data2 on the plotted PC plane
####                       gives both P values from sample permutations (sp.Pr(>F)) and  
#### panel 5: (bottom):    Table comparing directions of changes.
####                         Comparason1/2: Groups being compared
####                         Distance1/2:   The estimated distance between the centers of groups or per 1 unit of numeric change in each dataset.
####                         cos.theta:     cosine of the angel between the directions of changes in two datasets.
####                                        cos.theta = 1, exact same direction; = -1 exact opposite direction; 0, perpendicular
####                         P:             Bootstrap-determined P value of cos.theta != 0  
plot.starmaps <- function(starmap, pdf = "starmaps.pdf", pc.plane = c(1, 2),
                          col1 = NULL, pch1 = 21, cex1 = 2, col2 = NULL, pch2 = 18, cex2 = 3, alpha = 0.3){
  
  levels <- names(starmap)[-1]
  
  cond1 <- starmap$metadata$cond1
  cond2 <- starmap$metadata$cond2
  
  cond1.name <- starmap$metadata$cond1.name
  cond2.name <- starmap$metadata$cond2.name
  
  pdf(pdf, width = 7.5, height = 5.7)
  
  for(lv in levels){
    
    ilr1.PCs <- as.data.frame(starmap[[lv]]$ilr1.PCs)[, pc.plane]
    ilr2.Prj <- as.data.frame(starmap[[lv]]$ilr2.Prj)[, pc.plane]
    
    # panel 1 -- data1 PCs 
    xl <- min(ilr1.PCs[,1], ilr2.Prj[,1])
    xh <- max(ilr1.PCs[,1], ilr2.Prj[,1])
    yl <- min(ilr1.PCs[,2], ilr2.Prj[,2])
    yh <- max(ilr1.PCs[,2], ilr2.Prj[,2])
    
    p1 <- ggplot() + theme_bw() + xlim(xl, xh) + ylim(yl, yh) +
      geom_point(aes(x = PC1, y = PC2, fill = cond1), data = ilr1.PCs, shape = pch1, size = cex1, colour = "black")
    
    if(is.numeric(cond1)){
      if(is.null(col1))
        p1 <- p1 + scale_fill_viridis_c()
      else p1 <- p1+ scale_fill_gradientn(colors = col1)
    } else {
      if(is.null(col1))
        p1 <- p1 + scale_fill_brewer(palette = "Set1")
      else p1 <- p1 + scale_fill_manual(values = col1)
    }
    p1 <- p1 + labs(fill = cond1.name) + theme(legend.key.size = unit(0.15, "in"), aspect.ratio = 1)
    
    # panel 2 -- data2 projected on to data1 PCs
    p2 <- ggplot() + theme_bw() + xlim(xl, xh) + ylim(yl, yh) +
      geom_point(aes(x = PC1, y = PC2, fill = cond1), data=ilr1.PCs, shape = pch1, size = cex1, alpha = alpha, colour = "grey30")
    
    if(is.numeric(cond1)){
      if(is.null(col1))
        p2 <- p2 + scale_fill_viridis_c()
      else p2 <- p2 + scale_fill_gradientn(colors = col1)
    } else {
      if(is.null(col1))
        p2 <- p2 + scale_fill_brewer(palette = "Set1")
      else p2 <- p2 + scale_fill_manual(values = col1)
    }
      
    p2 <- p2 + geom_point(aes(x = PC1, y = PC2, colour = cond2), data=ilr2.Prj, shape = pch2, size = cex2)
    
    if(is.numeric(cond2)){
      if(is.null(col2))
        p2 <- p2 + scale_colour_viridis_c(option = "magma")
      else p2 <- p2 + scale_colour_gradientn(colors = col2)
    } else {
      if(is.null(col2))
        p2 <- p2 + scale_colour_brewer(palette = "Dark2")
      else p2 <- p2 + scale_colour_manual(values = col2)
    }
    
    p2 <- p2 + labs(colour = cond2.name, fill = cond1.name) + theme(legend.key.size = unit(0.15, "in"), aspect.ratio = 1)

    # panel 3 -- permanova table for data1
    tt <- gridExtra::ttheme_default(
      core = list(fg_params=list(cex = 0.5), padding=unit(c(2, 2), "mm")),
      colhead = list(fg_params=list(cex = 0.45), padding=unit(c(2, 2), "mm")),
      rowhead = list(fg_params=list(cex = 0.45)), padding=unit(c(2, 2), "mm"))
    
    p3 <- gridExtra::tableGrob(signif(starmap[[lv]]$permanova1, 3), theme = tt)
    
    # panel 4 -- permanova table for data2
    p4 <- gridExtra::tableGrob(signif(starmap[[lv]]$permanova2, 3), theme = tt)
    
    # panel 5 -- direction table
    direction <- starmap[[lv]]$direction
    direction[,c(2,4,5,6)] <- signif(direction[,c(2,4,5,6)], 3)
    
    p5 <- gridExtra::tableGrob(direction, theme = tt, rows = NULL)
    
    # plot
    gridExtra::grid.arrange(p1, p2, p3, p4, p5, heights = c(2, 1, 1),
                            layout_matrix = cbind(c(1, 3, 5), c(2, 4, 5)), 
                            top = lv)
    
  }# loop through taxa levels
  
  dev.off()
}


#### summary.starmaps: A function to get a summary table from the results of starmaps.
####
#### Arguments:
#### starmap:   The output object from the starmaps function.
####
#### Returns a table with each taxa level in each rows, columns include:
#### Level:            The taxonomic level at which the test was performed.
#### Overlap:          The proportion of features in data1 can be matched to data2
#### Data1.Group:      Groups being compared in dataset1.
#### Data1.F:          PERMANOVA F statistic for data1
#### Data1.Pr(>F):     PERMANOVA P value for data1 (sample permutation).
#### Data2.Group:      Groups being compared in dataset2.
#### Data2.F:          PERMANOVA F statistic for data2
#### Data2.sp.Pr(>F):  PERMANOVA P value for data2 (sample permutation).
#### Data2.fp.Pr(>F):  PERMANOVA P value for data2 (feature permutation).
#### cos(theta):       cosine of the angel between the directions of changes in two datasets.
####                   cos.theta = 1, exact same direction; = -1 exact opposite direction; 0, perpendicular.
#### Pr(cos(theta)=0): Bootstrap-determined P value of cos.theta = 0.
#### omnibus.Pval:     The overall test P value, taken as the max of all individual P values.
#### Similarity.call:  The group differences in two datasets are "alike", "opposite", or "n.s." not significantly similar.
summary.starmaps <- function(starmap){
  
  tbl <- list()
  levels <- names(starmap)[-1]
  for(lv in levels){
    overlap <- starmap[[lv]]$overlap
    grp1 <- as.character(starmap[[lv]]$direction[1,1])
    F1 <- starmap[[lv]]$permanova1[1, 4]
    pval1 <- starmap[[lv]]$permanova1[1, 6]
    grp2 <- as.character(starmap[[lv]]$direction[1,3])
    F2 <- starmap[[lv]]$permanova2[1, 4]
    pval2s <- starmap[[lv]]$permanova2[1, 6]
    pval2f <- starmap[[lv]]$permanova2[1, 7]
    pval2 <- max(pval2s, pval2f)
    cos.alpha <- starmap[[lv]]$direction[1, 5]
    pval3 <- starmap[[lv]]$direction[1, 6]
    omniP <- max(pval1, pval2, pval3)
    call0 <- ifelse(omniP >= 0.05, "n.s.", ifelse(cos.alpha > 0, "Alike", "Opposite"))
    
    tbl.lv <- data.frame(lv, overlap, grp1, F1, pval1, grp2, F2, pval2s, pval2f, pval2, cos.alpha, pval3, omniP, call0)
    names(tbl.lv) <- c("Level", "Overlap", "Data1.Group", "Data1.F", "Data1.Pr(>F)", 
                       "Data2.Group", "Data2.F", "Data2.sp.Pr(>F)", "Data2.fp.Pr(>F)", "Data2.Pr(>F)",
                       "cos(theta)", "Pr(cos(theta)=0)", "omnibus.Pval", "Similarity.call")
    tbl[[lv]] <- tbl.lv
  }
  tbl <- do.call(rbind, tbl)
  return(tbl)
}


#### get.taxa: A function to get and clean taxa names from row.names of a data matrix
####           Returns a matrix in which each column is a taxonimic level
####
#### Arguments:
#### counts:  A data matrix whose rownames are feature names, which begins at the Kingdom level
#### sep:     The character used to seperate each taxonomic level
get.taxa <- function(counts, sep=";"){
  
  taxa <- rownames(counts)
  lvs <- nchar(taxa) - nchar(gsub(sep, "", taxa))
  lowLv <- max(lvs) + 1
  st <- lowLv - lvs
  sps <- sapply(st, function(s){
    sp <- rep(sep, s)
    sp <- paste(sp, collapse = "")
  })
  
  taxa <- paste0(taxa, sps)
  
  taxa <- do.call(rbind, strsplit(taxa, sep))
  
  taxa <- gsub(" ", "", taxa)
  taxa <- gsub(".*__", "", taxa)
  taxa <- gsub("[[]", "", taxa)
  taxa <- gsub("]", "", taxa)
  taxa <- gsub("[.].*", "", taxa)
  taxa[which(taxa == "")] <- "Other"

  colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")[1:lowLv]
  
  return(taxa)
}


#### taxa.agg: a function to aggregate counts based on taxa
####
#### Arguments:
#### counts:  a matrix of counts with samples in columns and features in rows
#### taxa:    a matrix of taxa with each of the levels seperated in columns and rows matching features in counts
####          expecting output style from get.taxa()
#### lv:      a integer specifiying the taxonomic level to which the counts are aggregated
taxa.agg <- function(counts, taxa, lv = ncol(taxa)){
  counts <- as.data.frame(counts)
  counts$taxa <- apply(taxa[, 1:lv], 1, paste, collapse = ";")
  counts <- aggregate(.~ taxa, sum, data = counts)
  row.names(counts) <- counts$taxa
  counts$taxa <- NULL
  counts <- as.matrix(counts)
  
  return(counts)
}


#### match.taxa:  A function to match two sets of taxa
####              The function loops through each taxonomic level, from low to high;
####              At level i, if some taxa matches at the level i-1 but not i,
####              and if these taxa in one dataset contains "Other" at Level i,  
####              the taxa that matched at Level i-1 in the other dataset will be set to "Other",
####              so that they can be matached at level i.
####
#### Arguments:
#### taxa1/taxa2: Two matrices of taxa with the same set of levels (i.e., columns); 
####              Expecting output from get.taxa()
match.taxa <- function(taxa1, taxa2){
  
  if(!all(colnames(taxa1) == colnames(taxa2)))
    stop("taxonomic levels do not match!")
  
  lowLv <- ncol(taxa1)
  
  for(i in lowLv:2){
    
    ## comparing at the current level i
    taxa1L1 <- apply(taxa1[, 1:i], 1, paste, collapse=";")
    taxa2L1 <- apply(taxa2[, 1:i], 1, paste, collapse=";")
    commonL1 <- intersect(taxa1L1, taxa2L1)
    
    ## copmaring at Level i-1
    taxa1L0 <- apply(as.data.frame(taxa1[, 1:(i-1)]), 1, paste, collapse=";")
    taxa2L0 <- apply(as.data.frame(taxa2[, 1:(i-1)]), 1, paste, collapse=";")
    commonL0 <- intersect(taxa1L0, taxa2L0)
    
    ## indexing whether common (T) or unique (F)
    taxa1idx <- data.frame(L0 = taxa1L0 %in% commonL0, L1 = taxa1L1 %in% commonL1)
    taxa2idx <- data.frame(L0 = taxa2L0 %in% commonL0, L1 = taxa2L1 %in% commonL1)
    
    ## take the taxa at Level i-1 for those did not match at level i but matched at Level i-1 
    uniqueL0 <- c(taxa1L0[taxa1idx$L0 & (!taxa1idx$L1)], taxa2L0[taxa2idx$L0 & (!taxa2idx$L1)])
    
    ## loop through uniqueL0 taxa
    for (tx in uniqueL0){
      
      # convert if "Other" is found in the other dataset
      convert1 <- "Other" %in% taxa2[taxa2L0 == tx, i]
      convert2 <- "Other" %in% taxa1[taxa1L0 == tx, i]
      
      tx1L1 <- which(taxa1L0 == tx & taxa1idx$L1 == F)
      tx2L1 <- which(taxa2L0 == tx & taxa2idx$L1 == F)
      
      if (convert1) taxa1[tx1L1, i:lowLv] = "Other"
      if (convert2) taxa2[tx2L1, i:lowLv] = "Other"
    }
  }
  
  ## combine all taxa after matching
  matched <- rbind(taxa1, taxa2)
  
  return(matched)
}


#### PCilr: A function to perform isometric log-ratio transformation in perperation for PCA
####
#### Arguments:
#### x:    A vector of feature counts of a sample
#### dmc:  Number of Dirichlet Monte Carlo incidences
#### V:    A matrix, with columns giving the chosen basis of the clr-plane
####       By default, provided by compositions::ilrBase(D, "balanced"), where D is the number of features in x.
PCilr <- function(x, dmc = 1000, V = NULL){
  
  ### draw from Dirichlet distribution and get a point estimate, thus avioding constant 0s
  x <- x + 0.5
  dmc.x <- MCMCpack::rdirichlet(dmc, x)
  es <- colMeans(dmc.x)
  
  ### ilr
  if (is.null(V)) V <- compositions::ilrBase(D = length(es), "balanced")
  ilr.x <- compositions::ilr(es , V = V)
  
  return(ilr.x)
}


#### fperm.adonis2: A function to permute the ilr-transformed data2 by the features and compute PERMANOVA F value on permuted data
####
#### Arguments:
#### perm.i:  The i-th permutation
#### ilr:     A matrix of ilr-transformed data2
#### V:       The rotation matrix given by PCA of data1
#### sel:     The selection of samples in data2
#### PCs:     The selected PC space for data2 projections
#### cond:    A vector specifying selected experimental groups for ilr, can be either factor/character or integer/numeric.
fperm.adonis2 <- function(perm.i, ilr, V, sel, PCs, cond){
  fp.ilr <- ilr[, sample(ncol(ilr))]
  fp.Prj <- fp.ilr %*% V
  fp.Prj <- fp.Prj[sel, PCs]
  dat <- data.frame(cond = cond)
  fp.aov <- vegan::adonis2(fp.Prj ~ cond, data = dat, permutations = 0, method = "euclidean", 
                           contr.unordered = "contr.treatment")
  fval <- fp.aov["cond", "F"]
  names(fval) <- perm.i
  return(fval)
}


#### boot.angle: A function to perform a bootstrap of PC-projected data2 by the samples and compute the angle between data1 and data2 changes
####
#### Arguments:
#### boot.i:  The i-th bootstrap
#### prj:     A matrix of the projected values of ilr-transformed data2 in the selected PC space of data1
#### cond:    A vector specifying selected experimental groups for prj, can be either factor/character or integer/numeric.
#### vec1:    The vector of data1 group effect
boot.angle <- function(boot.i, prj, cond, vec1){
  
  if(is.numeric(cond))
    sl <- sample(1:length(cond), replace = T) else{
      cond <- factor(cond)
      sl1 <- sample(which(cond == levels(cond)[1]), summary(cond)[1], replace = T)
      sl2 <- sample(which(cond == levels(cond)[2]), summary(cond)[2], replace = T)
      sl <- c(sl1, sl2)
    }
  
  boot.data <- prj[sl, ]
  boot.cond <- cond[sl]
  
  boot.permanova <- vegan::adonis(boot.data ~ boot.cond, permutations = 0, method = "euclidean", 
                                  contr.unordered = "contr.treatment")
  
  boot.vec2 <- boot.permanova$coefficients[2, ]
  
  boot.cos.theta <- sum(vec1 * boot.vec2) / (sqrt(sum(vec1 * vec1)) * sqrt(sum(boot.vec2 * boot.vec2)))
  
  names(boot.cos.theta) <- boot.i
  return(boot.cos.theta)
}
  
