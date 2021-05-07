#' This R script can be applied to parse the GCT or GCTX files by using the cmapR package.
#'
#' Source of CMAP-LINCS GCTX file:
#'     GSE92742 - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742
#'     GSE70138 - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138
#'     GSE92743 - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92743
#'     
#' Date: 3/8/2018
#' Author: C.W.Weng (Jeff)

# install the ggplot2 package
install.packages("ggplot2")

# use the bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# install the cmpR package from bioconductor
BiocManager::install("cmapR")

# load packages
library(cmapR)
library(ggplot2)

# access the GCT class help
?`GCT-class`

## parse the entire file
# create a variable to store the path to the GCTX file
ds_path <- "GSE92742_Broad_LINCS_Level4_ZSPCINF_mlr12k_n1319138x12328.gctx"
my_ds <- parse.gctx(ds_path)

## parse metadata from a GCTX file
# extract the column metadata
col_meta_from_gctx <- read.gctx.meta(ds_path,dim="col")
# extract the row metadata
row_meta_from_gctx <- read.gctx.meta(ds_path,dim="row")

## parse a subset of a GCTX file
# read the column annotations as a data.frame
col_meta_path <- "GSE92742_Broad_LINCS_inst_info.txt"
col_meta <- read.delim(col_meta_path,sep="\t",stringsAsFactors=F)
# search the 'pert_iname' column
idx_1 <- which(col_meta$pert_iname=="proscillaridin")
idx_2 <- which(col_meta$pert_iname=="DMSO")
# get the corresponding sig_id
sig_ids <- col_meta$sig_id[idx_1]
sig_ids <- col_meta$sig_id[idx_2]
# get the corresponding distil_id
sig_id2distil_ids_1 <- col_meta$distil_id[idx_1]
distil_ids_1 <- vector(mode="character",length=0)
for(i in 1:length(sig_id2distil_ids_1)){
    print(paste0(sig_id2distil_ids_1[i]))
    sig_id2distil_ids_1_split <- as.vector(strsplit(sig_id2distil_ids_1[i],split='|', fixed=TRUE))
    distil_ids_1 <- unlist(c(distil_ids_1,sig_id2distil_ids_1_split),use.names=FALSE)
}
sig_id2distil_ids_2 <- col_meta$distil_id[idx_2]
distil_ids_2 <- vector(mode="character",length=0)
for(i in 1:length(sig_id2distil_ids_2)){
    print(paste0(sig_id2distil_ids_2[i]))
    sig_id2distil_ids_2_split <- as.vector(strsplit(sig_id2distil_ids_2[i],split='|', fixed=TRUE))
    distil_ids_2 <- unlist(c(distil_ids_2,sig_id2distil_ids_2_split),use.names=FALSE)
}
# get the corresponding inst_id
inst_ids_1 <- col_meta$inst_id[idx_1]
inst_ids_2 <- col_meta$inst_id[idx_2]
# read only those columns from the GCTX file by using the 'cid' parameter
instance_ds_1 <- parse.gctx(ds_path,cid=inst_ids_1)
instance_ds_2 <- parse.gctx(ds_path,cid=inst_ids_2)
# add annotations to a GCT object
(instance_ds_annotated_1 <- annotate.gct(instance_ds_1,col_meta,dim="col",keyfield="inst_id"))
(instance_ds_annotated_2 <- annotate.gct(instance_ds_2,col_meta,dim="col",keyfield="inst_id"))

## slice a GCTX object
idx1 <- which(instance_ds_annotated_1@cdesc$cell_id=="A549")
instance_ds_subset1 <- subset.gct(instance_ds_annotated_1,cid=idx1)
idx2 <- which(instance_ds_annotated_1@cdesc$cell_id=="HCC515")
instance_ds_subset2 <- subset.gct(instance_ds_annotated_1,cid=idx2)
idx3 <- which(instance_ds_annotated_2@cdesc$cell_id=="A549")
instance_ds_subset3 <- subset.gct(instance_ds_annotated_2,cid=idx3)
idx4 <- which(instance_ds_annotated_2@cdesc$cell_id=="HCC515")
instance_ds_subset4 <- subset.gct(instance_ds_annotated_2,cid=idx4)
instance_ds_first_25 <- subset.gct(instance_ds_annotated_1, rid=1:25)
# check the dimensionality of GCTX object
message("instance_ds_subset1:")
dim(instance_ds_subset1)
message("instance_ds_subset2:")
dim(instance_ds_subset2)
message("instance_ds_subset3:")
dim(instance_ds_subset3)
message("instance_ds_subset4:")
dim(instance_ds_subset4)
message("instance_ds_first_25:")
dim(instance_ds_first_25)

## merge GCTX objects together by columns
# because the two objects share common rows.
# this is equivalent to a 'cbind' operation on matrices
instance_ds_subset1_2 <- merge.gct(instance_ds_subset1,instance_ds_subset2,dim="column")
instance_ds_subset1_2_3 <- merge.gct(instance_ds_subset1_2,instance_ds_subset3,dim="column")
instance_ds_subset1_2_3_4 <- merge.gct(instance_ds_subset1_2_3,instance_ds_subset4,dim="column")
table(instance_ds_subset1_2_3_4@cdesc$cell_id)

## melt GCTX objects
# plot the z-scores of these 25 genes grouped by dose
(instance_ds_first_25_melted <- melt.gct(instance_ds_first_25))
ggplot(instance_ds_first_25_melted)+geom_boxplot(aes(x=pert_dose,y=value))+ylab("z-score")

## math operations on GCT objects
# compute the row and column means
row_means <- rowMeans(instance_ds_first_25@mat)
col_means <- colMeans(instance_ds_first_25@mat)
message("means:")
head(row_means)
head(col_means)
# using 'apply', compute the max of each row and column
row_max <- apply(instance_ds_first_25@mat,1,max)
col_max <- apply(instance_ds_first_25@mat,2,max)
message("maxes:")
head(row_max)
head(col_max)

## GCT-specific math functions
# transposing a GCT object - also swaps row and column annotations
(instance_ds_first_25_transpose <- transpose.gct(instance_ds_first_25))
# converting a GCT object's matrix to ranks
# the 'dim' option controls the direction along which the ranks are calculated
instance_ds_first_25_rank_by_column <- rank.gct(instance_ds_first_25,dim="col")
# plot z-score vs rank for genes (rows)
plot(instance_ds_first_25_rank_by_column@mat[1:25, ], instance_ds_first_25@mat[1:25, ],xlab="rank",ylab="z-score",main="z-score vs. rank")

## write the GCT or GCTX files
# output both GCT and GCTX format
write.gct(instance_ds_subset1_2_3_4,"instance_ds_subset1_2_3_4")
write.gctx(instance_ds_subset1_2_3_4, "instance_ds_subset1_2_3_4")
# write.gctx can also compress the dataset upon write,
# which can be controlled using the 'compression_level' option.
# the higher the value, the greater the compression, but the
# longer the read/write time
write.gctx(instance_ds_subset1_2_3_4, "instance_ds_subset1_2_3_4_compressed",compression_level=9)
