#!/usr/bin/env Rscript

library(tidyverse)
library(ComplexUpset)

args <- commandArgs(trailingOnly=TRUE)

nc_lists <- list(
  "RNAsamba (horse)" = read.table(args[1])$V1,
  "CPAT"             = read.table(args[2])$V1,
  "FEELnc"           = read.table(args[3])$V1
)

# get upset
nc_upset <- UpSetR::fromList(nc_lists)
nc_upset[colnames(nc_upset)] = nc_upset[colnames(nc_upset)] == 1
names = colnames(nc_upset)

# and plot
#png(file = tail(args,n = 1), height = 600, width = 600)
p <- ComplexUpset::upset(
  nc_upset, names, name = 'Coding potential', stripes = 'white',
  base_annotations=list(
      'Intersection size' = intersection_size(),
      'Intersection ratio' = intersection_ratio(text_mapping=aes(label=!!upset_text_percentage()))
    )
)
#dev.off()
ggsave(tail(args,n = 1), p, width = 6, height = 6, device = 'tiff', dpi = 600)
