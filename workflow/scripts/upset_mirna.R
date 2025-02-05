#!/usr/bin/env Rscript

library(tidyverse)
library(ComplexUpset)
library(scales)
library(patchwork)
library(ragg)

args <- commandArgs(trailingOnly=TRUE)

# read novel mirna dataframe and count number of overlaps
#df <- read.table(args[1],
#                 sep = ",", header = TRUE) %>%
#  group_by(mat_seq) %>%
#  mutate(overlap = ifelse(n() == 1, "1", ifelse(n() == 2, "2", ">2")))
num_tissues <- length(unique(read.csv(args[1]) %>% pull(sample)))

df <- read.table(args[1],
                 sep = ",", header = TRUE) %>%
  group_by(mat_seq) %>%
  mutate(overlap = n()) %>%
  ungroup() %>%
  mutate(overlap_cat = ifelse(overlap == 1, "Private", 
                              ifelse(overlap == 2, "Shared by two", 
                                     ifelse(overlap == num_tissues, "Shared by all", "Shared by many"))))

###########################
# UPSET
###########################

tissue_l <- list()

# prepare upset data format
for (i in unique(df$sample)) {
  sample <- i
  sample = df %>% filter(sample == i) %>% pull(mat_seq)
  tissue_l[[i]] <- sample
}

m <- UpSetR::fromList(tissue_l)
m[colnames(m)] <- m[colnames(m)] == 1
sets <- colnames(m)

# plot upset
query_by_degree = function(data, groups, params_by_degree, ...) {
  intersections = unique(upset_data(data, groups)$plot_intersections_subset)
  lapply(
    intersections,
    FUN=function(x) {
      members = strsplit(x, '-', fixed=TRUE)[[1]]
      args = c(
        list(intersect=members, ...),
        params_by_degree[[length(members)]]
      )
      do.call(upset_query, args)
    }
  )
}

# get length of sets (ie number of samples)
set_len <- length(sets)
# generate list of overlap degrees with colors for 1, 2, and full overlap
degrees <- list()

for (i in 1:set_len) {
  # print(i)
  if (i == 1){
    degrees[[i]] <- list(color="#c7e9b4", fill="#c7e9b4")
  } else if (i == 2){
    degrees[[i]] <- list(color="#7fcdbb", fill="#7fcdbb")
  } else if (i == set_len) {
    degrees[[i]] <- list(color="#2c7fb8", fill="#2c7fb8")
  } else {
    degrees[[i]] <- list(color="#41b6c4", fill="#41b6c4")
  }
}

p1 <- ComplexUpset::upset(
  m,
  sets,
  base_annotations = list(
    'Intersection size'=intersection_size(counts=FALSE)
  ),
  matrix = (
    intersection_matrix(geom=geom_point(size=5))
  ),
  sort_intersections_by='degree',
  wrap = TRUE,
  themes=upset_modify_themes(
    list(
      'intersections_matrix'=theme(
        text=element_text(size=15),
        axis.title = element_blank()),
      'overall_sizes'=theme(
        axis.text.x=element_text(size=10,angle=90)),
      'Intersection size'=theme(
        axis.title=element_blank(),
        axis.text=element_text(size=15))
    )
  ),
  queries = query_by_degree(
    m,
    sets,
    params_by_degree = degrees,
    only_components=c("intersections_matrix", "Intersection size")
  )
) +
  ggtitle(str_split(sets[1], "_(?=\\d)", n=2)[[1]]) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    aspect.ratio = 1
  )

       #axis.text=element_text(size=15))
###########################
# HISTOGRAM
###########################

# reorder overlap levels
#df$overlap <- factor(df$overlap, levels = c("1","2",">2"))
df$overlap_cat <- factor(df$overlap_cat, 
    levels = c("Private","Shared by two","Shared by many","Shared by all"))

p2 <- ggplot(df, 
       aes(x = read_count, fill = overlap_cat)) +
  geom_density(position = "stack") +
  scale_x_log10(labels = comma) +
  labs(x = expression("Read count"~(log["10"]))) +
 #scale_fill_manual(values = c("#bae4bc","#7bccc4","#43a2ca")) +
  scale_fill_manual(values = c("#c7e9b4","#7fcdbb","#41b6c4","#2c7fb8")) +
  theme_classic() +
  guides(fill = guide_legend(title = "miRNA overlap")) +
  theme(
    # legend.position = "none",
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.key.size = unit(1, "cm"),
    legend.position = c(0.9, 0.9),
    aspect.ratio = 1
  )

###########################
# COMBINE AND SAVE
###########################

ggsave(
  tail(args,n = 1),
  wrap_plots(p1,p2),
  device = agg_tiff,
  width = 10, height = 6, units = "in", res = 300,
  scaling = 0.4
)

#ggsave(tail(args,n = 1), wrap_plots(p1,p2), width = 400, height = 300, units = "mm", device = "tiff", dpi = 300)

