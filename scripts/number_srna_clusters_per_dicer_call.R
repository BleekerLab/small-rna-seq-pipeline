suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(svglite))

# capture the command-line arguments after --args (e.g. the shortstack results directory)
args <- commandArgs(trailingOnly = TRUE)
shortstack_result_directory = args[1]
output_png = args[2]
output_svg = args[3]

# get the sample names in a list
samples = list.dirs(path = shortstack_result_directory,recursive = F,full.names = F)


# read the shortstack directories in a named list
shortstack.dirs = sapply(X = samples,FUN = function(x) {file.path(shortstack_result_directory,x,"Results.txt")})

# get all ShortStack results file into a list
shortstacks = map(shortstack.dirs,function(x) {fread(x,data.table=FALSE,header = T,sep="\t",check.names = F)})

# combine all shortstack Results.txt files row-wise
all.clusters = plyr::ldply(shortstacks,data.frame)
colnames(all.clusters)[1]="sample"  # to replace the .id created by ldply

# filter to keep only relevant columns to produce the nclusters = f(Dicer Call)
all.clusters.subset = dplyr::select(all.clusters,"sample","Name","DicerCall")

# count number of clusters per dicer call
nclusters_per_dicer_call = all.clusters.subset %>% 
  group_by(sample,DicerCall) %>%
  summarise(number_of_clusters = length(Name))


# plot N clusters = f(DicerCall)
g <- ggplot(data = nclusters_per_dicer_call,aes(x=DicerCall,y=number_of_clusters,fill=DicerCall)) +
  geom_bar(stat="identity",color="black") +
  scale_fill_brewer(palette="Set3") +
  facet_wrap( ~ sample)

# save the plots
ggsave(filename = output_png,plot = g,width = 7,height = 5,dpi = 400)
ggsave(filename = output_svg,plot = g,width = 7,height = 5)

