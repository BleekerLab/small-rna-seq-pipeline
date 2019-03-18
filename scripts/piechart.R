suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(svglite))

# capture the command-line arguments after --args (e.g. the shortstack results directory)
args <- commandArgs(trailingOnly = TRUE)
shortstack_result_file = args[1]
piechart_plot_file = args[2]

# read the shortstack results 
df = fread(shortstack_result_file,data.table = F,header=T,sep="\t",check.names = F,stringsAsFactors = F)
colnames(df)[1]="sample"  # to replace the .id created by ldply

# filter to keep only relevant columns to produce the nclusters = f(Dicer Call)
df.subset = dplyr::select(df,"Name","MIRNA")
df.counts = as.data.frame(table(df.subset$MIRNA)) # count the occurences of MIRNA classes
colnames(df.counts)=c("MIRNA","Freq")
df.counts = df.counts %>% mutate(percentage = round(Freq / sum(Freq) * 100))
df.counts = arrange(df.counts,percentage)

# make the pie chart plot N clusters = f(DicerCall)
g <- ggplot(data = df.counts,aes(x="",y=percentage,fill=MIRNA)) +
   geom_bar(stat="identity",color="black",width=1) +
   scale_fill_brewer(palette="Set3") +
  coord_polar(theta = "y",start=0)

# save the plots
ggsave(filename = piechart_plot_file,plot = g,width = 7,height = 5,dpi = 400,device = "png")

