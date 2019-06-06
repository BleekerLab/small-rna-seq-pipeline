suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(svglite))

# capture the command-line arguments after --args (e.g. the shortstack results directory)
args <- commandArgs(trailingOnly = TRUE)
blast_file = args[1]
output_png = args[2]
output_svg = args[3]

# read the blast result file
df = read.delim(file = blast_file,header = T,stringsAsFactors = F)

# count occurences
MIR.freqs = as.data.frame(table(df$subject_id))
colnames(MIR.freqs)=c("MIR","Counts")

# plot N clusters = f(DicerCall)
g <- ggplot(data = MIR.freqs,aes(x=MIR,y=Counts,fill=MIR)) +
  geom_bar(stat="identity",color="black") +
  scale_fill_brewer(palette="Set3") +
  labs(x = "MIR gene family",y="Counts")

# save the plots
ggsave(filename = output_png,plot = g,width = 7,height = 5,dpi = 400,device = "png")
ggsave(filename = output_svg,plot = g,width = 7,height = 5,device = "svg")

