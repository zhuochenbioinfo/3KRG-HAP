library(ggplot2)
library(scales)
library(stringr)

args <- commandArgs(TRUE)

chrlist <- args[1]
input <- args[2]
output <- args[3]

data <- read.delim(chrlist,header=F)
bin <- read.delim(input,header=F)

plot.height <- 6
plot.width <- length(data$V1) * 0.5 + 4

pdf(output, width = plot.width, height = plot.height)

data$V1 <- factor(data$V1, levels = data$V1)
bin$V1 <- factor(bin$V1, levels = data$V1)

types <- length(levels(bin$V4))

ggplot(data=data) +
  geom_rect(aes(xmin = as.numeric(V1) - 0.2, 
				xmax = as.numeric(V1) + 0.2 , 
				ymax = V2, ymin = 0),
				colour="black", fill = "white") +
  #coord_flip() +
  geom_rect(data=bin, aes(xmin = as.numeric(V1) - 0.18, 
				xmax = as.numeric(V1) + 0.18, 
				ymax = V2, ymin = V3, fill=V4)) +
  #scale_fill_manual(values = c("aus"="#FFFF00", "ind"="#556B2F", "INDICA"="#9ACD32", "tej"="#00BFFF", "trj"="#8B008B", "JAPONICA"="#0000FF")) +
  scale_fill_manual(values = c("sp_aus"="#FF0000", "sp_ind"="#FFD700", "ss_INDx"="#FF8C00", "sp_tej"="#00FFFF", "sp_trj"="#800080", "ss_JAPx"="#4080C0", "admixed"="#C0C0C0")) +
  guides(fill=guide_legend(title="origin")) +
  theme(axis.text.x = element_text(colour = "black"),
		panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		panel.background = element_blank()) + 
  scale_x_discrete(position = "top", name = "chromosome", limits = data$V1) +
  scale_y_continuous(trans="reverse", labels = comma) +
  ylab("region (bp)")

dev.off()
