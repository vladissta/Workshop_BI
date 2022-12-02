library(ggplot2)

setwd('../pairend/jelly_res/')

histo <- read.table('jelly_histo_31')
histo_corrected <- read.table('jelly_histo_corrected_31')

mx <- histo[histo$V2 == max(histo[50:200,2]),]
mx_corrected <- histo_corrected[histo_corrected$V2 == max(histo_corrected[50:200,2]),]

ggplot() + 
  geom_line(data=histo_corrected, aes(x = V1, y = V2, col='red'), size=.5, alpha=.8) + 
  geom_line(data=histo, aes(x = V1, y = V2, col='blue'), size=.5, alpha=.8) +
  # geom_point(data = mx, aes(x=V1, y=V2), size=4, shape = 1) +
  # geom_text(data = mx, aes(x=V1, y=V2, label='125'), nudge_y = 2500, size = 5) +
  # geom_text(data = histo[12,], aes(x=V1, y=V2, label='12'), nudge_y = -2500, size = 5) +
  xlab('K-mer frequency') + ylab("Count") + 
  ylim(c(0,56000)) + 
  xlim(c(0,300)) +
  theme_bw() + 
  theme(axis.title = element_text(size = 12),
        axis.text= element_text(size = 12), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13), 
        legend.position = c(0.85, 0.8)) +
  scale_color_manual(values = c('blue', 'red'), 
                     labels = c('Raw', 'Corrected'), 
                     name = 'Library')

genome_estimation <- sum(as.numeric(full_histo[12:nrow(full_histo),1] * full_histo[12:nrow(full_histo),2]))/125
