library(ggplot2)

full_histo <- read.table('pairend/jelly_histo_31')

mx <- full_histo[full_histo$V2 == max(full_histo[50:200,2]),]

ggplot(full_histo, aes(x = V1, y = V2)) + geom_line(size=.8) +
  xlab('Number of times of occuring') + ylab("Count of occuring k-mers") + 
  ylim(c(0,54000)) + xlim(c(0,200)) + theme_bw() + 
  theme(axis.title = element_text(size = 16),
        axis.text= element_text(size = 12)) +
  geom_point(data = full_histo[12,], size = 4, shape=1) +
  geom_point(data = mx, aes(x=V1, y=V2), size=4, shape = 1) +
  geom_text(data = mx, aes(x=V1, y=V2, label='125'), nudge_y = 2500, size = 5) +
  geom_text(data = full_histo[12,], aes(label='13'), nudge_y = -2500, size = 5)

sum(as.numeric(full_histo[11:nrow(full_histo),1] * full_histo[12:nrow(full_histo),2]))/125
