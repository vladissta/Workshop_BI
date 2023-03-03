library(sleuth)
library(dplyr)
library(ggplot2)
library(latex2exp)

#<------------------READ-FILE--------------------------->

sample_id = dir(file.path("./kallisto_res/")) 
# %>% gsub('_out', '', .)

dirs <- sapply(sample_id, function(id) file.path('./kallisto_res', id))
conditions = rep(c('0_mins', '30_mins'), each = 2)

df = data.frame(sample = sample_id, path = dirs, conditions = conditions)

#<------------------DE--------------------------->

so = sleuth_prep(df, extra_bootstrap_summary = T, read_bootstrap_tpm = TRUE, ~conditions)
so = sleuth_fit(so, ~conditions, 'full')
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')

lfc<- kallisto_table(so) %>%
  group_by(target_id, conditions) %>%
  summarize(mean_ct=mean(est_counts)) %>%
  group_by(target_id) %>% 
  summarise(lFC = log2((mean_ct[conditions=='30_mins']+0.01)/(mean_ct[conditions=="0_mins"]+0.01)))

sleuth_table <- sleuth_results(so, test = "reduced:full",test_type = 'lrt')
sleuth_table = merge(sleuth_table, lfc, by = 'target_id') %>% na.omit() %>% arrange(qval)

#<------------------HEATMAP--------------------------->

pdf('heat.pdf')
plot_transcript_heatmap(so, transcripts = sleuth_table[1:50,]$target_id, 
                        color_high = 'orangered', color_low = 'white', 
                       fontsize_row = 5)
dev.off()

#<------------------SIGNIFICANT--------------------------->

de_sig_df = sleuth_table %>% filter(qval < 0.01, abs(lFC) > 2)

nrow(de_sig_df)
nrow(de_sig_down)
nrow(de_sig_up)

de_sig_up <- sleuth_table %>% filter(qval < 0.01, lFC > 2)

de_sig_down <- sleuth_table %>% filter(qval < 0.01, lFC < -2)

#<-------------------VOCLCANO-------------------------->

sleuth_table$is_sig = sleuth_table$target_id %in% de_sig_df$target_id

ggplot(sleuth_table) + geom_point(aes(x = lFC, y = -log10(qval), col = is_sig), size = 2) + 
  geom_hline(yintercept = -log10(0.01), col = 'black', linetype = 'dashed') +
  geom_vline(xintercept = 2, col = 'black', linetype = 'dashed') + 
  geom_vline(xintercept = -2, col = 'black', linetype = 'dashed') +
  scale_color_manual(values = c("TRUE" = "red", 'FALSE' = "skyblue"), name = "Significance") +
  xlab(TeX('$log_2$(fold change)')) + ylab(TeX('$log_{10}$(q-value)')) + 
  theme_classic() + theme(axis.title = element_text(size = 16),
                          legend.text = element_text(size = 14),
                          legend.title = element_text(size = 16), 
                          axis.text = element_text(size = 15)) +
  scale_x_continuous(breaks = c(-10,-2,0,2,10))


#<--------------------WRITE-TABLES----------------------->

write.table(de_sig_df[1:50,], 'de.tsv', sep = '\t', quote=FALSE)

write.table(de_sig_up[1:50,], 'de_up.tsv', sep = '\t', quote=FALSE)

write.table(de_sig_down[1:50,], 'de_down.tsv', sep = '\t', quote=FALSE)
