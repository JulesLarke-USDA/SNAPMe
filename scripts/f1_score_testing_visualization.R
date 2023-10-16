# Test difference in distributions of F1 scores for single and multi ingredient meal images
library(ggplot2)
library(dplyr)
dat <- read.csv('snapme_inv_cook_f1_score.csv')
dat[dat$Occ_Name == "", ]

## all F1 scores
dat$all <- 'Total'
all_plot <- ggplot(dat, aes(x = all, y= F1_score, fill = all)) +
  geom_violin(draw_quantiles = 0.5) +
  labs(x = '', y = 'F1 score', title = '') +
  scale_fill_manual(values = c('aquamarine')) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 8, color = 'black'),
        axis.title = element_text(size = 8))
all_plot

#create even distribution across groups for comparison
sum(is.na(dat$food_count))
xtabs(~ food_count, data = dat)

#new food count bins
dat <- mutate(dat, fc_bin = "2-3")
dat$fc_bin[as.numeric(dat$food_count) == 1 ] <- "1"
dat$fc_bin[as.numeric(dat$food_count) > 3 ] <- "4+"

xtabs(~ fc_bin, data = dat)

fc_group_plot <-
  ggplot(dat, aes(x = fc_bin, y = F1_score, fill = fc_bin)) +
  geom_violin(draw_quantiles = 0.5) +
  labs(x = 'Number of FoodCodes', y = '', title = '') +
  scale_fill_manual(values = c('darkorchid', 'forestgreen', 'orange')) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8, color = 'black'),
    axis.title = element_text(size = 8))
fc_group_plot

dat %>% group_by(fc_bin) %>% summarise(mean_f1_bin = mean(F1_score))

#kruskal test
kruskal.test(x = dat$F1_score, g = dat$fc_bin)
#post-hoc test using dunn's test
fb_dunn_fc <- dunn.test::dunn.test(x = dat$F1_score, g = dat$fc_bin, method='bh')
fb_dunn_fc

######## evaluate results for IM2RECIPE ########
im2r <- read.csv('snapme_im2r_f1_score.csv')
im2r[im2r$Occ_Name == "", ]

## all F1 scores
im2r$all <- 'Total'
im2r_all_plot <-
  ggplot(im2r, aes(x = all, y = F1_score, fill = all)) +
  geom_violin(draw_quantiles = 0.5) +
  labs(x = '', y = 'F1 score', title = '') +
  scale_fill_manual(values = c('aquamarine')) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8, color = 'black'),
    axis.title = element_text(size = 8)
  )
im2r_all_plot

mean(im2r$F1_score)

#create more even distro across groups for comparison
sum(is.na(im2r$food_count))
xtabs(~ food_count, data = im2r)

#new food count bins
im2r <- mutate(im2r, fc_bin = "2-3")
im2r$fc_bin[as.numeric(im2r$food_count) == 1 ] <- "1"
im2r$fc_bin[as.numeric(im2r$food_count) > 3 ] <- "4+"

xtabs(~ fc_bin, data = im2r)

### binned by food count groups:
im2r_fc_group_plot <-
  ggplot(im2r, aes(x = fc_bin, y = F1_score, fill = fc_bin)) +
  geom_violin(draw_quantiles = 0.5) +
  labs(x = 'Number of FoodCodes', y = '', title = '') +
  scale_fill_manual(values = c('darkorchid', 'forestgreen', 'orange')) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8, color = 'black'),
    axis.title = element_text(size = 8)
  )
im2r_fc_group_plot

im2r %>% group_by(fc_bin) %>% summarise(mean_f1_bin = mean(F1_score))

#kruskal test
kruskal.test(x = im2r$F1_score, g = im2r$fc_bin)
#post-hoc test using dunn's test
im2r_fc_dunn <- dunn.test::dunn.test(x = im2r$F1_score, g = im2r$fc_bin, method='bh')
im2r_fc_dunn

grid = cowplot::plot_grid(all_plot, fc_group_plot, im2r_all_plot, im2r_fc_group_plot, nrow = 2, rel_widths = c(1,3), labels = "AUTO", label_size = 10, greedy = TRUE)
grid

ggsave("R_output/snapme_f1_violin_fb_inv_and_im2r.png",
       plot = grid,
       width = 5,
       height = 3,
       dpi = 1000)
