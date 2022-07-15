# library(FSA)
# 
source("code/load_data.R")



alpha.stats <- table %>% 
  as.data.frame() %>% 
  pivot_longer(-c(Taxon, featureid),
               names_to = "sampleid",
               values_to = "counts") %>% 
  group_by(sampleid) %>% 
  summarize(obs_richness = specnumber(counts),
            shannon_e = diversity(counts, 
                                  index = "shannon",
                                  base = exp(1)),
            simpson = diversity(counts, 
                                index = "simpson")) %>% 
  left_join(metadata, .)



ll_data <- alpha.stats %>% 
  filter(project == "LL") %>% 
  filter(keep == "Y")
scfa_data <- alpha.stats %>% 
  filter(project == "SCFA") %>% 
  filter(keep == "Y")



ll_obs_test <- ll_data %>%
  drop_na() %>%
  t_test(data = .,
         formula = obs_richness ~ group) %>%
  add_significance() %>%
  add_xy_position(x = "group")
ll_simpson_test <- ll_data %>%
  drop_na() %>%
  t_test(data = .,
         formula = simpson ~ group) %>%
  add_significance() %>%
  add_xy_position(x = "group")
ll_shannon_test <- ll_data %>%
  drop_na() %>%
  t_test(data = .,
         formula = shannon_e ~ group) %>%
  add_significance() %>%
  add_xy_position(x = "group")

scfa_obs_test <- scfa_data %>%
  drop_na() %>%
  t_test(data = .,
         formula = obs_richness ~ group) %>%
  add_significance() %>%
  add_xy_position(x = "group") %>% 
  mutate(y.position = y.position +10)
scfa_simpson_test <- scfa_data %>%
  drop_na() %>%
  t_test(data = .,
         formula = simpson ~ group) %>%
  add_significance() %>%
  add_xy_position(x = "group")
scfa_shannon_test <- scfa_data %>%
  drop_na() %>%
  t_test(data = .,
         formula = shannon_e ~ group) %>%
  add_significance() %>%
  add_xy_position(x = "group") %>% 
  mutate(y.position = y.position + 0.10)


plot.alpha <- function(alpha.stats, alpha_stat, stat.test, y_axis) {
  alpha.stats %>% 
    drop_na() %>% 
    ggbarplot(., x = "group", y = alpha_stat, fill = "group",
              add = c("mean"),
              position = position_dodge(0.8)) +
    geom_point() +
    
    stat_pvalue_manual(stat.test, label = "p.signif",
                       hide.ns = TRUE,
                       linetype = 1,
                       tip.length = 0,
                       bracket.size = 0.75,
                       label.size = 5) +
    
    scale_fill_manual(values = c("white", "blue")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+  
    labs(y = y_axis) +
    stat_summary(fun = "mean",
                 geom = "point", 
                 color = "black") +
    stat_summary(aes(width = 0.25),
                 fun = mean,
                 geom = "errorbar",
                 color = "black",
                 fun.max = function(x) mean(x) + sd(x),
                 fun.min = function(x) mean(x)) +
    
    theme_cowplot() +
    theme(axis.title.x = element_blank(),
          legend.position = "none",
          strip.background = element_blank(),
          axis.text = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold")) 
  
}

ll_obs_richness.plot <- plot.alpha(ll_data, "obs_richness", 
                                   ll_obs_test, "Observed Richness")
ll_simpson.plot <- plot.alpha(ll_data,"simpson", 
                              ll_simpson_test, "Simpson Index")
ll_shannon.plot <- plot.alpha(ll_data, "shannon_e", 
                              ll_shannon_test, "Shannon Index")

scfa_obs_richness.plot <- plot.alpha(scfa_data, "obs_richness", 
                                     scfa_obs_test, "Observed Richness") 
scfa_simpson.plot <- plot.alpha(scfa_data,"simpson", 
                                scfa_simpson_test, "Simpson Index")
scfa_shannon.plot <- plot.alpha(scfa_data, "shannon_e", 
                                scfa_shannon_test, "Shannon Index")

plot_grid(ll_obs_richness.plot,
          ll_simpson.plot,
          ll_shannon.plot,
          nrow = 1)

# ggsave("plots/alpha_stats_LL.png",
#        dpi = 600,
#        width = 9,
#        height = 3,
#        units = c("in"),
#        bg = "white")


plot_grid(scfa_obs_richness.plot,
          scfa_simpson.plot,
          scfa_shannon.plot,
          nrow = 1)

# ggsave("plots/alpha_stats_SCFA.png",
#        dpi = 600,
#        width = 9,
#        height = 3,
#        units = c("in"),
#        bg = "white")





## Not sure as to why I can't input variable into kw and dunn test. 
## Need to change alpha stat within function before running for loop 
## below. ARRRRG.
# kw_sum<- function(value){
#   
#   sum_stats <- alpha.stats %>% 
#     drop_na() %>% 
#     group_by(group) %>% 
#     summarize(mean = mean({{value}}),
#               sd = sd({{value}}),
#               min = min({{value}}),
#               max = max({{value}}))
#   
#   kw <- alpha.stats %>%
#     drop_na() %>%
#     kruskal.test(shannon_e~group, data = .)
#   
#   wi.pairwise <- pairwise.wilcox.test(x = alpha.stats$shannon_e,
#                                       g = alpha.stats$group,
#                                       p.adjust.method = "BH",
#                                       paired = F)
#   
#   output <- list(sum_stats, kw, wi.pairwise)
#   return(output)
# }


# obs = kw_sum(obs_richness)
# simpson = kw_sum(simpson)
# shannon = kw_sum(shannon_e)

# capture.output(obs, file = "plots/stats/richness_kw.tsv",
#                append = T)
# capture.output(simpson, file = "plots/stats/simpson_kw.tsv",
#                append = T)
# capture.output(shannon, file = "plots/stats/shannon_kw.tsv",
#                append = T)
# 
# 



