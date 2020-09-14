library(tidyverse)
library(ggpubr)

# Import whitefly survival data
df.wf <- read.delim(file = "F2_KG_2015/wf_clipcage_assays.txt",
                   header= TRUE, sep = "\t")
df.wf$plant <- as.factor(df.wf$plant)

# Import metabolite data
df.metabolites <- read.delim(file = "F2_KG_2015/leafwash_ug_mg_tissue.txt",
                             header = TRUE, sep = "\t") %>%
  pivot_longer(cols = Zingiberene:epoxy_zingiberenol,
               names_to = "metabolite",
               values_to = "value")
df.metabolites$metabolite <- factor(df.metabolites$metabolite, 
                                    levels = c("Zingiberene", "Zingiberenol", "epoxy_zingiberenol"),
                                              ordered = TRUE)

####################
# Plot metabolites #
####################

df.metabolites %>%
  dplyr::group_by(genotype, metabolite) %>% 
  dplyr::summarise(mean_value = mean(value), se_metabolites = sd(value)/sqrt(n())) %>%
#  filter(!genotype %in% c("LA2167", "PI127826")) %>%
  
  ggplot()+
  geom_bar(aes(x = genotype, y = mean_value),
           stat = "identity") +
  geom_errorbar(aes(x = genotype, ymin = mean_value - se_metabolites, ymax = mean_value + se_metabolites), width = 0.3)+
  facet_wrap(~metabolite, scale = "free", ncol = 1) +
  theme_bw()

##########################
# plot whitefly survival #
##########################

df.wf$survival <- round((df.wf$alive/(df.wf$alive+df.wf$dead))*100)
df.wf %>% 
  dplyr::group_by(genotype) %>%
  dplyr::summarise(mean_survival = mean(survival), se_survival = sd(survival)/sqrt(n())) %>%
  arrange(mean_survival) %>%
  
  ggplot()+
  geom_bar(aes(x = reorder(reorder(genotype, mean_survival), mean_survival), y = mean_survival),
           stat = "identity") + 
  geom_errorbar(aes(x = genotype, ymin = mean_survival - se_survival, ymax = mean_survival + se_survival), width = 0.3)+
  theme_bw()

############################
# Metabolites vs suurvival #
############################

sum.survival <- df.wf %>% 
  dplyr::group_by(genotype) %>%
  dplyr::summarise(mean_survival = mean(survival), se_survival = sd(survival)/sqrt(n())) %>%
  select(-se_survival)

sum.metabolites <- df.metabolites %>%
  dplyr::group_by(genotype, metabolite) %>% 
  dplyr::summarise(mean_value = mean(value), se_metabolites = sd(value)/sqrt(n())) %>%
  select(-se_metabolites)

survival.metabolites <- inner_join(sum.survival, sum.metabolites, by = "genotype")

########
# Plot #
########

survival.metabolites %>%
ggplot(aes(x = mean_value,
         y = mean_survival)) +
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~metabolite, ncol = 1, scale = "free") +
  stat_cor(
    method = "pearson",
    label.x = 0.15,
    label.x.npc = 0, 
    label.y.npc = 0.9)+
  geom_text(label = survival.metabolites$genotype, vjust = 1, hjust = 1)
