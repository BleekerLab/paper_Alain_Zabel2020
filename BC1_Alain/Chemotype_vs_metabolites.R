library(tidyverse)

########################
# Load metabolite data #
########################
df.metabolites <- read.delim(file = "BC1_Alain/BC1_chemotypes.txt", sep = "\t",
                             check.names = FALSE) %>%
  select(genotype, chemotype)

######################
# Load whitefly data #
######################

df.wf <- read.delim(file = "BC1_Alain/BC1_wf_assay.txt", sep = "\t")
df.wf$survival <- round((df.wf$alive/(df.wf$alive+df.wf$dead))*100)
df.wf <- df.wf %>% select(genotype, survival)
wf.sum <- df.wf %>% group_by(genotype) %>% dplyr::summarise(mean_survival = mean(survival))


levels = c("7epiZ", "9HZ", "9H10epoZ")

df.metabolites.wf <- inner_join(df.metabolites, wf.sum, by = "genotype")

########
# Plot #
########
p.wf = 
df.metabolites.wf %>%
  ggplot(aes(x = chemotype, y = mean_survival))+
  geom_boxplot()+
  geom_point()+
  xlab("Chemotype")+
  ylab("B.tabaci survival (%)")+
  scale_x_discrete(labels = c("m" = "monoterpenes",
                              "z" = "7epiZ",
                              "z_d" = "7epiZ + derivatives")
                   )+
  theme_bw()

ggsave(filename = "BC1_Alain/wf_survival_per_chemotype.pdf", plot = p.wf, width = 3.5, height = 3)

##############
# Statistics #
##############

anova <- aov(data = df.metabolites.wf, mean_survival ~ chemotype)
summary(anova)


df.metabolites.wf %>%
  ggplot(aes(x = mean_survival)) +
  geom_density() +
  facet_wrap(~chemotype)
