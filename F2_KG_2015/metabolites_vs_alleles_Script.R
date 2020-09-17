library(tidyverse)

####################
# Leaf metabolites #
####################

# Import the avaialble genotyping data 
df.alleles <- read.delim(file = "F2_KG_2015/ShZO_allele_genotype.txt", header = TRUE, sep = "\t")


# Import metabolite data
df.metabolites <- read.delim(file = "F2_KG_2015/leafwash_ug_mg_tissue.txt",
                             header = TRUE, sep = "\t", check.names = FALSE) %>%
  pivot_longer(cols = '7epiZ':'9H10epoZ',
               names_to = "metabolite",
               values_to = "value")
df.metabolites$metabolite <- factor(df.metabolites$metabolite, 
                                    levels = c("7epiZ", "9HZ", "9H10epoZ"),
                                    ordered = TRUE)

# Summarise the metabolite data to get mean values and standard error
sum.metabolites <- df.metabolites %>%
  dplyr::group_by(genotype, metabolite) %>% 
  dplyr::summarise(mean_value = mean(value), se_metabolites = sd(value)/sqrt(n()))

# Join df's: keep the genotypes of which the genotypes are available
df.metabolites.alleles <- right_join(sum.metabolites, df.alleles, by = "genotype")

sum.metabolites$metabolite <- factor(sum.metabolites$metabolite, 
                                    levels = c("7epiZ", "9HZ", "9H10epoZ"),
                                    ordered = TRUE)



########
# Plot #
########
p.metabolites.alleles =
df.metabolites.alleles %>% #filter(!genotype %in% c("C32", "PI127826")) %>%
ggplot(aes(x = reorder(genotype, mean_value), y = mean_value, fill= ShZO_allele))+
  geom_bar(stat = "identity") +
  geom_errorbar(aes(x = genotype, ymin = mean_value - se_metabolites, ymax = mean_value + se_metabolites), width = 0.3)+
  geom_text(aes(label = round(mean_value, digits = 2), y = mean_value+0.04))+
  facet_wrap(~metabolite, scale = "free", ncol = 1) +
  xlab("F2-genotype")+
  ylab("Metabolite level (ug / mg fresh-leaf weight)")+
  scale_fill_manual("ShZO allele", values = c("homozygous_ShZO" = "black", "homozygous_Solyc01g008670" = "red"))+
  theme_bw()
  
  ggsave(file = "F2_KG_2015/metabolites_and_ShZO_alleles.pdf", plot = p.metabolites.alleles, height = 6, width = 6)
  

#####################
# Trichome extracts #
#####################

df.trichomes <- read.delim(file = "F2_KG_2015/terpene_terpenes_ng_mg_trichomes.txt", header = TRUE, sep = "\t", check.names = FALSE)
  sum.trichomes <- df.trichomes %>%
    dplyr::group_by(genotype, metabolite) %>% 
    dplyr::summarise(mean_value = mean(value), se_metabolites = sd(value)/sqrt(n()))
  
  sum.trichomes$metabolite <- factor(sum.trichomes$metabolite, 
                                       levels = c("7epiZ", "9HZ", "9H10epoZ"),
                                       ordered = TRUE)
  
  # Join df's: keep the genotypes of which the genotypes are available
  df.trichomes.alleles <- right_join(sum.trichomes, df.alleles, by = "genotype")
  
  p.trichomes.alleles =
    df.trichomes.alleles %>% #filter(!genotype %in% c("C32", "PI127826")) %>%
    ggplot(aes(x = factor(genotype, level = c("C32", "18", "35", "14", "34", "30", "PI127826")), 
               y = mean_value, 
               fill= P450_allele))+
    geom_bar(stat = "identity") +
    geom_errorbar(aes(x = genotype, ymin = mean_value - se_metabolites, ymax = mean_value + se_metabolites), width = 0.3)+
    geom_text(aes(label = round(mean_value, digits = 2), y = mean_value+0.5))+
    facet_wrap(~metabolite, scale = "free", ncol = 1) +
    xlab("F2-genotype")+
    ylab("Metabolite level (ng / mg trichomes)")+
    scale_fill_manual("P450 allele", values = c("homozygous_ShZO" = "black", "homozygous_Solyc01g008670" = "red"))+
    theme_bw()
  
  ggsave(file = "F2_KG_2015/trichomes_and_ShZO_alleles.pdf", plot = p.trichomes.alleles, height = 6, width = 6)
  
  

