library(tidyverse)

# Costum theme for plotting
my.theme = 
  theme(text = element_text(),
        axis.text.x = element_text(size = 8, angle = 0, colour = "black"), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.text.x = element_text(size=8, colour = "black")
  )+
  theme_bw()

#############################
# Data import & subsetting  #
#############################

df = read.table(file = "Fig_6_toxicity_assay/wf_assays.txt", header = TRUE, check.names = FALSE, sep = "\t")

df.sub = df %>% select(metabolite, dose_ug_cm2, relative_survival_to_mock) #%>% filter(dose_ug_cm2 != "mock")

##########################
# Summarise data and plot#
##########################
p1 = 
df.sub %>%
  dplyr::group_by(metabolite, dose_ug_cm2) %>% summarize(mean_survival = mean(relative_survival_to_mock),
                                                         se = sd(relative_survival_to_mock)/sqrt(5)) %>%
  mutate(dose_ug_cm2 = fct_relevel(dose_ug_cm2, 
                            "mock", "1", "10", "20", "50")) %>%
  
  ggplot(., aes(x = dose_ug_cm2, y = mean_survival))+
    geom_line(aes(x = dose_ug_cm2, y = mean_survival, group = metabolite, color = metabolite), size = 1)+
    geom_point(aes(x = dose_ug_cm2, y = mean_survival, group = metabolite, color = metabolite, shape = metabolite), size = 3)+
  scale_colour_manual(values = c("7epiZ" = "grey", "9HZ" = "red", "9H10epoZ" = "black"))+
  scale_shape_manual(values = c("7epiZ" = 17, "9HZ" = 19, "9H10epoZ" = 5))+
    geom_errorbar(aes(x = dose_ug_cm2, 
                      ymin = mean_survival - se, 
                      ymax = mean_survival + se, 
                      group = metabolite,
                      color = metabolite),
                  width = 0.2)+
  labs(x = "Applied Dose (ug/cm2)", y = "B. tabaci survival rate (relative to mock)")+
  my.theme

ggsave(file = "Fig_6_toxicity_assay/plots/dose_assay_wf.pdf", plot = p1, width = 5, height = 4)

##############
# Statistics #
##############

# Perform ANOVA on each treatment to see if concentrations give differnt survival
# Create a list by subsetting eacht treatment, then do the ANOVA on each of them
anova <- lapply(split(df.sub, 
             df.sub$metabolite), 
       function(d) {aov(relative_survival_to_mock ~ dose_ug_cm2, data = d) })

# Perform post-hoc test per metabolite and write to file
write.table(as.data.frame(TukeyHSD(anova$'7epiZ')$dose_ug_cm2),
            file = "Fig_6_toxicity_assay/statistics/Tukey_HSD_zingiberene.txt", sep = "\t", row.names = TRUE)

write.table(as.data.frame(TukeyHSD(anova$'9HZ')$dose_ug_cm2),
            file = "Fig_6_toxicity_assay/statistics/Tukey_HSD_9HZ.txt", sep = "\t", row.names = TRUE)

write.table(as.data.frame(TukeyHSD(anova$'7epiZ')$dose_ug_cm2),
            file = "Fig_6_toxicity_assay/statistics/Tukey_HSD_9H10epoZ.txt", sep = "\t", row.names = TRUE)

