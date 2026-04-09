library(tidyverse)
library(ggpubr)
library(readxl)
library(car)
library(dunn.test)
library(RRPP)
library(flexplot)


faba <- read_excel("../input/Biomass.xlsx", 
                   sheet = "fababean") %>%
  dplyr::rename(treatment = GBH_P) %>%
  mutate(treatment = str_replace(treatment, "KO",
                                 "C")) %>%
  mutate(treatment = str_replace(treatment, "KP",
                                 "P")) %>%
  mutate(treatment = str_replace(treatment, "GO",
                                 "G")) %>%
  mutate(treatment = str_replace(treatment, "GP",
                                 "GP")) %>%
  dplyr::rename("Treatment" = "treatment")



potatoes <- read_excel("../input/Biomass.xlsx", 
                       sheet = "potato") %>%
  dplyr::rename(treatment = GBH_P) %>%
  mutate(treatment = str_replace(treatment, "KO",
                                 "C")) %>%
  mutate(treatment = str_replace(treatment, "KP",
                                 "P")) %>%
  mutate(treatment = str_replace(treatment, "GO",
                                 "G")) %>%
  mutate(treatment = str_replace(treatment, "GP",
                                 "GP")) %>%
  rename("Weight" = "weight") %>%
  rename("Treatment" = "treatment")


oat <- read_excel("../input/Biomass.xlsx", 
                  sheet = "oat") %>%
  dplyr::rename(treatment = GBH_P) %>%
  mutate(treatment = str_replace(treatment, "KO",
                                 "C")) %>%
  mutate(treatment = str_replace(treatment, "KP",
                                 "P")) %>%
  mutate(treatment = str_replace(treatment, "GO",
                                 "G")) %>%
  mutate(treatment = str_replace(treatment, "GP",
                                 "GP")) %>%
  rename("Weight" = "weight") %>%
  rename("Treatment" = "treatment")


leveneTest(Weight ~ as.factor(Treatment), data = potatoes)
leveneTest(Weight ~ as.factor(Treatment), data = faba)
leveneTest(Weight ~ as.factor(Treatment), data = oat)

##### GLM #####

test_data <- potatoes %>%
  select(Ptrt, Gtrt, Weight) %>%
  mutate(Ptrt = ifelse(Ptrt == "P", TRUE, FALSE),
         Gtrt = ifelse(Gtrt == "G", TRUE, FALSE)) %>%
  rename("phosphate" = "Ptrt",
         "glyphosate" = "Gtrt")

test_pot <- glm(formula = Weight ~ glyphosate*phosphate,
             family = Gamma(link="log"),
             data = test_data)

p_sum <- as.data.frame(summary(test_pot)$coefficients) %>%
  rownames_to_column("Treatment") %>%
  mutate(plant = "Potato") %>%
  dplyr::relocate(plant)
  

test_data <- faba %>%
  select(Ptrt, Gtrt, Weight) %>%
  mutate(Ptrt = ifelse(Ptrt == "P", TRUE, FALSE),
         Gtrt = ifelse(Gtrt == "G", TRUE, FALSE)) %>%
  rename("phosphate" = "Ptrt",
         "glyphosate" = "Gtrt")

test_faba <- glm(formula = Weight ~ glyphosate*phosphate,
                family = Gamma(link="log"),
                data = test_data)

f_sum <- as.data.frame(summary(test_faba)$coefficients) %>%
  rownames_to_column("Treatment") %>%
  mutate(plant = "Faba bean") %>%
  dplyr::relocate(plant)


test_data <- oat %>%
  select(Ptrt, Gtrt, Weight) %>%
  mutate(Ptrt = ifelse(Ptrt == "P", TRUE, FALSE),
         Gtrt = ifelse(Gtrt == "G", TRUE, FALSE)) %>%
  rename("phosphate" = "Ptrt",
         "glyphosate" = "Gtrt")

test_oat <- glm(formula = Weight ~ glyphosate*phosphate,
                 family = Gamma(link="log"),
                 data = test_data)

o_sum <- as.data.frame(summary(test_oat)$coefficients) %>%
  rownames_to_column("Treatment") %>%
  mutate(plant = "Oat") %>%
  dplyr::relocate(plant)

# create a table

all_summ <- bind_rows(p_sum, f_sum, o_sum)

write_tsv(all_summ, "../output/biomass_GLM_coefficients.tsv")

##### wilcoxon test #####
pot_test <- ggpubr::compare_means(Weight ~ Treatment, data = potatoes)
faba_test <- ggpubr::compare_means(Weight ~ Treatment, data = faba)
oat_test <- ggpubr::compare_means(Weight ~ Treatment, data = oat)

test_p <- potatoes %>%
  group_by(Treatment) %>%
  summarise(n = n()) %>% ungroup()

test_f <- faba %>%
  group_by(Treatment) %>%
  summarise(n = n()) %>% ungroup()

test_o <-  oat %>%
  group_by(Treatment) %>%
  summarise(n = n()) %>% ungroup()

# adjust the significance line heights
pot_test <- pot_test %>% mutate(y.position = c(31, 34, 37, 40, 43, 46))


faba_test <- faba_test %>% mutate(y.position = c(77, 77, 87, 87, 95, 127))


oat_test <- oat_test %>% mutate(y.position = c(25, 25, 25, 25, 25, 25))


potatoes$Treatment <- factor(potatoes$Treatment, levels = c("C", "G", "GP", "P"))
faba$Treatment <- factor(faba$Treatment, levels = c("C", "G", "GP", "P"))
oat$Treatment <- factor(oat$Treatment, levels = c("C", "G", "GP", "P"))


a <- ggboxplot(potatoes, x = "Treatment", y = "Weight") + 
  stat_pvalue_manual(pot_test, label = "p.adj", hide.ns = T) +
  ylab("Dry biomass (g)")

b <- ggboxplot(faba, x = "Treatment", y = "Weight") + 
  stat_pvalue_manual(faba_test, label = "p.adj", hide.ns = T) +
  ylab("Dry biomass (g)")

c <- ggboxplot(oat, x = "Treatment", y = "Weight") + 
  stat_pvalue_manual(oat_test, label = "p.adj", hide.ns = T) +
  ylab("Dry biomass (g)")

ggarrange(a, b,c, ncol = 3, labels = c("a", "b", "c"))

ggsave(path = "../output", filename = "biomass_boxplot.pdf", width = 10, height = 5, device='pdf', dpi=300)

