library("readxl")
library("tidyr")
library("purrr")
library("tidyverse")
library("ggplot2")
library("ggbeeswarm")

#Import raw data file
setwd('C:\\Users\\omeacock\\OneDrive - Université de Lausanne\\Analysis\\Growth curves\\16_02_23_Ct_BetaLactamaseAssays\\')
resultBlock = read_excel('resultsWorkable.xlsx',range='A4:F44',sheet='R formatting')

lookup <- c(bioRep = 'Biological replicate', techRep = 'Technical replicate', amp = 'Ampicillin', betLac = 'Beta-lactamase activity (milliunits)', density='Culture density (CFUs/ml)')
resultBlock <- resultBlock %>% rename(lookup)

resultBlock %>% ggplot(aes(x=Timepoint,y=density,colour=as.factor(amp),group=interaction(bioRep,amp),shape=as.factor(bioRep))) +
  geom_point() +
  geom_line() +
  theme_bw() +
  scale_colour_brewer(palette = "Set1") +
  scale_y_continuous(trans='log10',limits=c(1000,1e10)) +
  labs(x = 'Time (hours)', y = 'Culture density (CFUs ml^-1)')

case_adjBet <- function(amp,betLac) {
  case_when(
    amp == '+' & betLac < 0.053 ~ 0.015,
    amp == '-' & betLac < 0.053 ~ 0.03,
    TRUE ~ betLac
  )
}

resultBlock <- resultBlock %>% mutate(standBeta = case_adjBet(amp,betLac))

resultBlock %>% ggplot(aes(x=Timepoint,y=standBeta,colour=as.factor(amp),group=interaction(bioRep,amp),shape=as.factor(bioRep))) +
  geom_beeswarm() +
  theme_bw() +
  scale_colour_brewer(palette = "Set1") + 
  annotate('segment',x=0, xend=80 ,y=0.053, yend=0.053) +
  labs(x = 'Time (hours)', y = 'Beta-lactamase activity (milliunits)') 
