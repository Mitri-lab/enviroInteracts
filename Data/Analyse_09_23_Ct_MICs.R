library("readxl")
library("tidyr")
library("purrr")
library("tidyverse")
library("ggplot2")

anchor = c(35,3) #Top left-hand corner of the growthcurve data block
noWs = 96 #Number of wells
noExpts = 3 #Number of replicates

sheetNames = c('Pro0 0.5','Pro0 1','Pro0 2','Pro0 5')
proConcs = c(0.5,1,2,5)
fileName = '07_10_23_Ct_MICs.xlsx'

for (i in 1:4){
  fulldat <- read_excel(fileName,sheet = sheetNames[i])
  datBlock <- fulldat[anchor[1]:anchor[1]+1,anchor[2]:(anchor[2] + noWs - 1)]
  colnames(datBlock) <- fulldat[anchor[1],anchor[2]:(anchor[2] + noWs - 1)]

  #Convert to numerical long table format
  datBlockSubLong <- pivot_longer(datBlock, cols=1:noWs, names_to = "WellID", values_to = "OD600")
  datBlockSubLong$OD600 <- as.numeric(datBlockSubLong$OD600)

  datBlockSubLong <- datBlockSubLong %>% mutate(rowID = substring(WellID,1,1)) %>%
    mutate(colID = substring(WellID,2,nchar(WellID))) %>%
    mutate(techRep = recode(colID,'1'=1,'4'=1,'7'=1,'10'=1,'2'=2,'5'=2,'8'=2,'11'=2,'3'=3,'6'=3,'9'=3,'12'=3)) %>%
    mutate(startPro = proConcs[i]) %>%
    mutate(startAmp = recode(colID,'1'=0,'2'=0,'3'=0,'4'=10,'5'=10,'6'=10,'7'=20,'8'=20,'9'=20,'10'=30,'11'=30,'12'=30)) %>%
    mutate(MICamp = recode(rowID,'A'=0,'B'=400,'C'=800,'D'=1200,'E'=1600,'F'=2000,'G'=2400,'H'=NULL))

  if (i == 1) {
    datBlockLong <- datBlockSubLong 
  } else {
    datBlockLong <- datBlockLong %>% bind_rows(datBlockSubLong)
  }
}

#Correct ODs by subtracting mean blank well OD
meanBlk <- datBlockLong %>% filter(rowID == 'H') %>% summarise(mean(OD600))
meanBlk = meanBlk[[1]]
datBlockLong <- datBlockLong %>%
  mutate(correctOD = OD600 - meanBlk)

#For each technical replicate, find the first ampicillin concentration where the OD600 falls below 20% of the maximum
datBlockLongNest <- datBlockLong %>% group_by(techRep,startAmp,startPro) %>% 
  nest() %>% mutate(MICpos = map_dbl(data, ~ which.max(.x$correctOD < .x$correctOD[1]/5))) %>%
  mutate(MIC = map_dbl(data, ~ .x$MICamp[MICpos]))

#Format for heatmap presentation
datBlockLongNest <- datBlockLongNest %>%
  mutate(ampInd = recode(startAmp,'0'=1,'10'=2,'20'=3,'30'=4)) %>%
  mutate(proInd = recode(startPro,'0.5'=1,'1'=2,'2'=3,'5'=4)) %>%
  mutate(yCoord = (ampInd-1) * noExpts + techRep)

datBlockLongNest %>% ggplot(aes(y=yCoord,x=proInd)) +
  geom_tile(aes(fill=MIC), color = 'grey50') +
  scale_fill_distiller(palette='RdPu',limits=c(0,2400)) +
  labs(x='Pro',y='Amp') + 
  scale_y_discrete(breaks = c(1,2,3,4), labels = c(0, 10, 20, 30)) + 
  scale_x_discrete(breaks = c(1,2,3,4), labels = c(0.5, 1, 2, 5))
