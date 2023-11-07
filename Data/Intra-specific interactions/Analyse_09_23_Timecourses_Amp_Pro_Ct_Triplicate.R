library("readxl")
library("tidyr")
library("purrr")
library("tidyverse")
library("ggplot2")
library("zoo")
library("RColorBrewer")

anchor = c(40,4) #Top left-hand corner of the growthcurve data block
noTs = 220 #Number of timepoints
maxT = 100 #Maximum timepoint to show
dt = 0.5
noWs = 96 #Number of wells
startRat = 4 #Ratio of high-density to low-density ODs

#Import raw data file
fileList = c('14_09_23_Timecourses_Amp_Pro_Ct.xlsx','25_09_23_Timecourses_Amp_Pro_Ct.xlsx','01_10_23_Timecourses_Amp_Pro_Ct.xlsx')
noExpts = length(fileList)

#Generate concatenated table containing all data
for (i in 1:noExpts) {
  fulldat = read_excel(fileList[[i]])
  curveBlock = fulldat[anchor[1]:(anchor[1] + noTs - 1),anchor[2]:(anchor[2] + noWs - 1)]
  colnames(curveBlock) = fulldat[anchor[1]-1,anchor[2]:(anchor[2] + noWs - 1)]
  curveBlock$Time = seq(noTs) * dt
  
  #Convert to numerical long table format
  curveBlockSubLong <- pivot_longer(curveBlock, cols=1:noWs, names_to = "WellID", values_to = "OD600")
  curveBlockSubLong$OD600 <- as.numeric(curveBlockSubLong$OD600)
  
  #Annotate with experiment metadata
  if (i == 1) {
    curveBlockSubLong <- curveBlockSubLong %>% mutate(rowID = substring(WellID,1,1)) %>%
      mutate(colID = substring(WellID,2,nchar(WellID))) %>%
      mutate(techRep = recode(colID,'1'=1,'4'=1,'7'=1,'10'=1,'2'=2,'5'=2,'8'=2,'11'=2,'3'=3,'6'=3,'9'=3,'12'=3)) %>%
      mutate(bioRep = 1) %>%
      mutate(startOD = if_else(rowID %in% c('A','C','E','G'),0.04,0.01)) %>%
      mutate(startPro = recode(colID,'1'=0.5,'2'=0.5,'3'=0.5,'4'=1,'5'=1,'6'=1,'7'=2,'8'=2,'9'=2,'10'=5,'11'=5,'12'=5)) %>%
      mutate(startAmp = recode(rowID,'A'=0,'B'=0,'C'=10,'D'=10,'E'=20,'F'=20,'G'=30,'H'=30))
    
    curveBlockLong <- curveBlockSubLong
  } else if (i == 2) {
    curveBlockSubLong <- curveBlockSubLong %>% mutate(rowID = substring(WellID,1,1)) %>%
      mutate(colID = substring(WellID,2,nchar(WellID))) %>%
      mutate(techRep = recode(colID,'1'=1,'4'=1,'7'=1,'10'=1,'2'=2,'5'=2,'8'=2,'11'=2,'3'=3,'6'=3,'9'=3,'12'=3)) %>%
      mutate(bioRep = 2) %>%
      mutate(startOD = if_else(rowID %in% c('A','C','E','G'),0.01,0.04)) %>%
      mutate(startPro = recode(colID,'1'=0.5,'2'=0.5,'3'=0.5,'4'=1,'5'=1,'6'=1,'7'=2,'8'=2,'9'=2,'10'=5,'11'=5,'12'=5)) %>%
      mutate(startAmp = recode(rowID,'A'=0,'B'=0,'C'=10,'D'=10,'E'=20,'F'=20,'G'=30,'H'=30))
    
    curveBlockLong <- curveBlockLong %>% bind_rows(curveBlockSubLong)
  } else if (i ==3) {
    curveBlockSubLong <- curveBlockSubLong %>% mutate(rowID = substring(WellID,1,1)) %>%
      mutate(colID = substring(WellID,2,nchar(WellID))) %>%
      mutate(techRep = recode(colID,'1'=1,'4'=1,'7'=1,'10'=1,'2'=2,'5'=2,'8'=2,'11'=2,'3'=3,'6'=3,'9'=3,'12'=3)) %>%
      mutate(bioRep = 3) %>%
      mutate(startOD = if_else(rowID %in% c('A','C','E','G'),0.01,0.04)) %>%
      mutate(startPro = recode(rowID,'A'=0.5,'B'=0.5,'C'=1,'D'=1,'E'=2,'F'=2,'G'=5,'H'=5)) %>%
      mutate(startAmp = recode(colID,'1'=0,'2'=0,'3'=0,'4'=10,'5'=10,'6'=10,'7'=20,'8'=20,'9'=20,'10'=30,'11'=30,'12'=30))
    
    curveBlockLong <- curveBlockLong %>% bind_rows(curveBlockSubLong)
  }
}
  
#Group together equivalent conditions
curveBlockNest <- curveBlockLong %>% group_by(techRep,bioRep,startOD,startPro,startAmp) %>%
  nest()

#Remove background signal
ODpt04Starts <- curveBlockLong$startOD == 0.04 & curveBlockLong$Time == 0.5
ODpt01Starts <- curveBlockLong$startOD == 0.01 & curveBlockLong$Time == 0.5

#OD contributions from a 0.04 starting culture
bkgd <- (mean(curveBlockLong$OD600[ODpt04Starts]) - mean(curveBlockLong$OD600[ODpt01Starts]))/(1-1/startRat)

curveBlockLong <- curveBlockNest %>% 
  mutate(tgtOD = case_when(startOD == 0.04 ~ bkgd,startOD == 0.01 ~ bkgd/startRat)) %>%
  unnest(cols = c(data,tgtOD)) %>% group_by(techRep,bioRep,startOD,startPro,startAmp,WellID) %>% nest() %>%
  mutate(correctOD = map(data, ~ .x$OD600 - mean(.x$OD600[2:4]) + .x$tgtOD[1])) %>%
  unnest(cols = c(data,correctOD))

my_blu = brewer.pal(n = 9, "Blues")[c(3,5,7,8,9)]
my_grn = brewer.pal(n = 9, "Greens")[c(3,5,7,9)]

#Plot growthcurves for highest proline conc.
curveBlockLong %>% filter(startPro == 5, bioRep == 3) %>% 
  group_by(startAmp,startOD) %>%
  ggplot(aes(x=Time,y=correctOD,colour=factor(startAmp),linetype=factor(startOD))) +
  geom_line(aes(group=interaction(factor(startAmp),factor(startOD),factor(techRep))),size=1,alpha=0.4) +
  geom_smooth(span=0.25,size=2,aes(group=interaction(factor(startAmp),factor(startOD))),se=FALSE) +
  scale_colour_manual(values=my_blu) +
  scale_x_continuous(limits=c(0,maxT)) +
  #scale_y_continuous(limits=c(0,0.6))+
  theme_bw() +
  theme(legend.justification=c(0,1),legend.position=c(0.05, 0.95)) +
  labs(title='Ct, [Pro] = 5 \u3BCM',color = '[Amp]',linetype=bquote(OD[0]),x='Time (hr)', y=bquote(OD[600]))

#Plot growthcurves for lowest Amp conc.
curveBlockLong %>% filter(startAmp == 0, bioRep == 2) %>% 
  group_by(startPro,startOD) %>%
  ggplot(aes(x=Time,y=correctOD,colour=factor(startPro),linetype=factor(startOD))) +
  geom_line(aes(group=interaction(factor(startPro),factor(startOD),factor(techRep))),size=1,alpha=0.4) +
  geom_smooth(span=0.25,size=2,aes(group=interaction(factor(startPro),factor(startOD))),se=FALSE) +
  scale_colour_manual(values=my_grn) +
  scale_x_continuous(limits=c(0,maxT)) +
  theme_bw() +
  theme(legend.justification=c(0,1),legend.position=c(0.05, 0.95)) +
  labs(title='Ct, [Amp] = 0 \u3BCM',color = '[Pro]',linetype=bquote(OD[0]),x='Time (hr)', y=bquote(OD[600]))

#Calculate median population sizes across replicates
meanedData <- curveBlockLong %>%
  group_by(startAmp,startOD,startPro,Time,bioRep) %>% nest() %>%
  mutate(meanedOD = map_dbl(data,~ mean(na.omit(.x$correctOD)))) %>%
  select(startOD,startPro,startAmp,Time,meanedOD,bioRep)

#Calculate interactions between high and low starting ODs
Interacts <- meanedData %>% filter(startOD == 0.04) %>% 
  rename(OD_HiStart = meanedOD) %>%
  ungroup() %>% 
  select(startPro,startAmp,Time,OD_HiStart,bioRep) %>%
  mutate(OD_LoStart = pull(filter(ungroup(meanedData),startOD == 0.01),meanedOD)) %>%
  mutate(InteractionEst = OD_HiStart/startRat - OD_LoStart) %>%
  mutate(InteractRatEst = (OD_HiStart/startRat)/OD_LoStart)

#Create maximum positive and final interaction measurement heatmaps
Interacts <- Interacts %>% mutate(reliableInteractRats = ifelse(OD_HiStart>0.005 & OD_LoStart>0.005,InteractRatEst,NA)) #Eliminate noisy, early abundance ratios

interactsNest <- Interacts %>% group_by(startPro,startAmp,bioRep) %>% nest() %>%
  mutate(maxPos = map_dbl(data, ~ max(.x$InteractionEst,na.rm=TRUE))) %>% 
  mutate(maxPosWhen = map_dbl(data, ~ .x$Time[which.max(.x$InteractionEst)])) %>%
  mutate(finNeg = map_dbl(data, ~ .x$InteractionEst[.x$Time == maxT]))

#Format interactions for plotting in heatmaps
interactsNest <-interactsNest %>% 
  mutate(ampInd = recode(startAmp,'0'=1,'10'=2,'20'=3,'30'=4)) %>%
  mutate(proInd = recode(startPro,'0.5'=1,'1'=2,'2'=3,'5'=4)) %>%
  mutate(yCoord = (ampInd-1) * noExpts + bioRep)

sampVals = c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.425,0.45,0.475,0.4875,0.5,0.5125,0.525,0.55,0.575,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)

interactsNest %>% ggplot(aes(y=yCoord,x=proInd)) +
  geom_tile(aes(fill=-maxPos), color = 'grey50') +
  scale_fill_distiller(palette='BrBG',limits=c(-0.05,0.05),values=sampVals) +
  labs(x='Pro',y='Amp') + 
  scale_y_discrete(breaks = c(1,2,3,4), labels = c(0, 10, 20, 40)) + 
  scale_x_discrete(breaks = c(1,2,3,4), labels = c(0.5, 1, 2, 5))

interactsNest %>% ggplot(aes(y=yCoord,x=proInd)) +
  geom_tile(aes(fill=-finNeg), color = 'grey50') +
  scale_fill_distiller(palette='BrBG',limits=c(-0.4,0.4),values=sampVals, guide='colourbar') + 
  labs(x='Pro',y='Amp') + 
  scale_y_discrete(breaks = c(1,2,3,4), labels = c(0, 10, 20, 40)) + 
  scale_x_discrete(breaks = c(1,2,3,4), labels = c(0.5, 1, 2, 5))

#Plot interaction estimate timecourses, along with extreme values

interactsNest %>% unnest(cols = c(data,maxPos,maxPosWhen,finNeg)) %>% 
  filter(startPro == 5) %>% 
  ggplot(aes(x=Time, y=InteractionEst, colour=factor(startAmp), lineStyle=factor(bioRep))) +
  geom_line(size=1,alpha=0.4) +
  geom_smooth(span=0.25,size=1.5,aes(group=factor(startAmp)),se=FALSE) +
  scale_colour_manual(values=my_blu) +
  annotate("segment", x = 0, xend = maxT, y = 0, yend = 0, colour = "black", linetype = "dotted") +
  scale_y_continuous(limits=c(-0.4,0.05)) +
  scale_x_continuous(limits=c(0,maxT)) + 
  geom_point(aes(x=maxPosWhen,y=maxPos),size = 3,colour="purple",stat='unique') +
  geom_point(aes(x=maxPosWhen,y=maxPos,colour=factor(startAmp)),size=1.5,stat='unique') +
  geom_point(aes(x=maxT,y=finNeg),size = 3,colour="orange",stat='unique') +
  geom_point(aes(x=maxT,y=finNeg,colour=factor(startAmp)),size=1.5,stat='unique') +
  theme_bw() +
  theme(legend.position='none') +
  labs(title='Ct, [Pro] = 5 mM',x='Time',y='Abundance difference')

interactsNest %>% unnest(cols = c(data,maxPos,maxPosWhen,finNeg)) %>% 
  filter(startAmp == 0) %>% 
  ggplot(aes(x=Time, y=InteractionEst, colour=factor(startPro), lineStyle=factor(bioRep))) +
  geom_line(size=1,alpha=0.4) +
  geom_smooth(span=0.25,size=1.5,aes(group=factor(startPro)),se=FALSE) +
  scale_colour_manual(values=my_grn) +
  annotate("segment", x = 0, xend = maxT, y = 0, yend = 0, colour = "black", linetype = "dotted") +
  scale_y_continuous(limits=c(-0.4,0.05)) +
  scale_x_continuous(limits=c(0,maxT)) + 
  geom_point(aes(x=maxPosWhen,y=maxPos),size = 3,colour="purple",stat='unique') +
  geom_point(aes(x=maxPosWhen,y=maxPos,colour=factor(startPro)),size=1.5,stat='unique') +
  geom_point(aes(x=maxT,y=finNeg),size = 3,colour="orange",stat='unique') +
  geom_point(aes(x=maxT,y=finNeg,colour=factor(startPro)),size=1.5,stat='unique') +
  theme_bw() +
  theme(legend.position='none') +
  labs(title='Ct, [Amp] = 0 \u3BCM',x='Time',y='Abundance difference')

# Plot the idea between the abundance difference as an interaction measure
loOD = meanedData %>% filter(bioRep == 3, startPro == 5, startAmp == 30, startOD == 0.01)
hiOD = meanedData %>% filter(bioRep == 3, startPro == 5, startAmp == 30, startOD == 0.04)

plot(loOD$Time,loOD$meanedOD,type='l',xlim=c(0,maxT),lwd=2)
lines(hiOD$Time,hiOD$meanedOD,lwd=2,lty=2)
lines(hiOD$Time,hiOD$meanedOD/4,lwd=2,lty=2,col='grey')

plot(hiOD$Time,(hiOD$meanedOD/4)-loOD$meanedOD,type='l',xlim=c(0,maxT),lwd=2)
lines(c(0,maxT),c(0,0), lty=3)

# Plot endpoint abundances (at 100 hours) as a function of proline concentration
finalODs <- curveBlockLong %>% filter(Time==maxT) %>%
  group_by(startPro,startAmp) %>% nest() %>% 
  mutate(meanOD = map_dbl(data, ~mean(.x$correctOD))) %>%
  mutate(sdOD = map_dbl(data, ~sd(.x$correctOD)))

origFit <- lm(meanOD ~ 0 + startPro,data = finalODs)

finalODs %>% ggplot(aes(x=startPro,y=meanOD,color=startAmp,group=startAmp)) +
  geom_pointrange(aes(ymin=meanOD-sdOD, ymax=meanOD+sdOD), size=0.3)  +
  scale_x_continuous(limits=c(0,5)) +
  scale_color_gradient(low='gray80',high='blue') +
  theme_bw() +
  annotate('segment',x=0,xend=5,y=0,yend=coef(origFit)*5,color='black',linetype='dashed') +
  labs(title='Ct yields at experiment end',x='Pro_0',y='Yield')