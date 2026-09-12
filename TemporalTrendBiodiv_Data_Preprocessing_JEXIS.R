# Code for the paper:
#Meng, B., et al.
# Temporal decoupling of biodiversity effects on grassland productivity and stability
#
# This script preprocesses the data and quantifies the net biodiversity effect .
#
# Raw data downloaded from: https://jexis.idiv.de/
# Includes:
#   - Aboveground plant biomass from main plots in the Jena Experiment (2003–2021)
#      ID:38, 40,41,124,127,176,320
#   - Aboveground plant biomass from the dBEF Experiment (2017–2021)
#      ID:12,71,142,143,380,381,512
#   - Precipitation, temperature, and SPEI data for the Jena Experiment
#      ID:484,507

#
# Reproducibility information
# This script was run successfully under the following R session:
#
#R version 4.6.0 (2026-04-24 ucrt)
#Platform: x86_64-w64-mingw32/x64
#Running under: Windows 11 x64 (build 26200)

##contact: Bo Meng
##bomeng@cau.edu.cn




#The function used to calculate the biodiversity effect of interspecific interaction on community variability
#matrix input
#X_mix is the observed community
#X_mono is the hypothetical community; note that, such as in grassland experiments, the hypothetical species biomass should be miu_mono = 1/n*miu_mix；n is Species richness
#columns is time series, rows is species
variability_partition<- function(X_mono,X_mix){
  
  X_mono<- as.matrix(X_mono)
  X_mix<- as.matrix(X_mix)
  
  CV.mix = sd(rowSums(X_mix))/mean(rowSums(X_mix))
  CV.mono = sd(rowSums(X_mono))/mean(rowSums(X_mono))
  # average biomass of community 
  miu.mix = mean(rowSums(X_mix))
  miu.mono = mean(rowSums(X_mono))
  
  tmp<- data.frame()
  for (i in 1:ncol(X_mix)) {
    CV.mix.i = sd(X_mix[,i])/mean(X_mix[,i])
    CV.mono.i = sd(X_mono[,i])/mean(X_mono[,i])
    P.mix.i = mean(X_mix[,i])/miu.mix
    P.mono.i = mean(X_mono[,i])/miu.mono 
    delta.CV = CV.mix.i-CV.mono.i
    delta.P = P.mix.i-P.mono.i
    tmp<- rbind(tmp, c(CV.mix.i,CV.mono.i,P.mix.i,P.mono.i,delta.CV,delta.P))
  }
  #statistic Base line
  CV.base = mean(tmp[,2])
  # partitioning in species variability  
  CVS.mix = sum(tmp[,1]*tmp[,3], na.rm = T)#species variability in mixture (observed community)
  CVS.mono = sum(tmp[,2]*tmp[,4], na.rm = T)#species variability in monoculture (expected community)
  
  AE.SpVar = sum(tmp[,5]*tmp[,3], na.rm = T)/CVS.mono
  SE.SpVar = sum(tmp[,2]*tmp[,6], na.rm = T)/CVS.mono
  cov.mix <- cov(X_mix)
  cov.mono <- cov(X_mono)
  cor.mix <- cor(X_mix)
  cor.mono <- cor(X_mono)
  delta.cor<- cor.mix-cor.mono
  
  # partitioning in species synchrony   
  phi.mix = sum(cov.mix)/(sum(sqrt(diag(cov.mix))))^2 #species synchrony in mixture (observed community); there is the φ^2
  phi.mono = sum(cov.mono)/(sum(sqrt(diag(cov.mono))))^2 #species synchrony in monoculture (expected community); there is the φ^2
  
  omega.mix<- (cov.mix/cor.mix)/(sum(sqrt(diag(cov.mix))))^2
  omega.mono<- (cov.mono/cor.mono)/(sum(sqrt(diag(cov.mono))))^2
  omega.mix[is.na(omega.mix)]<-0
  omega.mono[is.na(omega.mono)]<-0
  delta.omega<- omega.mix - omega.mono
  
  SE.Syn = sum(cor.mono * delta.omega, na.rm =T) / phi.mono
  AE.Syn = sum(delta.cor * omega.mix, na.rm =T) / phi.mono
  
  return(c(CV.o=CV.mix, CV.e=CV.mono, CV.b= CV.base, CVS.o=CVS.mix, CVS.e=CVS.mono, phi.o=sqrt(phi.mix), phi.e= sqrt(phi.mono),
           AE.SpVar = AE.SpVar, SE.SpVar = SE.SpVar, AE.Syn = AE.Syn, SE.Syn = SE.Syn))
  # CV.o is the variability of observed community
  # CV.e is the variability of expected community without interspecific interaction; ecological reference
  # CV.b or CV.null is the average monoculture variability across component species; Statistical baseline
  # CVS.o is the species variability in observed community
  # CVS.e is the species variability in expected community
  # phi.o is the species synchrony in observed community
  # phi.e is the species synchrony in expected community
  # AE.SpVar is the average effect of interspecific interaction on species variability
  # SE.SpVar is the selection effect of interspecific interaction on species variability
  # AE.Syn is the average effect of interspecific interaction on species synchrony
  # SE.Syn is the selection effect of interspecific interaction on species synchrony
}


### The function used to calculate the biodiversity effect on community biomass
### we used an alternative definition of net biodiversity effect by the ratio of mixture and monoculture biomass, see Method
### n is the Species richness
BEF<- function(X_mono,X_mix,n){
  NBE.F= sum( (colMeans(X_mix)/colMeans(X_mono)-1) * colMeans(X_mono))/sum(colMeans(X_mono))+1
  SE =   cov( (colMeans(X_mix)/colMeans(X_mono)-1) , colMeans(X_mono))*(n-1)/sum(colMeans(X_mono))
  CE =   sum(colMeans(X_mix)/colMeans(X_mono)-1) * sum(colMeans(X_mono))/n/sum(colMeans(X_mono))
  return(c(NBE.F=NBE.F,SE=SE,CE=CE))
}

# For details on data filtering, see the Methods section of:
#   Meng, B., et al. (2026). Temporal decoupling of biodiversity effects on grassland productivity and stability.
#   - Aboveground plant biomass from main plots in the Jena Experiment (2003–2021)
#      ID:38, 40,41,124,127,176,320
jena_data<- read.csv("...JenaMainExp(filtered).csv")
#   - Precipitation, temperature, and SPEI data for the Jena Experiment
#      ID:484,507
climate <- read.csv("...JenaClimateData.csv")
#   - Aboveground plant biomass from the dBEF Experiment (2017–2021)
#      ID:12,71,142,143,380,381,512
dBEF_data<-  read.csv("...csv")
# ----------------------------
# Quantify the impacts of biodiversity on ecosystem stability and functioning, and partition the contributions of interspecific interaction
# main plots in the Jena Experiment 
# ----------------------------
PlantedN<- jena_data%>%
  group_by(plotcode)%>%
  summarise(Richness = mean(Richness))

Jena_time<- data.frame()
plotcode<- unique(jena_data$plotcode)
for (i in 1:length(plotcode)) {
  for (j in 0:14) {
    mix.i<- jena_data[jena_data$plotcode ==plotcode[i],]%>%
      ungroup()%>%
      filter(year %in% c((j+2003):(j+2007)))%>%
      select(year,species_other,AGB.mix)%>%
      spread(species_other,AGB.mix)%>%
      select(-year)
    mix.i[is.na(mix.i)]<- 0
    mono.i<- jena_data[jena_data$plotcode == plotcode[i],]%>%
      ungroup()%>%
      filter(year %in% c((j+2003):(j+2007)))%>%
      select(year,species_other,AGB.mono)%>%
      spread(species_other,AGB.mono)%>%
      select(-year)
    
    cov.mono<- cov(colMeans(mono.i), apply(mono.i, 2, function(x){sd(x)/mean(x)}))
    miu_mix<- sum(colMeans(mix.i))
    n <- ncol(mix.i)
    miu_mono<- sum(colMeans(mono.i))
    block<- str_split_fixed(plotcode[i],  "A", n =2)[,1]
    ## These functions are implemented in `PartitioningFramework.R`
    a<- variability_partition(mono.i, mix.i)
    CESE<- BEF(mono.i, mix.i,n)
    
    evenness_mix<- diversity(colMeans(mix.i),index = "shannon")/log(n)
    evenness_mono<- diversity(colMeans(mono.i),index = "shannon")/log(n)
    
    CorS_mix<- mean(cor(mix.i)[upper.tri(cor(mix.i))], na.rm = T)
    CorS_mono<- mean(cor(mono.i)[upper.tri(cor(mono.i))], na.rm = T)
    
    climate.j<- climate%>%
      filter(year %in% c((j+2003):(j+2007)))%>%
      summarise(T_mean = mean(T_air),
                SPEI_mean = mean(spei_p),
                P_mean = mean(rain),
                SPEI_Lowest = min(spei_p) )
    climate.j<- as.matrix(climate.j)[1,]
    Order<- j+1
    plot<- plotcode[i]
    year_begin<- (j*2+2003)
    year_end<- (j*2+2008)
    Jena_time<- rbind(Jena_time, c(block,plot,n,Order,climate.j,miu_mix,miu_mono, evenness_mix,evenness_mono,CorS_mono,cov.mono, CESE, a)) 
    
  }}
colnames(Jena_time)<-c("block","plot","n","Order","T_mean","SPEI_mean","P_mean","Driest","miu_mix","miu_mono", "evenness_mix","evenness_mono", "CorS_mono","cov.mono","NBE.F","SE","CE", "CV.mix", "CV.mono","CV.avg", "CVS.mix", "CVS.mono", "phi.mix", "phi.mono","AE.SpVar", "SE.SpVar", "AE.Syn", "SE.Syn")
Jena_time[,-c(1:2)]<-apply(Jena_time[,-c(1:2)],2,as.numeric)
leaveplot<- Jena_time%>%
  group_by(plot)%>%
  summarise(n = mean (n),
            NBE.F = mean(NBE.F))%>%
  filter(n>1)%>%
  na.omit()

Jena_time<-Jena_time%>%
  filter(plot %in% leaveplot$plot)%>%
  left_join(PlantedN, by = c("plot" = "plotcode"))

Jena_time$NIIE <-  Jena_time$CV.mono/Jena_time$CV.avg
Jena_time$NIDE <- Jena_time$CV.mix/Jena_time$CV.mono
Jena_time$NBE <-  Jena_time$CV.mix/Jena_time$CV.avg
Jena_time$NIDE.SpVar <- Jena_time$CVS.mix/Jena_time$CVS.mono
Jena_time$NIDE.Syn <- Jena_time$phi.mix/Jena_time$phi.mono
Jena_time$NIIE.Syn<- Jena_time$phi.mono
Jena_time$NIIE.SpVar<- Jena_time$CVS.mono/Jena_time$CV.avg



# ----------------------------
# Quantify the impacts of biodiversity on ecosystem stability and functioning, and partition the contributions of interspecific interaction
# dBEF Experiment 
# ----------------------------
dBEF_time<- data.frame()  
plotcode<- unique(Jena_time$plot)
for (i in 1:length(plotcode)) {
  for (j in c("D1","D3")) {
    
    mix.i<- dBEF_data[dBEF_data$plotcode ==plotcode[i],]%>%
      filter(treatment == j)%>%
      ungroup()%>%
      select(year,species_other,AGB.mix)%>%
      spread(species_other,AGB.mix)%>%
      select(-year)
    mono.i<- dBEF_data[dBEF_data$plotcode == plotcode[i],]%>%
      filter(treatment == j)%>%
      ungroup()%>%
      select(year,species_other,AGB.mono)%>%
      spread(species_other,AGB.mono)%>%
      select(-year)
    mix.i<- mix.i[,colMeans(mono.i)>0]
    mono.i<- mono.i[,colMeans(mono.i)>0]
    miu_mix<- sum(colMeans(mix.i))
    n <- ncol(mix.i)
    miu_mono<- sum(colMeans(mono.i))
    #block<- str_split_fixed(plotcode[i],  "A", n =2)[,1]
    a<- variability_partition(mono.i, mix.i)
    CESE<- BEF(mono.i, mix.i,n)
    evenness_mix<- diversity(colMeans(mix.i),index = "shannon")/log(n)
    evenness_mono<- diversity(colMeans(mono.i),index = "shannon")/log(n)
    
    CorS_mix<- mean(cor(mix.i)[upper.tri(cor(mix.i))], na.rm = T)
    CorS_mono<- mean(cor(mono.i)[upper.tri(cor(mono.i))], na.rm = T)
    
    plot<- plotcode[i]
    N <- dBEF_data[dBEF_data$plotcode ==plotcode[i],]$sowndiv[1]
    dBEF_time<- rbind(dBEF_time, c(plot,j,n,N,miu_mix,miu_mono, evenness_mix,evenness_mono, CorS_mono, CESE, a)) 
    
  }  }
colnames(dBEF_time)<-c("plot","Order","n","Richness","miu_mix","miu_mono", "evenness_mix","evenness_mono","CorS_mono","NBE.F","SE","CE", "CV.mix", "CV.mono","CV.avg", "CVS.mix", "CVS.mono", "phi.mix", "phi.mono","AE.SpVar", "SE.SpVar", "AE.Syn", "SE.Syn")
dBEF_time[,-c(1,2)]<-apply(dBEF_time[,-c(1,2)],2,as.numeric)

dBEF_time$NIIE <-  dBEF_time$CV.mono/dBEF_time$CV.avg
dBEF_time$NIDE <- dBEF_time$CV.mix/dBEF_time$CV.mono
dBEF_time$NBE <- dBEF_time$CV.mix/dBEF_time$CV.avg
dBEF_time$NIDE.SpVar <- dBEF_time$CVS.mix/dBEF_time$CVS.mono
dBEF_time$NIDE.Syn <- dBEF_time$phi.mix/dBEF_time$phi.mono
dBEF_time$NIIE.Syn<- dBEF_time$phi.mono
dBEF_time$NIIE.SpVar<- dBEF_time$CVS.mono/dBEF_time$CV.avg
dBEF_time$Exp<-"D"
dBEF_time[dBEF_time$Order== "D1",]$Order<- "New"
dBEF_time[dBEF_time$Order== "D3",]$Order<- "Old"



# ----------------------------
###The effect of SPEI on NIIE.CV and NIDE,CV
# ----------------------------

jena_detrend<- data.frame()
for (j in 1: length(unique(Jena_time$Richness))) {
  Jena_time.j<- Jena_time[Jena_time$Richness == unique(Jena_time$Richness)[j], ]
  jena_detrend.j<- Jena_time.j%>% select(Richness, plot,  Order,SPEI_mean,Driest)
  for (i in 4:10) {
    Jena_time.i<- Jena_time.j%>%
      select(plot,Richness,Order,NBE,NIIE,NIIE.Syn, NIIE.SpVar,NIDE,NIDE.SpVar,NIDE.Syn)
    Jena_time.i<- Jena_time.i[,c(1,2,3,i)]
    colnames(Jena_time.i)[4]<- "variable"
    model<- gls(variable~log(Order), correlation=corAR1(form=~1|plot), Jena_time.i)
    res<- model$residuals
    jena_detrend.j<- cbind(jena_detrend.j, res)
    
  }
  jena_detrend<- rbind(jena_detrend, jena_detrend.j)
}
colnames(jena_detrend)[6:12]<- c("NBE","NIIE","NIIE.Syn","NIIE.SpVar","NIDE","NIDE.SpVar","NIDE.Syn")

