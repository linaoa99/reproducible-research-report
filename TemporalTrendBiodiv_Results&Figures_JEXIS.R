# Code for the paper:
#Meng, B., et al.
# Temporal decoupling of biodiversity effects on grassland productivity and stability
#
# This script performs the data analysis and generates the figures.
#
# Raw data downloaded from: https://jexis.idiv.de/

# Includes:
#   - Aboveground plant biomass from main plots in the Jena Experiment (2003–2021)
#      ID:38, 40,41,124,127,176,320
#   - Aboveground plant biomass from the dBEF Experiment (2017–2021)
#      ID:12,71,142,143,380,381,512
#   - Precipitation, temperature, and SPEI data for the Jena Experiment
#      ID:484,507

# For details on data filtering, see the Methods section of:
#   Meng, B., et al. (2026). Temporal decoupling of biodiversity effects on grassland productivity and stability.

#
# Reproducibility information
# This script was run successfully under the following R session:
#
#R version 4.6.0 (2026-04-24 ucrt)
#Platform: x86_64-w64-mingw32/x64
#Running under: Windows 11 x64 (build 26200)

##contact: Bo Meng
##bomeng@cau.edu.cn


# ----------------------------
# Load required packages
# ----------------------------

library("nlme")
#nlme_3.1-169
library("MuMIn")
#MuMIn_1.48.19
library("tidyverse")
#tidyverse_2.0.0
library("vegan")
#vegan_2.7-3
library("ggpubr")
#ggpubr_0.6.3
library("RColorBrewer")
#RColorBrewer_1.1-3

# ----------------------------
# Read the data that have been published together with this R script
jena_detrend<-read.csv("jena_detrend_SPEI.csv")
Jena_time<-read.csv("Jena_time_NBE.csv")
dBEF_time<-read.csv("dBEF_time_NBE.csv")


# ----------------------------
### temporal trends analysis
# ----------------------------
slope.df<- data.frame()
for (i in 1:4) {
  Jena_time.i<- Jena_time[Jena_time$Richness == unique(Jena_time$Richness)[i],-29]
  #Jena_time.i[,29:35]<- apply(Jena_time.i[,29:35], 2, function(x){log(x)})
  slope.c<-  c()
  slopeSe.c<-  c()
  intercept.c<-  c()
  interceptSe.c<-  c()
  for (j in 11:35) {
    Jena_time.ij<-   as.data.frame(cbind(Jena_time.i[,c(1,2,4)], Jena_time.i[,j])) 
    colnames(Jena_time.ij)<-c("block","plot","Order","variable") 
    model<- lme(variable ~ log(Order),random = ~1|plot,correlation=corAR1(form=~1|plot), Jena_time.ij )
    intercept<- model$coefficients$fixed[1]
    
    intercept.se<- sqrt(diag(vcov(model)))[1]
    slope<- model$coefficients$fixed[2]
    slope.se<- sqrt(diag(vcov(model)))[2]
    slope.c<- c(slope.c, slope)
    slopeSe.c<- c(slopeSe.c, slope.se)
    intercept.c<-  c(intercept.c, intercept)
    interceptSe.c<-  c(interceptSe.c, intercept)
  }
  Var  = colnames(Jena_time.i[,11:35])
  slope.df<- rbind(slope.df, data.frame( Richness = unique(Jena_time$Richness)[i], VarName  = colnames(Jena_time.i[,11:35]), intercept = intercept.c,slope =  slope.c, interceptSe = interceptSe.c, slopeSe =  slopeSe.c))
  
}
subtraction<- function(x){
  d = x[2]-x[1]
  return(d)
}
dBEF_time_NBECV<- dBEF_time%>%
  select(plot,Richness, Order, NBE.F, CE,SE,CV.mix,CVS.mix,phi.mix, NBE, NIIE,NIIE.SpVar,NIIE.Syn, NIDE, NIDE.SpVar,NIDE.Syn)%>%
  group_by(plot,Richness)%>%
  summarise(dNBE.F = subtraction(NBE.F),
            dCE = subtraction(CE),
            dSE = subtraction(SE),
            dCV.mix = subtraction( log(CV.mix)),
            dCVS.mix = subtraction(log(CVS.mix)),
            dphi.mix = subtraction(log(phi.mix)),
            dNBE = subtraction(log(NBE)),
            dNIIE = subtraction(log(NIIE)),
            dNIIE.SpVar= subtraction(log(NIIE.SpVar)),
            dNIIE.Syn= subtraction(log(NIIE.Syn)),
            dNIDE = subtraction(log(NIDE)),
            dNIDE.SpVar= subtraction(log(NIDE.SpVar)),
            dNIDE.Syn= subtraction(log(NIDE.Syn)))


color_palette.M<- c(brewer.pal(9, "YlOrRd")[c(3,5,7,9)])
color_palette.D<- c(brewer.pal(9, "GnBu")[c(3,5,7,9)])
color_palette.dom<- c(brewer.pal(9, "YlOrRd")[c(3:9)])
color_palette.main<- c(brewer.pal(9, "GnBu")[c(3:9)])

Jena_time.f<- na.omit( Jena_time)%>%mutate(Richness = as.factor(Richness))
# ----------------------------
###Fig.1
# ----------------------------
anova(lme(NBE~ log(Order)*Richness,random = ~1|plot,  correlation=corAR1(form=~1|plot), Jena_time))
Fig.1a<-ggplot(Jena_time.f)+
  theme_bw()+
  geom_smooth(aes(x =Order, y=NBE, color = Richness, fill = Richness),alpha=0.2,span=1, size=0.8, linetype = 1, method = 'lm',formula=y~log(x),show.legend = F,fullrange = T)+              
  scale_linetype_manual(values = c(1,1,1,1))+
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.5)+
  scale_color_manual(values =  color_palette.M)+
  scale_fill_manual(values =  color_palette.M)+
  xlab("time order")+
  ylab(expression(NBE[CV]))+
  theme(axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        panel.border = element_rect(fill=NA,color="black", size=0.4, linetype="solid"))

anova(lme(NIIE~ log(Order)*Richness,random = ~1|plot,  correlation=corAR1(form=~1|plot), Jena_time))
Fig.1b<- ggplot(Jena_time.f)+
  theme_bw()+
  geom_smooth(aes(x =Order, y=NIIE, color = Richness, fill = Richness),alpha=0.2,span=1, size=0.8, linetype = 1, method = 'lm',formula=y~log(x),show.legend = F,fullrange = T)+              
  scale_linetype_manual(values = c(1,1,1,1))+
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.5)+
  scale_color_manual(values =  color_palette.M)+
  scale_fill_manual(values =  color_palette.M)+
  xlab("time order")+
  ylab(expression(NIIE[CV]))+
  theme(axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        panel.border = element_rect(fill=NA,color="black", size=0.4, linetype="solid"))

anova(lme(NIDE~ log(Order)*Richness,random = ~1|plot,  correlation=corAR1(form=~1|plot), Jena_time))
Fig.1c<- ggplot(Jena_time.f)+
  theme_bw()+
  geom_smooth(aes(x =Order, y=NIDE, color = Richness, fill = Richness),alpha=0.2,span=1, size=0.8, linetype = 1, method = 'lm',formula=y~log(x),show.legend = F,fullrange = T)+              
  scale_linetype_manual(values = c(1,1,1,1))+
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.5)+
  scale_color_manual(values =  color_palette.M)+
  scale_fill_manual(values =  color_palette.M)+
  xlab("time order")+
  ylab(expression(NIDE[CV]))+
  theme(axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        panel.border = element_rect(fill=NA,color="black", size=0.4, linetype="solid"))

model_NBE<- lme(NBE~ Order:as.factor(Richness)+as.factor(Richness)-1,random = ~1|plot,  correlation=corAR1(form=~1|plot), Jena_time)
summary(model_NBE)
Fig.1d<- ggplot( slope.df[slope.df$VarName == "NBE",])+
  geom_bar(aes(x = as.factor(Richness), y = slope, fill = as.factor(Richness)),color = "black", position = "dodge",stat = "identity", size = 0.4,alpha= 0.8,width=0.8)+
  geom_errorbar(aes(x = as.factor(Richness), ymin = slope - slopeSe, ymax = slope + slopeSe), size = 0.4, width = 0.3)+
  scale_fill_manual(values =  color_palette.M)+
  theme_bw()+
  ylim(-0.4,0.4)+
  geom_hline(yintercept = 0,linetype = 2)+
  xlab("Sown diversity")+
  ylab(expression(Temporal~trend~on~NBE[CV]))+
  theme(axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        plot.margin = margin(t = 5, r = 5, b = 20, l = 5),
        panel.border = element_rect(fill=NA,color="black", size=0.4, linetype="solid"))

model_NIIE<- lme(log(NIIE)~ Order:as.factor(Richness)+as.factor(Richness)-1,random = ~1|plot,  correlation=corAR1(form=~1|plot), Jena_time)
summary(model_NIIE)
Fig.1e<- ggplot( slope.df[slope.df$VarName == "NIIE",])+
  geom_bar(aes(x = as.factor(Richness), y = slope, fill = as.factor(Richness)),color = "black", position = "dodge",stat = "identity", size = 0.4, alpha= 0.8,width=0.8)+
  geom_errorbar(aes(x = as.factor(Richness), ymin = slope - slopeSe, ymax = slope + slopeSe), size = 0.4,width = 0.3)+
  scale_fill_manual(values =color_palette.M)+
  ylim(-0.4,0.4)+
  theme_bw()+
  geom_hline(yintercept = 0,linetype = 2)+
  xlab("Sown diversity")+
  ylab(expression(Temporal~trend~on~NIIE[CV]))+
  theme(axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        plot.margin = margin(t = 5, r = 5, b = 20, l = 5),
        panel.border = element_rect(fill=NA,color="black", size=0.4, linetype="solid"))

model_NIDE<- lme(log(NIDE)~ Order:as.factor(Richness)+as.factor(Richness)-1,random = ~1|plot,  correlation=corAR1(form=~1|plot), Jena_time)
summary(model_NIDE)
Fig.1f<- ggplot( slope.df[slope.df$VarName == "NIDE",])+
  geom_bar(aes(x = as.factor(Richness), y = slope, fill = as.factor(Richness)),color = "black", position = "dodge",stat = "identity", size = 0.4, alpha= 0.8,width=0.8)+
  geom_errorbar(aes(x = as.factor(Richness), ymin = slope - slopeSe, ymax = slope + slopeSe), size = 0.4, width = 0.3)+
  scale_fill_manual(values = color_palette.M)+
  ylim(-0.4,0.4)+
  theme_bw()+
  geom_hline(yintercept = 0,linetype = 2, linewidth = 0.5)+
  xlab("Sown diversity")+
  ylab(expression(Temporal~trend~on~NIDE[CV]))+
  theme(axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        plot.margin = margin(t = 5, r = 5, b = 20, l = 5),
        panel.border = element_rect(fill=NA,color="black", size=0.4, linetype="solid"))

summary(lm(dNBE~1+as.factor(Richness)-1, dBEF_time_NBECV))
Fig.1g<- ggplot(dBEF_time_NBECV)+geom_boxplot(aes(x= as.factor(Richness), y= dNBE, fill = as.factor(Richness)),shape = 21, size = 0.4,  alpha = 0.8,width=0.8)+
  theme_bw()+
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.5)+
  scale_fill_manual(values = color_palette.D)+
  ylim(-3.3, 1.8)+
  xlab("Sown diversity")+
  ylab(expression(Delta~NBE[CV]))+
  theme(axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        panel.border = element_rect(fill=NA,color="black", size=0.4, linetype="solid"))

summary(lm(dNIIE~1+as.factor(Richness)-1, dBEF_time_NBECV))
Fig.1h<- ggplot(dBEF_time_NBECV)+geom_boxplot(aes(x= as.factor(Richness), y= dNIIE, fill = as.factor(Richness)),shape = 21, size = 0.4,  alpha = 0.8,width=0.8)+
  theme_bw()+
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.5)+
  scale_fill_manual(values = color_palette.D)+
  ylim(-3.3, 1.8)+
  xlab("Sown diversity")+
  ylab(expression(Delta~NIIE[CV]))+
  theme(axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        panel.border = element_rect(fill=NA,color="black", size=0.4, linetype="solid"))

summary(lm(dNIDE~1+as.factor(Richness)-1, dBEF_time_NBECV))
Fig.1i<- ggplot(dBEF_time_NBECV)+geom_boxplot(aes(x= as.factor(Richness), y= dNIDE, fill = as.factor(Richness)),shape = 21, size = 0.4,  alpha = 0.8,width=0.8)+
  theme_bw()+
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.5)+
  scale_fill_manual(values = color_palette.D)+
  ylim(-3.3, 1.8)+
  xlab("Sown diversity")+
  ylab(expression(Delta~NIDE[CV]))+
  theme(axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        panel.border = element_rect(fill=NA,color="black", size=0.4, linetype="solid"))
Fig1<- ggarrange(Fig.1a, Fig.1b, Fig.1c, Fig.1d, Fig.1e, Fig.1f,Fig.1g,Fig.1h,Fig.1i,
                 widths = c(1, 1, 1),heights = c(1, 0.8, 0.8),
                 align = "hv", common.legend = T)

# ----------------------------
####Fig2
# ----------------------------
anova(lme(NIIE.SpVar~ log(Order)*Richness,random = ~1|plot, correlation=corAR1(form=~1|plot), Jena_time))
Fig2.a<- ggplot(Jena_time.f)+
  theme_bw()+
  geom_smooth(aes(x =Order, y=NIIE.SpVar, color = Richness, fill = Richness),alpha=0.2,span=1, size=0.7, linetype = 1, method = 'lm',formula=y~log(x),show.legend = F,fullrange = T)+              
  scale_linetype_manual(values = c(1,1,1,1))+
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.5)+
  scale_color_manual(values =  color_palette.M)+
  scale_fill_manual(values =  color_palette.M)+
  xlab("time order")+
  ylab(expression(NIIE[SpVar]))+
  theme(axis.title.x=element_text(vjust=0, size=14),
        axis.title.y=element_text(hjust=0.5, size=14),
        axis.text.x=element_text(vjust=0,size=12),
        axis.text.y=element_text(hjust=0,size=12),
        panel.border = element_rect(fill=NA,color="black", size=0.4, linetype="solid"))

anova(lme(NIIE.Syn~ log(Order)*Richness,random = ~1|plot,  correlation=corAR1(form=~1|plot), Jena_time))
Fig2.b<- ggplot(Jena_time.f)+
  theme_bw()+
  geom_smooth(aes(x =Order, y=NIIE.Syn, color = Richness, fill = Richness),alpha=0.2,span=1, size=0.7, linetype = 1, method = 'lm',formula=y~log(x),show.legend = F,fullrange = T)+              
  scale_linetype_manual(values = c(1,1,1,1))+
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.5)+
  scale_color_manual(values =  color_palette.M)+
  scale_fill_manual(values =  color_palette.M)+
  xlab("time order")+
  ylab(expression(NIIE[Syn]))+
  theme(axis.title.x=element_text(vjust=0, size=14),
        axis.title.y=element_text(hjust=0.5, size=14),
        axis.text.x=element_text(vjust=0,size=12),
        axis.text.y=element_text(hjust=0,size=12),
        panel.border = element_rect(fill=NA,color="black", size=0.4, linetype="solid"))

anova(lme(evenness_mono~ log(Order)*Richness,random = ~1|plot,  correlation=corAR1(form=~1|plot), Jena_time))
Fig2.c<- ggplot(Jena_time.f)+
  theme_bw()+
  geom_smooth(aes(x =Order, y=evenness_mono, color = Richness, fill = Richness),alpha=0.2,span=1, size=0.7, linetype = 1, method = 'lm',formula=y~log(x),show.legend = F,fullrange = T)+              
  scale_linetype_manual(values = c(1,1,1,1))+
  #geom_hline(yintercept = 1, linetype = 2, linewidth = 0.5)+
  scale_color_manual(values =  color_palette.M)+
  scale_fill_manual(values =  color_palette.M)+
  xlab("time order")+
  ylab(expression(NIIE[CV]))+
  theme(axis.title.x=element_text(vjust=0, size=14),
        axis.title.y=element_text(hjust=0.5, size=14),
        axis.text.x=element_text(vjust=0,size=12),
        axis.text.y=element_text(hjust=0,size=12),
        panel.border = element_rect(fill=NA,color="black", size=0.4, linetype="solid"))

anova(lme(CorS_mono~ log(Order)*Richness,random = ~1|plot,  correlation=corAR1(form=~1|plot), Jena_time))
Fig2.d<- ggplot(Jena_time.f)+
  theme_bw()+
  geom_smooth(aes(x =Order, y= CorS_mono, color = Richness, fill = Richness),alpha=0.2,span=1, size=0.7, linetype = 1, method = 'lm',formula=y~log(x),show.legend = F,fullrange = T)+              
  scale_linetype_manual(values = c(1,1,1,1))+
  
  scale_color_manual(values =  color_palette.M)+
  scale_fill_manual(values =  color_palette.M)+
  xlab("time order")+
  ylab(expression(NIIE[CV]))+
  theme(axis.title.x=element_text(vjust=0, size=14),
        axis.title.y=element_text(hjust=0.5, size=14),
        axis.text.x=element_text(vjust=0,size=12),
        axis.text.y=element_text(hjust=0,size=12),
        panel.border = element_rect(fill=NA,color="black", size=0.4, linetype="solid"))


Fig2<- ggarrange(Fig2.a, Fig2.b, Fig2.c, Fig2.d,
                 ncol=2,nrow = 2,
                 align = "hv", common.legend = T)


# ----------------------------
###Fig.3
# ----------------------------
Fig3a<- ggplot(Jena_time.f)+
  theme_bw()+
  geom_smooth(aes(x =Order, y=NIDE.SpVar, color = Richness, fill = Richness),alpha=0.2,span=1, size=0.7, linetype = 1, method = 'lm',formula=y~log(x),show.legend = F,fullrange = T)+              
  scale_linetype_manual(values = c(1,1,1,1))+
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.3)+
  scale_color_manual(values =  color_palette.M)+
  scale_fill_manual(values =  color_palette.M)+
  xlab("time order")+
  ylab(expression(NIDE[SpVar]))+
  theme(axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        panel.border = element_rect(fill=NA,color="black", size=0.4, linetype="solid"))
Fig3d<- ggplot(Jena_time.f)+
  theme_bw()+
  geom_smooth(aes(x =Order, y=NIDE.Syn, color = Richness, fill = Richness),alpha=0.2,span=1, size=0.7, linetype = 1, method = 'lm',formula=y~log(x),show.legend = F,fullrange = T)+              
  scale_linetype_manual(values = c(1,1,1,1))+
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.3)+
  scale_color_manual(values =  color_palette.M)+
  scale_fill_manual(values =  color_palette.M)+
  xlab("Time order")+
  ylab(expression(NIDE[Syn]))+
  theme(axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        panel.border = element_rect(fill=NA,color="black", size=0.4, linetype="solid"))

anova(lme(AE.SpVar~ log(Order)*Richness,random = ~1|plot,  correlation=corAR1(form=~1|plot), Jena_time))
model_AE.SpVar<- lme(AE.SpVar~ log(Order):as.factor(Richness)+as.factor(Richness)-1,random = ~1|plot,  correlation=corAR1(form=~1|plot), Jena_time)
summary(model_AE.SpVar)
Fig3b<- ggplot( slope.df[slope.df$VarName == "AE.SpVar",])+
  geom_hline(yintercept = 0,linetype = 2, linewidth = 0.3)+
  geom_bar(aes(x = as.factor(Richness), y = slope, fill = as.factor(Richness)),color = "black", position = "dodge",stat = "identity", size = 0.4,alpha= 0.8,width=0.8)+
  geom_errorbar(aes(x = as.factor(Richness), ymin = slope - slopeSe, ymax = slope + slopeSe), size = 0.4, width = 0.3)+
  scale_fill_manual(values =  color_palette.M)+
  theme_bw()+
  ylim(-0.3,0.4)+
  xlab("Sown diversity")+
  ylab(expression(Temporal~trend~on~AE[SpVar]))+
  theme(axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        panel.border = element_rect(fill=NA,color="black", size=0.4, linetype="solid"))
anova(lme(SE.SpVar~ log(Order)*Richness,random = ~1|plot,  correlation=corAR1(form=~1|plot), Jena_time))
model_SE.SpVar<- lme(SE.SpVar~ log(Order):as.factor(Richness)+as.factor(Richness)-1,random = ~1|plot,  correlation=corAR1(form=~1|plot), Jena_time)
summary(model_SE.SpVar)
Fig3c<- ggplot( slope.df[slope.df$VarName == "SE.SpVar",])+
  geom_hline(yintercept = 0,linetype = 2, linewidth = 0.3)+
  geom_bar(aes(x = as.factor(Richness), y = slope, fill = as.factor(Richness)),color = "black", position = "dodge",stat = "identity", size = 0.4, alpha= 0.8,width=0.8)+
  geom_errorbar(aes(x = as.factor(Richness), ymin = slope - slopeSe, ymax = slope + slopeSe), size = 0.4,width = 0.3)+
  scale_fill_manual(values =color_palette.M)+
  ylim(-0.3,0.4)+
  theme_bw()+
  xlab("Sown diversity")+
  ylab(expression(Temporal~trend~on~SE[SpVar]))+
  theme(axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        panel.border = element_rect(fill=NA,color="black", size=0.4, linetype="solid"))

model_AE.Syn<- lme(AE.Syn~ log(Order):as.factor(Richness)+as.factor(Richness)-1,random = ~1|plot,  correlation=corAR1(form=~1|plot), Jena_time)
anova(lme(AE.Syn~ log(Order)*Richness,random = ~1|plot,  correlation=corAR1(form=~1|plot), Jena_time))
summary(model_AE.Syn)
Fig3e<- ggplot( slope.df[slope.df$VarName == "AE.Syn",])+
  geom_hline(yintercept = 0,linetype = 2, linewidth = 0.3)+
  geom_bar(aes(x = as.factor(Richness), y = slope, fill = as.factor(Richness)),color = "black", position = "dodge",stat = "identity", size = 0.4, alpha= 0.8,width=0.8)+
  geom_errorbar(aes(x = as.factor(Richness), ymin = slope - slopeSe, ymax = slope + slopeSe), size = 0.4, width = 0.3)+
  scale_fill_manual(values = color_palette.M)+
  ylim(-0.3,1)+
  theme_bw()+
  xlab("Sown diversity")+
  ylab(expression(Temporal~trend~on~AE[Syn]))+
  theme(axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        panel.border = element_rect(fill=NA,color="black", size=0.4, linetype="solid"))

model_SE.Syn<- lme(SE.Syn~ log(Order):as.factor(Richness)+as.factor(Richness)-1,random = ~1|plot,  correlation=corAR1(form=~1|plot), Jena_time)
anova(lme(SE.Syn~ log(Order)*Richness,random = ~1|plot,  correlation=corAR1(form=~1|plot), Jena_time))
summary(model_SE.Syn)
Fig3f<- ggplot( slope.df[slope.df$VarName == "SE.Syn",])+
  geom_hline(yintercept = 0,linetype = 2, linewidth = 0.3)+
  geom_bar(aes(x = as.factor(Richness), y = slope, fill = as.factor(Richness)),color = "black", position = "dodge",stat = "identity", size = 0.4, alpha= 0.8,width=0.8)+
  geom_errorbar(aes(x = as.factor(Richness), ymin = slope - slopeSe, ymax = slope + slopeSe), size = 0.4, width = 0.3)+
  scale_fill_manual(values = color_palette.M)+
  ylim(-0.3,1)+
  theme_bw()+
  xlab("Sown diversity")+
  ylab(expression(Temporal~trend~on~SE[Syn]))+
  theme(axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        panel.border = element_rect(fill=NA,color="black", size=0.4, linetype="solid"))
Fig3<- ggarrange(Fig3a, Fig3b, Fig3c, Fig3d,Fig3e, Fig3f,
                 ncol=3,nrow = 2,
                 align = "hv", common.legend = T)
# ----------------------------
###Fig.4
# ----------------------------
Jena_time.f2<-mutate(Jena_time, Order = as.factor(Order))

anova(lme(NBE~ Order*NBE.F, random = ~1|Richness/plot, correlation=corAR1(form=~1|Richness/plot), Jena_time[Jena_time$Order %in% c(1,1:3*5),]))
MNBE.NBE<- lme(NBE~ NBE.F:Order+Order-1, random = ~1|Richness/plot, correlation=corAR1(form=~1|Richness/plot), Jena_time.f2[Jena_time.f2$Order %in% c(1,1:3*5),])
pred_MNBE.NBE<- predict(MNBE.NBE,  Jena_time.f2[Jena_time.f2$Order %in% c(1,1:3*5),], level = 0, se.fit=T)
summary(MNBE.NBE)

MNBE.NBE_D<- lme(NBE~ NBE.F, random = ~1|Richness, correlation=corAR1(form=~1|Richness/plot), dBEF_time[dBEF_time$Order == "New",])
pred_MNBE.NBE_D<- predict(MNBE.NBE_D, dBEF_time[dBEF_time$Order == "New",], level = 0, se.fit=T)
Fig.4a<- ggplot(Jena_time[Jena_time$Order %in% c(1,1:3*5),])+geom_point( aes(x= NBE.F, y= NBE), shape = 21, color = "grey", alpha = 0.6, size = 1.3)+
  geom_smooth(data = Jena_time[Jena_time$Order %in% c(1,1:3*5),], aes(x =NBE.F, y= pred_MNBE.NBE$fit, group = rev(Order), color = rev(Order) ), size=1.2,method = 'lm',formula=y~x, fill = NA, alpha = 0.1)+
  geom_smooth(data = dBEF_time[dBEF_time$Order == "New",], aes(x =NBE.F, y= pred_MNBE.NBE_D$fit),color = color_palette.D[3], size=1.2,method = 'lm',formula=y~x, fill = NA, alpha = 0.1)+
  scale_color_gradientn(values = c (0,0.4,0.6,0.8,1.0),colors =alpha(color_palette.M, 0.9))+
  #xlim(0,9.3)+
  #ylim(-0,2)+
  theme_bw()+
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.5)+
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.5)+
  xlab(expression(NBE[Bio]))+
  ylab(expression(NBE[CV]))+
  theme(axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        panel.border = element_rect(fill=NA,color="black", size=0.4, linetype="solid"))

anova(lme(NIIE~ Order*NBE.F, random = ~1|Richness/plot, correlation=corAR1(form=~1|Richness/plot), Jena_time[Jena_time$Order %in% c(1,1:3*5),]))
MNIIE.NBE<- lme(NIIE~ NBE.F:Order+Order-1, random = ~1|Richness/plot, correlation=corAR1(form=~1|Richness/plot), Jena_time.f2[Jena_time.f2$Order %in% c(1,1:3*5),])
summary(MNIIE.NBE)
pred_MNIIE.NBE<- predict(MNIIE.NBE,  Jena_time.f2[Jena_time.f2$Order %in% c(1,1:3*5),], level = 0, se.fit=T)
MNIIE.NBE_D<- lme(NIIE~ NBE.F, random = ~1|Richness, correlation=corAR1(form=~1|Richness/plot), dBEF_time[dBEF_time$Order == "New",])
pred_MNIIE.NBE_D<- predict(MNIIE.NBE_D, dBEF_time[dBEF_time$Order == "New",], level = 0, se.fit=T)
Fig.4b<- ggplot(Jena_time[Jena_time$Order %in% c(1,1:3*5),])+geom_point( aes(x= NBE.F, y= NIIE), shape = 21, color = "grey", alpha = 0.6, size = 1.3)+
  geom_smooth(data = Jena_time[Jena_time$Order %in% c(1,1:3*5),], aes(x =NBE.F, y= pred_MNIIE.NBE$fit, group = rev(Order), color = rev(Order) ), size=1.2,method = 'lm',formula=y~x, fill = NA, alpha = 0.1)+
  geom_smooth(data = dBEF_time[dBEF_time$Order == "New",], aes(x =NBE.F, y= pred_MNIIE.NBE_D$fit),color = color_palette.D[3], size=1.2,method = 'lm',formula=y~x, fill = NA, alpha = 0.1)+
  scale_color_gradientn(values = c (0,0.4,0.6,0.8,1.0),colors =alpha(color_palette.M, 0.9))+
  #lim(0,9.3)+
  ylim(0,2.5)+
  theme_bw()+
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.5)+
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.5)+
  #ylim(-3, 1.8)+
  xlab(expression(NBE[Bio]))+
  ylab(expression(NIDE[CV]))+
  theme(axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        panel.border = element_rect(fill=NA,color="black", size=0.4, linetype="solid"))

anova(lme(NIDE~ Order*NBE.F, random = ~1|Richness/plot, correlation=corAR1(form=~1|Richness/plot), Jena_time[Jena_time$Order %in% c(1,1:3*5),]))
MNIDE.NBE<- lme(NIDE~ NBE.F:Order+Order-1, random = ~1|Richness/plot, correlation=corAR1(form=~1|Richness/plot), Jena_time.f2[Jena_time.f2$Order %in% c(1,1:3*5),])
pred_MNIDE.NBE<- predict(MNIDE.NBE,  Jena_time.f2[Jena_time.f2$Order %in% c(1,1:3*5),], level = 0, se.fit=T)
MNIDE.NBE_D<- lme(NIDE~ NBE.F, random = ~1|Richness, correlation=corAR1(form=~1|Richness/plot), dBEF_time[dBEF_time$Order == "New",])
pred_MNIDE.NBE_D<- predict(MNIDE.NBE_D, dBEF_time[dBEF_time$Order == "New",], level = 0, se.fit=T)
Fig.4c<- ggplot(Jena_time[Jena_time$Order %in% c(1,1:3*5),])+geom_point( aes(x= NBE.F, y= NIDE), shape = 21, color = "grey", alpha = 0.6, size = 1.3)+
  geom_smooth(data = Jena_time[Jena_time$Order %in% c(1,1:3*5),], aes(x =NBE.F, y= pred_MNIDE.NBE$fit, group = rev(Order), color = rev(Order) ), size=1.2,method = 'lm',formula=y~x, fill = NA, alpha = 0.1)+
  geom_smooth(data = dBEF_time[dBEF_time$Order == "New",], aes(x =NBE.F, y= pred_MNIDE.NBE_D$fit),color = color_palette.D[3], size=1.2,method = 'lm',formula=y~x, fill = NA, alpha = 0.1)+
  scale_color_gradientn(values = c (0,0.4,0.6,0.8,1.0),colors =alpha(color_palette.M, 0.9))+
  #xlim(-0,9.3)+
  #ylim(-0,4)+
  theme_bw()+
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.3)+
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3)+
  xlab(expression(NBE[Bio]))+
  ylab(expression(NIDE[CV]))+
  theme(axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        panel.border = element_rect(fill=NA,color="black", size=0.4, linetype="solid"))
Fig4<- ggarrange(Fig.4a, Fig.4b, Fig.4c, ncol=3,align = "hv", common.legend = T)





anova(lme(NIIE~I(SPEI_mean)+I((SPEI_mean)^2), random = ~1|Richness,  jena_detrend))

Pre.SPEI.NIIE<- predict((lme(NIIE~I(SPEI_mean)+I(SPEI_mean^2), random = ~1|Richness,  jena_detrend)),  jena_detrend , se.fit = T,level =0)
Pre.SPEI.NIIE16<- predict( lm(NIIE~I(SPEI_mean)+I(SPEI_mean^2),  jena_detrend[jena_detrend$Richness == 16,]), jena_detrend[jena_detrend$Richness == 16,] , se.fit = T, level =0)
Pre.SPEI.NIIE8<- predict( lm(NIIE~I(SPEI_mean)+I(SPEI_mean^2),  jena_detrend[jena_detrend$Richness == 8,]), jena_detrend[jena_detrend$Richness == 8,] , se.fit = T, level =0)
Pre.SPEI.NIIE4<- predict( lm(NIIE~I(SPEI_mean)+I(SPEI_mean^2),  jena_detrend[jena_detrend$Richness == 4,]), jena_detrend[jena_detrend$Richness == 4,] , se.fit = T, level =0)
Pre.SPEI.NIIE2<- predict( lm(NIIE~I(SPEI_mean)+I(SPEI_mean^2),  jena_detrend[jena_detrend$Richness == 2,]), jena_detrend[jena_detrend$Richness == 2,] , se.fit = T, level =0)

# ----------------------------
###Fig.5
# ----------------------------
jena_detrend$group<- paste(jena_detrend$Order, jena_detrend$Richness, sep = "_")
jena_detrend$group<- factor(jena_detrend$group, levels = c(unique(jena_detrend$group)))


model<- (lme(NIDE~I(SPEI_mean), random = ~1|Richness,  jena_detrend))
model.2<- update(model, .~. +I((SPEI_mean)^2))
AICc(model, model.2)

anova(model.2)
model.P<- (lm(NIDE~SPEI_mean:as.factor(Richness)-1+I(SPEI_mean^2),  jena_detrend))
summary(model.2)

se<- function(x){sd(x)/sqrt(length(x)-1)}

jena_detrend_Richness<- jena_detrend%>%
  group_by(Richness,Order)%>%
  summarise(SPEI_mean = mean(SPEI_mean),
            Driest = mean(Driest),
            NIDE.mean = mean(NIDE),
            NIDE.SE = se(NIDE),
            NIIE.mean = mean(NIIE),
            NIIE.SE = se(NIIE))

Fig.5a <- ggplot(jena_detrend_Richness)+  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.5)+ geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5)+ 
  geom_point(aes(x = SPEI_mean, y = NIIE.mean,  color = as.factor(Richness)), size = 1.5, alpha = 0.6)+
  geom_errorbar(aes(x = SPEI_mean, y = NIIE.mean, ymin =NIIE.mean- NIIE.SE, ymax = NIIE.mean+ NIIE.SE,  color = as.factor(Richness)), size = 0.3)+
  #geom_smooth(data = jena_detrend, aes(x = SPEI_mean, y = Pre.SPEI.NIIE2$fit, group = as.factor(Richness),color = as.factor(Richness)))+
  geom_smooth(data = jena_detrend[jena_detrend$Richness == 16,] , aes(x = SPEI_mean, y = Pre.SPEI.NIIE16$fit),color = color_palette.M[4], size = 0.5)+
  geom_smooth(data = jena_detrend[jena_detrend$Richness == 8,] , aes(x = SPEI_mean, y = Pre.SPEI.NIIE8$fit),color = color_palette.M[3], size = 0.5)+
  geom_smooth(data = jena_detrend[jena_detrend$Richness == 4,] , aes(x = SPEI_mean, y = Pre.SPEI.NIIE4$fit),color = color_palette.M[2], size = 0.5)+
  geom_smooth(data = jena_detrend[jena_detrend$Richness == 2,] , aes(x = SPEI_mean, y = Pre.SPEI.NIIE2$fit),color = color_palette.M[1], size = 0.5)+
  geom_smooth(data = jena_detrend, aes(x = SPEI_mean, y = Pre.SPEI.NIIE$fit),color = "black", size = 0.8)+
  geom_ribbon(data = jena_detrend, aes(x =SPEI_mean, ymin= Pre.SPEI.NIIE$fit-1.96*Pre.SPEI.NIIE$se.fit, ymax=Pre.SPEI.NIIE$fit+1.96*Pre.SPEI.NIIE$se.fit), fill= "grey",alpha = 0.3,linetype=0)+
  scale_color_manual(values = color_palette.M)+
  #geom_point(aes(x = 15, y = zD1), size = 5, shape =21, fill= color_palette.D[3], alpha= 0.3)+
  theme_bw()+
  xlim(-0.82, 0.58)+
  ylim(-0.7, 0.7)+
  xlab('The average SPEI')+
  ylab(expression(The~detrended~NIIE[CV]))+
  theme(axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        panel.border = element_rect(fill=NA,color="black", size=0.4, linetype="solid"))


#jena_detrend.f<- mutate(jena_detrend, Richness = as.factor(Richness)) 
Pre.SPEI<- predict((lme(NIDE~I(SPEI_mean)+I(SPEI_mean^2), random = ~1|Richness,  jena_detrend)),  jena_detrend , se.fit = T,level =0)
Pre.SPEI16<- predict( lm(NIDE~I(SPEI_mean)+I(SPEI_mean^2),  jena_detrend[jena_detrend$Richness == 16,]), jena_detrend[jena_detrend$Richness == 16,] , se.fit = T, level =0)
Pre.SPEI8<- predict( lm(NIDE~I(SPEI_mean)+I(SPEI_mean^2),  jena_detrend[jena_detrend$Richness == 8,]), jena_detrend[jena_detrend$Richness == 8,] , se.fit = T, level =0)
Pre.SPEI4<- predict( lm(NIDE~I(SPEI_mean)+I(SPEI_mean^2),  jena_detrend[jena_detrend$Richness == 4,]), jena_detrend[jena_detrend$Richness == 4,] , se.fit = T, level =0)
Pre.SPEI2<- predict( lm(NIDE~I(SPEI_mean)+I(SPEI_mean^2),  jena_detrend[jena_detrend$Richness == 2,]), jena_detrend[jena_detrend$Richness == 2,] , se.fit = T, level =0)

Fig.5b <- ggplot(jena_detrend_Richness)+  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.5)+ geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5)+
  geom_point(aes(x = SPEI_mean, y = NIDE.mean,  color = as.factor(Richness)), size = 1.5, alpha = 0.6)+
  geom_errorbar(aes(x = SPEI_mean, y = NIDE.mean, ymin =NIDE.mean- NIDE.SE, ymax = NIDE.mean+ NIDE.SE,  color = as.factor(Richness)), size = 0.3)+
  #geom_smooth(data = jena_detrend, aes(x = SPEI_mean, y = Pre.SPEI2$fit, group = as.factor(Richness),color = as.factor(Richness)))+
  geom_smooth(data = jena_detrend[jena_detrend$Richness == 16,] , aes(x = SPEI_mean, y = Pre.SPEI16$fit),color = color_palette.M[4], size = 0.5)+
  geom_smooth(data = jena_detrend[jena_detrend$Richness == 8,] , aes(x = SPEI_mean, y = Pre.SPEI8$fit),color = color_palette.M[3], size = 0.5)+
  geom_smooth(data = jena_detrend[jena_detrend$Richness == 4,] , aes(x = SPEI_mean, y = Pre.SPEI4$fit),color = color_palette.M[2], size = 0.5)+
  geom_smooth(data = jena_detrend[jena_detrend$Richness == 2,] , aes(x = SPEI_mean, y = Pre.SPEI2$fit),color = color_palette.M[1], size = 0.5)+
  geom_smooth(data = jena_detrend, aes(x = SPEI_mean, y = Pre.SPEI$fit),color = "black", size = 0.8)+
  geom_ribbon(data = jena_detrend, aes(x =SPEI_mean, ymin= Pre.SPEI$fit-1.96*Pre.SPEI$se.fit, ymax=Pre.SPEI$fit+1.96*Pre.SPEI$se.fit), fill= "grey",alpha = 0.3,linetype=0)+
  scale_color_manual(values = color_palette.M)+
  #geom_point(aes(x = 15, y = zD1), size = 5, shape =21, fill= color_palette.D[3], alpha= 0.3)+
  theme_bw()+
  xlim(-0.82, 0.58)+
  ylim(-0.7, 0.7)+
  xlab('The average SPEI')+
  ylab(expression(The~detrended~NIDE[CV]))+
  theme(axis.title.x=element_text(vjust=0, size=12),
        axis.title.y=element_text(hjust=0.5, size=12),
        axis.text.x=element_text(vjust=0,size=10),
        axis.text.y=element_text(hjust=0,size=10),
        panel.border = element_rect(fill=NA,color="black", size=0.4, linetype="solid"))


Fig5<- ggarrange(Fig.5a, Fig.5b, 
                 widths = c(1, 1),nrow = 1,
                 align = "hv", common.legend = T)

