library(rgrass)
library(terra)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggcorrplot)

#### Dataset ####
bver50<-read_RAST(c("DEM","Slope_50","Walk_NF","DEM_easterness","DEM_northerness","Slope_50_Curvature","Topidx_50","TPI_9"))
siti<-read_VECT('Siti_pixel')
random<-read_VECT("BVER_random")
sitestab<-data.frame(extract(bver50,siti),name=rep(1,length(siti)))
nositestab<-data.frame(extract(bver50,random),name=rep(0,length(random)))
tab<-rbind(sitestab,nositestab)
tab$name<-as.factor(tab$name)
sitestab<-subset(tab,tab$name==1)
nositestab<-subset(tab,tab$name==0)

#### Collinearity ####

collin<-round(cor(select(tab, -ID, -name),method="pearson"),4)
collin

#### Empirical Cumulative Distribution Plot ####

jpeg("Variables_ecd.jpg",width=4000,height=3000, units="px",res=300)
ggarrange(
  ggplot(tab, aes(x=DEM,group = name,colour= name))+
    labs(color = "Site")+
    stat_ecdf(linewidth=1)+
    ylab(NULL),
  ggplot(tab, aes(x=Slope_50,group = name,colour= name))+
    labs(color = "Site")+
    stat_ecdf(linewidth=1)+
    ylab(NULL),
  ggplot(tab, aes(x=Walk_NF,group = name,colour= name))+
    labs(color = "Site")+
    stat_ecdf(linewidth=1)+
    ylab(NULL),
  ggplot(tab, aes(x=DEM_easterness,group = name,colour= name))+
    labs(color = "Site")+
    stat_ecdf(linewidth=1)+
    ylab(NULL),
  ggplot(tab, aes(x=DEM_northerness,group = name,colour= name))+
    labs(color = "Site")+
    stat_ecdf(linewidth=1)+
    ylab(NULL),
  ggplot(tab, aes(x=Slope_50_Curvature,group = name,colour= name))+
    labs(color = "Site")+
    stat_ecdf(linewidth=1)+
    ylab(NULL),  
  ggplot(tab, aes(x=Topidx_50,group = name,colour= name))+
    labs(color = "Site")+
    stat_ecdf(linewidth=1)+
    ylab(NULL),
  ggplot(tab, aes(x=TPI_9,group = name,colour= name))+
    labs(color = "Site")+
    stat_ecdf(linewidth=1)+
    ylab(NULL),
  ncol=4,nrow=2,align='hv',common.legend = TRUE, legend="bottom")
dev.off()

#### Kolmogorov-Smirnov ####

sites_num<-select(sitestab, -ID, -name)
nosites_num<-select(nositestab, -ID, -name)

library(Matching)

#Function to analyse all the record at once
#x = sites
#y = no-sites
#nsim = number of repetitions

ksb.col<-function(x,y,nsim){
  ksb_tab<-list()
  for(i in 1:ncol(x)){
    ksb_tab[[i]]<-ks.boot(x[,i],y[,i],nboots=nsim)
  }
  names(ksb_tab)<-colnames(x)
  return(lapply(ksb_tab,summary))
}

ks_file<-capture.output(ksb.col(sites_num,nosites_num,10000))
ks_file

#### Univariate Logistic Regression ####

univar<-function(x,y){
  lrunivar<-list()
  for(i in x){
    f<-formula(paste("name","~",i))
    lrunivar[[i]]<-glm(f,data=y,family=binomial(logit))
  }
  return(lrunivar)
}

indep<-c("DEM", "Slope_50", "Walk_NF", "DEM_easterness", "DEM_northerness",
         "Slope_50_Curvature", "Topidx_50", "TPI_9")

glm_uni<-univar(indep,tab)
glm_uni

#### Multivariate Logistic Regression ####

library(MASS)

glm_mul<-glm(name~DEM+Slope_50+Walk_NF+DEM_northerness+Slope_50_Curvature+TPI_9,
             data=tab,family=binomial(logit))

glmfile<-capture.output(summary(glm_mul))
glmfile

#### Akaike Information Criterion ####

AIC_glm<-stepAIC(glm_mul)
summary(AIC_glm)
aicfile<-capture.output(summary(AIC_glm))

#### Bayesian Information Criterion ####

BIC_glm<-stepAIC(glm_mul,k=log(length(tab$name)))
summary(BIC_glm)
bicfile<-capture.output(summary(BIC_glm))

#### Coefficients Standardization ####
#crf. Vaughn & Crawford 2009: 550

library(fmsb)

summary(BIC_glm)$coef[-1,1]
(BIC_glm)$coef

bx<-summary(BIC_glm)$coef[-1,1]
bx

R<-NagelkerkeR2(BIC_glm)
R<-R$R2
R

dem<-sd(as.numeric(unlist(BIC_glm$model[2])))
dem
walk<-sd(as.numeric(unlist(BIC_glm$model[3])))
walk
north<-sd(as.numeric(unlist(BIC_glm$model[4])))
north
tpi<-sd(as.numeric(unlist(BIC_glm$model[5])))
tpi
sx<-cbind(dem,walk,north,tpi)
sx

sy<-sd(as.numeric(unlist(BIC_glm$model[1])))
sy

(bx*sx*R)/sy

#### Variance inflation factor ####

library(car)

vif_MUL<-vif(glm_mul)
vifile<-capture.output(vif_MUL)

vif_AIC<-vif(AIC_glm)
vifileAIC<-capture.output(vif_AIC)

vif_BIC<-vif(BIC_glm)
vifileBIC<-capture.output(vif_BIC)

#### AUC - ROC ####

library(pROC)

roc<-roc(BIC_glm$y,BIC_glm$fitted.values)
capture.output(roc)

#### Residuals Analysis (BIC) ####

BIC_resid<-residuals(BIC_glm)
bicresfile<-capture.output(summary(BIC_resid))

Pearson_resid<-residuals(BIC_glm, "pearson")
pearresfile<-capture.output(summary(Pearson_resid))

site_coor<-crds(siti, df=TRUE, list=FALSE)
nosite_coor<-crds(random, df=TRUE, list=FALSE)
coor<-as.matrix(rbind(site_coor,nosite_coor))

library(pgirmess)

resid_correl<-correlog(coor,BIC_resid,method="Moran")
plot(resid_correl)

#### Predictive Surface ####

capture.output(BIC_glm$coefficients)
BIC_glm$coefficients[1]
BIC_glm$coefficients[2]
BIC_glm$coefficients[3]
BIC_glm$coefficients[4]
BIC_glm$coefficients[5]

#Implemented in GRASS:
#r.mapcalc expression="Regression =  0.2975856742  + (DEM@BVER_50 * -0.0038020775) + ( Walk_NF@BVER_50 * -0.0001349823 ) + ( DEM_northerness@BVER_50 * -0.2269732731 ) + ( TPI_9@BVER_50 * 0.4136050599  )" --overwrite
#r.mapcalc expression="Probability_Surface = ( exp ( Regression@BVER_50 ) ) / ( 1 + ( exp ( Regression@BVER_50 ) ) )""

#### Predictive Surface Data ####

predsurf<-read_RAST("Probability_Surface")

sitespred<-extract(predsurf, siti)
sitespred<-select(sitespred,-ID)
colnames(sitespred)<-"pred_mod"

nositespred<-extract(predsurf, random)
nositespred<-select(nositespred,-ID)
colnames(nositespred)<-"pred_mod"

#### Kolmogorov-Smirnov Test (Pred. Values)####

library(Matching)

kspredfile<-capture.output(summary(ks.boot(sitespred,nositespred,10000)))

#### Upload raster with predictive steps ####
#Created in GRASS with Map Algebra
#e.g. r.mapcalc expression="Prob01 = Probability_Surface@BVER_50 >= 0.1"

Cat_Pred<-read_RAST(c("Prob01","Prob02","Prob03",
                      "Prob04","Prob05","Prob06",
                      "Prob07","Prob08","Prob09"))

#### Kvamme's Gain Calculation ####

#sites in each predictive step

sitespredcat<-extract(Cat_Pred, siti)
sitespredcat<-select(sitespredcat, -ID)
totsites<-colSums(sitespredcat, na.rm=TRUE)

#no_sites in sites area

nositespredcat<-extract(Cat_Pred, random)
nositespredcat
nositespredcat<-select(nositespredcat, -ID)
wrong_ns<-colSums(nositespredcat,na.rm=TRUE)

#CORRECT SITES 

#sites presence % for each predictive step

percsit<-round((totsites/length(siti)),digits=2)
percsit

##steps area

areapredcat1<-freq(Cat_Pred,bylayer=T,usenames=TRUE, value=1)
areapredcat0<-freq(Cat_Pred,bylayer=T,usenames=TRUE, value=0)
areapredcat<-rbind(areapredcat1, areapredcat0)
areapredcat<-areapredcat %>% 
  group_by(layer) %>% 
  mutate(area = sum(count)) %>%
  filter (value > 0) %>%
  mutate (percarea = round(count/area,2))
percarea<-areapredcat$percarea
percarea

#Gain calculation

Sites_KG<-(1-(percarea/percsit))

KG_sites<-rbind(percarea,percsit,Sites_KG)
colnames(KG_sites)<-c(">0.1",">0.2",">0.3",">0.4",">0.5",">0.6",">0.7",
                      ">0.8",">0.9")
row.names(KG_sites)<-c("% Area","% Siti","Kvamme's Gain")
round(KG_sites, digits=2)

#WRONG NO-SITES 

#wrongly predicted no-sites %

percnosit<-round((wrong_ns/length(random)),digits=2)  
percnosit 

#Gain calculation 

NSS_KG<-(1-(percarea/percnosit))

KG_ns_sites<-rbind(percarea,percnosit,NSS_KG)
colnames(KG_ns_sites)<-c(">0.1",">0.2",">0.3",">0.4",">0.5",">0.6",">0.7",
                         ">0.8",">0.9")
row.names(KG_ns_sites)<-c("% Area","% No Sites","Kvamme's Gain")
round(KG_ns_sites, digits=2)

#CORRECT NO-SITES

right_ns<-length(random)-colSums(nositespredcat, na.rm=TRUE)
right_ns

(wrong_ns+right_ns)/length(random)

percnono<-round((right_ns/length(random)),digits=2)  

##un-predicted area

areanopredcat1<-freq(Cat_Pred,bylayer=T,usenames=TRUE,value=1)
areanopredcat0<-freq(Cat_Pred,bylayer=T,usenames=TRUE,value=0)
areanopredcat<-rbind(areanopredcat1,areanopredcat0)
areanopredcat<-areanopredcat %>% 
  group_by(layer) %>% 
  mutate(area = sum(count)) %>%
  filter (value == 0) %>%
  mutate (percnoarea = round(count/area,2))
percnoarea<-areanopredcat$percnoarea

#Gain Calculation

NSNS_KG<-(1-(percnoarea/percnono))

KG_ns_ns<-rbind(percnoarea,percnono,NSNS_KG)
colnames(KG_ns_ns)<-c(">0.1",">0.2",">0.3",">0.4",">0.5",">0.6",">0.7",
                      ">0.8",">0.9")
row.names(KG_ns_ns)<-c("% Area Unpred","% No Sites","Kvamme's Gain")
round(KG_ns_ns, digits=2)

####CUMULATIVE PREDICTION PLOT####

#sites
percsit["Prob0"]<-1
percsit["Prob10"]<-0
percfull<-as.data.frame(percsit)
percprob<-rbind(percfull[10,],percfull[1,],percfull[2,],percfull[3,],
                percfull[4,],percfull[5,],percfull[6,],percfull[7,],
                percfull[8,],percfull[9,],percfull[11,])

#no sites
percnono["Prob0"]<-1
percnono["Prob10"]<-0
percfullns<-as.data.frame(percnono)
percprobns<-rbind(percfullns[11,],percfullns[1,],percfullns[2,],percfullns[3,],
                  percfullns[4,],percfullns[5,],percfullns[6,],percfullns[7,],
                  percfullns[8,],percfullns[9,],percfullns[10,])

#table
cutoff<-c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)

tabperc<-cbind(cutoff,(percprob*100),(percprobns*100))
colnames(tabperc)<-c("cutoff","%sites","%random")

#plot

jpeg("Correct_Prediction.jpg",width=2000,height=2000,units="px",
     res=300)
par(oma= c(1, 0, 1, 0))
plot(tabperc[,2]~tabperc[,1],xlab="predicted site probability",ylab="percent
correct predictions",main="Correct Prediction %",sub="Sites in red & Non-Sites in black",
     col="red",pch = 19,axes=F,font.sub=3)
lines(tabperc[,2]~tabperc[,1],col="red")
par(new=TRUE)
plot(tabperc[,3]~tabperc[,1],xlab="",ylab="",pch = 19,axes=F)
lines(tabperc[,3]~tabperc[,1])
abline(v=0.3,lty=2,col="forestgreen")
abline(v=0.6,lty=2,col="firebrick")
Axis(side = 1, at = seq(0,1,by = 0.1), labels=seq(0,1,by = 0.1), cex=0.2)
Axis(side = 1, at = seq(0,1,by = 0.05),labels=F)
Axis(side = 2, at = seq(0,100,by = 10), labels=seq(0,100,10),las=1)
Axis(side = 2, at = seq(0,100,by = 5), labels=F)
dev.off()

#### EXTERNAL CONTROL SITES ####

Pred_control<-read_RAST('Probability_Surface_control')
predsurf<-read_RAST("Probability_Surface")

siti_como<-read_VECT('Siti_Como')
siti_lecco<-read_VECT('Siti_Lecco')
siti_garda<-read_VECT('Siti_Garda')
siti_bver<-read_VECT('Siti_pixel')

pred_como<-data.frame(extract(Pred_control,siti_como),name=rep('como',length(siti_como)))
pred_lecco<-data.frame(extract(Pred_control,siti_lecco),name=rep('lecco',length(siti_lecco)))
pred_garda<-data.frame(extract(Pred_control,siti_garda),name=rep('garda',length(siti_garda)))
pred_bver<-data.frame(extract(predsurf,siti_bver),name=rep('bver',length(siti_bver)))
colnames(pred_lecco)<-c("ID","pred_mod","area")
colnames(pred_garda)<-c("ID","pred_mod","area")
colnames(pred_bver)<-c("ID","pred_mod","area")
colnames(pred_como)<-c("ID","pred_mod","area")
control_sites<-rbind(pred_como,pred_lecco,pred_garda,pred_bver)
control_sites$area<-factor(control_sites$area, levels = c('bver','como','lecco','garda'))
  
  jpeg("Predictive_Control_Sites.jpg",width=3000,height=3000, units="px",res=300)
ggplot(control_sites, aes(y=pred_mod, group=area, fill=area))+
  stat_boxplot(geom ='errorbar') + 
  geom_boxplot()+
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.1), expand = c(0, 0))+
  labs(y='Predictive Value')+
  scale_fill_manual(values =c("cornflowerblue","gold3","lightcoral","olivedrab3"))+
  ggtitle('Valori Predittivi Siti GNA')+
  geom_hline(yintercept = 0.3,linetype='longdash',colour='red4', linewidth=1,5)+
  geom_hline(yintercept = 0.6,linetype='longdash',colour='red4', linewidth=1,5)+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()

#Kolmogorov-Smirnov test (Preditcive Values)
#Each control dataset against original BVER dataset

library(Matching)

bver_garda_ks<-capture.output(summary(ks.boot(pred_bver$pred_mod,na.omit(pred_garda$pred_mod),10000)))
bver_como_ks<-capture.output(summary(ks.boot(pred_bver$pred_mod,pred_como$pred_mod,10000)))
bver_lecco_ks<-capture.output(summary(ks.boot(pred_bver$pred_mod,na.omit(pred_lecco$pred_mod),10000)))
como_lecco_ks<-capture.output(summary(ks.boot(pred_como$pred_mod,na.omit(pred_lecco$pred_mod),10000)))
como_garda_ks<-capture.output(summary(ks.boot(pred_como$pred_mod,na.omit(pred_garda$pred_mod),10000)))
lecco_garda_ks<-capture.output(summary(ks.boot(na.omit(pred_lecco$pred_mod),na.omit(pred_garda$pred_mod),10000)))

