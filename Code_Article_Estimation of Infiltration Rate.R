#Preparing libraries
library(spdep)
library(ape)
library(sp)
library(corrplot)
library(spatialreg)
library(MVA)
library(Hmisc)
library(nortest)
library(readxl)
library(psych) 

#Data reading
XD<- read_excel("C:/Users/ASUS/Documents/Maestria en Geomatica/TESIS 2023/Datos/BD Segunda Version Articulo- Poligonos URH/Matriz Modelo Final.xlsx",sheet="X")

#Extracting X,Y coordinates
XYdata<- read_excel("C:/Users/ASUS/Documents/Maestria en Geomatica/TESIS 2023/Datos/BD Segunda Version Articulo- Poligonos URH/Matriz Modelo Final.xlsx",sheet="XY")

#Creating Data Frame
datx=data.frame(XD)
names(datx)

#Response variable
LNINF=datx$INFILTRACION_LOG

#Creating distance matrix
LI.d=as.matrix(dist(XYdata, diag=T, upper=T))
LI.d.inv <-as.matrix(1/LI.d**2)
diag(LI.d.inv) <- 0
W=as.matrix(LI.d.inv)
SUMAS=apply(W,1,sum)
We=W/SUMAS
apply(We,1,sum)

#Creating distance matrix required for the Spatial Reg library
contnb=dnearneigh(coordinates(XYdata),0,380000,longlat = F)
dlist <- nbdists(contnb, XYdata)
dlist <- lapply(dlist, function(x) 1/x**2)
Wve=nb2listw(contnb,glist=dlist,style = "W")
Wve

#Final Model-SARAR
SARAR1=sacsarlm(formula=INFILTRACION_LOG~CO+DAMC+RVI_6Marzo_Media+
                  SMI_1Noviembre_Media+SMI_13Noviembre_Media+DEM+A+COBFRSUP
                ,data=datx,listw=Wve)

#Summary of the model
summary(SARAR1)

#Normality in the results
resi5=SARAR1$residuals
shapiro.test(resi5)
ad.test(resi5)
cvm.test(resi5)
lillie.test(resi5)

#Estimated values
SARAR1$fitted.values

#Plot of observed vs. estimated values
LIE5=as.data.frame(SARAR1$fitted.values)
DF5=data.frame(datx$INFILTRACION_LOG,LIE5)
colnames(DF5) <- c("LI","LIE5") 
plot(DF5$LI,DF5$LIE5,cex=0.5,pch=19, xlab="Valores observados Ln Infiltración",ylab="Valores estimados Ln Infiltración ",
     col = c("green", "blue"))
grid(5,5)
legend("bottomright", legend = c("Observado", "Estimado"), col = c("green", "blue"), pch = c(16, 16))

#Pearson correlation
cor(DF5$LI,DF5$LIE5)

#Spatial dependence in residuals
moran.plot(resi5,Wve,pch=19,cex=0.5)
Moranresi4=(Moran.I(resi5,We))
Moranresi4

#Mean square error
RMSEm5 <- sqrt(mean(((SARAR1$fitted.values) - LNINF)^2))
RMSEm5

#Descriptive statistics
#Ln Infiltracion
describe(datx[,1])
#CO
describe(datx[,3])
#DAMC
describe(datx[,5])
#RVI 6 Marzo
describe(datx[,6])
#SMI 1 Noviembre
describe(datx[,7])
#SMI 13 Noviembre
describe(datx[,8])
#DEM
describe(datx[,9])
#A
describe(datx[,2])
#COBFRSUP
describe(datx[,4])

#Exploration of spatial dependence
#Moran's Index
#Ln Infiltracion
MoranMO=(Moran.I(datx[,1], We));MoranMO
#CO
MoranMO1=(Moran.I(datx[,3], We));MoranMO1
#DAMC
MoranMO2=(Moran.I(datx[,5], We));MoranMO2
#RVI 6 Marzo
MoranMO3=(Moran.I(datx[,6], We));MoranMO3
#SMI 1 Noviembre
MoranMO4=(Moran.I(datx[,7], We));MoranMO4
#SMI 13 Noviembre
MoranMO5=(Moran.I(datx[,8], We));MoranMO5
#DEM
MoranMO6=(Moran.I(datx[,9], We));MoranMO6
#A
MoranMO7=(Moran.I(datx[,2], We));MoranMO7
#COBFRSUP
MoranMO8=(Moran.I(datx[,4], We));MoranMO8

#Impacts
W <- as(Wve, "CsparseMatrix")
trMatb <- trW(W, type="mult")
trMC <- trW(W, type="MC")
summary(impacts(SARAR1, tr=trMatb, R=200), zstats=TRUE, short=TRUE)

#SARAR model changing weights

#Creating distance matrix
LI.d=as.matrix(dist(XYdata, diag=T, upper=T))
LI.d.inv <-as.matrix(1/LI.d**2)
diag(LI.d.inv) <- 0
W=as.matrix(LI.d.inv)
SUMAS=apply(W,1,sum)
We=W/SUMAS
apply(We,1,sum)
#Function to generate weight-dependent models
sararm <- function(b) {
  #Creating distance matrix required for the Spatial Reg library
  contnb=dnearneigh(coordinates(XYdata),0,380000,longlat = F)
  dlist <- nbdists(contnb, XYdata)
  dlist <- lapply(dlist, function(x) 1/x**b)
  Wve=nb2listw(contnb,glist=dlist,style = "W")
  #Final Model-SARAR
  SARAR1=sacsarlm(formula=INFILTRACION_LOG~CO+DAMC+RVI_6Marzo_Media+
                    SMI_1Noviembre_Media+SMI_13Noviembre_Media+DEM+A+COBFRSUP
                  ,data=datx,listw=Wve)
  #Summary of the model
  summary(SARAR1)
}

#Loop to extract parameters of interest from models

modelos=list(); rho=c();p_rho=c();lambda=c();p_lambda=c();coef=list();
p_CO=c();p_DAMC=c();p_RVI_6Marzo_Media=c();p_SMI_1Noviembre_Media=c();
p_SMI_13Noviembre_Media=c();p_DEM=c();p_A=c();p_COBFRSUP=c();corr=c();
Morani=c();Normalidad=c();pesos=c()
for(i in 1:200){
  modelos[[i]]=sararm(0.2+i/100)
  rho[i]=modelos[[i]]$rho
  lambda[i]=modelos[[i]]$lambda
  p_rho[i]=(1-pnorm(modelos[[i]]$rho/modelos[[i]]$rho.se))*2
  p_lambda[i]=(1-pnorm(modelos[[i]]$lambda/modelos[[i]]$lambda.se))*2
  coef[[i]]=as.numeric(modelos[[i]]$coefficients)
  corr[i]=cor(modelos[[i]]$fitted.values,XD$INFILTRACION_LOG)
  Morani[i]=Moran.I(modelos[[i]]$residuals,We**(0.2+i/100))$p.value
  pesos[i]=(We**((0.2+i/100)))[1,which.max((We**(0.2+i/100))[1,])]
  Normalidad[i]=shapiro.test(modelos[[i]]$residuals)$p.value
  p_CO[i]=(pnorm(coef[[i]][2]/modelos[[i]]$rest.se[2]))
  p_DAMC[i]=(pnorm(coef[[i]][3]/modelos[[i]]$rest.se[3]))
  p_RVI_6Marzo_Media[i]=(1-pnorm(coef[[i]][4]/modelos[[i]]$rest.se[4]))*2
  p_SMI_1Noviembre_Media[i]=(1-pnorm(coef[[i]][5]/modelos[[i]]$rest.se[5]))*2
  p_SMI_13Noviembre_Media[i]=(pnorm(coef[[i]][6]/modelos[[i]]$rest.se[6]))
  p_DEM[i]=(1-pnorm(coef[[i]][7]/modelos[[i]]$rest.se[7]))*2
  p_A[i]=(1-pnorm(coef[[i]][8]/modelos[[i]]$rest.se[8]))*2
  p_COBFRSUP[i]=(pnorm(coef[[i]][9]/modelos[[i]]$rest.se[9]))
}

# Defining colors for graphics

colorr=ifelse(p_rho<0.1,"red","blue")
colorl=ifelse(p_lambda<0.1,"red","blue")
colorCO=ifelse(p_CO<0.1,"red","blue")
colorDAMC=ifelse(p_DAMC<0.1,"red","blue")
colorRVI_6Marzo_Media=ifelse(p_RVI_6Marzo_Media<0.1,"red","blue")
colorSMI_1Noviembre_Media=ifelse(p_SMI_1Noviembre_Media<0.1,"red","blue")
colorSMI_13Noviembre_Media=ifelse(p_SMI_13Noviembre_Media<0.1,"red","blue")
colorDEM=ifelse(p_DEM<0.1,"red","blue")
colorA=ifelse(p_A<0.1,"red","blue")
colorCOBFRSUP=ifelse(p_COBFRSUP<0.1,"red","blue")
colorMORAN=ifelse(Morani<0.1,"red","blue")
colorNormalidad=ifelse(Normalidad<0.1,"red","blue")
dff=data.frame(rho=rho,lambda=lambda, b=seq(0.21,2.2,0.01),pv_rho=p_rho, pv_lambda=p_lambda,
               pv_CO=p_CO,pv_DAMC=p_DAMC,pv_RVI_6Marzo_Media=p_RVI_6Marzo_Media,
               pv_SMI_1Noviembre_Media=p_SMI_1Noviembre_Media,
               pv_SMI_13Noviembre_Media=p_SMI_13Noviembre_Media,
               pv_DEM=p_DEM,pv_A=p_A,pv_COBFRSUP=p_COBFRSUP,r=corr, pv_MORAN=Morani,
               pv_Shapiro=Normalidad, peso_max=pesos)

# Graphics

ggplot(dff, aes(x=b, y=rho, color=colorr)) +
  geom_point()
ggplot(dff, aes(x=b, y=lambda, color=colorl)) +
  geom_point()
ggplot(dff, aes(x=b, y=-log10(pv_CO), color=colorCO)) +
  geom_point()
ggplot(dff, aes(x=b, y=-log(pv_DAMC), color=colorDAMC)) +
  geom_point()
ggplot(dff, aes(x=b, y=-log10(pv_RVI_6Marzo_Media), color=colorRVI_6Marzo_Media)) +
  geom_point()
ggplot(dff, aes(x=b, y=-log10(pv_SMI_1Noviembre_Media), color=colorSMI_1Noviembre_Media)) +
  geom_point()
ggplot(dff, aes(x=b, y=-log10(pv_SMI_13Noviembre_Media), color=colorSMI_13Noviembre_Media)) +
  geom_point()
ggplot(dff, aes(x=b, y=-log10(pv_DEM), color=colorDEM)) +
  geom_point()
ggplot(dff, aes(x=b, y=-log10(pv_A), color=colorA)) +
  geom_point()
ggplot(dff, aes(x=b, y=-log10(pv_COBFRSUP), color=colorCOBFRSUP)) +
  geom_point()
ggplot(dff, aes(x=b, y=100*corr**2)) +
  geom_point(colour="green")+
  labs(y = "Correlation between observed and estimated values")
ggplot(dff, aes(x=b, y=pv_MORAN)) +
  geom_point(colour=colorr)+
  labs(y = "p value (Moran Index)")
ggplot(dff, aes(x=b, y=pv_Shapiro)) +
  geom_point()
ggplot(dff, aes(x=b, y=peso_max)) +
  geom_point(colour="red")+
  labs(y = "Maximum weight")

  
  