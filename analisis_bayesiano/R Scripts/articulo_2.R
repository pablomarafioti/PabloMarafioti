setwd("C:/Users/Juan/Desktop/PABLO/doctorado/tesis_de_doctorado")
load("regresion_2.RData")

#devtools::install_github("glmmTMB/glmmTMB/glmmTMB")
#devtools::install_github("strengejacke/ggeffects")

source("func_aux.R")
library(broom)
library(broom.mixed)
library(lsmeans)
library(WriteXLS)
library(readxl)
#library(stringdist)
#library(quanteda)
library(corrplot)
library(languageR)
library(car)
#library(lmtest)
library("sjPlot")
library("sjstats")
library(Hmisc)
library("PerformanceAnalytics")
library(geepack)
library(MuMIn)
library("ggplot2")  
library("aods3")
library(reshape)
library(lattice)
library("gridExtra")
source("allFit.r")
library(caret)
library(clickR)
library("xtable")
#library(koRpus)

#load("datos_impute_total.RData")

#library(nnet)


# datos para el modelo

load("datos_articulo.RData")


############# multinomial con nnet (selección de variables)

library(nnet)
library(MuMIn)

multi.fit.1 = multinom(RES_CAT ~ Fabs.SC.f+MORF.f+STEM.f+MOD+ES+ANIM+GRAMS+FAM.LEX.f+IMA.CONC.f+CUMRES.f+LDA
                       , data=datos)

multi.fit.2 = multinom(RES_CAT ~ POS + DIS + EST1 + EST2 + EST3 + EST4 + EST5 + EST6 + EST7 + GRUPO6, data=datos)




# las siguientes tres líneas realizan la selección de modelos. 
# Los resultados están gardados en:
# "seleccion_grupo_1.RData" y "seleccion_grupo_2.RData"

load("seleccion_grupo_1.RData")
load("seleccion_grupo_2.RData")

#### corrida de seleccón de modelos

options(na.action = "na.fail")

modelos.m.1  <- dredge(multi.fit.1, rank = "AIC")
modelos.m.2  <- dredge(multi.fit.2, rank = "AIC")

####


sub.m.1   <- subset(modelos.m.1, 1/8 < weight/max(modelos.m.1$weight))[, c("logLik","AIC", "delta", "weight")]
importance.1 <- importance(modelos.m.1)
conf.av.1 <- summary(model.avg(subset(modelos.m.1, 1/8 < weight/max(modelos.m.1$weight))))

# importancia de las variables: grupo 1
#
#                     Fabs.SC.f ANIM MORF.f CUMRES.f FAM.LEX.f MOD  ES  
#Sum of weights:      1.00      1.00 1.00   1.00     1.00      0.99 0.90
#N containing models: 1024      1024 1024   1024     1024      1024 1024
#                    IMA.CONC.f LDA  STEM.f GRAMS
#Sum of weights:      0.72       0.44 0.22   0.11 
#N containing models: 1024       1024 1024   1024 

#Promedio de los coeficientes con FULL AVERAGE: grupo 1

# Model-averaged coefficients:  
#   (full average) 
# Estimate Std. Error z value Pr(>|z|)    
# 1((Intercept))  -3.798675   0.691426   5.494  < 2e-16 ***
# 1(ANIM1)         0.186352   0.340951   0.547 0.584677    
# 1(CUMRES.f1)     0.483107   0.341845   1.413 0.157586    
# 1(CUMRES.f2)     0.348968   0.356727   0.978 0.327951    
# 1(ES1)           0.474374   0.609295   0.779 0.436237    
# 1(ES2)         -14.654605   7.663950   1.912 0.055857 .  
# 1(FAM.LEX.f1)   -0.030567   0.274903   0.111 0.911464    
# 1(Fabs.SC.f1)   -1.171358   0.319983   3.661 0.000251 ***
# 1(IMA.CONC.f1)   0.546806   0.397111   1.377 0.168525    
# 1(MOD1)         -0.739842   1.056345   0.700 0.483691    
# 1(MOD2)          0.645929   0.309275   2.089 0.036751 *  
# 1(MOD3)         -0.628208   0.429849   1.461 0.143889    
# 1(MORF.f1)       0.522214   0.534416   0.977 0.328486    
# 1(MORF.f2)      -0.056713   0.688804   0.082 0.934380    
# 2((Intercept))  -3.543079   0.923904   3.835 0.000126 ***
# 2(ANIM1)         1.832519   0.431428   4.248 2.16e-05 ***
# 2(CUMRES.f1)     0.506814   0.445414   1.138 0.255184    
# 2(CUMRES.f2)    -0.348771   0.555664   0.628 0.530223    
# 2(ES1)           0.232550   0.790903   0.294 0.768735    
# 2(ES2)           0.815147   1.343208   0.607 0.543940    
# 2(FAM.LEX.f1)   -1.316818   0.424472   3.102 0.001921 ** 
# 2(Fabs.SC.f1)    0.856733   0.531008   1.613 0.106656    
# 2(IMA.CONC.f1)  -0.588362   0.483260   1.217 0.223419    
# 2(MOD1)         -8.792169 237.281485   0.037 0.970442    
# 2(MOD2)          0.345775   0.517900   0.668 0.504358    
# 2(MOD3)          0.208514   0.684397   0.305 0.760619    
# 2(MORF.f1)      -1.441987   0.760852   1.895 0.058062 .  
# 2(MORF.f2)      -0.982811   0.722167   1.361 0.173539    
# 3((Intercept))  -0.724484   0.294401   2.461 0.013860 *  
# 3(ANIM1)         0.102687   0.164525   0.624 0.532534    
# 3(CUMRES.f1)     0.261252   0.174742   1.495 0.134895    
# 3(CUMRES.f2)     0.469879   0.173292   2.711 0.006698 ** 
# 3(ES1)          -0.493425   0.298111   1.655 0.097889 .  
# 3(ES2)         -18.804264   8.831732   2.129 0.033240 *  
# 3(FAM.LEX.f1)   -0.330112   0.137103   2.408 0.016050 *  
# 3(Fabs.SC.f1)   -0.402001   0.167858   2.395 0.016625 *  
# 3(IMA.CONC.f1)  -0.034270   0.134761   0.254 0.799263    
# 3(MOD1)         -0.349807   0.557954   0.627 0.530694    
# 3(MOD2)          0.353311   0.175064   2.018 0.043572 *  
# 3(MOD3)          0.481862   0.198703   2.425 0.015307 *  
# 3(MORF.f1)      -1.089739   0.230864   4.720 2.40e-06 ***
# 3(MORF.f2)      -0.690523   0.303956   2.272 0.023099 *  
# 4((Intercept))  -3.840682   0.598226   6.420  < 2e-16 ***
# 4(ANIM1)         0.329694   0.267250   1.234 0.217333    
# 4(CUMRES.f1)     1.114311   0.363315   3.067 0.002162 ** 
# 4(CUMRES.f2)     1.296583   0.364776   3.554 0.000379 ***
# 4(ES1)           0.640538   0.445816   1.437 0.150782    
# 4(ES2)           0.577832   0.765561   0.755 0.450379    
# 4(FAM.LEX.f1)   -0.464110   0.229848   2.019 0.043466 *  
# 4(Fabs.SC.f1)   -0.391212   0.287804   1.359 0.174052    
# 4(IMA.CONC.f1)  -0.088896   0.225493   0.394 0.693412    
# 4(MOD1)         -0.334067   1.063878   0.314 0.753514    
# 4(MOD2)          1.073710   0.293452   3.659 0.000253 ***
# 4(MOD3)          0.641313   0.348480   1.840 0.065722 .  
# 4(MORF.f1)      -0.259490   0.428594   0.605 0.544884    
# 4(MORF.f2)      -0.346109   0.453848   0.763 0.445696    
# 1(LDA1)         -0.421387   0.849390   0.496 0.619820    
# 2(LDA1)         -0.124367   0.761539   0.163 0.870274    
# 3(LDA1)          0.203100   0.301416   0.674 0.500426    
# 4(LDA1)         -0.477541   0.739245   0.646 0.518290    
# 1(STEM.f1)       0.041864   0.164849   0.254 0.799531    
# 2(STEM.f1)       0.139224   0.349838   0.398 0.690655    
# 3(STEM.f1)       0.006307   0.073628   0.086 0.931734    
# 4(STEM.f1)       0.069017   0.188854   0.365 0.714774    
# 1(GRAMS1)       -0.001542   0.071557   0.022 0.982805    
# 2(GRAMS1)       -0.014626   0.144079   0.102 0.919144    
# 3(GRAMS1)       -0.001894   0.036999   0.051 0.959171    
# 4(GRAMS1)       -0.024954   0.131935   0.189 0.849986   
# 


sub.m.2   <- subset(modelos.m.2, 1/8 < weight/max(modelos.m.2$weight))[, c("logLik","AIC", "delta", "weight")]
importance.2 <- importance(modelos.m.2)
conf.av.2 <- summary(model.avg(subset(modelos.m.2, 1/8 < weight/max(modelos.m.2$weight))))

# importancia de las variables: grupo 1
#
#                     EST5 EST1 GRUPO6 EST2 EST6 EST7 EST3 EST4 DIS  POS 
#Sum of weights:      1.00 0.98 0.96   0.93 0.82 0.75 0.62 0.48 0.27 0.06
#N containing models:  512  512  512    512  512  512  512  512  512  512


# Model-averaged coefficients:  
#   (full average) 
# Estimate Std. Error z value Pr(>|z|)    
# 1((Intercept))  -2.90363    0.32448   8.949  < 2e-16 ***
# 1(EST11)        -1.69340    0.53435   3.169  0.00153 ** 
# 1(EST21)         0.15674    0.35839   0.437  0.66185    
# 1(EST31)        -7.20105    5.14786   1.399  0.16186    
# 1(EST41)         0.41010    0.60983   0.672  0.50128    
# 1(EST51)        -0.19183    0.52880   0.363  0.71678    
# 1(EST61)         0.68508    0.53260   1.286  0.19834    
# 1(EST71)         0.87443    0.66117   1.323  0.18598    
# 1(GRUPO62)      -0.82030    0.54830   1.496  0.13463    
# 1(GRUPO63)       0.43508    0.32206   1.351  0.17672    
# 1(GRUPO64)      -0.03435    0.41953   0.082  0.93475    
# 1(GRUPO65)      -3.11112    4.13368   0.753  0.45167    
# 1(GRUPO66)      -0.60006    0.78751   0.762  0.44608    
# 2((Intercept))  -2.94567    0.37032   7.954  < 2e-16 ***
# 2(EST11)        -0.64825    0.47094   1.376  0.16867    
# 2(EST21)        -2.62878    1.05326   2.496  0.01257 *  
# 2(EST31)        -5.98897  282.54167   0.021  0.98309    
# 2(EST41)         0.03399    0.41985   0.081  0.93548    
# 2(EST51)       -24.21786   47.11784   0.514  0.60726    
# 2(EST61)       -16.56482  157.95212   0.105  0.91648    
# 2(EST71)        -0.17700    0.61361   0.288  0.77300    
# 2(GRUPO62)       0.02489    0.62884   0.040  0.96843    
# 2(GRUPO63)      -0.82864    0.64747   1.280  0.20061    
# 2(GRUPO64)      -0.53474    0.64869   0.824  0.40974    
# 2(GRUPO65)      -2.25393   37.91483   0.059  0.95260    
# 2(GRUPO66)       1.41916    0.54584   2.600  0.00932 ** 
# 3((Intercept))  -1.84529    0.17564  10.506  < 2e-16 ***
# 3(EST11)        -0.08324    0.19464   0.428  0.66889    
# 3(EST21)         0.18217    0.19818   0.919  0.35798    
# 3(EST31)         1.54847    1.31764   1.175  0.23992    
# 3(EST41)        -0.30223    0.47310   0.639  0.52294    
# 3(EST51)         0.03291    0.24880   0.132  0.89478    
# 3(EST61)         0.45273    0.33298   1.360  0.17394    
# 3(EST71)        -0.24058    0.38954   0.618  0.53684    
# 3(GRUPO62)       0.15625    0.25095   0.623  0.53351    
# 3(GRUPO63)       0.54432    0.16716   3.256  0.00113 ** 
# 3(GRUPO64)       0.19519    0.21233   0.919  0.35795    
# 3(GRUPO65)       0.97393    0.30012   3.245  0.00117 ** 
# 3(GRUPO66)      -0.04683    0.43715   0.107  0.91469    
# 4((Intercept))  -2.91423    0.31587   9.226  < 2e-16 ***
# 4(EST11)        -0.68363    0.40040   1.707  0.08775 .  
# 4(EST21)         0.24025    0.35177   0.683  0.49462    
# 4(EST31)         1.27501    1.37289   0.929  0.35304    
# 4(EST41)         0.42365    0.56362   0.752  0.45226    
# 4(EST51)        -0.97053    0.45981   2.111  0.03479 *  
# 4(EST61)         0.97582    0.52002   1.877  0.06059 .  
# 4(EST71)         0.88708    0.57227   1.550  0.12111    
# 4(GRUPO62)       0.92972    0.36932   2.517  0.01182 *  
# 4(GRUPO63)       0.19720    0.30780   0.641  0.52173    
# 4(GRUPO64)       0.06694    0.36490   0.183  0.85444    
# 4(GRUPO65)      -2.49226    4.06373   0.613  0.53968    
# 4(GRUPO66)       0.39790    0.56618   0.703  0.48219    
# 1(DIS1)         -2.72609    3.96193   0.688  0.49141    
# 1(DIS2)          2.43663    4.75543   0.512  0.60838    
# 1(DIS3)         -1.98565    3.56355   0.557  0.57738    
# 2(DIS1)         -2.75502    7.55272   0.365  0.71528    
# 2(DIS2)          2.45189   37.99968   0.065  0.94855    
# 2(DIS3)         -1.74712   73.88021   0.024  0.98113    
# 3(DIS1)          0.06778    0.20115   0.337  0.73614    
# 3(DIS2)          0.05752    0.28078   0.205  0.83770    
# 3(DIS3)          0.03746    0.26804   0.140  0.88886    
# 4(DIS1)         -2.53086    4.23569   0.598  0.55017    
# 4(DIS2)          2.45208    4.77291   0.514  0.60743    
# 4(DIS3)         -2.05883    4.03105   0.511  0.60953  


# tablas de latex

xtbl.importance.1 <- xtable(tidy(importance.1), 
                            caption = "Importancia relativa de las predictoras: grupo 1", digits = 2,label = "xtable:example")

xtbl.seleccion.1 <- xtable(data.frame(sub.m.1), 
              caption = "Seleccion de modelos: grupo 1. 
             logLik: verosimilitud del modelo, AICc: AKAIKE corregido, delta: delta de AKAIKE, weight: pesos de 
            AKAIKE", digits = 3,label = "xtable:example")

xtbl.full.average.1 <- xtable(conf.av.1$coefmat.full, 
                              caption = "Promedio de los coeficientes con FULL AVERAGE : grupo 1", digits = 3,label = "xtable:example")


xtbl.full.average.1.1 <- xtable(conf.av.1$coefmat.full[conf.av.1$coefmat.full[,4]< 0.05,], 
                                caption = "Promedio de los coeficientes con FULL AVERAGE (p < 0.05): grupo 1", digits = 3,label = "xtable:example")


xtbl.importance.2 <- xtable(tidy(importance.2), 
                            caption = "Importancia relativa de las predictoras: grupo 2", digits = 2,label = "xtable:example")

xtbl.seleccion.2 <- xtable(data.frame(sub.m.2), 
                           caption = "Seleccion de modelos: grupo 2. 
    logLik: verosimilitud del modelo, AICc: AKAIKE corregido, delta: delta de AKAIKE, weight: pesos de 
    AKAIKE", digits = 3,label = "xtable:example")

xtbl.full.average.2 <- xtable(conf.av.2$coefmat.full, 
                              caption = "Promedio de los coeficientes con FULL AVERAGE : grupo 2", digits = 3,label = "xtable:example")

xtbl.full.average.2.1 <- xtable(conf.av.2$coefmat.full[conf.av.2$coefmat.full[,4]< 0.05,], 
                                caption = "Promedio de los coeficientes con FULL AVERAGE (p < 0.05): grupo 2", digits = 3,label = "xtable:example")

print(xtbl.full.average.1, include.rownames=TRUE, tabular.environment="longtable", floating=FALSE)
print(xtbl.full.average.2, include.rownames=TRUE, tabular.environment="longtable", floating=FALSE)


########################## Modelo multinomial bayesiano

library("MCMCglmm")

## cada modelo tarda alrededor de 1 día, están guardados en 
## "modelo_1.RData", "modelo_2.RData", "modelo_3.RData"

load("modelo_1.RData")
load("modelo_2.RData")
load("modelo_3.RData")

### corrida de modelos
### las siguientes líneas corren lo modelos (¡¡¡ tardan un día !!!)

model.update <- updateable(MCMCglmm)


k <- length(levels(datos$RES_CAT))
IJ <- (1/k) * (diag(k-1) + matrix(1,k-1,k-1)) # (1/k)*(I+J), J=unit matrix, ambas de dim J-1

#model.update <- updateable(MCMCglmm)

v <- sum(diag(IJ))
prior <- list(R = list(fix=1, V=IJ, nu=1), 
              B = list(mu = rep(0, 88), V = kronecker(IJ, diag(22)) * (v + pi^2/3)),
              G = list(G1 = list(V = diag(k-1), nu= 1000, alpha.mu = rep(0,k-1), alpha.V = diag(k-1))))


global.model.2 <- MCMCglmm(RES_CAT ~ -1 + trait +  trait:(TIME + MOD + Fabs.SC.f + ANIM + ES + MORF.f + FAM.LEX.f + CUMRES.f + EST1 + EST2 + EST5 + GRUPO6),
                           random = ~ idh(trait):ID.SESION, rcov = ~idh(trait):units,
                           prior = prior,
                           burnin = 100000,
                           nitt = 4100000,
                           thin = 2000,
                           Pr = TRUE,
                           family = "categorical",
                           data = datos)

v <- sum(diag(IJ))
pred_k <- length(levels(datos$FAM.LEX.f))
prior <- list(R = list(fix=1, V=IJ, n=k-1),
              B = list(mu = rep(0, 88), V = kronecker(IJ, diag(22)) * (v + pi^2/3)),
              G = list(G1 = list(V = diag(k-1), nu= 1000, alpha.mu = rep(0,k-1), alpha.V = diag(k-1)),
                       G2 = list(V = diag(pred_k), nu= 1000, alpha.mu = rep(0,pred_k), alpha.V = diag(pred_k))
                       ))

global.model.3 <- MCMCglmm(RES_CAT ~ -1 + trait + trait:(TIME + MOD + Fabs.SC.f + ANIM + ES + MORF.f + FAM.LEX.f + CUMRES.f + EST1 + EST2 + EST5 + GRUPO6),
                           random = ~ idh(trait):ID.SESION + idh(FAM.LEX.f):ID.SESION, rcov = ~idh(trait):units,
                           prior = prior,
                           burnin = 100000,
                           nitt = 4100000,
                           thin = 2000,
                           pr = TRUE,
                           family = "categorical",
                           data = datos)


v <- sum(diag(IJ))
prior <- list(R = list(fix=1, V=IJ, nu=1), 
              B = list(mu = rep(0, 88), V = kronecker(IJ, diag(22)) * (v + pi^2/3)),
              G = list(G1 = list(V = diag(5), nu= 1000, alpha.mu = rep(0,5), alpha.V = diag(5))
              ))


global.model.4 <- MCMCglmm(RES_CAT ~ -1 + trait + trait:(TIME + Fabs.SC.f + MOD + ANIM + ES + MORF.f + CUMRES.f +  FAM.LEX.f + EST1 + EST2 + EST5 + GRUPO6),
                           random = ~ idh(trait + TIME):ID.SESION, rcov = ~idh(trait):units,
                           prior = prior,
                           burnin = 100000,
                           nitt = 4100000,
                           thin = 2000,
                           pr = TRUE,
                           family = "categorical",
                           data = datos)

###############

# Cuadro 5
# DIC - ICC (seleccón de modelos, elige modelo 2: "global.model.3")

DIC(global.model.2)

ICC.1.1 <- global.model.2$VCV[,1]/(global.model.2$VCV[,1] + pi^2/3)
posterior.mode(ICC.1.1)
effectiveSize(ICC.1.1)
HPDinterval(ICC.1.1)
ICC.1.2 <- global.model.2$VCV[,2]/(global.model.2$VCV[,2] + pi^2/3)
posterior.mode(ICC.1.2)
effectiveSize(ICC.1.2)
HPDinterval(ICC.1.2)
ICC.1.3 <- global.model.2$VCV[,3]/(global.model.2$VCV[,3] + pi^2/3)
posterior.mode(ICC.1.3)
effectiveSize(ICC.1.3)
HPDinterval(ICC.1.3)
ICC.1.4 <- global.model.2$VCV[,4]/(global.model.2$VCV[,4] + pi^2/3)
posterior.mode(ICC.1.4)
effectiveSize(ICC.1.4)
HPDinterval(ICC.1.4)

DIC(global.model.3)

ICC.2.1 <- global.model.3$VCV[,1]/(global.model.3$VCV[,1] + global.model.3$VCV[,5] + global.model.3$VCV[,6] + pi^2/3)
posterior.mode(ICC.2.1)
effectiveSize(ICC.2.1)
HPDinterval(ICC.2.1)
ICC.2.2 <- global.model.3$VCV[,2]/(global.model.3$VCV[,2] + global.model.3$VCV[,5] + global.model.3$VCV[,6] + pi^2/3)
posterior.mode(ICC.2.2)
effectiveSize(ICC.2.2)
HPDinterval(ICC.2.2)
ICC.2.3 <- global.model.3$VCV[,3]/(global.model.3$VCV[,3] + global.model.3$VCV[,5] + global.model.3$VCV[,6] + pi^2/3)
posterior.mode(ICC.2.3)
effectiveSize(ICC.2.3)
HPDinterval(ICC.2.3)
ICC.2.4 <- global.model.3$VCV[,4]/(global.model.3$VCV[,4] + global.model.3$VCV[,5] + global.model.3$VCV[,6] + pi^2/3)
posterior.mode(ICC.2.4)
effectiveSize(ICC.2.4)
HPDinterval(ICC.2.4)

DIC(global.model.4)

ICC.3.1 <- global.model.4$VCV[,1]/(global.model.4$VCV[,1] + global.model.4$VCV[,5] + pi^2/3)
posterior.mode(ICC.3.1)
effectiveSize(ICC.3.1)
HPDinterval(ICC.3.1)
ICC.3.2 <- global.model.4$VCV[,2]/(global.model.4$VCV[,2] + global.model.4$VCV[,5] + pi^2/3)
posterior.mode(ICC.3.2)
effectiveSize(ICC.3.2)
HPDinterval(ICC.3.2)
ICC.3.3 <- global.model.4$VCV[,3]/(global.model.4$VCV[,3] + global.model.4$VCV[,5] + pi^2/3)
posterior.mode(ICC.3.3)
effectiveSize(ICC.3.3)
HPDinterval(ICC.3.3)
ICC.3.4 <- global.model.4$VCV[,4]/(global.model.4$VCV[,4] + global.model.4$VCV[,5]  + pi^2/3)
posterior.mode(ICC.3.4)
effectiveSize(ICC.3.4)
HPDinterval(ICC.3.4)


############# Análisis de supuestos y resultados (modelo 2)

summary(global.model.3)
heidel.diag(global.model.3$VCV)[2,]
heidel.diag(global.model.3$Sol)

plot(global.model.3$VCV[,4])

autocorr.diag(global.model.3$VCV)
autocorr.diag(global.model.3$Sol[,1:88])


sqrt(summary(global.model.3)$Gcovariances[,1])

# Posterior de las varianzas y sus modas y  SD

apply(global.model.3$VCV[,1:6],2,posterior.mode)

##modas
# traitRES_CAT.1.ID.SESION traitRES_CAT.2.ID.SESION traitRES_CAT.3.ID.SESION traitRES_CAT.4.ID.SESION     FAM.LEX.f0.ID.SESION 
# 1.139763699              0.015831973              0.002720140              0.004658977              0.001918262 
# FAM.LEX.f1.ID.SESION 
# 0.562581843 

# sqrt(moda)
# traitRES_CAT.1.ID.SESION traitRES_CAT.2.ID.SESION traitRES_CAT.3.ID.SESION traitRES_CAT.4.ID.SESION     FAM.LEX.f0.ID.SESION 
# 1.06759716               0.12582517               0.05215496               0.06825670               0.04379797 
# FAM.LEX.f1.ID.SESION 
# 0.75005456 

sqrt(apply(global.model.3$VCV[,1:6],2,posterior.mode))

ref <- data.frame(global.model.3$VCV[,1:6])
rdf <- melt(ref, value.name="variance")
ggplot(rdf, aes(x=variance,color=variable))+geom_density()

#efectos aletorios: moda y media de la posterior 

rmode <- posterior.mode(global.model.3$Sol[,-c(1:88)]) 
rmode.cat.1 <- rmode[1:52]
rmode.cat.2 <- rmode[53:104]
rmode.cat.3 <- rmode[105:156]
rmode.cat.4 <- rmode[157:208]

rmode.fam.lex.0 <- rmode[209:260]
rmode.fam.lex.1 <- rmode[261:312]

rmean <- colMeans(global.model.3$Sol[,-c(1:88)])

rmean.cat.1 <- rmean[1:52]
rmean.cat.2 <- rmean[53:104]
rmean.cat.3 <- rmean[105:156]
rmean.cat.4 <- rmean[157:208]

rmean.fam.lex.0 <- rmean[209:260]
rmean.fam.lex.1 <- rmean[261:312]


rmedian <- apply(global.model.3$Sol[,-c(1:88)],2,median)

rmedian.cat.1 <- rmedian[1:52]
rmedian.cat.2 <- rmedian[53:104]
rmedian.cat.3 <- rmedian[105:156]
rmedian.cat.4 <- rmedian[157:208]

rmedian.fam.lex.0 <- rmedian[209:260]
rmedian.fam.lex.1 <- rmedian[261:312]


cor(rmode.cat.1, rmean.cat.1) # 0.9482757
cor(rmode.cat.2, rmean.cat.2) #  0.5997795
cor(rmode.cat.3, rmean.cat.3) #  0.3292984
cor(rmode.cat.4, rmean.cat.4) #  0.2530744
cor(rmode.fam.lex.0, rmean.fam.lex.0) #  0.04633252
cor(rmode.fam.lex.1, rmean.fam.lex.1) #  0.966191


cor(rmedian.cat.1, rmean.cat.1) #  0.9996711
cor(rmedian.cat.2, rmean.cat.2) #  0.9923571
cor(rmedian.cat.3, rmean.cat.3) #  0.9950963
cor(rmedian.cat.4, rmean.cat.4) #  0.9913889
cor(rmedian.fam.lex.0, rmean.fam.lex.0) #  0.9848327
cor(rmedian.fam.lex.1, rmean.fam.lex.1) #  0.9996111

shapiro.cat.1  <- shapiro.test(rmean.cat.1) # W = 0.91185, p-value = 0.0009489
shapiro.cat.2  <- shapiro.test(rmean.cat.2) # W = 0.83788, p-value = 5.208e-06
shapiro.cat.3  <- shapiro.test(rmean.cat.3) # W = 0.98548, p-value = 0.7733
shapiro.cat.4  <- shapiro.test(rmean.cat.4) # W = 0.97815, p-value = 0.4506

shapiro.fam.lex.0  <- shapiro.test(rmean.fam.lex.0) # W = 0.97654, p-value = 0.391
shapiro.fam.lex.1  <- shapiro.test(rmean.fam.lex.1) # W = 0.9563, p-value = 0.05415


# qqplots de factores aletorios.

library(lattice)
library(plyr)
library(ggplot2)

cat_1_u <- rmean.cat.1
reStack <- ldply(rmean.cat.1)

png("qqplot_cat_1.png")
print(qqmath( ~ cat_1_u, data=reStack, scales=list(relation="free"),
              prepanel = prepanel.qqmathline,
              panel = function(x, ...) {
                panel.qqmathline(x, ...)
                panel.qqmath(x, ...)
              },
              layout=c(1,1)))
dev.off()

cat_2_u <- rmean.cat.2
reStack <- ldply(rmean.cat.2)

png("qqplot_cat_2.png")
print(qqmath( ~ cat_2_u, data=reStack, scales=list(relation="free"),
              prepanel = prepanel.qqmathline,
              panel = function(x, ...) {
                panel.qqmathline(x, ...)
                panel.qqmath(x, ...)
              },
              layout=c(1,1)))
dev.off()

cat_3_u <- rmean.cat.3
reStack <- ldply(rmean.cat.3)

png("qqplot_cat_3.png")
print(qqmath( ~ cat_3_u, data=reStack, scales=list(relation="free"),
              prepanel = prepanel.qqmathline,
              panel = function(x, ...) {
                panel.qqmathline(x, ...)
                panel.qqmath(x, ...)
              },
              layout=c(1,1)))
dev.off()

cat_4_u <- rmean.cat.4
reStack <- ldply(rmean.cat.4)

png("qqplot_cat_4.png")
print(qqmath( ~ cat_4_u, data=reStack, scales=list(relation="free"),
              prepanel = prepanel.qqmathline,
              panel = function(x, ...) {
                panel.qqmathline(x, ...)
                panel.qqmath(x, ...)
              },
              layout=c(1,1)))
dev.off()

fam_lex_0_u <- rmean.fam.lex.0
reStack <- ldply( rmean.fam.lex.0)

png("qqplot_fam_lex_0.png")
print(qqmath( ~ fam_lex_0_u, data=reStack, scales=list(relation="free"),
              prepanel = prepanel.qqmathline,
              panel = function(x, ...) {
                panel.qqmathline(x, ...)
                panel.qqmath(x, ...)
              },
              layout=c(1,1)))
dev.off()

fam_lex_1_u <- rmean.fam.lex.1
reStack <- ldply( rmean.fam.lex.1)

png("qqplot_fam_lex_1.png")
print(qqmath( ~ fam_lex_1_u, data=reStack, scales=list(relation="free"),
              prepanel = prepanel.qqmathline,
              panel = function(x, ...) {
                panel.qqmathline(x, ...)
                panel.qqmath(x, ...)
              },
              layout=c(1,1)))
dev.off()

# trayectorias de efectos aleatorios

# cat 1

mean_frame <- data.frame(order = c(rep(1:4,12),rep(2:3,2)), rmean.cat.1)

t_v0_1    <- data.frame(v0i = c(mean_frame$rmean.cat.1[which(mean_frame$order==1)], 
mean_frame$rmean.cat.1[which(mean_frame$order==2)], 
mean_frame$rmean.cat.1[which(mean_frame$order==3)], 
mean_frame$rmean.cat.1[which(mean_frame$order==4)]), 
id = as.numeric(sapply(1:52, function(x){strsplit(as.character(unique(datos$ID.SESION)), 
                                         split="[.]")[[x]][1]})), 
sesion = as.numeric(sapply(1:52, function(x){strsplit(as.character(unique(datos$ID.SESION)), 
                                        split="[.]")[[x]][2]})) 
)            


t_v0_1[t_v0_1$id==1,"id"] <-  "SONIA"
t_v0_1[t_v0_1$id==2,"id"] <-  "NATI"
t_v0_1[t_v0_1$id==3,"id"] <-  "JAKO"
t_v0_1[t_v0_1$id==4,"id"] <-  "MIRKA"

png("trayectory_random_1.png")
ggplot(data=t_v0_1, aes(x=sesion,y=v0i)) +
  geom_point() + 
  xlab("Sesion") + 
  ylab("Random Effects") + 
  geom_hline(yintercept = 0.00, linetype="dotted", color = "blue", size=1.1) +
  facet_wrap( ~ id)
dev.off()

# cat 2

mean_frame <- data.frame(order = c(rep(1:4,12),rep(2:3,2)), rmean.cat.2)

t_v0_2    <- data.frame(v0i = c(mean_frame$rmean.cat.2[which(mean_frame$order==1)], 
                                mean_frame$rmean.cat.2[which(mean_frame$order==2)], 
                                mean_frame$rmean.cat.2[which(mean_frame$order==3)], 
                                mean_frame$rmean.cat.2[which(mean_frame$order==4)]), 
                        id = as.numeric(sapply(1:52, function(x){strsplit(as.character(unique(datos$ID.SESION)), 
                                                                          split="[.]")[[x]][1]})), 
                        sesion = as.numeric(sapply(1:52, function(x){strsplit(as.character(unique(datos$ID.SESION)), 
                                                                              split="[.]")[[x]][2]})) 
)            


t_v0_2[t_v0_2$id==1,"id"] <-  "SONIA"
t_v0_2[t_v0_2$id==2,"id"] <-  "NATI"
t_v0_2[t_v0_2$id==3,"id"] <-  "JAKO"
t_v0_2[t_v0_2$id==4,"id"] <-  "MIRKA"

png("trayectory_random_2.png")
ggplot(data=t_v0_2, aes(x=sesion,y=v0i)) +
  geom_point() + 
  xlab("Sesion") + 
  ylab("Random Effects") + 
  geom_hline(yintercept = 0.00, linetype="dotted", color = "blue", size=1.1) +
  facet_wrap( ~ id)
dev.off()

# cat 3

mean_frame <- data.frame(order = c(rep(1:4,12),rep(2:3,2)), rmean.cat.3)

t_v0_3    <- data.frame(v0i = c(mean_frame$rmean.cat.3[which(mean_frame$order==1)], 
                                mean_frame$rmean.cat.3[which(mean_frame$order==2)], 
                                mean_frame$rmean.cat.3[which(mean_frame$order==3)], 
                                mean_frame$rmean.cat.3[which(mean_frame$order==4)]), 
                        id = as.numeric(sapply(1:52, function(x){strsplit(as.character(unique(datos$ID.SESION)), 
                                                                          split="[.]")[[x]][1]})), 
                        sesion = as.numeric(sapply(1:52, function(x){strsplit(as.character(unique(datos$ID.SESION)), 
                                                                              split="[.]")[[x]][2]})) 
)            


t_v0_3[t_v0_3$id==1,"id"] <-  "SONIA"
t_v0_3[t_v0_3$id==2,"id"] <-  "NATI"
t_v0_3[t_v0_3$id==3,"id"] <-  "JAKO"
t_v0_3[t_v0_3$id==4,"id"] <-  "MIRKA"

png("trayectory_random_3.png")
ggplot(data=t_v0_3, aes(x=sesion,y=v0i)) +
  geom_point() + 
  xlab("Sesion") + 
  ylab("Random Effects") + 
  geom_hline(yintercept = 0.00, linetype="dotted", color = "blue", size=1.1) +
  facet_wrap( ~ id)
dev.off()

# cat 4

mean_frame <- data.frame(order = c(rep(1:4,12),rep(2:3,2)), rmean.cat.4)

t_v0_4    <- data.frame(v0i = c(mean_frame$rmean.cat.4[which(mean_frame$order==1)], 
                                mean_frame$rmean.cat.4[which(mean_frame$order==2)], 
                                mean_frame$rmean.cat.4[which(mean_frame$order==3)], 
                                mean_frame$rmean.cat.4[which(mean_frame$order==4)]), 
id = as.numeric(sapply(1:52, function(x){strsplit(as.character(unique(datos$ID.SESION)), 
                                         split="[.]")[[x]][1]})), 
sesion = as.numeric(sapply(1:52, function(x){strsplit(as.character(unique(datos$ID.SESION)), 
                                         split="[.]")[[x]][2]})) 
)            


t_v0_4[t_v0_4$id==1,"id"] <-  "SONIA"
t_v0_4[t_v0_4$id==2,"id"] <-  "NATI"
t_v0_4[t_v0_4$id==3,"id"] <-  "JAKO"
t_v0_4[t_v0_4$id==4,"id"] <-  "MIRKA"

png("trayectory_random_4.png")
ggplot(data=t_v0_4, aes(x=sesion,y=v0i)) +
  geom_point() + 
  xlab("Sesion") + 
  ylab("Random Effects") + 
  geom_hline(yintercept = 0.00, linetype="dotted", color = "blue", size=1.1) +
  facet_wrap( ~ id)
dev.off()


# efectos fijos


# posterior.OR <-  apply(global.model.3$Sol[,1:88],2,exp)
# posterior.mean.OR <- colMeans(posterior.OR)
# eff.posterior.OR  <- sapply(1:88, function(x){effectiveSize(as.mcmc(posterior.OR[,x]))})                 
# IC.posterior.OR  <- sapply(1:88, function(x){HPDinterval(as.mcmc(posterior.OR[,x]))})                    

cor(apply(global.model.3$Sol[,1:88],2,mean), apply(global.model.3$Sol[,1:88],2,posterior.mode))
# 0.9888705

SOL.3 <- summary(global.model.3)$solutions

SOL <- data.frame(SOL.3, OR = exp(SOL.3[,1]), l.CI.OR = exp(SOL.3[,2]), u.CI.OR = exp(SOL.3[,3]))

colnames(SOL) <- c("mean", "l.CI","u.CI", "eff", "pMCMC", "OR", "l.CI.OR", "u.CI.OR")

apply(SOL,2, function(x){round(x,digits=3)})

SOL.sig <- SOL[SOL$pMCMC < 0.05,]

SOL.sig[,5] <- round(SOL.sig[,5], digits=3)
SOL.sig[,5] <-    ifelse(SOL.sig$pMCMC < 0.01, "p<0.001",
                  ifelse(SOL.sig$pMCMC<0.0001, "p<0.0001", SOL.sig$pMCMC)                                   
                  )
                  

xtbl.seleccion.0 <- xtable(SOL.sig, 
                           caption = "ajuste del modelo (efectos con p<0.05). 
                           mean = media de la posterior, l CI = Int. de credibilidad del 95 por ciento, 
                           extremo inferior, u CI = Int. de credibilidad del 95 por ciento, extremo superior;
                           eff = tama??o muestral efectivo; pMCMC = p valor;
                           OR = exp(mean), l CI OR = exp(l CI), u CI OR = exp(u CI)"
                           ,digits = 3,label = "xtable:example")

xtbl.seleccion.full <- xtable(SOL, 
                           caption = "ajuste del modelo (todos los predictores). 
                           mean = media de la posterior, l CI = Int. de credibilidad del 95 por ciento, 
                           extremo inferior, u CI = Int. de credibilidad del 95 por ciento, extremo superior;
                           eff = tama??o muestral efectivo; pMCMC = p valor;
                           OR = exp(mean), l CI OR = exp(l CI), u CI OR = exp(u CI)"
                           ,digits = 3,label = "xtable:example")


print(xtbl.seleccion.full, include.rownames=TRUE, tabular.environment="longtable", floating=FALSE)


## Forestplot

library("MCMCvis")

objeto.mcmc <- coda::as.mcmc.list(global.model.3$Sol)[,row.names(SOL.sig)]


MCMCplot(objeto.mcmc, 
         rank = TRUE,
         xlab = 'ESTIMATE', col = 'red')



### efectos fijos (ajuste)

#                             mean   l.CI   u.CI      eff pMCMC    OR l.CI.OR u.CI.OR
# traitRES_CAT.1            -2.281 -3.443 -1.081 2000.000 0.000 0.102   0.032   0.339
# traitRES_CAT.2            -3.071 -4.420 -1.890 1572.478 0.000 0.046   0.012   0.151
# traitRES_CAT.3            -0.724 -1.432 -0.085 2000.000 0.037 0.485   0.239   0.919
# traitRES_CAT.4            -2.772 -3.801 -1.786 1594.264 0.000 0.063   0.022   0.168
# traitRES_CAT.1:TIME       -0.012 -0.044  0.019 1332.941 0.458 0.988   0.957   1.020
# traitRES_CAT.2:TIME       -0.007 -0.043  0.029 1093.957 0.734 0.993   0.958   1.029
# traitRES_CAT.3:TIME        0.009 -0.004  0.024 1870.334 0.184 1.009   0.996   1.024
# traitRES_CAT.4:TIME        0.009 -0.012  0.029 1874.825 0.384 1.009   0.989   1.030
# traitRES_CAT.1:MOD1       -0.169 -1.743  1.276 2000.000 0.847 0.845   0.175   3.582
# traitRES_CAT.2:MOD1       -0.583 -2.443  1.412 1713.370 0.569 0.558   0.087   4.104
# traitRES_CAT.3:MOD1       -0.489 -1.645  0.532 1687.632 0.382 0.613   0.193   1.702
# traitRES_CAT.4:MOD1       -0.556 -1.936  1.012 1609.775 0.475 0.574   0.144   2.750
# traitRES_CAT.1:MOD2        0.347 -0.287  0.935 1812.128 0.293 1.415   0.751   2.546
# traitRES_CAT.2:MOD2        0.194 -0.755  1.124 1276.952 0.644 1.215   0.470   3.077
# traitRES_CAT.3:MOD2        0.328 -0.032  0.706 2038.619 0.089 1.388   0.968   2.025
# traitRES_CAT.4:MOD2        0.883  0.315  1.500 1596.517 0.001 2.419   1.370   4.480
# traitRES_CAT.1:MOD3       -0.898 -1.814 -0.163 1173.143 0.032 0.407   0.163   0.850
# traitRES_CAT.2:MOD3       -0.259 -1.443  0.758 1208.766 0.683 0.772   0.236   2.135
# traitRES_CAT.3:MOD3        0.431  0.007  0.840 1738.849 0.045 1.538   1.007   2.316
# traitRES_CAT.4:MOD3        0.472 -0.157  1.163 1655.423 0.149 1.603   0.855   3.201
# traitRES_CAT.1:Fabs.SC.f1 -1.277 -2.036 -0.491 1596.650 0.001 0.279   0.130   0.612
# traitRES_CAT.2:Fabs.SC.f1  0.071 -0.799  1.099 1278.310 0.900 1.073   0.450   3.002
# traitRES_CAT.3:Fabs.SC.f1 -0.424 -0.799 -0.012 1851.367 0.035 0.654   0.450   0.988
# traitRES_CAT.4:Fabs.SC.f1 -0.927 -1.485 -0.250 1530.948 0.003 0.396   0.226   0.778
# traitRES_CAT.1:ANIM1      -0.121 -0.772  0.543 1447.283 0.738 0.886   0.462   1.721
# traitRES_CAT.2:ANIM1       1.254  0.521  2.008 1190.987 0.003 3.503   1.684   7.445
# traitRES_CAT.3:ANIM1       0.132 -0.192  0.467 1809.376 0.436 1.141   0.826   1.596
# traitRES_CAT.4:ANIM1       0.519 -0.026  1.031 1673.390 0.063 1.680   0.974   2.803
# traitRES_CAT.1:ES1         0.541 -0.647  1.677 1797.495 0.382 1.718   0.524   5.351
# traitRES_CAT.2:ES1         0.066 -1.274  1.383 1680.962 0.933 1.068   0.280   3.989
# traitRES_CAT.3:ES1        -0.308 -1.049  0.379 2000.000 0.426 0.735   0.350   1.461
# traitRES_CAT.4:ES1         0.482 -0.519  1.454 1799.522 0.355 1.620   0.595   4.281
# traitRES_CAT.1:ES2        -0.672 -2.717  1.237 1867.843 0.525 0.511   0.066   3.445
# traitRES_CAT.2:ES2         0.111 -1.610  1.923 1505.873 0.872 1.117   0.200   6.842
# traitRES_CAT.3:ES2        -2.007 -3.532 -0.624 1601.730 0.004 0.134   0.029   0.536
# traitRES_CAT.4:ES2         0.046 -1.304  1.484 1548.543 0.952 1.047   0.272   4.409
# traitRES_CAT.1:MORF.f1    -0.191 -1.150  0.779 2000.000 0.667 0.826   0.317   2.178
# traitRES_CAT.2:MORF.f1    -0.962 -2.058  0.274 1537.687 0.100 0.382   0.128   1.315
# traitRES_CAT.3:MORF.f1    -1.233 -1.714 -0.746 2000.000 0.000 0.291   0.180   0.474
# traitRES_CAT.4:MORF.f1    -0.630 -1.342  0.097 1733.746 0.091 0.533   0.261   1.102
# traitRES_CAT.1:MORF.f2    -0.342 -1.359  0.757 1867.880 0.540 0.710   0.257   2.132
# traitRES_CAT.2:MORF.f2    -0.086 -1.293  1.039 1725.062 0.870 0.917   0.274   2.826
# traitRES_CAT.3:MORF.f2    -0.949 -1.557 -0.287 1647.699 0.005 0.387   0.211   0.751
# traitRES_CAT.4:MORF.f2    -0.497 -1.299  0.374 1724.647 0.264 0.608   0.273   1.454
# traitRES_CAT.1:FAM.LEX.f1 -0.464 -1.105  0.173 1544.699 0.163 0.629   0.331   1.189
# traitRES_CAT.2:FAM.LEX.f1 -1.362 -2.201 -0.475 1098.049 0.001 0.256   0.111   0.622
# traitRES_CAT.3:FAM.LEX.f1 -0.563 -0.953 -0.141 2000.000 0.004 0.570   0.386   0.868
# traitRES_CAT.4:FAM.LEX.f1 -0.897 -1.432 -0.360 1533.184 0.001 0.408   0.239   0.698
# traitRES_CAT.1:CUMRES.f1   0.137 -0.573  0.862 1863.065 0.736 1.147   0.564   2.369
# traitRES_CAT.2:CUMRES.f1   0.404 -0.443  1.269 1052.338 0.351 1.498   0.642   3.558
# traitRES_CAT.3:CUMRES.f1  -0.038 -0.449  0.373 2000.000 0.855 0.962   0.638   1.452
# traitRES_CAT.4:CUMRES.f1   0.704  0.025  1.322 1653.229 0.026 2.021   1.026   3.749
# traitRES_CAT.1:CUMRES.f2   0.046 -1.058  1.130 2000.000 0.951 1.047   0.347   3.097
# traitRES_CAT.2:CUMRES.f2  -0.510 -1.849  0.874 1403.991 0.468 0.601   0.157   2.396
# traitRES_CAT.3:CUMRES.f2  -0.222 -0.834  0.433 2000.000 0.467 0.801   0.434   1.542
# traitRES_CAT.4:CUMRES.f2   0.503 -0.379  1.401 1616.208 0.296 1.653   0.685   4.060
# traitRES_CAT.1:EST11      -1.773 -2.720 -0.924  914.828 0.000 0.170   0.066   0.397
# traitRES_CAT.2:EST11      -0.444 -1.288  0.497 1306.433 0.347 0.641   0.276   1.644
# traitRES_CAT.3:EST11      -0.077 -0.476  0.314 2000.000 0.691 0.926   0.621   1.368
# traitRES_CAT.4:EST11      -0.893 -1.539 -0.229 1130.814 0.008 0.410   0.215   0.796
# traitRES_CAT.1:EST21       0.118 -0.641  0.818 1760.448 0.755 1.125   0.527   2.267
# traitRES_CAT.2:EST21      -1.510 -2.773 -0.198  672.114 0.021 0.221   0.062   0.820
# traitRES_CAT.3:EST21       0.272 -0.124  0.672 2000.000 0.205 1.313   0.884   1.959
# traitRES_CAT.4:EST21       0.102 -0.481  0.722 1711.588 0.740 1.107   0.618   2.058
# traitRES_CAT.1:EST51      -0.865 -2.059  0.175 1556.623 0.120 0.421   0.128   1.191
# traitRES_CAT.2:EST51      -2.139 -3.540 -0.729  856.514 0.000 0.118   0.029   0.483
# traitRES_CAT.3:EST51      -0.107 -0.715  0.427 2000.000 0.730 0.899   0.489   1.533
# traitRES_CAT.4:EST51      -1.855 -2.756 -1.088 1372.963 0.000 0.156   0.064   0.337
# traitRES_CAT.1:GRUPO62    -0.603 -1.757  0.401 1510.496 0.275 0.547   0.173   1.493
# traitRES_CAT.2:GRUPO62    -0.520 -1.809  0.640 1549.728 0.396 0.594   0.164   1.897
# traitRES_CAT.3:GRUPO62     0.170 -0.412  0.781 2000.000 0.565 1.186   0.663   2.183
# traitRES_CAT.4:GRUPO62     0.664 -0.083  1.472 1576.769 0.100 1.943   0.920   4.357
# traitRES_CAT.1:GRUPO63    -0.139 -0.987  0.661 1723.035 0.756 0.871   0.373   1.937
# traitRES_CAT.2:GRUPO63    -0.504 -1.695  0.654 1029.938 0.405 0.604   0.184   1.923
# traitRES_CAT.3:GRUPO63     0.074 -0.342  0.488 1647.641 0.747 1.077   0.711   1.629
# traitRES_CAT.4:GRUPO63    -0.488 -1.150  0.204 1587.010 0.163 0.614   0.317   1.226
# traitRES_CAT.1:GRUPO64    -0.218 -1.119  0.687 1710.223 0.661 0.804   0.327   1.987
# traitRES_CAT.2:GRUPO64    -0.278 -1.364  0.945 1299.757 0.659 0.757   0.256   2.573
# traitRES_CAT.3:GRUPO64    -0.008 -0.448  0.499 1813.817 0.970 0.992   0.639   1.647
# traitRES_CAT.4:GRUPO64    -0.067 -0.854  0.638 1650.738 0.864 0.935   0.426   1.893
# traitRES_CAT.1:GRUPO65    -0.586 -2.041  0.842 1737.677 0.435 0.557   0.130   2.321
# traitRES_CAT.2:GRUPO65    -0.591 -2.236  0.893 1723.644 0.466 0.554   0.107   2.443
# traitRES_CAT.3:GRUPO65     0.519 -0.169  1.146 2000.000 0.117 1.681   0.845   3.145
# traitRES_CAT.4:GRUPO65    -1.207 -2.366  0.050 1367.219 0.036 0.299   0.094   1.051
# traitRES_CAT.1:GRUPO66    -0.242 -1.545  1.177 1369.080 0.763 0.785   0.213   3.246
# traitRES_CAT.2:GRUPO66     0.752 -0.351  2.051 1864.481 0.237 2.121   0.704   7.777
# traitRES_CAT.3:GRUPO66    -0.251 -1.139  0.675 2000.000 0.605 0.778   0.320   1.964
# traitRES_CAT.4:GRUPO66     0.138 -0.959  1.388 1836.897 0.797 1.148   0.383   4.005


##### TYPES más frecuentes

# lista de casos
datos[datos$Fabs.SC.f==1 & datos$RES_CAT!=0 & datos$RES_CAT!=2, c("RES_CAT","TYPE", "INSTANCIA", "ID.SESION", "LINEA")])

# tabla de cantidad de TYPES por tipo de error
as.matrix(table(datos[datos$Fabs.SC.f==1 & datos$RES_CAT!=0 & datos$RES_CAT!=2, c("TYPE", "RES_CAT")]))[,-c(1,3)]

# ejemplos por TYPE

datos[datos$Fabs.SC.f==1 & datos$RES_CAT!=0 & datos$RES_CAT!=2 & datos$TYPE=="d-n-as-*es", c("RES_CAT","TYPE", "INSTANCIA")]

datos[datos$Fabs.SC.f==1 & datos$RES_CAT!=0 & datos$RES_CAT!=2 & datos$TYPE=="d-n-as-*es", c("RES_CAT","TYPE", "INSTANCIA","ID.SESION")]

datos[datos$Fabs.SC.f==1 & datos$RES_CAT!=0 & datos$RES_CAT!=2 & datos$TYPE=="d-n-as-*es", c("RES_CAT","TYPE", "INSTANCIA","ID.SESION", "LINEA")]

datos[datos$Fabs.SC.f==1 & datos$RES_CAT!=0 & datos$RES_CAT!=2 & datos$TYPE=="d-n-as-as", c("RES_CAT","TYPE", "INSTANCIA","ID.SESION", "LINEA")]

datos[datos$Fabs.SC.f==1 & datos$RES_CAT!=0 & datos$RES_CAT!=2 & datos$TYPE=="d-n-os-os", c("RES_CAT","TYPE", "INSTANCIA","ID.SESION", "LINEA")]

datos[datos$Fabs.SC.f==1 & datos$RES_CAT!=0 & datos$RES_CAT!=2 & datos$TYPE=="l-n-as-*es", c("RES_CAT","TYPE", "INSTANCIA","ID.SESION", "LINEA")]

datos[datos$Fabs.SC.f==1 & datos$RES_CAT!=0 & datos$RES_CAT!=2 & datos$TYPE=="l-n-as-as", c("RES_CAT","TYPE", "INSTANCIA","ID.SESION", "LINEA")]

datos[datos$Fabs.SC.f==1 & datos$RES_CAT!=0 & datos$RES_CAT!=2 & datos$TYPE=="l-n-j-os-os-os", c("RES_CAT","TYPE", "INSTANCIA","ID.SESION", "LINEA")]

datos[datos$Fabs.SC.f==1 & datos$RES_CAT!=0 & datos$RES_CAT!=2 & datos$TYPE=="l-n-os-*es", c("RES_CAT","TYPE", "INSTANCIA","ID.SESION", "LINEA")]

datos[datos$Fabs.SC.f==1 & datos$RES_CAT!=0 & datos$RES_CAT!=2 & datos$TYPE=="l-n-os-es", c("RES_CAT","TYPE", "INSTANCIA","ID.SESION", "LINEA")]

datos[datos$Fabs.SC.f==1 & datos$RES_CAT!=0 & datos$RES_CAT!=2 & datos$TYPE=="l-n-os-os", c("RES_CAT","TYPE", "INSTANCIA","ID.SESION", "LINEA")]

