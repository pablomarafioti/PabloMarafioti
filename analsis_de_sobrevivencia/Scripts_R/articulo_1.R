library(survival)
library(survminer)
library(MuMIn)
library("ggplot2") 
library("xtable")
library(broom)
library(broom.mixed)
library(tidyverse)
library(forestplot)
library("mstate")

load("datos_suv.RData") ## datos para el modelo de eventos multiples

############# modelo de eventos multiples


## exploratorio

fit <- survfit(Surv(tstop-tstart ,status) ~ ID, data = datos.suv)
ggsurvplot(fit, pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, risk.table.y.text.col = TRUE)


res <- pairwise_survdiff(Surv(tstop-tstart, status) ~ ID,
                         data = datos.suv)

# Pairwise comparisons using Log-Rank test 
# 
# data:  datos.suv and ID 
# 
#         1       2       3      
# 2 0.04985 -       -      
# 3 0.00045 4.5e-10 -      
# 4 0.00045 1.5e-10 0.27017
# 
# P value adjustment method: BH 

datos.suv$time <- datos.suv$tstop - datos.suv$tstart
ggsurvevents(Surv(datos.suv$time, datos.suv$status), type = "cumulative")

## modelos

# AG model

model.AG=coxph(Surv(tstart, tstop, status) ~ Fabs.SC.f + MORF.f + STEM.f + MOD + ES + ANIM + 
                 GRAMS + FAM.LEX.f + IMA.CONC.f + LDA + EST1 + EST2 + EST3 + EST4 + 
                 EST5 + EST6 + EST7 + cluster(id), method="breslow", data=datos.suv)


model.AG.strata=coxph(Surv(tstart, tstop, status) ~ Fabs.SC.f + MORF.f + STEM.f + MOD + ES + ANIM + 
                        GRAMS + FAM.LEX.f + IMA.CONC.f + LDA + EST1 + EST2 + EST3 + EST4 + 
                        EST5 + EST6 + EST7 + cluster(id) + strata(ID), method="breslow", data=datos.suv)


# Condictional (PWP)

model.PWP=coxph(Surv(tstart, tstop, status) ~ Fabs.SC.f + MORF.f + STEM.f + MOD + ES + ANIM + 
                  GRAMS + FAM.LEX.f + IMA.CONC.f + LDA + EST1 + EST2 + EST3 + EST4 + 
                  EST5 + EST6 + EST7 + cluster(id) + strata(event.2), method="breslow", data=datos.suv)


# Fraility 

model.frailty=coxph(Surv(tstart, tstop, status) ~ Fabs.SC.f + MORF.f + STEM.f + MOD + ES + ANIM + 
                      GRAMS + FAM.LEX.f + IMA.CONC.f + LDA + EST1 + EST2 + EST3 + EST4 + 
                      EST5 + EST6 + EST7 + frailty(id), data=datos.suv)



model.frailty.strata=coxph(Surv(tstart, tstop, status) ~ Fabs.SC.f + MORF.f + STEM.f + MOD + ES + ANIM + 
                             GRAMS + FAM.LEX.f + IMA.CONC.f + LDA + EST1 + EST2 + EST3 + EST4 + 
                             EST5 + EST6 + EST7 + frailty(id) + strata(ID), data=datos.suv)

# comparacion de modelos

seleccion <- data.frame(model.sel(model.AG, model.AG.strata,  model.PWP,model.frailty,
                                  model.frailty.strata,rank="AIC"))
model_sel  <- seleccion[,c("logLik","AIC","delta","weight")]


xtable(model_sel, caption = "Seleccion de modelos. 
logLik: verosimilitud del modelo, AIC: AKAIKE, delta: delta de AKAIKE, weight: pesos de 
       AKAIKE", digits = 3,label = "xtable:example")


#                         logLik      AIC       delta        weight
# model.frailty.strata -2727.688 5516.660    0.000000  6.666949e-01
# model.AG.strata      -2738.023 5518.046    1.386548  3.333051e-01
# model.PWP            -2993.258 6028.516  511.855783 4.740849e-112
# model.frailty        -3245.910 6556.913 1040.253417 8.625170e-227
# model.AG             -3258.391 6558.782 1042.122507 3.387661e-227

(summary(model.PWP)$coefficients)[which(summary(model.PWP)$coefficients[,6]<0.05),]
(summary(model.WLW)$coefficients)[which(summary(model.WLW)$coefficients[,6]<0.05),]
(summary(model.frailty)$coefficients)[which(summary(model.frailty)$coefficients[,6]<0.05),]
(summary(model.frailty.strata)$coefficients)[which(summary(model.frailty.strata)$coefficients[,6]<0.05),]
(summary(model.AG.strata)$coefficients)[which(summary(model.AG.strata)$coefficients[,6]<0.05),]

# seleccion de modelos (model.AG.strata)

# lo siguiente hace la seleccion de modelos, cuyo resultado es el objeto 
# sel.dredge.WLW, que se ha grabado y se lo recupera con "load". 

#options(na.action = "na.fail")
#sel.dredge.WLW <- dredge(model.AG.strata,fixed=c("strata(ID)") ,rank="AIC")

load(sel.dredge.WLW)

sub.m   <- subset(sel.dredge.WLW[-131072,], 1/8 < weight/max(sel.dredge.WLW[-131072,]$weight))
importance.m <- importance(sel.dredge.WLW[-131072,])
average.m <- summary(model.avg(subset(sel.dredge.WLW[-131072,], 1/8 < weight/max(sel.dredge.WLW[-131072,]$weight))))

xtbl.importance <- xtable(tidy(importance.m)[-1,], 
                          caption = "Importancia Relativa de las predictoras.", digits = 2,label = "xtable:example")

xtbl.full.average <- xtable(average.m$coefmat.full, 
                          caption = "Promedio de los coeficientes con FULL AVERAGE.", digits = 3,label = "xtable:example")

#Importancia Relativa de las predictoras.
#
#                        strata(ID) ANIM   EST1   EST5   MORF.f ES    
# Sum of weights:        1.00       1.00   0.99   0.98   0.89   0.74
# N containing models: 131071      65536  65536  65536  65536  65536
#                        IMA.CONC.f FAM.LEX.f EST4   MOD    LDA    STEM.f EST3  
# Sum of weights:        0.59       0.57      0.50   0.45   0.44   0.39   0.39
# N containing models:  65536      65536     65536  65536  65536  65536  65536
#                       Fabs.SC.f EST7   EST2   EST6   GRAMS 
# Sum of weights:        0.37      0.34   0.32   0.30   0.29
# N containing models:  65536     65536  65536  65536  65536


# Model-averaged coefficients:  
#   (full average) 
#              Estimate Std. Error z value Pr(>|z|)    
# ANIM1        0.4119295  0.1145429   3.596 0.000323 ***
# ES1          0.2436273  0.2495078   0.976 0.328851    
# ES2         -0.6037448  0.5679692   1.063 0.287787    
# EST11       -0.4945333  0.1177934   4.198 2.69e-05 ***
# EST41       -0.3321252  0.3548875   0.936 0.349345    
# EST51       -0.5976854  0.1854098   3.224 0.001266 ** 
# FAM.LEX.f1  -0.1067561  0.1079167   0.989 0.322543    
# IMA.CONC.f1 -0.1143326  0.1168591   0.978 0.327886    
# MORF.f1     -0.5185097  0.1735834   2.987 0.002816 ** 
# MORF.f2     -0.4905616  0.2786522   1.760 0.078326 .  
# STEM.f1      0.0367694  0.0862336   0.426 0.669821    
# LDA1        -0.0871724  0.1741489   0.501 0.616679    
# EST71        0.0334122  0.1213880   0.275 0.783123    
# MOD1        -0.1941258  0.3724266   0.521 0.602196    
# MOD2         0.0456849  0.0956244   0.478 0.632826    
# MOD3         0.0002968  0.0740295   0.004 0.996801    
# EST31        0.1088265  0.3019582   0.360 0.718546    
# EST21        0.0091748  0.0517235   0.177 0.859208    
# Fabs.SC.f1  -0.0153529  0.0581048   0.264 0.791605    
# GRAMS1       0.0031205  0.0322558   0.097 0.922930    
# EST61        0.0065680  0.0629282   0.104 0.916874   


# modelo reducido

model.AG.strata.2=coxph(Surv(tstart, tstop, status) ~ MORF.f + ANIM + EST1 + EST5  + 
                          cluster(id) + strata(ID), data=datos.suv)

# analisis de supuestos

ggcoxdiagnostics(model.AG.strata.2, type = "dfbetas",
                 linear.predictions = FALSE, ggtheme = theme_bw())

ggcoxdiagnostics(model.AG.strata.2, type = "deviance",
                 linear.predictions = FALSE, ggtheme = theme_bw())

res.deviance <- residuals(model.AG.strata.2, type = "deviance")

sum(abs(res.deviance)>2)/length(res.deviance) # 0.05681191, 103


shon.1 <- cox.zph(model.AG.strata.2, transform="rank")
shon.2 <- cox.zph(model.AG.strata.2, transform="km")
shon.3 <- cox.zph(model.AG.strata.2, transform=function(t){t})
shon.4 <- cox.zph(model.AG.strata.2, transform=function(t){log(t)})

shoenfeld <- data.frame(df= c(rep(1,5),5),rank.chi2 = shon.1$table[,2], rank.p = shon.1$table[,3],
                        km.chi2 = shon.2$table[,2], km.p = shon.2$table[,3],
                        time.chi2 = shon.3$table[,2], time.p = shon.3$table[,3],
                        log.time.chi2 = shon.4$table[,2], log.time.p = shon.4$table[,3])

xtbl.test.shoenfeld <- xtable(shoenfeld, 
                              caption = "Test para theta = 0, según diferentes funciones del tiempo: estadístico de chi cuadrado y p valor", 
                              digits = 4,label = "xtable:example")

ggcoxzph(shon.3, var = "ANIM1")

## pruebas para supuesto de proporcionalidad
#
#         df  rank.chi2       rank.p      km.chi2         km.p    time.chi2
# MORF.f1  1  0.1429597 0.7053562266  0.008689879 0.9257291671 2.237024e-04
# MORF.f2  1  0.7135987 0.3982517118  0.508751648 0.4756797579 5.284034e-01
# ANIM1    1 11.4895321 0.0006998925 11.002649614 0.0009098173 1.105479e+01
# EST11    1  0.4290532 0.5124538834  0.264519138 0.6070323073 3.078444e-01
# EST51    1  1.5289449 0.2162708648  1.344321541 0.2462732912 1.374359e+00
# GLOBAL   5 11.8178969 0.0373700319 11.269061043 0.0462986951 1.128094e+01
#             time.p  log.time.chi2  log.time.p
# MORF.f1 0.988066737      0.201104 0.653831197
# MORF.f2 0.467278866      1.155523 0.282395979
# ANIM1   0.000884583      7.243048 0.007117607
# EST11   0.579005966      1.180771 0.277199124
# EST51   0.241064391      1.435552 0.230860403
# GLOBAL  0.046085617      7.753002 0.170383802


### modelo con coeficiente de ANIM dependiendo del tiempo linealmente: 
# beta_anim = beta_anim_0 + beta_anim_1*(x*t)  


datos.suv$ANIM1 <- as.numeric(as.character(datos.suv$ANIM))
model.AG.time=coxph(Surv(tstart, tstop ,status) ~ MORF.f + ANIM1 + tt(ANIM1) + 
                      EST1 + EST5  + cluster(id)  + strata(ID), 
                    tt=function(x,t, ...){x*t},   data=datos.suv)

resid <- cox.zph(model.AG.strata.2, transform=function(t){t})
plot(resid[3])
abline(coef(model.AG.time)[3:4], col=2)

### modelo final (ajuste)

sumario <- data.frame(summary(model.AG.time)$coefficients, summary(model.AG.time)$conf.int[,3:4])

xtable(sumario, 
       caption = "Modelo de eventos múltiples. 
       coef: betas estimados, exp(coef): hazard ratios, se(coef): error típico de betas estimados, 
       robust se: error típico de beta con estimador sandwhich, 
       z: coef / robust se, Pr(z): p-valor,
       lower 95: extremo izquierdo de intervalo de confianza de 95 por ciento para hazard ratio,
       upper 95: extremo derecho de intervalo de confianza de 95 por ciento para hazard ratio
       ", digits = 3,label = "xtable:example")

#ajuste del modelo final
#
#                   coef exp.coef.     se.coef.   robust.se         z
# MORF.f1   -0.397677539 0.6718786 0.1350412858 0.102667645 -3.873446
# MORF.f2   -0.256488574 0.7737638 0.1633310694 0.144973790 -1.769207
# ANIM1      0.657709690 1.9303661 0.1764824565 0.195910373  3.357197
# tt(ANIM1) -0.001456267 0.9985448 0.0006192538 0.000603753 -2.412024
# EST11     -0.400986420 0.6696592 0.1125910764 0.098361225 -4.076672
# EST51     -0.396407044 0.6727328 0.1689806255 0.155444900 -2.550145
#               Pr...z.. lower..95 upper..95
# MORF.f1   1.073074e-04 0.5494143 0.8216402
# MORF.f2   7.685942e-02 0.5823801 1.0280408
# ANIM1     7.873703e-04 1.3148629 2.8339938
# tt(ANIM1) 1.586425e-02 0.9973639 0.9997271
# EST11     4.568493e-05 0.5522409 0.8120430
# EST51     1.076781e-02 0.4960525 0.9123417



############# modelo de riesgos competitivos

load("suv_cr.RData") # datos 


# exploratorio

fit1 <- survfit(Surv(tiempos, status, type="mstate") ~ 1, data=suv.cr)
fit2 <- survfit(Surv(tiempos, status, type="mstate") ~ 1 + strata(ID), data=suv.cr)
ggcompetingrisks(fit1)
ggcompetingrisks(fit2)


# datos para modelo de riesgos competitivos en formato "long"

tmat <- trans.comprisk(4, names = c("0", "1", "2","3","4"))

suv.cr$stat1 <- as.numeric(suv.cr$status == 1)
suv.cr$stat2 <- as.numeric(suv.cr$status == 2)
suv.cr$stat3 <- as.numeric(suv.cr$status == 3)
suv.cr$stat4 <- as.numeric(suv.cr$status == 4)

cr.long <- msprep(time = c(NA, "tiempos", "tiempos", "tiempos", "tiempos"), status = c(NA,
                                                                                       "stat1", "stat2", "stat3", "stat4"), data = suv.cr, keep = colnames(suv.cr)[1:23], trans = tmat)

formula <- paste("Surv(time, status)", "~", "Fabs.SC.f + MORF.f + STEM.f + MOD + ES + ANIM + 
  GRAMS + FAM.LEX.f + IMA.CONC.f + LDA + EST1 + EST2 + EST3 + EST4 + 
  EST5 + EST6 + EST7 + strata(trans)")

model.WLW.cr=coxph(as.formula(paste(formula, "+cluster(ID.SESION)")), method="breslow", data=cr.long)

model.frailty.cr=coxph(as.formula(paste(formula, "+frailty(ID.SESION)")), method="breslow", data=cr.long)

AIC(model.WLW.cr,  model.frailty.cr)

#                        df      AIC
# model.WLW.cr     21.00000 6703.988
# model.frailty.cr 55.39555 6628.656

(summary(model.frailty.cr)$coefficients)[which(summary(model.frailty.cr)$coefficients[,6]<0.05),]

# seleccion de variables

# lo siguiente hace la seleccion de modelos, cuyo resultado es el objeto 
# sel.dredge.frailty.c, que se ha grabado y se lo recupera con "load". 

#options(na.action = "na.fail")
#sel.dredge.frailty.cr <- dredge(model.frailty.cr,fixed=c("frailty(ID.SESION)", "strata(trans)") ,rank="AIC")

load("sel_dredge_frailty.RData")

sub.cr   <- subset(sel.dredge.frailty.cr[-which(sel.dredge.frailty.cr$df==34),], 
                   1/8 < weight/max(sel.dredge.frailty.cr[-which(sel.dredge.frailty.cr$df==34),]$weight))
importance.cr <- importance(sel.dredge.frailty.cr[-which(sel.dredge.frailty.cr$df==34),])
average.cr <- summary(model.avg(subset(sel.dredge.frailty.cr[-which(sel.dredge.frailty.cr$df==34),],
                                       1/8 < weight/max(sel.dredge.frailty.cr[-which(sel.dredge.frailty.cr$df==34),]$weight))))

xtbl.importance <- xtable(tidy(importance.cr), 
                          caption = "Importancia Relativa de las predictoras.", digits = 2,label = "xtable:example")

xtbl.full.average <- xtable(average.cr$coefmat.full, 
                            caption = "Promedio de los coeficientes con FULL AVERAGE.", digits = 3,label = "xtable:example")

#importancia relativa de las predictoras
#
#                        frailty(ID.SESION) strata(trans) MOD    EST5   MORF.f
# Sum of weights:        1.00               1.00          1.00   1.00   0.99
# N containing models: 131071             131071         65536  65536  65536
#                        FAM.LEX.f EST1   Fabs.SC.f ANIM   ES     EST4   STEM.f
# Sum of weights:        0.99      0.97   0.93      0.93   0.83   0.70   0.56
# N containing models:  65536     65536  65536     65536  65536  65536  65536
#                        EST3   EST7   EST6   IMA.CONC.f EST2   GRAMS  LDA   
# Sum of weights:        0.52   0.44   0.43   0.36       0.36   0.35   0.28
# N containing models:  65536  65536  65536  65536      65536  65536  65536


# Model-averaged coefficients:  
#   (full average) 
#                Estimate Std. Error z value Pr(>|z|)    
# ANIM1        0.3522849  0.1170828   3.009 0.002622 ** 
# ES1          0.3164449  0.2621420   1.207 0.227374    
# ES2         -0.6942047  0.5707971   1.216 0.223908    
# EST11       -0.4328352  0.1269894   3.408 0.000653 ***
# EST41       -0.5336625  0.3689767   1.446 0.148084    
# EST51       -0.7856536  0.2003748   3.921 8.82e-05 ***
# FAM.LEX.f1  -0.3472709  0.1002007   3.466 0.000529 ***
# Fabs.SC.f1  -0.3320068  0.1293046   2.568 0.010240 *  
# MOD1        -0.4192980  0.4293435   0.977 0.328766    
# MOD2         0.5182303  0.1255434   4.128 3.66e-05 ***
# MOD3         0.3696102  0.1427668   2.589 0.009628 ** 
# MORF.f1     -0.7189908  0.1828845   3.931 8.45e-05 ***
# MORF.f2     -0.6672735  0.2933506   2.275 0.022926 *  
# STEM.f1      0.1080345  0.1290175   0.837 0.402388    
# EST31        0.3115539  0.4754402   0.655 0.512277    
# EST71        0.1176634  0.2139193   0.550 0.582294    
# EST61        0.0707710  0.1590072   0.445 0.656262    
# GRAMS1      -0.0267939  0.0763420   0.351 0.725609    
# IMA.CONC.f1 -0.0215752  0.0690600   0.312 0.754727    
# EST21        0.0170118  0.0741187   0.230 0.818464    
# LDA1         0.0004987  0.0844858   0.006 0.995290 


# modelo reducido


which(colnames(cr.long)=="MOD1.2"|colnames(cr.long)=="EST51.2") # 49, 105

betas.cr  <- colnames(cr.long)[c(32:43,48:59,68:71,76:79, 88:91,104:107)]
which(betas.cr=="MOD1.2"|betas.cr=="EST51.2") # 14, 38

formula.strata <- paste("Surv(time, status)", "~",
                        paste( betas.cr[-c(14, 38)], collapse = " + "), 
                        "+ frailty(ID.SESION)")

# modelo 1:  una hazard baseline para cada estrato.

model.frailty.cr.s=coxph(as.formula(paste(formula.strata, "+ strata(trans)")), data=cr.long)

# modelo 2:  una hazard baseline comun para todos los estratos.

model.frailty.cr.0=coxph(Surv(time, status) ~ Fabs.SC.f + MORF.f + MOD + ANIM + FAM.LEX.f + 
                           EST1 + EST5 + strata(trans) + frailty(ID.SESION), data=cr.long)

# comparacion: elige modelo 1

AIC(model.frailty.cr.0, model.frailty.cr.s)

#                         df      AIC
# model.frailty.cr.0 43.64588 6626.958
# model.frailty.cr.s 71.43354 6595.025

sumario.cr.s <- data.frame(summary(model.frailty.cr.s)$coefficients[-39,-3], 
                           summary(model.frailty.cr.s)$conf.int[,-2])

xtable(sumario.cr.s, 
       caption = "Modelo de riesgos competitivos. 
       coef: betas estimados, exp(coef): hazard ratios, se(coef): error típico de betas estimados, 
       robust se: error típico de beta con estimador sandwhich, 
       z: coef / robust se, Pr(z): p-valor,
       lower 95: extremo izquierdo de intervalo de confianza de 95 por ciento para hazard ratio,
       upper 95: extremo derecho de intervalo de confianza de 95 por ciento para hazard ratio
       ", digits = 3,label = "xtable:example")

#                     coef  se.coef.       Chisq DF            p exp.coef.
# Fabs.SC.f1.1 -0.97157826 0.3034958 10.24825887  1 1.368137e-03 0.3784852
# Fabs.SC.f1.2  1.05646109 0.5239136  4.06619038  1 4.374984e-02 2.8761744
# Fabs.SC.f1.3 -0.31544983 0.1516801  4.32517444  1 3.755256e-02 0.7294607
# Fabs.SC.f1.4 -0.25819313 0.2685459  0.92438387  1 3.363265e-01 0.7724460
# MORF.f1.1     0.27779833 0.4802269  0.33463115  1 5.629447e-01 1.3202199
# MORF.f1.2    -1.37913589 0.6122134  5.07467751  1 2.427779e-02 0.2517960
# MORF.f1.3    -0.83654440 0.1725287 23.51013564  1 1.242574e-06 0.4332049
# MORF.f1.4    -0.43954595 0.3554626  1.52904625  1 2.162556e-01 0.6443289
# MORF.f2.1     0.24490554 0.5659510  0.18725766  1 6.652089e-01 1.2775006
# MORF.f2.2    -0.50205407 0.6343142  0.62645820  1 4.286575e-01 0.6052861
# MORF.f2.3    -0.90728579 0.2184656 17.24734564  1 3.281557e-05 0.4036182
# MORF.f2.4     0.36951835 0.3887627  0.90344725  1 3.418591e-01 1.4470375
# MOD1.1       -0.50475239 1.0362067  0.23728154  1 6.261765e-01 0.6036550
# MOD1.3       -0.32208971 0.5225334  0.37994928  1 5.376304e-01 0.7246332
# MOD1.4       -0.29409515 1.0388722  0.08014040  1 7.771072e-01 0.7452056
# MOD2.1        0.55067970 0.2937844  3.51350525  1 6.087054e-02 1.7344315
# MOD2.2        0.37252399 0.5070285  0.53981359  1 4.625100e-01 1.4513933
# MOD2.3        0.40909291 0.1596225  6.56833707  1 1.038089e-02 1.5054516
# MOD2.4        1.06181100 0.2807450 14.30442102  1 1.554992e-04 2.8916030
# MOD3.1       -0.69056688 0.4040632  2.92087495  1 8.744044e-02 0.5012918
# MOD3.2        0.48123942 0.6469819  0.55327128  1 4.569841e-01 1.6180786
# MOD3.3        0.48674255 0.1709423  8.10773571  1 4.407675e-03 1.6270077
# MOD3.4        0.73024663 0.3244982  5.06424440  1 2.442435e-02 2.0755924
# ANIM1.1       0.47026794 0.2766798  2.88892250  1 8.919056e-02 1.6004230
# ANIM1.2       1.65423454 0.3916083 17.84392438  1 2.397836e-05 5.2290757
# ANIM1.3       0.08909168 0.1404036  0.40264083  1 5.257286e-01 1.0931809
# ANIM1.4       0.49063012 0.2423089  4.09986682  1 4.288659e-02 1.6333451
# FAM.LEX.f1.1 -0.11008672 0.2634346  0.17463216  1 6.760275e-01 0.8957565
# FAM.LEX.f1.2 -1.29487116 0.4251148  9.27770513  1 2.319599e-03 0.2739332
# FAM.LEX.f1.3 -0.27122998 0.1250447  4.70483864  1 3.007783e-02 0.7624411
# FAM.LEX.f1.4 -0.50976101 0.2219679  5.27415650  1 2.164427e-02 0.6006391
# EST11.1      -1.72073974 0.4740907 13.17369001  1 2.839071e-04 0.1789337
# EST11.2      -0.11825775 0.4543958  0.06773149  1 7.946688e-01 0.8884670
# EST11.3      -0.09674760 0.1383341  0.48912748  1 4.843168e-01 0.9077851
# EST11.4      -0.99596597 0.3098990 10.32877321  1 1.309722e-03 0.3693665
# EST51.1      -0.79202285 0.5072010  2.43845767  1 1.183926e-01 0.4529277
# EST51.3      -0.15545582 0.2165162  0.51550546  1 4.727652e-01 0.8560249
# EST51.4      -1.38824284 0.4094643 11.49473032  1 6.979378e-04 0.2495134
#               lower..95  upper..95
# Fabs.SC.f1.1 0.20879149  0.6860963
# Fabs.SC.f1.2 1.03005207  8.0310303
# Fabs.SC.f1.3 0.54186560  0.9820016
# Fabs.SC.f1.4 0.45633257  1.3075395
# MORF.f1.1    0.51508207  3.3838891
# MORF.f1.2    0.07584586  0.8359223
# MORF.f1.3    0.30891334  0.6075053
# MORF.f1.4    0.32102382  1.2932366
# MORF.f2.1    0.42133039  3.8734635
# MORF.f2.2    0.17459486  2.0984079
# MORF.f2.3    0.26303441  0.6193398
# MORF.f2.4    0.67540513  3.1002392
# MOD1.1       0.07920782  4.6005483
# MOD1.3       0.26021780  2.0178990
# MOD1.4       0.09727170  5.7090745
# MOD2.1       0.97518587  3.0847992
# MOD2.2       0.53728128  3.9207443
# MOD2.3       1.10102159  2.0584378
# MOD2.4       1.66789250  5.0131334
# MOD3.1       0.22706564  1.1066997
# MOD3.2       0.45528950  5.7505796
# MOD3.3       1.16381315  2.2745524
# MOD3.4       1.09882555  3.9206259
# ANIM1.1      0.93051743  2.7526122
# ANIM1.2      2.42709813 11.2658127
# ANIM1.3      0.83019528  1.4394739
# ANIM1.4      1.01583755  2.6262233
# FAM.LEX.f1.1 0.53450775  1.5011562
# FAM.LEX.f1.2 0.11906557  0.6302357
# FAM.LEX.f1.3 0.59671652  0.9741920
# FAM.LEX.f1.4 0.38875342  0.9280107
# EST11.1      0.07065542  0.4531469
# EST11.2      0.36463548  2.1648295
# EST11.3      0.69220203  1.1905105
# EST11.4      0.20121993  0.6780223
# EST51.1      0.16760948  1.2239372
# EST51.3      0.55999931  1.3085349
# EST51.4      0.11182970  0.5567118


# supuesto de proporcionalidad


shon.cr.1 <- cox.zph(model.frailty.cr.s.3, transform="rank")
shon.cr.2 <- cox.zph(model.frailty.cr.s.3, transform="km")
shon.cr.3 <- cox.zph(model.frailty.cr.s.3, transform=function(t){t})
shon.cr.4 <- cox.zph(model.frailty.cr.s.3, transform=function(t){log(t)})


shoenfeld.s.cr <- data.frame(df= c(rep(1,dim(shon.cr.1$table)[1]-1),dim(shon.cr.1$table)[1]-1),
                             rank.chi2 = shon.cr.1$table[,2], 
                             rank.p = shon.cr.1$table[,3],
                             km.chi2 = shon.cr.2$table[,2], km.p = shon.cr.2$table[,3],
                             t.chi2 = shon.cr.3$table[,2], t.p = shon.cr.3$table[,3],
                             log.t.chi2 = shon.cr.4$table[,2], log.t.p = shon.cr.4$table[,3])

xtbl.test.shoenfeld.s <- xtable(shoenfeld.s.cr, 
                                caption = "Test para theta = 0, según diferentes funciones del tiempo: estadístico de chi cuadrado y p valor", 
                                digits = 4,label = "xtable:example")


# df    rank.chi2     rank.p      km.chi2       km.p       t.chi2
# Fabs.SC.f1.1  1 1.424228e-01 0.70588424 1.718710e-01 0.67845433 2.095386e-01
# Fabs.SC.f1.2  1 8.737320e-02 0.76754353 4.866124e-02 0.82540941 5.896771e-02
# Fabs.SC.f1.3  1 3.481355e+00 0.06206390 4.714294e+00 0.02991285 4.994608e+00
# MORF.f1.2     1 1.430764e-01 0.70524162 6.048276e-02 0.80573454 5.582867e-02
# MORF.f1.3     1 2.120626e-02 0.88421848 2.166599e-01 0.64159679 2.567657e-01
# MORF.f2.3     1 2.631384e-01 0.60797221 3.643825e-01 0.54608232 3.385452e-01
# MOD2.3        1 9.175209e-01 0.33812674 9.863113e-01 0.32064560 1.278955e+00
# MOD2.4        1 1.413652e-01 0.70692760 1.876548e-02 0.89104095 4.364471e-02
# MOD3.3        1 1.297793e+00 0.25461684 1.041882e+00 0.30738410 9.623506e-01
# MOD3.4        1 1.132622e+00 0.28721682 8.330631e-01 0.36138830 8.268537e-01
# ANIM1.2       1 4.003649e-04 0.98403610 2.468846e-05 0.99603553 1.970500e-04
# ANIM1.4       1 1.097743e+00 0.29476196 1.470494e+00 0.22526781 1.364085e+00
# FAM.LEX.f1.2  1 1.446396e-01 0.70371147 5.152336e-02 0.82043341 5.869085e-02
# FAM.LEX.f1.3  1 6.575213e+00 0.01034086 4.846589e+00 0.02770086 4.266089e+00
# FAM.LEX.f1.4  1 1.441976e-01 0.70414320 5.083765e-03 0.94315857 8.732413e-06
# EST11.1       1 5.495928e+00 0.01906081 4.213090e+00 0.04011320 4.063396e+00
# EST11.4       1 4.203280e-03 0.94830720 1.775932e-01 0.67344981 1.269162e-01
# EST51.4       1 6.369048e-01 0.42483385 5.994693e-01 0.43878057 5.602237e-01
# GLOBAL       18 2.150997e+01 0.25446917 1.954361e+01 0.35909109 1.934260e+01
# t.p  log.t.chi2     log.t.p
# Fabs.SC.f1.1 0.64712926  0.51220851 0.474184367
# Fabs.SC.f1.2 0.80813502  0.59777154 0.439429551
# Fabs.SC.f1.3 0.02542642  1.85753179 0.172909938
# MORF.f1.2    0.81321470  0.77265294 0.379397594
# MORF.f1.3    0.61235090  0.02193186 0.882268507
# MORF.f2.3    0.56067020  0.16322087 0.686208285
# MOD2.3       0.25809340  0.77577196 0.378437353
# MOD2.4       0.83451603  0.45721397 0.498928917
# MOD3.3       0.32659535  0.91997728 0.337480938
# MOD3.4       0.36318389  2.03191138 0.154026791
# ANIM1.2      0.98880010  0.42883068 0.512563274
# ANIM1.4      0.24283075  1.07138917 0.300631875
# FAM.LEX.f1.2 0.80857721  1.25572754 0.262461362
# FAM.LEX.f1.3 0.03888033  7.15568507 0.007472675
# FAM.LEX.f1.4 0.99764220  0.86789493 0.351538741
# EST11.1      0.04382228  4.41858472 0.035549510
# EST11.4      0.72165107  0.26696459 0.605375257
# EST51.4      0.45417011  0.06974954 0.791701793
# GLOBAL       0.37101765 22.93069975 0.193267246

## modelo final

vet2 <- survSplit(Surv(time, status) ~ ., data= cr.long, cut=c(300), episode= "tgroup", id="ind")

model.frailty.cr.s.4=coxph(Surv(time, status) ~ Fabs.SC.f1.1 + Fabs.SC.f1.2 + Fabs.SC.f1.3:strata(tgroup) + 
                             MORF.f1.2 + MORF.f1.3 + MORF.f2.3 + MOD2.3 + MOD2.4 + MOD3.3 + MOD3.4 + ANIM1.2 + ANIM1.4 + 
                             FAM.LEX.f1.2 + FAM.LEX.f1.3:strata(tgroup) + FAM.LEX.f1.4 + EST11.1:strata(tgroup) + EST11.4 + EST51.4 + 
                             frailty(ID.SESION) + strata(trans), data=vet2) 

cox.zph(model.frailty.cr.s.4)

# chequeo de supuestos


ggcoxdiagnostics(model.frailty.cr.s.4, type = "dfbetas",
                 linear.predictions = FALSE,  ggtheme = theme_bw())

ggcoxdiagnostics(model.frailty.cr.s.4, type = "deviance",
                 linear.predictions = FALSE, ggtheme = theme_bw())

res.deviance <- residuals(model.frailty.cr.s.4, type = "deviance")

sum(abs(res.deviance)>2)/length(res.deviance) # 0.03758792, 513

# ojo: hay 1857*4 = 7428 observaciones !


shon.cr.1 <- cox.zph(model.frailty.cr.s.4, transform="rank")
shon.cr.2 <- cox.zph(model.frailty.cr.s.4, transform="km")
shon.cr.3 <- cox.zph(model.frailty.cr.s.4, transform=function(t){t})
shon.cr.4 <- cox.zph(model.frailty.cr.s.4, transform=function(t){log(t)})

shoenfeld.cr <- data.frame(df= c(rep(1,dim(shon.cr.1$table)[1]-1),dim(shon.cr.1$table)[1]-1),
                           rank.chi2 = shon.cr.1$table[,2], 
                           rank.p = shon.cr.1$table[,3],
                           km.chi2 = shon.cr.2$table[,2], km.p = shon.cr.2$table[,3],
                           t.chi2 = shon.cr.3$table[,2], t.p = shon.cr.3$table[,3],
                           log.t.chi2 = shon.cr.4$table[,2], log.t.p = shon.cr.4$table[,3])

xtbl.test.shoenfeld <- xtable(shoenfeld.cr, 
                              caption = "Test para theta = 0, según diferentes funciones del tiempo: estad??stico de chi cuadrado y p valor", 
                              digits = 4,label = "xtable:example")

# ajuste modelo final

size <- abs(1 - summary(model.frailty.cr.s.4)$conf.int[,1])*100
dir <- ifelse(size/100>1, "UP","DOWN")

sumario.cr <- data.frame(summary(model.frailty.cr.s.4)$conf.int[,c(1,3:4)],size, dir)

labels  <- rownames(summary(model.frailty.cr.s.4)$conf.int)
labels[16:21] <- c("Fabs.SC.f1.3:group=1", "Fabs.SC.f1.3:group=2", "FAM.LEX.f1.3:group=1",
                   "FAM.LEX.f1.3:group=2", "EST11.1:group=1", "EST11.1:group=2")
rownames(sumario.cr) <- labels


xtable(sumario.cr, 
       caption = "Modelo de riesgos competitivos. 
       exp(coef): hazard ratios, size: tama??o del efecto en porcentaje, dir: direcci??n del efecto  
       lower 95: extremo izquierdo de intervalo de confianza de 95 por ciento para hazard ratio,
       upper 95: extremo derecho de intervalo de confianza de 95 por ciento para hazard ratio.
       ", digits = 3,label = "xtable:example")


#Ajuste modelo final
#                       exp.coef.   lower..95 upper..95       size  dir
# Fabs.SC.f1.1         0.50818823 0.292835223 0.8819133  49.181177 DOWN
# Fabs.SC.f1.2         2.09255574 0.997689403 4.3889306 109.255574   UP
# MORF.f1.2            0.36089078 0.178402173 0.7300480  63.910922 DOWN
# MORF.f1.3            0.43318107 0.309637442 0.6060179  56.681893 DOWN
# MORF.f2.3            0.38375598 0.257680206 0.5715171  61.624402 DOWN
# MOD2.3               1.49772562 1.104046867 2.0317815  49.772562 DOWN
# MOD2.4               3.04283124 1.803576485 5.1335899 204.283124   UP
# MOD3.3               1.65179405 1.192863549 2.2872889  65.179405 DOWN
# MOD3.4               2.69673606 1.587915705 4.5798309 169.673606   UP
# ANIM1.2              4.65698315 2.244025882 9.6645463 365.698315   UP
# ANIM1.4              1.52008055 0.948475050 2.4361683  52.008055 DOWN
# FAM.LEX.f1.2         0.31141350 0.142229267 0.6818454  68.858650 DOWN
# FAM.LEX.f1.4         0.60295938 0.390659552 0.9306313  39.704062 DOWN
# EST11.4              0.32081179 0.176242619 0.5839689  67.918821 DOWN
# EST51.4              0.41308131 0.196070650 0.8702790  58.691869 DOWN
# Fabs.SC.f1.3:group=1 0.94878913 0.671267520 1.3410463   5.121087 DOWN
# Fabs.SC.f1.3:group=2 0.49090603 0.314336313 0.7666589  50.909397 DOWN
# FAM.LEX.f1.3:group=1 0.89319753 0.657726506 1.2129689  10.680247 DOWN
# FAM.LEX.f1.3:group=2 0.62132531 0.426999949 0.9040871  37.867469 DOWN
# EST11.1:group=1      0.04931303 0.006796346 0.3578063  95.068697 DOWN
# EST11.1:group=2      0.97884292 0.300665176 3.1867125   2.115708 DOWN



# hazard ratio: comparacion por tipo de error y variable 

forestplot(labeltext=labels, 
           mean=summary(model.frailty.cr.s.4)$conf.int[,1],
           lower=summary(model.frailty.cr.s.4)$conf.int[,3], 
           upper=summary(model.frailty.cr.s.4)$conf.int[,4],
           title="Hazard Ratio", zero=1)




