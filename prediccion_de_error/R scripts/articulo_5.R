library(mlr)

## Predicción de estatus de error.

## para reproducir los cuadros 8 a 11 simplemente 
## hay que cargar los resultados de los clasificadores.

load("RES_TODOS.RData")

for (i in 1:4){print(RES_TODOS[[1]][[i]])}

# para LATEX

for (i in 1:4){
  
  print(  
    xtable::xtable(RES_TODOS[[1]][[i]], digits = 3, 
                   caption = " (1) Clasifcadores: support vector machine (svm), random forest (rf),
rpart, Gradient Boosting Machine (gbm), eXtreme Gradient Boosting (xgb),
logistic regresion (logreg), ensemble; 
(2) tipos de pesos: undersampling (under), oversampling (over), smote, 
class weight (cw), base (sin pesos)")
  )
  
}

# inportancia de predictores (cuadro 12)

load("res80SONIA.RData")
imp_SONIA  <- mlr::getFeatureImportance(RES_80[[2]][[2]][[5]])$res
(imp_SONIA[order(imp_SONIA, decreasing = T)])[1:10]


load("res80NATI.RData")
imp_NATI  <- mlr::getFeatureImportance(RES_80[[2]][[3]][[2]])$res
(imp_NATI[order(imp_NATI, decreasing = T)])[1:10]

 
load("res80JAKO.RData")
imp_JAKO  <- mlr::getFeatureImportance(RES_80[[2]][[4]][[2]])$res
(imp_JAKO[order(imp_JAKO, decreasing = T)])[1:10]


load("res80MIRKA.RData")
imp_MIRKA  <- mlr::getFeatureImportance(RES_80[[2]][[3]][[5]])$res
(imp_MIRKA[order(imp_MIRKA, decreasing = T)])[1:10]


#################################################################
### lo siguiente muestra la manera de llegar a estos resultados.

### atributos dinámicos
## el resultado está en el siguiente archivo

load("DATA.LIST.RData")


# a continuación se muestra cómo llegar a ésto.
##################

library(tsfeatures)
library(wmtsa)
library(statcomp)
library(data.table)
library(zoo)
library(e1071)
library(crqa)

setwd("C:/Users/Juan/Desktop/PABLO/doctorado/tesis_de_doctorado")
load("datos_TESIS.RData")

## funciones


WT <- function(y){
  y <- as.numeric(y)
  WT <- wavMODWT(y, wavelet="haar",  keep.series=T)
  E <- summary(WT)[[1]][,10]
  VWT <- wavVar(y, wavelet="haar")
  V <- -log10(summary(VWT)$vmat[2,])
  res <- cbind(t(data.frame(E)), t(data.frame(V)))
  return(res)
}


## atributos a partir de la respuesta acumulada hasta la observacion anterior.

load("datos_articulo.RData")

datos <- datos_articulo

datos$RES_CUM <- unlist(tapply(as.numeric(as.character(datos$RES_BIN)), datos$ID, 
                               function(x){ cumsum(shift(x, fill=0))}))

DATA.LIST <- list()

DATA.TABLE <- data.table(Close=datos$RES_CUM, RES=datos$RESBIN ,ID=datos$ID, ID.SESION=datos$ID.SESION)

n <- 45

for (k in 1:4){
  
  DT <- DATA.TABLE[DATA.TABLE$ID==k,]
  DT[, Median:=rollapply(Close,n,median,fill=NA, align="right")]
  DT[, Burst:=rollapply(Close,n, function(x){(sd(x)-mean(x))/(sd(x)+mean(x))},fill=NA, align="right")]
  DT[, Max_High:=rollapply(Close,n, max,fill=NA, align="right")]
  DT[, Min_Low:=rollapply(Close,n,min,fill=NA, align="right")]
  DT[, Mean:=rollapply(Close,n,mean,fill=NA, align="right")]
  DT[, Mean_T:=rollapply(Close,n,function(x){mean(x,trim=0.05)},fill=NA, align="right")]
  DT[, SD_Diff:=rollapply(Close,n,function(x){sd(diff(x))},fill=NA, align="right")]
  DT[, SD:=rollapply(Close,n,sd,fill=NA, align="right")]
  DT[, CV:=rollapply(Close,n,function(x){sd(x)/mean(x)},fill=NA, align="right")]
  DT[, CV_Diff:=rollapply(Close,n,function(x){sd(diff(x))/mean(diff(x))},fill=NA, align="right")]
  DT[, Skew:=rollapply(Close,n,skewness,fill=NA, align="right")]
  DT[, Kurt:=rollapply(Close,n,kurtosis,fill=NA, align="right")]
  DT[, Ben_dist:=rollapply(Close,n,function(x){mean(x[which(x>mean(x))])/mean(x[which(x<mean(x))])},fill=NA, align="right")]
  DT[, Hs:=rollapply(Close,n,function(x){global_complexity(x,ndemb=2)[[1]]},fill=NA, align="right")]
  DT[, C:=rollapply(Close,n,function(x){global_complexity(x,ndemb=2)[[2]]},fill=NA, align="right")]
  DT[, FE:=rollapply(Close,n,function(x){fe <- which(diff(x)>0)[1]+1; 
  if (is.na(fe)) {return(0)} else {return(fe)}},fill=NA, align="right")]
  DT[, CP:=rollapply(Close,n,function(x){crossing_points(x)},fill=NA, align="right")]
  DT[, FS:=rollapply(Close,n,function(x){flat_spots(x)},fill=NA, align="right")]
  DT[, STAB:=rollapply(Close,n,function(x){stability(x)},fill=NA, align="right")]
  DT[, LUMP:=rollapply(Close,n,function(x){lumpiness(x)},fill=NA, align="right")]
  DT[, HURST:=rollapply(Close,n,function(x){hurst(x)},fill=NA, align="right")]
  DT[, ACF_1:=rollapply(Close,n,function(x){acf_features(x)[1]},fill=NA, align="right")]
  DT[, ACF_2:=rollapply(Close,n,function(x){firstmin_ac(x)},fill=NA, align="right")]
  DT[, ACF_3:=rollapply(Close,n,function(x){firstzero_ac(x)},fill=NA, align="right")]
  DT[, MOTIV:=rollapply(Close,n,function(x){motiftwo_entro3(x)},fill=NA, align="right")]
  DT[, SHIFT:=rollapply(Close,n,function(x){max_level_shift(x)[1]},fill=NA, align="right")]
  
  DT[, RR:= rollapply(DT$Close,n,function(x){
    crqa::crqa(x, x, delay = 1, embed = 2, radius = 0.0001, rescale = 0, 
               normalize = 0, mindiagline = 2, minvertline = 2)$RR},
    fill=NA, align="right")]
  
  DT[, DET:= rollapply(DT$Close,n,function(x){
    crqa::crqa(x, x, delay = 1, embed = 2, radius = 0.0001, rescale = 0, 
               normalize = 0, mindiagline = 2, minvertline = 2)$DET},
    fill=NA, align="right")]
  
  DT[, TT:= rollapply(DT$Close,n,function(x){
    crqa::crqa(x, x, delay = 1, embed = 2, radius = 0.0001, rescale = 0, 
               normalize = 0, mindiagline = 2, minvertline = 2)$TT},
    fill=NA, align="right")]
  
  DT[, LAM:= rollapply(DT$Close,n,function(x){
    crqa::crqa(x, x, delay = 1, embed = 2, radius = 0.0001, rescale = 0, 
               normalize = 0, mindiagline = 2, minvertline = 2)$LAM},
    fill=NA, align="right")]
  
  
  WAV <- rollapply(DT$Close,n,function(x){WT(x)},fill=NA, align="right")
  DT <- cbind(DT,WAV)
  
  DATA.LIST[[k]] <- DT
  
}

#save(DATA.LIST, file = "DATA.LIST.RData")

## atributos de las redes y de las cocnordancias

library(dplyr, warn.conflicts=FALSE)
library(tidyr, warn.conflicts=FALSE)
library(igraph)
load("DATA.LIST.RData")
load("ES_list.RData")

giant.component <- function(graph) { 
  cl <- clusters(graph) 
  induced.subgraph(graph, which(cl$membership == which.max(cl$csize)))} 

################################################################
#repetir para k = 1, 2, 3, 4
#--> el resultado final está en el archivo: lp_list.RData



featureS <- c("SESION","MOD","LDA","GRAMS","ES","MORF.f", "STEM.f", "CUMRES",
              "CONTR_WHERE","Fabs_C", "EST1", "EST2", "EST3", "EST4", "EST5", "EST6",
              "EST7", "IMA.CONC.f", "FAM.LEX.f", "ANIM","TYPE", "RES_BIN",
              
              colnames(DATA.LIST[[k]])[-c(1:4)])



lp.data  <- data.frame(token.matrix.list[[k]], datos[datos$ID==k, featureS[1:22]], 
                       DATA.LIST[[k]][,-c(1:4)])
colnames(lp.data) <- c("from_n", "to_n", featureS)

lp.red <- graph_from_edgelist(as.matrix(lp.data[,1:2]), directed = F)

E(lp.red)$sesion <- as.numeric(as.character(lp.data$SESION)) 
E(lp.red)$mod <- lp.data$MOD
E(lp.red)$es <-  lp.data$ES 
E(lp.red)$morf <- lp.data$MORF.f
E(lp.red)$stem <- lp.data$STEM.f
E(lp.red)$ima.conc <- lp.data$IMA.CONC.f
E(lp.red)$fam.lex <- lp.data$FAM.LEX.f 
E(lp.red)$anim <- lp.data$ANIM
E(lp.red)$est1 <- lp.data$EST1
E(lp.red)$est2 <- lp.data$EST2
E(lp.red)$est3 <- lp.data$EST3
E(lp.red)$est4 <- lp.data$EST4
E(lp.red)$est5 <- lp.data$EST5
E(lp.red)$est6 <- lp.data$EST6
E(lp.red)$est7 <- lp.data$EST7
E(lp.red)$gram <- lp.data$GRAMS 
E(lp.red)$type <- lp.data$TYPE
E(lp.red)$res <- as.numeric(as.character(lp.data$RES_BIN)) 
E(lp.red)$lda <-  lp.data$LDA 
E(lp.red)$cumres <-  lp.data$CUMRES


ESP_1 <- as.vector(ifelse(ES.df.list[[k]]$ES1==1, "*es", ES.df.list[[k]]$ESP1))
ESP_2 <- as.vector(ifelse(ES.df.list[[k]]$ES2==1, "*es", ES.df.list[[k]]$ESP2))

ESP_1[ESP_1=="*es"] <- 6
ESP_1[ESP_1=="es"] <- 5
ESP_1[ESP_1=="as"] <- 4
ESP_1[ESP_1=="os"] <- 3
ESP_1[ESP_1=="is"] <- 2
ESP_1[ESP_1=="us"] <- 1
ESP_1 <- as.integer(ESP_1)

ESP_2[ESP_2=="*es"] <- 6
ESP_2[ESP_2=="es"] <- 5
ESP_2[ESP_2=="as"] <- 4
ESP_2[ESP_2=="os"] <- 3
ESP_2[ESP_2=="is"] <- 2
ESP_2[ESP_2=="us"] <- 1
ESP_2 <- as.integer(ESP_2)

E(lp.red)$esp1 <- ESP_1
E(lp.red)$esp2 <- ESP_2

DATA.LIST.K  <- as.data.frame(DATA.LIST[[k]])

E(lp.red)$Median <- DATA.LIST.K[,"Median"]
E(lp.red)$Burst <- DATA.LIST.K[,"Burst"]
E(lp.red)$Max_High <- DATA.LIST.K[,"Max_High"]
E(lp.red)$Min_Low <- DATA.LIST.K[,"Min_Low"]
E(lp.red)$Mean <- DATA.LIST.K[,"Mean"]
E(lp.red)$Mean_T  <- DATA.LIST.K[,"Mean_T"]
E(lp.red)$SD_Diff  <- DATA.LIST.K[,"SD_Diff"]
E(lp.red)$SD  <- DATA.LIST.K[,"SD"]
E(lp.red)$CV  <- DATA.LIST.K[,"CV"]
E(lp.red)$CV_Diff  <- DATA.LIST.K[,"CV_Diff"]
E(lp.red)$Skew  <- DATA.LIST.K[,"Skew"]
E(lp.red)$Kurt  <- DATA.LIST.K[,"Kurt"]
E(lp.red)$Ben_dist  <- DATA.LIST.K[,"Ben_dist"]
E(lp.red)$Hs  <- DATA.LIST.K[,"Hs"]
E(lp.red)$C  <- DATA.LIST.K[,"C"]
E(lp.red)$FE  <- DATA.LIST.K[,"FE"]
E(lp.red)$CP  <- DATA.LIST.K[,"CP"]
E(lp.red)$STAB  <- DATA.LIST.K[,"STAB"]
E(lp.red)$LUMP  <- DATA.LIST.K[,"LUMP"]
E(lp.red)$HURST <- DATA.LIST.K[,"HURST"]
E(lp.red)$ACF_1 <- DATA.LIST.K[,"ACF_1"]
E(lp.red)$ACF_2 <- DATA.LIST.K[,"ACF_2"]
E(lp.red)$ACF_3 <- DATA.LIST.K[,"ACF_3"]
E(lp.red)$MOTIV <- DATA.LIST.K[,"MOTIV"]
E(lp.red)$SHIFT <- DATA.LIST.K[,"SHIFT"]
E(lp.red)$RR <- DATA.LIST.K[,"RR"]
E(lp.red)$DET <- DATA.LIST.K[,"DET"]
E(lp.red)$TT <- DATA.LIST.K[,"TT"]
E(lp.red)$LAM <- DATA.LIST.K[,"LAM"]
E(lp.red)$d1 <- DATA.LIST.K[,"d1"]
E(lp.red)$d2 <- DATA.LIST.K[,"d2"]
E(lp.red)$d3 <- DATA.LIST.K[,"d3"]
E(lp.red)$d4 <- DATA.LIST.K[,"d4"]
E(lp.red)$d5 <- DATA.LIST.K[,"d5"]
E(lp.red)$s5 <- DATA.LIST.K[,"s5"]
E(lp.red)$D1 <- DATA.LIST.K[,41]
E(lp.red)$D2 <- DATA.LIST.K[,42]
E(lp.red)$D3 <- DATA.LIST.K[,43]
E(lp.red)$D4 <- DATA.LIST.K[,44]
E(lp.red)$D5 <- DATA.LIST.K[,45]


lp.giant <- giant.component(lp.red) 

lp.giant.edge <- as.data.frame(cbind(get.edgelist(lp.giant), get.edgelist(lp.giant, names=FALSE)))
colnames(lp.giant.edge) <- c("from_n", "to_n", "from", "to")
lp.giant.edge$from <- as.numeric(as.character(lp.giant.edge$from)) 
lp.giant.edge$to <- as.numeric(as.character(lp.giant.edge$to)) 

lp.cosi <- linkprediction::proxfun(lp.giant, method="cos_l", value="edgelist") %>%
  filter(from < to) %>%
  rename(cosi=value)
lp.lp <- linkprediction::proxfun(lp.giant, method="lp", value="edgelist") %>%
  filter(from < to) %>%
  rename(lp=value)
lp.pa <- linkprediction::proxfun(lp.giant, method="pa", value="edgelist") %>%
  filter(from < to) %>%
  rename(pa=value)
lp.katz <- linkprediction::proxfun(lp.giant, method="katz", value="edgelist") %>%
  filter(from < to) %>%
  rename(ka=value)
lp.act_n <- linkprediction::proxfun(lp.giant, method="act_n", value="edgelist") %>%
  filter(from < to) %>%
  rename(act_n=value)
lp.rwr <- linkprediction::proxfun(lp.giant, method="rwr", value="edgelist") %>%
  filter(from < to) %>%
  rename(rwr=value)
lp.l <- linkprediction::proxfun(lp.giant, method="l", value="edgelist") %>%
  filter(from < to) %>%
  rename(l=value)
lp.mf <- linkprediction::proxfun(lp.giant, method="mf", value="edgelist") %>%
  filter(from < to) %>%
  rename(mf=value)
lp.act <- linkprediction::proxfun(lp.giant, method="act", value="edgelist") %>%
  filter(from < to) %>%
  rename(act=value)
# lp.hpi <- linkprediction::proxfun(lp.giant, method="hpi", value="edgelist") %>%
#   filter(from < to) %>%
#   rename(hpi=value)
# lp.hdi <- linkprediction::proxfun(lp.giant, method="hdi", value="edgelist") %>%
#   filter(from < to) %>%
#   rename(hdi=value)
# lp.lhn_local <- linkprediction::proxfun(lp.giant, method="lhn_local", value="edgelist") %>%
#   filter(from < to) %>%
#   rename(lhn_local=value)

lp <- lp.giant.edge %>%
  left_join(lp.pa, by=c("from", "to")) %>%
  left_join(lp.cosi, by=c("from", "to")) %>%
  left_join(lp.katz, by=c("from", "to")) %>%
  left_join(lp.act_n, by=c("from", "to")) %>%
  left_join(lp.act, by=c("from", "to")) %>%
  left_join(lp.rwr, by=c("from", "to")) %>%
  left_join(lp.l, by=c("from", "to")) %>%
  left_join(lp.mf, by=c("from", "to")) %>%
  #left_join(lp.lhn, by=c("from", "to")) %>%
  #left_join(lp.hpi, by=c("from", "to")) %>%
  #left_join(lp.hdi, by=c("from", "to")) %>%
  #left_join(lp.lhn_local, by=c("from", "to")) %>%
  mutate_at(
    c("pa", "cosi", "ka", "act_n", "act","rwr", "l", "mf"),
    funs(ifelse(is.na(.), 0, .))
  )

lp$sesion <- E(lp.giant)$sesion
lp$mod <- factor(E(lp.giant)$mod-1)
lp$es <- factor(E(lp.giant)$es-1)
lp$morf <- factor(E(lp.giant)$morf-1)
lp$stem <- factor(E(lp.giant)$stem-1)
lp$ima.conc <- factor(E(lp.giant)$ima.conc-1)
lp$fam.lex <- factor(E(lp.giant)$fam.lex-1)
lp$anim <- factor(E(lp.giant)$anim-1)
lp$est1 <- factor(E(lp.giant)$est1-1)
lp$est2 <- factor(E(lp.giant)$est2-1)
lp$est3 <- factor(E(lp.giant)$est3-1)
lp$est4 <- factor(E(lp.giant)$est4-1)
lp$est5 <- factor(E(lp.giant)$est5-1)
lp$est6 <- factor(E(lp.giant)$est6-1)
lp$est7 <- factor(E(lp.giant)$est7-1)
lp$target <- factor(E(lp.giant)$res)
lp$lda <- factor(E(lp.giant)$lda-1)
lp$gram <- factor(E(lp.giant)$gram-1)
lp$esp1 <- factor(E(lp.giant)$esp1)
lp$esp2 <- factor(E(lp.giant)$esp2)
lp$cumres <- E(lp.giant)$cumres

lp$Median <- E(lp.giant)$Median
lp$Burst <- E(lp.giant)$Burst
lp$Max_High <- E(lp.giant)$Max_High
lp$Min_Low <- E(lp.giant)$Min_Low
lp$Mean <- E(lp.giant)$Mean
lp$Mean_T <- E(lp.giant)$Mean_T
lp$SD_Diff <- E(lp.giant)$SD_Diff
lp$SD <- E(lp.giant)$SD
lp$CV <- E(lp.giant)$CV
lp$CV_Diff <- E(lp.giant)$CV_Diff
lp$Skew <- E(lp.giant)$Skew
lp$Kurt <- E(lp.giant)$Kurt
lp$Hs <- E(lp.giant)$Hs
lp$C <- E(lp.giant)$C
lp$FE <- E(lp.giant)$FE
lp$CP <- E(lp.giant)$CP
lp$FS <- E(lp.giant)$FS
lp$STAB <- E(lp.giant)$STAB
lp$LUMP <- E(lp.giant)$LUMP
lp$HURST <- E(lp.giant)$HURST
lp$ACF_1 <- E(lp.giant)$ACF_1
lp$ACF_2 <- E(lp.giant)$ACF_2
lp$ACF_3 <- E(lp.giant)$ACF_3
lp$MOTIV <- E(lp.giant)$MOTIV
lp$SHIFT <- E(lp.giant)$SHIFT
lp$RR <- E(lp.giant)$RR
lp$DET <- E(lp.giant)$DET
lp$TT <- E(lp.giant)$TT
lp$LAM <- E(lp.giant)$LAM
lp$d1 <- E(lp.giant)$d1
lp$d2 <- E(lp.giant)$d2
lp$d3 <- E(lp.giant)$d3
lp$d4 <- E(lp.giant)$d4
lp$d5 <- E(lp.giant)$d5
lp$s5 <- E(lp.giant)$s5
lp$D1 <- E(lp.giant)$D1
lp$D2 <- E(lp.giant)$D2
lp$D3 <- E(lp.giant)$D3
lp$D4 <- E(lp.giant)$D4
lp$D5 <- E(lp.giant)$D5


nombre_enlace  <- apply(lp[,c("from_n","to_n")],1, function(x){paste(x, collapse=" ")})
contar_enlace  <- table(nombre_enlace)
freq_enlace <- sapply(1:length(nombre_enlace), function(x){ contar_enlace[match(nombre_enlace[x], names(contar_enlace))]})
#freq_error  <- tapply(as.numeric(as.character(lp$target)), nombre_enlace, sum) 

lp$freq_enlace <- as.vector(freq_enlace)


## features: 1 grams

dtm <- quanteda::dfm(as.vector(nombre_enlace), verbose = TRUE, toLower = TRUE,
                     removeNumbers = FALSE, removePunct = TRUE, removeSeparators = TRUE,
                     removeTwitter = FALSE, stem = FALSE,
                     keptFeatures = NULL, language = "spanish", ngrams=1) #con mono- y bi-gramas!
features.freq <- quanteda::topfeatures(dtm, n=50) #los 300 features m??s frecuentes
dtm.freq <- dtm[,names(features.freq)] # selecci??n de columnas
dtm.tfidf <- quanteda::dfm_tfidf(dtm.freq) 

lp <- data.frame(lp, quanteda::convert(dtm.tfidf, to = "data.frame"))

# lista de bases de datos
# lp_SONIA <- lp
# lp_NATI <- lp
# lp_JAKO <- lp
# lp_MIRKA <- lp
# lp_list <- list(lp_SONIA, lp_NATI, lp_JAKO, lp_MIRKA)

# save(lp_list, file = "lp_list.RData")

##########################################################
##conjuntos de entrenamiento y validación 
#
# --> el resultado final está en el archivo: tasks.list.RData

col_out <- match(c("from_n","to_n","from","to","sesion","est3", "document"),
                 colnames(lp))

#col_out <- match(c("from_n","to_n","from","to","sesion","est3"),
#                                 colnames(lp))

if (k==1|k==4){barrier <- 9}else{barrier <- 11}

lp.train <- lp[lp$sesion<barrier,-col_out]
#lp.train.res <- lp[lp$sesion<9,"target"]
lp.test <- data.frame(lp[lp$sesion>(barrier-1),-col_out])
#target_test <- lp[lp$sesion>8, "target"]

labels <- lp.train$target
ts_label <- lp.test$target

target.id  <- which(colnames(lp.train)=="target")

new_tr <- model.matrix(~.+0,
                       model.frame(~.+0,data = lp.train[,-target.id],na.action=NULL))            
new_ts <- model.matrix(~.+0,
                       model.frame(~.+0,data = lp.test[,-target.id],na.action=NULL))  

new_tr <- apply(new_tr,2,as.numeric)
new_ts <- apply(new_ts,2,as.numeric)

labels <- as.numeric(labels)-1
ts_label <- as.numeric(ts_label)-1



library(mlr)

library(parallel)
library(parallelMap) 
parallelStartSocket(cpus = 10)
parallelLibrary("mlr")

set.seed(123, "L'Ecuyer") 


train_set  <-  apply(new_tr, 2, function(x){ifelse(is.na(x),0,x)})
test_set  <-  apply(new_ts, 2, function(x){ifelse(is.na(x),0,x)})


train_set1 <- train_set[,-c(9:35, 77:dim(train_set)[2])]
test_set1 <-  test_set[,-c(9:35, 77:dim(test_set)[2])]
train_set1  <- apply(train_set1,2, function(x){BBmisc::normalize(x, method = "standardize")} ) 
test_set1  <- apply(test_set1,2,function(x){BBmisc::normalize(x, method = "standardize")})

train_set_task <- cbind(train_set[,c(9:35, 77:126)], train_set1)
test_set_task <- cbind(test_set[, c(9:35, 77:126)], test_set1)

traintask <- makeClassifTask(data = data.frame(train_set_task, y = factor(labels)),
                             target = "y", positive = "1")
testtask <- makeClassifTask(data = data.frame(test_set_task, y = factor(ts_label)),
                            target = "y", positive = "1")


## primer filtrado de datos

filtered.train.task = filterFeatures(traintask, abs = 80,
                                     base.methods = c("FSelector_relief", "auc","praznik_JMI", 
                                                      "randomForest_importance"))                           

traintask.features.filtered   <- getTaskFeatureNames(filtered.train.task)    

testtask.features  <- getTaskFeatureNames(testtask) 

features.drop <-   testtask.features[-match(traintask.features.filtered, testtask.features)]                       
filtered.test.task   <- dropFeatures(testtask, features.drop)

#traintask <- removeConstantFeatures(traintask)
#testtask <- removeConstantFeatures(testtask)
#trainTask < - normalizeFeatures(traintask, method = "standardize", 
#                                on.constant = "quiet")
#testTask <- normalizeFeatures(testtask, method = "standardize", 
#                              on.constant = "quiet")



#imbalance::imbalanceRatio(data.frame(labels), classAttr="labels")

tr.task <- filtered.train.task
ts.task <- filtered.test.task

#tr.task <-   traintask
#ts.task <-   testtask

tr.task.list <- list()
ts.task.list <- list()

tr.task.list[[1]] <-  tr.task
ts.task.list[[1]] <-  ts.task

tr.task.list[[2]] <-  tr.task
ts.task.list[[2]] <-  ts.task

tr.task.list[[3]] <-  tr.task
ts.task.list[[3]] <-  ts.task

tr.task.list[[4]] <-  tr.task
ts.task.list[[4]] <-  ts.task

tasks.list <- list(tr.task.list, ts.task.list)

save(tasks.list, file = "tasks.list.RData")

######################################################
##imbalance
#
# --> la matriz de costos final está en el archivo: cost.measure.list.RData

IR <- tr.task$task.desc$class.distribution[1] / tr.task$task.desc$class.distribution[2]


#sonia: 0-1 = 1; 1-0 = 6
#jako: 0-1 =  1; 1-0 = 10
#nati: 0-1 =  2; 1-0 = 6; 1-1 = -2
#mirka: 0-1 = 1; 1-0 = 2; 1-1 = -1

cost.measure.list <- list()
cero_uno <- c(1,1,2,1)
uno_cero <- c(6,10,6,2)
uno_uno <- c(0,0,-2,-1)

for (i in 1:4){
  
  costo_pos  <- uno_cero[i]
  costo_neg <-  cero_uno[i]
  costo_uno_uno <- uno_uno[i]
  costo_max <- max(costo_pos, costo_neg, costo_uno_uno)
  
  costs = matrix(c(0, costo_pos, costo_neg, costo_uno_uno), 2)
  colnames(costs) = rownames(costs) = getTaskClassLevels(tr.task)
  costs.measure = makeCostMeasure(id = "costs", name = "costs",
                                  costs = costs, best = 0, worst = costo_max)
  
  cost.measure.list[[i]] <- costs.measure
  
}

#save(cost.measure.list, file = "cost.measure.list.RData")

### clasificadores


lrn_svm <- makeLearner("classif.ksvm", predict.type = "prob",
                       scaled = F, fix.factors.prediction = TRUE )

lrn_rpart <- makeLearner("classif.rpart", predict.type = "prob",
                         fix.factors.prediction = TRUE )

lrn_rf <- makeLearner("classif.randomForest", predict.type = "prob",
                      fix.factors.prediction = T)
lrn_rf$par.vals <- list(importance = TRUE)

lrn_logreg <- makeLearner("classif.logreg", predict.type = "prob",
                          fix.factors.prediction = TRUE )

lrn_gbm <- makeLearner("classif.gbm", predict.type = "prob",
                       fix.factors.prediction = TRUE )
lrn_gbm$par.vals <- list(train.fraction = 0.7)
# lrn_gbm$par.set = c(lrn_gbm$par.set, makeParamSet(
#   makeNumericLearnerParam("nTrain"))) 
# lrn_gbm$par.set$pars$nTrain$default <- tr.task$task.desc$size
# lrn_gbm$par.vals <- list(nTrain = tr.task$task.desc$size)


#lrn_gbm <- makeLearner("classif.h2o.gbm", predict.type = "prob",
#                      fix.factors.prediction = TRUE )


lrn_xgb <- makeLearner("classif.xgboost", predict.type = "prob",
                       fix.factors.prediction = TRUE )



lrn_dl <- makeLearner("classif.h2o.deeplearning", predict.type = "prob",
                      fix.factors.prediction = TRUE, 
                      par.vals = list(shuffle_training_data = T, fast_mode = T, l1 = 1e-6)
)



pars.svm   <- 
  list(
    #makeNumericParam("epsilon", lower = 0.01, upper = 0.05),
    makeNumericParam("C", lower = -8, upper = 15, trafo = function(x) 2^x),
    makeNumericParam("sigma", lower = -15, upper = 3, trafo = function(x) 2^x),
    makeNumericParam("wcw.weight", lower = 1, upper = 10),
    makeNumericParam("usw.rate", lower = 0.6*(1/IR), upper = 1),
    makeNumericParam("osw.rate", lower = 1, upper = 1.5*IR),
    makeNumericParam("sw.rate", lower = 1, upper = 1.5*IR),
    makeDiscreteParam("fw.abs", values = seq_len(getTaskNFeats(tr.task)))
    
  ) 


pars <- pars.svm

pars_svm_list <- list(makeParamSet(pars[[1]], pars[[2]],pars[[7]]),
                      makeParamSet(pars[[1]], pars[[2]], pars[[3]], pars[[7]]),
                      makeParamSet(pars[[1]], pars[[2]], pars[[4]], pars[[7]]),
                      makeParamSet(pars[[1]], pars[[2]], pars[[5]], pars[[7]]),
                      makeParamSet(pars[[1]], pars[[2]], pars[[6]], pars[[7]])
)

# pars_svm_list <- list(makeParamSet(pars[[1]], pars[[2]]),
#                       makeParamSet(pars[[1]], pars[[2]], pars[[3]]),
#                       makeParamSet(pars[[1]], pars[[2]], pars[[4]]),
#                       makeParamSet(pars[[1]], pars[[2]], pars[[5]]),
#                       makeParamSet(pars[[1]], pars[[2]], pars[[6]])
# )




pars_rpart   <- 
  list(
    makeNumericParam("cp", lower = 0.0001, upper = 0.1),
    makeIntegerParam("minsplit", lower = 1, upper = 10),
    makeIntegerParam("minbucket", lower = 5, upper = 50),
    makeNumericParam("wcw.weight", lower = 1, upper = 10),
    makeNumericParam("usw.rate", lower = 0.6*(1/IR), upper = 1),
    makeNumericParam("osw.rate", lower = 1, upper = 1.5*IR),
    makeNumericParam("sw.rate", lower = 1, upper = 1.5*IR),
    makeDiscreteParam("fw.abs", values = seq_len(getTaskNFeats(tr.task)))
  ) 

pars <- pars_rpart 

pars_rpart_list <- list(makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[8]]),
                        makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[4]],  pars[[8]]),
                        makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[5]],  pars[[8]]),
                        makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[6]],  pars[[8]]),
                        makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[7]],  pars[[8]])
)

# pars_rpart_list <- list(makeParamSet(pars[[1]], pars[[2]], pars[[3]]),
# makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[4]]),
# makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[5]]),
# makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[6]]),
# makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[7]])
# )

pars_rf   <- 
  list(
    makeIntegerParam("ntree",lower = 50, upper = 500),
    makeIntegerParam("mtry", lower = 3, upper = 10),
    makeIntegerParam("nodesize", lower = 10, upper = 50),
    makeNumericParam("wcw.weight", lower = 1, upper = 10),
    makeNumericParam("usw.rate", lower = 0.6*(1/IR), upper = 1),
    makeNumericParam("osw.rate", lower = 1, upper = 1.5*IR),
    makeNumericParam("sw.rate", lower = 1, upper = 1.5*IR),
    makeDiscreteParam("fw.abs", values = seq_len(getTaskNFeats(tr.task)))
  ) 

pars <- pars_rf 

pars_rf_list <- list(makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[8]]),
                     makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[4]],  pars[[8]]),
                     makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[5]],  pars[[8]]),
                     makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[6]],  pars[[8]]),
                     makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[7]],  pars[[8]])
)

# pars_rf_list <- list(makeParamSet(pars[[1]], pars[[2]], pars[[3]]),
# makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[4]]),
# makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[5]]),
# makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[6]]),
# makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[7]])
# )


pars_xgb   <- 
  list(
    makeIntegerParam("nrounds",lower=200,upper=600),
    makeIntegerParam("max_depth",lower=3,upper=20),
    makeNumericParam("lambda",lower = 2^-10,upper = 2^10),
    makeNumericParam("wcw.weight", lower = 1, upper = 10),
    makeNumericParam("usw.rate", lower = 0.6*(1/IR), upper = 1),
    makeNumericParam("osw.rate", lower = 1, upper = 1.5*IR),
    makeNumericParam("sw.rate", lower = 1, upper = 1.5*IR),
    makeDiscreteParam("fw.abs", values = seq_len(getTaskNFeats(tr.task)))
  ) 

pars <- pars_xgb

pars_xgb_list <- list(makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[8]]),
                      makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[4]],  pars[[8]]),
                      makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[5]],  pars[[8]]),
                      makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[6]],  pars[[8]]),
                      makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[7]],  pars[[8]])
)

# pars_xgb_list <- list(makeParamSet(pars[[1]], pars[[2]], pars[[3]]),
# makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[4]]),
# makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[5]]),
# makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[6]]),
# makeParamSet(pars[[1]], pars[[2]], pars[[3]],  pars[[7]])
# )


pars_gbm   <- 
  list(
    makeIntegerParam("n.trees", lower = 100, upper = 5000), #number of trees
    makeIntegerParam("interaction.depth", lower = 2, upper = 10), #depth of tree
    makeNumericParam("bag.fraction", lower = 0.7, upper = 1),
    makeIntegerParam("n.minobsinnode", lower = 5, upper = 15),
    makeNumericParam("shrinkage",lower = 0.001, upper = 0.5),
    makeNumericParam("wcw.weight", lower = 1, upper = 10),
    makeNumericParam("usw.rate", lower = 0.6*(1/IR), upper = 1),
    makeNumericParam("osw.rate", lower = 1, upper = 1.5*IR),
    makeNumericParam("sw.rate", lower = 1, upper = 1.5*IR),
    makeDiscreteParam("distribution", values = "bernoulli"),
    makeDiscreteParam("fw.abs", values = seq_len(getTaskNFeats(tr.task)))
  ) 

pars <- pars_gbm

pars_gbm_list <- list(makeParamSet(pars[[1]], pars[[2]], pars[[3]], pars[[4]], pars[[5]],pars[[10]], pars[[11]]),
                      makeParamSet(pars[[1]], pars[[2]], pars[[3]],pars[[4]],  pars[[5]],pars[[6]],  pars[[10]], pars[[11]]),
                      makeParamSet(pars[[1]], pars[[2]], pars[[3]], pars[[4]], pars[[5]],pars[[7]],  pars[[10]], pars[[11]]),
                      makeParamSet(pars[[1]], pars[[2]], pars[[3]],pars[[4]],  pars[[5]],pars[[8]],  pars[[10]], pars[[11]]),
                      makeParamSet(pars[[1]], pars[[2]], pars[[3]], pars[[4]], pars[[5]],pars[[9]],  pars[[10]], pars[[11]])
)

# pars_gbm_list <- list(makeParamSet(pars[[1]], pars[[2]], pars[[3]], pars[[4]], pars[[5]],pars[[10]]),
# makeParamSet(pars[[1]], pars[[2]], pars[[3]],pars[[4]],  pars[[5]],pars[[6]],  pars[[10]]),
# makeParamSet(pars[[1]], pars[[2]], pars[[3]], pars[[4]], pars[[5]],pars[[7]],  pars[[10]]),
# makeParamSet(pars[[1]], pars[[2]], pars[[3]],pars[[4]],  pars[[5]],pars[[8]],  pars[[10]]),
# makeParamSet(pars[[1]], pars[[2]], pars[[3]], pars[[4]], pars[[5]],pars[[9]],  pars[[10]])
# )


# pars_gbm   <- 
#   list(
#     makeIntegerParam("n.trees", lower = 100, upper = 5000), #number of trees
#     makeIntegerParam("max_depth", lower = 2, upper = 10), #depth of tree
#     makeIntegerParam("min_rows", lower = 5, upper = 20),
#     makeNumericParam("learn_rate",lower = 0.001, upper = 0.5),
#     makeNumericParam("wcw.weight", lower = 1, upper = 10),
#     makeNumericParam("usw.rate", lower = 0.6*(1/IR), upper = 1),
#     makeNumericParam("osw.rate", lower = 1, upper = 1.5*IR),
#     makeNumericParam("sw.rate", lower = 1, upper = 1.5*IR),
#     makeDiscreteParam("fw.abs", values = seq_len(getTaskNFeats(tr.task))),
#   ) 
# 
# pars <- pars_gbm
# 
# pars_gbm_list <- list(makeParamSet(pars[[1]], pars[[2]], pars[[3]], pars[[4]],pars[[9]]),
#                       makeParamSet(pars[[1]], pars[[2]], pars[[3]],pars[[4]], pars[[5]],  pars[[9]]),
#                       makeParamSet(pars[[1]], pars[[2]], pars[[3]], pars[[4]],pars[[6]],  pars[[9]]),
#                       makeParamSet(pars[[1]], pars[[2]], pars[[3]],pars[[4]], pars[[7]],  pars[[9]]),
#                       makeParamSet(pars[[1]], pars[[2]], pars[[3]], pars[[4]],pars[[8]],  pars[[9]])
# )



pars_reglog   <- 
  list(
    makeNumericParam("wcw.weight", lower = 1, upper = 10),
    makeNumericParam("usw.rate", lower = 0.6*(1/IR), upper = 1),
    makeNumericParam("osw.rate", lower = 1, upper = 1.5*IR),
    makeNumericParam("sw.rate", lower = 1, upper = 1.5*IR),
    makeDiscreteParam("fw.abs", values = seq_len(getTaskNFeats(tr.task)))
  ) 

pars <- pars_reglog

pars_reglog_list <- list(makeParamSet(pars[[5]]),
                         makeParamSet(pars[[1]],  pars[[5]]),
                         makeParamSet(pars[[2]],  pars[[5]]),
                         makeParamSet(pars[[3]],  pars[[5]]),
                         makeParamSet(pars[[4]],  pars[[5]])
)

# pars_reglog_list <- list(NULL,
#                          makeParamSet(pars[[1]]),
#                          makeParamSet(pars[[2]]),
#                          makeParamSet(pars[[3]]),
#                          makeParamSet(pars[[4]])
# )






pars_dl   <- 
  list(
    makeDiscreteParam("hidden", values = list(a = c(50,50), b = c(80,80), 
                                              c = c(128,64,32), d = c(32,32,32,32))),
    #makeIntegerVectorParam("hidden", len = 3, lower = 10, upper = 300),
    makeIntegerParam("epochs", lower = 30, upper = 300),
    makeNumericParam("rho", lower = 0.95, upper = 0.999),
    makeNumericParam("epsilon", lower = 1e-8, upper = 1e-6),
    #makeDiscreteParam("activation", values = c("Rectifier", "Tanh", "Maxout")),
    #makeDiscreteParam("activation", values = c("Rectifier","Tanh","Maxout",
    #          "RectifierWithDropout","TanhWithDropout","MaxoutWithDropout")),
    #makeNumericParam("l1", lower = 0.0001, upper = 1),
    #makeNumericParam("l2", lower = 0.0001, upper = 1),
    #makeNumericParam("input_dropout_ratio", lower = 0, upper = 0.9),
    #makeIntegerVectorParam("hidden", len = 3, lower = 10, upper = 300)
    #makeDiscreteParam("rate", values = list(0.01,0.02)),
    #makeDiscreteParam("rate_annealing", values = list(1e-8,1e-7,1e-6))                
    #makeNumericParam("rate", lower = 0.01, upper = 0.03),
    #makeNumericParam("rate_annealing", lower = 1e-8, upper = 1e-6)
    makeNumericParam("wcw.weight", lower = 1, upper = 10),
    makeNumericParam("usw.rate", lower = 0.6*(1/IR), upper = 1),
    makeNumericParam("osw.rate", lower = 1, upper = 1.5*IR),
    makeNumericParam("sw.rate", lower = 1, upper = 1.5*IR),
    makeDiscreteParam("fw.abs", values = seq_len(getTaskNFeats(tr.task)))
  ) 

pars <- pars_dl

pars_dl_list <- list(makeParamSet(pars[[1]], pars[[2]], pars[[3]], pars[[4]],pars[[9]]),
                     makeParamSet(pars[[1]], pars[[2]], pars[[3]],pars[[4]], pars[[5]],  pars[[9]]),
                     makeParamSet(pars[[1]], pars[[2]], pars[[3]], pars[[4]],pars[[6]],  pars[[9]]),
                     makeParamSet(pars[[1]], pars[[2]], pars[[3]],pars[[4]], pars[[7]],  pars[[9]]),
                     makeParamSet(pars[[1]], pars[[2]], pars[[3]], pars[[4]],pars[[8]],  pars[[9]])
)

# pars_dl_list <- list(makeParamSet(pars[[1]], pars[[2]], pars[[3]], pars[[4]]),
# makeParamSet(pars[[1]], pars[[2]], pars[[3]],pars[[4]], pars[[5]]),
# makeParamSet(pars[[1]], pars[[2]], pars[[3]], pars[[4]],pars[[6]]),
# makeParamSet(pars[[1]], pars[[2]], pars[[3]],pars[[4]], pars[[7]]),
# makeParamSet(pars[[1]], pars[[2]], pars[[3]], pars[[4]],pars[[8]])
# )




clasif.list <- c(lrn_svm$id, lrn_rf$id, lrn_rpart$id,
                 lrn_xgb$id, lrn_logreg$id, lrn_gbm$id, lrn_dl$id)

ps.list <- list(pars_svm_list, pars_rf_list, pars_rpart_list, 
                pars_xgb_list, pars_reglog_list, pars_gbm_list, pars_dl_list )

ctrl = makeTuneMultiCritControlRandom(maxit = 100L)
rdesc = makeResampleDesc("CV", iters = 5L)

n_clasif <- length(clasif.list)

res.list <- list(list(), list(), list(), list(), list(), list(), list())
model.list <- list(list(), list(), list(), list(), list(), list(), list())
pred.list <- list(list(), list(), list(), list(), list(), list(), list())
CM.LIST <- list(list(), list(), list(), list(), list(), list(), list())
chart.eval.list <- vector("list", length = 7)

parallelStop()
parallelStartSocket(10)



for (i in 1:(n_clasif-1)){
  
  lrn.raw <- clasif.list[[i]]
  lrn.under = makeUndersampleWrapper(lrn.raw, usw.cl = "0")
  lrn.over = makeOversampleWrapper(lrn.raw, osw.cl = "1")
  lrn.smote = makeSMOTEWrapper(lrn.raw, sw.nn = 5)
  
  if (i == 1){
    lrn.cw = makeWeightedClassesWrapper(lrn.raw,  wcw.param = "class.weights")}else{
      lrn.cw = makeWeightedClassesWrapper(lrn.raw)}
  
  lrn.list <- list(lrn.raw, lrn.cw, lrn.under, lrn.over,  lrn.smote)
  
  
  for (k in 1:5){
    
    lrn.filter = makeFilterWrapper(lrn.list[[k]], fw.method = "praznik_JMI")
    #lrn.filter = makeFilterWrapper(lrn.list[[k]],  fw.method = "randomForest_importance")
    #lrn.filter = makeFilterWrapper(lrn.list[[k]], fw.method = "E-Borda",
    #fw.base.methods = c("FSelector_relief", "auc","praznik_JMI", 
    #"randomForest_importance")) 
    
    
    #lrn.filter = makeFeatSelWrapper(lrn.list[[k]], resampling = rdesc,
    #control = makeFeatSelControlRandom(maxit = 5), show.info = FALSE)
    
    
    
    lrn.filter = setPredictType(lrn.filter, "prob")
    
    res = tuneParamsMultiCrit(lrn.filter, task = tr.task, control = ctrl,
                              resampling = rdesc, par.set = ps.list[[i]][[k]], show.info = FALSE,
                              #measures = list(bac, f1, costs.measure))
                              measures = list(fpr, fnr, costs.measure))
    #measures = list(fpr, fnr))
    
    meas.vector <- vector()
    meas.model <- list()
    meas.pred <- list()
    
    for (h in 1:length(res$x)){
      
      lrn_tune = setHyperPars(lrn.filter, par.vals = res$x[[h]])
      
      lrn_model <- train(learner = lrn_tune, task = tr.task)
      
      #opt.feat  <- lrn_model$learner.model$opt.result$x
      #task.feat <- getTaskFeatureNames(tr.task)
      #opt.tr.task <- dropFeatures(tr.task, task.feat[-match(opt.feat,task.feat)])
      #lrn_model <- train(learner = lrn_tune, task = opt.tr.task)                           
      
      meas.model[[h]] <-  lrn_model
      
      lrn_pred <- predict(lrn_model, ts.task)
      tune.res = tuneThreshold(pred = lrn_pred, measure = costs.measure)
      #tune.res = tuneThreshold(pred = lrn_pred, measure = bac)
      lrn_pred <- setThreshold(lrn_pred, tune.res$th)
      
      meas.pred[[h]] <- lrn_pred
      meas.vector[h] <- performance(lrn_pred, costs.measure)
      #meas.vector[h] <- performance(lrn_pred, bac)
    }
    
    meas.best <- which.min(meas.vector)
    
    # lrn_tune = setHyperPars(lrn.filter, par.vals = res$x[[meas.best]])
    # lrn_model <- train(learner = lrn_tune, task = tr.task)
    # lrn_pred <- predict(lrn_model, ts.task)
    # 
    # tune.res = tuneThreshold(pred = lrn_pred, measure = costs.measure)
    # #tune.res = tuneThreshold(pred = lrn_pred, measure = bac)
    # 
    # lrn_pred <- setThreshold(lrn_pred, tune.res$th)
    
    
    res.list[[i]][[k]] <- res
    model.list[[i]][[k]] <- meas.model[[meas.best]]
    pred.list[[i]][[k]] <- meas.pred[[meas.best]]
    
  }
  
  
  EVAL.TABLE.LIST <- list()
  
  for (k in 1:5){
    
    CM <-  caret::confusionMatrix (data = pred.list[[i]][[k]]$data$response, 
                                   reference = pred.list[[i]][[k]]$data$truth, positive = "1", mode = "prec_recall")
    
    CM.LIST[[i]][[k]] <- CM
    
    BAC <- CM$byClass["Balanced Accuracy"]
    PRE <- CM$byClass["Precision"]
    REC <- CM$byClass["Recall"]
    F1 <-  CM$byClass["F1"]
    AUC <- performance(pred.list[[i]][[k]], measures = list(auc))
    C <- performance(pred.list[[i]][[k]], measures = list(costs.measure))
    
    
    EVAL.TABLE.LIST[[k]] <- data.frame(bac = BAC, pr = PRE, rec = REC, f1 = F1, auc = AUC, c = C)
    print(paste0("tune: ", i, k))
    print(CM$table)
    
    
  }
  
  chart.eval <- rbind(EVAL.TABLE.LIST[[1]], EVAL.TABLE.LIST[[2]], EVAL.TABLE.LIST[[3]],
                      EVAL.TABLE.LIST[[4]], EVAL.TABLE.LIST[[5]])
  
  rownames(chart.eval)   <- c("base", "cw", "over", "under", "smote")
  
  chart.eval.list[[i]] <- chart.eval
  
  print(chart.eval.list[[i]])
  
}


#RES_80 <- list(res.list, model.list, pred.list, CM.LIST, chart.eval.list) 

#save(RES_80, file = "res80.RData")
#load("res80.RData")

#getFilteredFeatures(RES_80[[2]][[1]][[3]])

### ENSEMBLE
# library(mlr)
# library(caret)

load("cost.measure.list.RData")
load("tasks.list.RData")

setwd("C:/Users/Juan/Desktop/PABLO/doctorado/articulo complejidad/capitulo_redes")
load("res80SONIA.RData")
res_sonia <- RES_80
load("res80NATI.RData")
res_nati <- RES_80
load("res80JAKO.RData")
res_jako <- RES_80
load("res80MIRKA.RData")
res_mirka <- RES_80

RES_CLASIF <- list(res_sonia, res_nati, res_jako, res_mirka)

STACK.LIST <- list()
FEATURES.LIST <- list()
ELECCION <- list()

clasif.names <- c("ensemble" ,"svm", "rf", "rpart", "xgb", "logreg", "gbm")


for (i in 2:4){
  
  mejores.clasif  <- unlist(lapply(RES_CLASIF[[i]][[5]], function(x){which.max(x$bac)}))
  mejores.modelos  <- sapply(1:6, function(x){RES_CLASIF[[i]][[2]][[x]][[mejores.clasif[x]]]})
  mejores.features  <- sapply(1:6, function(x){getFilteredFeatures(RES_CLASIF[[i]][[2]][[x]][[mejores.clasif[x]]])})
  mejores.performance  <- t(sapply(1:6, function(x){RES_CLASIF[[i]][[5]][[x]][mejores.clasif[x],]}))
  mejores.nombres <- sapply(1:6, function(x){row.names(RES_CLASIF[[i]][[5]][[x]][mejores.clasif[x],])})
  base.learners = apply(mejores.modelos,2,function(x){x$learner})
  STACK = makeStackedLearner(base.learner = base.learners, 
                             super.learner = NULL, method = "hill.climb", predict.type = "prob")   
  fitted = mlr::train(STACK, tasks.list[[1]][[i]])
  pred.stack = predict(fitted, tasks.list[[2]][[i]])
  tune.pred.th = tuneThreshold(pred = pred.stack, measure = cost.measure.list[[i]])
  pred.stack.th <- setThreshold(pred.stack, tune.pred.th$th)
  
  
  CM <- caret::confusionMatrix (data = pred.stack.th$data$response, 
                                reference = pred.stack.th$data$truth, 
                                positive = "1", mode = "prec_recall")
  
  BAC <- CM$byClass["Balanced Accuracy"]
  PRE <- CM$byClass["Precision"]
  REC <- CM$byClass["Recall"]
  F1 <-  CM$byClass["F1"]
  AUC <- performance(pred.stack.th , measures = list(auc))
  C <- performance(pred.stack.th , measures = list(costs.measure))
  
  STACK.LIST[[i]] <- data.frame(bac = BAC, pr = PRE, rec = REC, f1 = F1, auc = AUC, c = C)
  FEATURES.LIST[[i]] <- mejores.features
  
  res.data.frame <- as.data.frame(rbind(STACK.LIST[[i]], mejores.performance))
  res.data.frame <- cbind(clasif = clasif.names, tipo = c("", mejores.nombres), res.data.frame)
  
  rownames(res.data.frame) <- as.character(1:7)
  
  res.data.frame <- res.data.frame[order(unlist(res.data.frame$bac), decreasing =T),]
  
  ELECCION[[i]] <- res.data.frame
  
}

RES_TODOS <- list(ELECCION, STACK.LIST, FEATURES.LIST)
#save(RES_TODOS, file = "RES_TODOS.RData")

#load("RES_TODOS.RData")

#################################################################
## Clustering por mixtura de bernoulli, LASSO mixto  y GLMM.

library(flexmix)
library(ggplot2)
library(dplyr)

load("lp_list.RData")

mix_data <- rbind(lp_list[[1]][,1:73],lp_list[[2]][,1:73],lp_list[[3]][,1:73],lp_list[[4]][,1:73])  
mix_data[is.na(mix_data)] <- 0 
col_out <- match(c("from_n","to_n","from","to","est3"),
                 colnames(mix_data))
mix_data <- mix_data[,-col_out]
id <- c(rep("s",dim(lp_list[[1]])[1]), rep("n",dim(lp_list[[2]])[1]), rep("j",dim(lp_list[[3]])[1]), rep("m",dim(lp_list[[4]])[1]))
id_sesion <- apply(data.frame(id,as.character(mix_data$sesion)),1,function(x){paste(x, collapse = "")})
sesion <- mix_data$sesion
mix_data <- mix_data[,-which(colnames(mix_data)=="sesion")]

mix_data <- mix_data %>% mutate_if(is.numeric, scale)
mix_data$target <-  as.numeric(as.character(mix_data$target))

id_sesion_cum <- tapply(rep(1,length(id_sesion)), id_sesion, cumsum)
TIME <- as.vector(unlist(id_sesion_cum[unique(id_sesion)]))

mix_data <- data.frame(mix_data,id_sesion, sesion, TIME, id)

mix_data$esp1 <- ifelse(mix_data$esp1==2, 1, mix_data$esp1)
mix_data$esp1 <- factor(mix_data$esp1)

mix_data$esp2 <- ifelse(mix_data$esp2==2, 1, mix_data$esp2)
mix_data$esp2 <- factor(mix_data$esp2)


betaMixFix <- stepFlexmix(target ~ 1 | id_sesion,
                          model = FLXMCmvbinary(),
                          k = 2, nrep = 5, data = mix_data,
                          control = list(tolerance = 1e-15, iter.max = 1000, minprior = 0.005))


media.cluster <- round(c(as.numeric(betaMixFix@components$Comp.1[[1]]@parameters),
                         as.numeric(betaMixFix@components$Comp.2[[1]]@parameters)), digits = 3)


medias_sesiones <- tapply(mix_data$target, mix_data$id_sesion, mean)

mix_clusters <- betaMixFix@cluster

sesiones_1 <- unique(mix_data[which(mix_clusters==1),"id_sesion"])
sesiones_2 <- unique(mix_data[which(mix_clusters==2),"id_sesion"])

plot.1  <- ggplot(data.frame(x=sesiones_1, y=medias_sesiones[sesiones_1]), aes(x = x, y = y)) +
  geom_point(color="green") +
  geom_hline(yintercept=media.cluster[1], color = "red", linetype = "dashed") +           
  labs(x = "sesiones", y = "media") + theme_bw() + ggtitle("Cluster 1") +
  annotate(geom="text", x=3, y=0.1, color="red",
           label=paste0("mu: ",media.cluster[1]))  

plot.2  <- ggplot(data.frame(x=sesiones_2, y=medias_sesiones[sesiones_2]), aes(x = x, y = y)) +
  geom_point(color="green")  +
  geom_hline(yintercept=media.cluster[2], color = "red", linetype = "dashed") +           
  labs(x = "sesiones", y = "media") + theme_bw() + ggtitle("Cluster 2") +
  annotate(geom="text", x=3, y=0.1, color="red",
           label=paste0("mu: ",media.cluster[2])) 


## Figura 1

png(paste0("Clusters_mix.png"))
cowplot::plot_grid(plot.1, plot.2,  nrow = 2, ncol = 1)
dev.off()


which_target <- which(colnames(mix_data)=="target")
which_id_sesion <- which(colnames(mix_data)=="id_sesion")
which_sesion <- which(colnames(mix_data)=="sesion")
which_time <- which(colnames(mix_data)=="TIME")

mix_1  <- mix_data[which(mix_clusters==1),] 
mix_2  <- mix_data[which(mix_clusters==2),] 
mix_1$sesion <- factor(mix_1$sesion)
mix_2$sesion <- factor(mix_2$sesion)

mix_1$es <- as.factor(as.character(mix_1$es))
mix_1$esp1 <- as.factor(as.character(mix_1$esp1))
mix_1$esp2 <- as.factor(as.character(mix_1$esp2))
mix_1$mod <- as.factor(as.character(mix_1$mod))
mix_1$morf <- as.factor(as.character(mix_1$morf))


mix_2$es <- as.factor(as.character(mix_2$es))
mix_2$esp1 <- as.factor(as.character(mix_2$esp1))
mix_2$esp2 <- as.factor(as.character(mix_2$esp2))
mix_2$mod <- as.factor(as.character(mix_2$mod))
mix_2$morf <- as.factor(as.character(mix_2$morf))


costantes_1 <- which(apply(mix_1,2,sd)==0)
costantes_2 <- which(apply(mix_2,2,sd)==0)

factores <-  colnames(mix_data)[9:26]      
factores.fla <- sapply(1:length(factores), function(x){paste0("as.factor(",factores[x],")")})
fla.nombres <- colnames(mix_data)
fla.nombres[9:26] <- factores.fla 

formula_1 <- as.formula(paste0("target ~ 1 + ", paste(fla.nombres[-c(which_target,which_sesion,which_id_sesion, which_time,costantes_1)], collapse = " + ")))
formula_2 <- as.formula(paste0("target ~ 1 + ", paste(fla.nombres[-c(which_target,which_sesion,which_id_sesion, which_time,costantes_2)], collapse = " + ")))

## CV para lambda y seleccion de variables (cuadro 13)

# las variables seleccionadas están en: LASSO_seleccion.RData

load("LASSO_seleccion.RData")
LASSO_seleccion[[1]]
LASSO_seleccion[[2]]

##### LASSO mixto

source("cv.glmmLasso-master/R/calc.lossFunction.r ")
source("cv.glmmLasso-master/R/computeLambda.R ")
source("cv.glmmLasso-master/R/cv.glmmLasso.R ")
source("cv.glmmLasso-master/R/glmmLasso_MultLambdas.r ")

lambda <- seq(0,500,by=5)

mod.cv.1 <- cv.glmmLasso(fix = formula_1, rnd = list(id_sesion=~1), data = mix_1, 
                         family = binomial(link = "logit"), kfold = 5, lambdas = lambda,
                         nlambdas = NULL, lambda.final = 'lambda.1se', loss = calc_deviance)

mod.cv.2 <- cv.glmmLasso(fix = formula_2, rnd = list(id_sesion=~1), data = mix_2, 
                         family = binomial(link = "logit"), kfold = 5, lambdas = lambda,
                         nlambdas = NULL, lambda.final = 'lambda.1se', loss = calc_deviance)

sel.1 <- sort(abs(mod.cv.1$glmmLasso.final$coefficients), decreasing = T)

var.sel.1 <- unique(ifelse(substr(names(sel.1), start = 1, stop = 9)=="as.factor", 
                           substr(names(sel.1), start = 11, stop = nchar(names(sel.1))-2), names(sel.1)))
var.sel.1 <- var.sel.1[-which(var.sel.1=="(Intercept)")] 


sel.2 <- sort(abs(mod.cv.2$glmmLasso.final$coefficients), decreasing = T)

var.sel.2 <- unique(ifelse(substr(names(sel.2), start = 1, stop = 9)=="as.factor", 
                           substr(names(sel.2), start = 11, stop = nchar(names(sel.2))-2), names(sel.2)))
var.sel.2 <- var.sel.2[-which(var.sel.2=="(Intercept)")]


# Lambda y selección de primeras doce variables

# mod.cv.1$lambda.min; mod.cv.2$lambda.min

# var.sel.1[1:12]
# var.sel.2[1:12]

# LASSO_seleccion <- list()
# LASSO_seleccion[[1]] <-  var.sel.1[1:12]
# LASSO_seleccion[[2]] <-  var.sel.2[1:12]
# save(LASSO_seleccion, file = "LASSO_seleccion.RData")




## GLMM en cada cluster con variables seleccionadas.

library(glmmTMB) 
library(MuMIn)

## lo siguiente hace la selección de modelos. El resultado está
## en el archivo: models.dredge.RData  


m.2 <- glmmTMB(target ~ esp1 + morf + mod + est5 + lda + est2 + C + anim + est6 + esp2 + Skew + stem +
                 (1|id_sesion), family = binomial(link="logit"), data = mix_1, verbose = FALSE, 
               control = glmmTMBControl(optCtrl=list(iter.max=1e4,eval.max=1e4),collect=TRUE))

modelos.m.2  <- dredge(m.2, rank = "AIC")


m.3 <- glmmTMB(target ~ es + esp1 + morf + est6 + fam.lex + lda + stem + est4 + est1 + est7 + mod + esp2 +
                 (1|id_sesion), family = binomial(link="logit"), data = mix_2, verbose = FALSE, 
               control = glmmTMBControl(optCtrl=list(iter.max=1e4,eval.max=1e4),collect=TRUE))

modelos.m.3  <- dredge(m.3, rank = "AIC")


models.dredge  <- list(modelos.m.2, modelos.m.3)

#save(models.dredge,  file = "models.dredge.RData")

#####

load("models.dredge.RData")  
modelos.m.2  <- models.dredge[[1]] 
modelos.m.3  <- models.dredge[[2]] 


importance.m.2 <- importance(modelos.m.2)

# coeficientes promediados:
conf.av.m.2 <- summary(model.avg(subset(modelos.m.2, 1/8 < weight/max(modelos.m.2$weight))))
# mejor modelo:
modelos.m.2.best <- get.models(modelos.m.2, subset = 1)[[1]] 
summary(modelos.m.2.best)

importance.m.3 <- importance(modelos.m.3)

# coeficientes promediados:
conf.av.m.3 <- summary(model.avg(subset(modelos.m.3, 1/8 < weight/max(modelos.m.3$weight))))
# mejor modelo:
modelos.m.3.best <- get.models(modelos.m.3, subset = 1)[[1]] 
summary(modelos.m.3.best)



### figura de Odds Ratio


png("OR_cluster1.png")
sjPlot::plot_model(modelos.m.2.best, type = "est",  show.values = T)
dev.off()

png("OR_cluster2.png")
sjPlot::plot_model(modelos.m.3.best, type = "est",  show.values = T)
dev.off()

# LATEX:

xtable::xtable(conf.av.m.2$coefmat.full, digits = 4, 
               caption = "Cluster 1: Coeficientes promediados.")
xtable::xtable(conf.av.m.3$coefmat.full, digits = 4, 
               caption = "Cluster 2: Coeficientes promediados.")

#t1 <- broom.mixed::tidy(modelos.m.3.best, conf.int = TRUE)

# comparaciones multiples

library(multcomp)

glht_glmmTMB <- function (model, ..., component="cond") {
  glht(model, ...,
       coef. = function(x) fixef(x)[[component]],
       vcov. = function(x) vcov(x)[[component]],
       df = NULL)
}
modelparm.glmmTMB <- function (model, coef. = function(x) fixef(x)[[component]],
                               vcov. = function(x) vcov(x)[[component]],
                               df = NULL, component="cond", ...) {
  multcomp:::modelparm.default(model, coef. = coef., vcov. = vcov.,
                               df = df, ...)
}


g1_esp1 <- glht(modelos.m.2.best, linfct = mcp(esp1 = "Tukey"))
g2_esp1 <- glht(modelos.m.3.best, linfct = mcp(esp1 = "Tukey"))

g2_esp2 <- glht(modelos.m.3.best, linfct = mcp(esp2 = "Tukey"))

g1_mod <- glht(modelos.m.2.best, linfct = mcp(mod = "Tukey"))
g2_mod <- glht(modelos.m.3.best, linfct = mcp(mod = "Tukey"))


# Figuras de comparaciones multiples

plot(g1_esp1)
plot(g2_esp1)
plot(g2_esp2)
plot(g1_mod)
plot(g2_mod)


# LATEX:

xtable::xtable(broom::tidy(summary(g1_esp1)), 
               caption = "Cluster 1. Comparaciones múltiples: ESP1")

xtable::xtable(broom::tidy(summary(g1_mod)), 
               caption = "Cluster 1. Comparaciones múltiples: MOD")

xtable::xtable(broom::tidy(summary(g2_esp1)), 
               caption = "Cluster 2. Comparaciones múltiples: ESP1")

xtable::xtable(broom::tidy(summary(g2_esp2)), 
               caption = "Cluster 2. Comparaciones múltiples: ESP2")

xtable::xtable(broom::tidy(summary(g2_mod)), 
               caption = "Cluster 2. Comparaciones múltiples: MOD")


