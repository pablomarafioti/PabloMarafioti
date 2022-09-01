

## Código en R del articulo: Regímenes latentes de error 
##    en el aprendizaje de la concordancia plural en ELE


####### MOTIVOS

load("DATA.RData")


library(jmotif)
library(zoo)
library(xtable)


n <- 3
p <- 2
s <- 3


DATA.LIST <- list()

for (i in 1:51){
  
  RESPONSE   <- ts(as.numeric(as.character(DATA.DF[DATA.DF$ID.SESION==unique(DATA.DF$ID.SESION)[i],"B"])))
  
  
  SAX <- rollapply(RESPONSE,n,function(x){series_to_string(paa(x, s),p)},fill=NA, align="right")
  
  STRINGE <- rollapply(RESPONSE,n,function(x){paste(as.character(x),collapse="")},fill=NA, align="right")
  
  DATA.LIST[[i]] <- data.frame(STRINGE,SAX, DATA.DF[DATA.DF$ID.SESION==unique(DATA.DF$ID.SESION)[i],])
  DATA.LIST[[i]] <- (DATA.LIST[[i]])[-c(1:(n-1)),] 
}

data.motif <- DATA.LIST[[1]]
for (i in 2:51){data.motif <- rbind(data.motif, DATA.LIST[[i]])}


MOTIVOS <- tapply(data.motif$STRINGE, factor(data.motif$SAX), function(x){unique(x)})


xtable.motif  <- xtable(data.frame(motif = c(names(MOTIVOS)), 
                                   secuencia = unlist(sapply(1:length(names(MOTIVOS)), function(x){paste(as.character(MOTIVOS[[x]]), collapse=" ")}))), 
                        caption = "Motivos y secuencias que representan",label = "xtable:example")

names(MOTIVOS)  


library(quanteda)
#library(lsa)


BAG <- unlist(tapply(data.motif$SAX, factor(data.motif$ID.SESION), function(x){paste(x, collapse=" ")}))

BAG <- BAG[as.character(unique(data.motif$ID.SESION))]

dtm <- dfm(as.vector(BAG), verbose = TRUE) 
#dtm <- dfm_trim(dtm, min_termfreq = 2)
#dtm <- dfm_weight(dtm, scheme = "logcount")
rownames(dtm) <- unique(data.motif$ID.SESION)


dtm.tfidf <- dfm_tfidf(dtm, scheme_tf="count", scheme_df = "inverse", 
                       k=0, smoothing = 1, base=2) # pasar de freq a tfidf 

library(text2vec)

set.seed(1234)
n_topics  = 3


G <- expand.grid(alfa = seq(1 / n_topics, 51 / n_topics,1), beta = seq(0.01, 1 / n_topics, 0.01))

freq.rel.sesion <- as.matrix(tapply(as.numeric(as.character(datos$RES_BIN)),
                                    datos$ID.SESION, function(x){sum(x)/length(x)}))
freq.sesion <- as.matrix(tapply(as.numeric(as.character(datos$RES_BIN)),
                                datos$ID.SESION, function(x){sum(x)}))
n.sesion <- as.matrix(tapply(as.numeric(as.character(datos$RES_BIN)),
                             datos$ID.SESION, function(x){length(x)}))

DFs <- list()
KAPPA <- vector() 

for (i in 1:nrow(G)){
set.seed(1234)
lda_model = LDA$new(n_topics =  n_topics, doc_topic_prior =  G[i,1], topic_word_prior =  G[i,2])
doc_topic_distr = 
  lda_model$fit_transform(as(dtm.tfidf, "dgCMatrix"), n_iter = 1000, 
                          convergence_tol = 0.001, n_check_convergence = 25, 
                          progressbar = FALSE)

topics <- unlist(apply(doc_topic_distr,1, function(x){which(x==max(x))[1]}))

v  <- ifelse(freq.rel.sesion<0.21, 1, 
                        ifelse(freq.rel.sesion>0.21 & freq.rel.sesion<0.35,2,3))
v <- v[-1,]

C <- merge(x = data.frame(t = topics, r = names(topics)), 
           y = data.frame(v = v, r = names(v)),
           by = "r"
            )

conf.matrix <- caret::confusionMatrix(factor(C$t), factor(C$v))
  
KAPPA[[i]] <- conf.matrix$overall[2]

  }

ganador <- which.max(KAPPA) # 279, 
KAPPA[ganador] # kappa = 0.8120393
G[ganador,] 
#        alfa beta
#278 6.333333 0.17


# modelo elegido

set.seed(1234)
n_topics  = 3

lda_model = LDA$new(n_topics =  n_topics, doc_topic_prior =  G[ganador,1], topic_word_prior =  G[ganador,2])
doc_topic_distr = 
  lda_model$fit_transform(as(dtm.tfidf, "dgCMatrix"), n_iter = 1000, 
                          convergence_tol = 0.001, n_check_convergence = 25, 
                          progressbar = FALSE)

###### tablas con resultados
res_lda <- list(lda_model, doc_topic_distr)
save(res_lda, file ="res_lda.RData")

load("res_lda.RData")
doc_topic_distr <- res_lda[[2]]

topics <- unlist(apply(doc_topic_distr,1, function(x){which(x==max(x))[1]}))
topics <- data.frame(topics)
rownames(topics) <- unique(data.motif$ID.SESION)

df.rel.par   <- merge(data.frame(freq.rel.sesion, freq.sesion, n.sesion), topics, by=0)
colnames(df.rel.par) <- c("sesiones","freq.rel","freq.abs","n","topics")

for (i in 1:3){
  print(range(df.rel.par[df.rel.par$topics==i,2]))
}


df.rel.par$v  <- ifelse(df.rel.par$freq.rel<0.21,1, 
                        ifelse(df.rel.par$freq.rel>0.21 & df.rel.par$freq.rel<0.35,2,3))

df.topicos <-  df.rel.par[match(as.character(unique(data.motif$ID.SESION)),as.character(df.rel.par$sesiones)),]
colnames(df.topicos) <- c("sesiones","freq.rel","freq.abs","n","topicos","ref.")
df.topicos$sesiones <- as.character(df.topicos$sesiones)

xtable.topicos  <- xtable(df.topicos,  
                          caption = "Clasificación de tópicos. sesiones: sesión por cada alumno; 
freq.rel: frecuencia relativa de error; freq. abs: frecuencia absoluta de error;
n: número de instancias por sesión; tópicos: tópico asignado a cada sesión;
ref.: tópico de referencia",
                          label = "xtable:example")


conf.matrix <- caret::confusionMatrix(factor(df.rel.par$topics), factor(df.rel.par$v))


xtable.conf  <- xtable(conf.matrix$table, 
                       caption = "Matriz de confusión. Filas: predicción, Columnas: referencia",
                       label = "xtable:example")



xtable.conf.med  <- xtable(t(conf.matrix$byClass[,c(1:4,7,11)]), 
                           caption = "Matriz de confusión. Medidas de evaluación de la clasificación",
                           label = "xtable:example")

# medida AUC
caTools::colAUC(df.rel.par$topics, df.rel.par$v, plot=F)


# Graficos

library(LDAvis)
library(servr)

lda_model$plot()


######### Markov Chain con covariables

head(data.motif)

Y = data.motif$STRINGE
Y = ifelse(Y == "000", 0,  
           ifelse(Y == "100", 1,
                  ifelse(Y == "010",2, 
                         ifelse(Y == "001", 3, 
                                ifelse(Y == "111", 4,
                                       ifelse(Y == "110",5,
                                              ifelse(Y == "011", 6, 7))))))
           )

Names <- c("000", "100","010", "001", "111", "110", "011", "101")

data.long <- data.motif[,-c(1:3,7:9,20)]
data.long$Y <- Y

##################################

library(depmixS4)

fit_modelos <- list()
proba_emision <- list()
matriz_M0 <- list()
matriz_M1 <- list()
n_estados = 3


for (i in 1:4){
print(i)
set.seed(1)
  
mod <- depmix(Y ~ 1, data = data.long[data.long$ID==i,], nstates = n_estados, family = multinomial("identity"),
       transition = ~ LDA + EP + GRAM + MOD + IMA_CONC + FAM_LEX + AN + EST1 + EST2 + EST5 + FREQ_S, 
       ntimes = nrow(data.long[data.long$ID==i,]))
#fm <- fit(mod, verbose = FALSE,  em=em.control(maxit=2000, tol=1e-08, crit="relative"))
#fm <- fit(mod, em=em.control(maxit=2000))
#summary(fm, which = "transition")

# ajustar el modelo 20 veces y calcular BIC
fit_mod <- list()
ll <- rep(NA,20)
for(j in 1:20){
  
  fit_mod[[j]] <- fit(mod, verbose = FALSE,  em=em.control(maxit=2000))
  ll[j] <- BIC(fit_mod[[j]])
}

# el modelo con BIC más bajo
fit.mod <- fit_mod[[which.min(ll)]]

# matriz de probabilidades de emisión

a<-lapply(1:length(fit.mod@response[[1]]), function(y) {
  sapply(fit.mod@response, function(x) x[[y]]@parameters$coefficients)
})

#A <- apply(a[[1]], 2, function(x){exp(x)/sum(exp(x))})
A <- a[[1]]
row.names(A) <- Names
colnames(A) <- c("E1", "E2", "E3")
round(A,4)

# estados por individuo según algoritmo viterbi
#pst <- posterior(fit.mod)
#colMeans(post[,2:4])


# Matriz de transición entre estados latentes cuando las covariables valen cero


 eta1 <- as.matrix(attributes(fit.mod@transition[[1]])$parameters$coefficients)
 eta2 <- as.matrix(attributes(fit.mod@transition[[2]])$parameters$coefficients)
 eta3 <- as.matrix(attributes(fit.mod@transition[[3]])$parameters$coefficients)
 
  TM_0  <- matrix(c(exp(eta1[1,])/sum(exp(eta1[1,])),
                  exp(eta2[1,])/sum(exp(eta2[1,])),
                  exp(eta3[1,])/sum(exp(eta3[1,]))), 3, 3, byrow = T) 
  stateNames <- c("E1", "E2", "E3")
  row.names(TM_0) <- stateNames; colnames(TM_0) <- stateNames
  round(TM_0, 2)
  rowSums(TM_0)
  
  # Matriz de transición entre estados latentes cuando las covariables valen uno
  
  
  TM_1  <- matrix(c(exp( c(sum(eta1[,1]),  sum(eta1[,2]), sum(eta1[,3])))/ sum(exp( c(sum(eta1[,1]),  sum(eta1[,2]), sum(eta1[,3])))),
                    exp( c(sum(eta2[,1]),  sum(eta2[,2]), sum(eta2[,3])))/ sum(exp( c(sum(eta2[,1]),  sum(eta2[,2]), sum(eta2[,3])))),
                    exp( c(sum(eta3[,1]),  sum(eta3[,2]), sum(eta3[,3])))/ sum(exp( c(sum(eta3[,1]),  sum(eta3[,2]), sum(eta3[,3]))))
                    ), 3, 3, byrow = T) 
  row.names(TM_1) <- stateNames; colnames(TM_1) <- stateNames 
  round(TM_1, 2)
  rowSums(TM_1)

fit_modelos[[i]] <- fit.mod
proba_emision[[i]] <- A
matriz_M0[[i]] <- TM_0
matriz_M1[[i]] <- TM_1 
  
}

resultados <- list()
resultados[[1]] <- list(fit_modelos[[1]], proba_emision[[1]], matriz_M0[[1]], matriz_M1[[1]]) # SONIA
resultados[[2]] <- list(fit_modelos[[2]], proba_emision[[2]], matriz_M0[[2]], matriz_M1[[2]]) # NATI
resultados[[3]] <- list(fit_modelos[[3]], proba_emision[[3]], matriz_M0[[3]], matriz_M1[[3]]) # JAKO
resultados[[4]] <- list(fit_modelos[[4]], proba_emision[[4]], matriz_M0[[4]], matriz_M1[[4]]) # MIRKA
 
save(resultados, file = "resultados.RData")

# Matriz de transición entre estados latentes cuando cada covariable vale 1 y las demás cero

load("resultados.RData")

m_m_m <- list()  
M_M <- list()
A <- c(2,1,2,3) # estado oculto con menos errores 


for (i in 1:4) {
  m_m <- list()
  
  C_C <- colnames(as.matrix(attributes(resultados[[i]][[1]]@transition[[1]])$x))
  C_C[1] <- "Intercept"
  n_n <- length(C_C)
  X_X <- diag(n_n)
  X_X[,1] <- rep(1, n_n)
  
  
  Eta1 <- as.matrix(attributes(resultados[[i]][[1]]@transition[[1]])$parameters$coefficients)
  Eta2 <- as.matrix(attributes(resultados[[i]][[1]]@transition[[2]])$parameters$coefficients)
  Eta3 <- as.matrix(attributes(resultados[[i]][[1]]@transition[[3]])$parameters$coefficients)
  
  for (j in 1:n_n){
    
    x_B_1_1 <-   t(matrix(X_X[j,])) %*%  matrix(Eta1[,1])
    x_B_1_2 <-   t(matrix(X_X[j,])) %*%  matrix(Eta1[,2])
    x_B_1_3 <-   t(matrix(X_X[j,])) %*%  matrix(Eta1[,3])
    
    x_B_2_1 <-   t(matrix(X_X[j,])) %*%  matrix(Eta2[,1])
    x_B_2_2 <-   t(matrix(X_X[j,])) %*%  matrix(Eta2[,2])
    x_B_2_3 <-   t(matrix(X_X[j,])) %*%  matrix(Eta2[,3])
    
    x_B_3_1 <-   t(matrix(X_X[j,])) %*%  matrix(Eta3[,1])
    x_B_3_2 <-   t(matrix(X_X[j,])) %*%  matrix(Eta3[,2])
    x_B_3_3 <-   t(matrix(X_X[j,])) %*%  matrix(Eta3[,3])
    
    
    
    T_M  <-   matrix(c(exp( c(x_B_1_1,  x_B_1_2, x_B_1_3))/sum(exp( c(x_B_1_1,   x_B_1_2, x_B_1_3))),
                      exp( c(x_B_2_1,  x_B_2_2, x_B_2_3))/ sum( exp( c(x_B_2_1,  x_B_2_2, x_B_2_3))),
                      exp( c(x_B_3_1,  x_B_3_2, x_B_3_3))/ sum( exp( c(x_B_3_1,  x_B_3_2, x_B_3_3)))
    ), 3, 3, byrow = T) 
    stateNames <- c("E1", "E2", "E3")
    row.names(T_M) <- stateNames; colnames(T_M) <- stateNames
    
    m_m[[j]] <- T_M

        
  }
  m_m_m[[i]] <- m_m
  M  <- matrix(unlist(lapply(m_m_m[[i]], function(x){x[A[[i]],]})),length(m_m_m[[i]]), 3, byrow=T)
  row.names(M) <- C_C; colnames(M) <- stateNames
  M_M[[i]] <- M
  }



res_MM <- list(m_m_m, M_M)

save(res_MM, file = "res_MM.RData")

lapply(m_m_m[[1]], function(x){round(x,2)})

# Tablas con probabilidades de de emisión 

load("resultados.RData")


suj <- c("SONIA", "NATI", "JAKO", "MIRKA")

for (i in 1:4) {
  
  df1 <- data.frame(E = (resultados[[i]][[2]])[,1], S = row.names(resultados[[i]][[2]]))
  df2 <- data.frame(E = (resultados[[i]][[2]])[,2], S = row.names(resultados[[i]][[2]]))
  df3 <- data.frame(E = (resultados[[i]][[2]])[,3], S = row.names(resultados[[i]][[2]]))
 
  ta <- data.frame(df1[,1], df2[,1], df3[,1])  
  colnames(ta) <- c("E1", "E2", "E3") 
  rownames(ta) <- df1[,2]
  
  print(xtable::xtable(ta, digits = 4, caption = paste(suj[i], ": Probabilidades de transición")))
}




# Tabla con betas estimados

for (i in 1:4) {

Eta1 <- as.matrix(attributes(resultados[[i]][[1]]@transition[[1]])$parameters$coefficients)
Eta2 <- as.matrix(attributes(resultados[[i]][[1]]@transition[[2]])$parameters$coefficients)
Eta3 <- as.matrix(attributes(resultados[[i]][[1]]@transition[[3]])$parameters$coefficients)


print(xtable::xtable(Eta1, digits = 4, caption = paste(suj[i], ": Betas estimados (modelo desde estado 1)")))
print(xtable::xtable(Eta2, digits = 4, caption = paste(suj[i], ": Betas estimados (modelo desde estado 2)")))
print(xtable::xtable(Eta3, digits = 4, caption = paste(suj[i], ": Betas estimados (modelo desde estado 3)")))
}


# Tabla con probabilidades de transición (para cada predictora)

load("res_MM.RData")

tabla_mat <- list()

for (i in 1:4){
  Matrices <- res_MM[[1]][[i]][-1]
  col <- rownames(res_MM[[2]][[i]])[-1]
  
  tabla_m <- matrix(unlist(lapply(Matrices, function(x) {as.vector(t(x))})), length(col), 9, byrow = T)
  rownames(tabla_m) <- col
  colnames(tabla_m) <- c("1a1", "1a2", "1a3", "2a1", "2a2", "2a3", 
                         "3a1", "3a2", "3a3")
  tabla_mat[[i]] <- round(tabla_m,2)
  write.csv(tabla_mat[[i]], paste(i,"_tablas_proba_trans",".csv"), row.names = T)
}
  



# graficos

load("resultados.RData")


# graficar matriz de transición (cuando todas las covariables valen cero)  
  library("qgraph")
  

for (i in 1:4){

 
  
GgrTM_0 <- qgraph(resultados[[i]][[3]], edge.labels = TRUE, edge.color = "black", trans = 1, 
                    mode = "direct", vsize = 10, label.cex = 2,
                    border.color = "black", mar = c(10,10,10,10), label.color = "black",
                    edge.label.cex = 1.5, repulsion = 0.1, asize = 2, curveAll = T) 
  
qgraph(GgrTM_0, filename = paste0("TM_0_", i) , filetype = "png", height = 5, width = 5)

  # graficar matriz de transición  (cuando todas las covariables valen 1)
  #library("qgraph")
  GgrTM_1 <- qgraph(resultados[[i]][[4]], edge.labels = TRUE, edge.color = "black", trans = 1, 
                    mode = "direct", vsize = 10, label.cex = 2,
                    border.color = "black", mar = c(10,10,10,10), label.color = "black",
                    edge.label.cex = 1.5, repulsion = 0.1, asize = 2, curveAll = T)
  
  qgraph(GgrTM_1, filename = paste0("TM_1_", i) , filetype = "png", height = 5, width = 5) 

  }    

# graficar matrices de transición  (cuando CADA covariable valen 1 y las demás valen 0)
#library("qgraph")

load("res_MM.RData")


for (i in 1:4){
  Matrices <- res_MM[[i]][[i]][-1]
  col <- rownames(res_MM[[2]][[i]])[-1]
  layout(matrix(c(1:15),5,3, byrow = T))  
   
  for (j in 1:length(Matrices)){

GgrTM_2 <- qgraph(Matrices[[j]], edge.labels = TRUE, edge.color = "black", trans = 1, 
                  mode = "direct", vsize = 10, label.cex = 2, title = col[j],
                  border.color = "black", mar = c(10,10,10,10), label.color = "black",
                  edge.label.cex = 1.5, repulsion = 0.1, asize = 2, curveAll = T)

   
QG <- qgraph(c(GgrTM_2, GgrTM_2_2), filename = paste0("TM_2_", i) , filetype = "png", height = 5, width = 5) 




}    
}




  
# graficar probabilidades de emision  
library(gridExtra)
library(ggplot2)

for (i in 1:4) {

df1 <- data.frame(E = (resultados[[i]][[2]])[,1], S = row.names(resultados[[i]][[2]]))
df2 <- data.frame(E = (resultados[[i]][[2]])[,2], S = row.names(resultados[[i]][[2]]))
df3 <- data.frame(E = (resultados[[i]][[2]])[,3], S = row.names(resultados[[i]][[2]]))

P1   <- ggplot(df1, aes(x = S, y = E)) + geom_bar(stat = "identity") + 
ylab("Probabilidad") + xlab("Símbolo") + ggtitle("Estado 1")
P2   <- ggplot(df2, aes(x = S, y = E)) + geom_bar(stat = "identity") + 
  ylab("Probabilidad") + xlab("Símbolo") + ggtitle("Estado 2")
P3   <- ggplot(df3, aes(x = S, y = E)) + geom_bar(stat = "identity") + 
  ylab("Probabilidad") + xlab("Símbolo") + ggtitle("Estado 3")

P <- grid.arrange(P1, P2, P3)

ggsave(filename = paste0("emision_",i,".png"), plot = P,  width = 5,
       height = 5)

}

# graficar Betas 

for (i in 1:4){

Eta1 <- as.matrix(attributes(resultados[[i]][[1]]@transition[[1]])$parameters$coefficients)
Eta2 <- as.matrix(attributes(resultados[[i]][[1]]@transition[[2]])$parameters$coefficients)
Eta3 <- as.matrix(attributes(resultados[[i]][[1]]@transition[[3]])$parameters$coefficients)

ETA <- list(Eta1, Eta2, Eta3)

for (j in 1:3) {
   for (k in 2:3) {

dat <- data.frame(Index = 1:nrow(ETA[[j]]), label = rownames(ETA[[j]]), B = ETA[[j]][,k])

inf <- round(min(dat$B))
sup <- round(max(dat$B))
each <- round(length(inf:sup)/10)

plot1 <- ggplot(dat, aes(y = Index, x = B)) +
  geom_point(shape = 18, size = 5) +  
  #geom_errorbarh(aes(xmin = LL, xmax = UL), height = 0.25) +
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed", cex = 1, alpha = 0.5) +
  scale_y_continuous(name = "", breaks=1:length(dat$label), labels = dat$label, trans = "reverse") +
  scale_x_continuous(name = "Betas", breaks= seq(inf, sup, each)) +
  #xlab("Betas (escala del logit)") + 
  ggtitle(paste0("E",j, " a ", "E",k)) +
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x.bottom = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"))


ggsave(filename = paste0("OR_",i,"_estado_",j, "_a_estado_",k,".png"), plot = plot1,  width = 8,
       height = 6)

}
}
}

# simulación

# matrices de transición para CADA concordancia


mmm <- list()  

for (i in 1:4) {
mm <- list()
nn <- resultados[[i]][[1]]@ntimes
XX <- as.matrix(attributes(resultados[[i]][[1]]@transition[[1]])$x)
colnames(XX)[1] <- "Intercept"

Eta1 <- as.matrix(attributes(resultados[[i]][[1]]@transition[[1]])$parameters$coefficients)
Eta2 <- as.matrix(attributes(resultados[[i]][[1]]@transition[[2]])$parameters$coefficients)
Eta3 <- as.matrix(attributes(resultados[[i]][[1]]@transition[[3]])$parameters$coefficients)

for (j in 1:nn){

  xB_1_1 <-   t(matrix(XX[j,])) %*%  matrix(Eta1[,1])
  xB_1_2 <-   t(matrix(XX[j,])) %*%  matrix(Eta1[,2])
  xB_1_3 <-   t(matrix(XX[j,])) %*%  matrix(Eta1[,3])
  
  xB_2_1 <-   t(matrix(XX[j,])) %*%  matrix(Eta2[,1])
  xB_2_2 <-   t(matrix(XX[j,])) %*%  matrix(Eta2[,2])
  xB_2_3 <-   t(matrix(XX[j,])) %*%  matrix(Eta2[,3])
  
  xB_3_1 <-   t(matrix(XX[j,])) %*%  matrix(Eta3[,1])
  xB_3_2 <-   t(matrix(XX[j,])) %*%  matrix(Eta3[,2])
  xB_3_3 <-   t(matrix(XX[j,])) %*%  matrix(Eta3[,3])
  
  
  
  TM  <-   matrix(c(exp( c(xB_1_1,  xB_1_2, xB_1_3))/ sum(exp( c(xB_1_1,  xB_1_2, xB_1_3))),
                    exp( c(xB_2_1,  xB_2_2, xB_2_3))/ sum(exp( c(xB_2_1,  xB_2_2, xB_2_3))),
                    exp( c(xB_3_1,  xB_3_2, xB_3_3))/ sum(exp( c(xB_3_1,  xB_3_2, xB_3_3)))
  ), 3, 3, byrow = T) 
stateNames <- c("E1", "E2", "E3")
row.names(TM) <- stateNames; colnames(TM) <- stateNames

mm[[j]] <- TM

}
mmm[[i]] <- mm
}

# simular secuencia de estados 
library(markovchain)

set.seed(1)
regMC_S <- list()
for (i in 1:4){
nn <- resultados[[i]][[1]]@ntimes
regMC <- vector()
for (j in 1:nn){
MC <- new("markovchain",transitionMatrix=mmm[[i]][[j]],
             states=c("E1","E2","E3"),
             name= "MarkovChain") 

regMC[j] <- rmarkovchain(n=1,object=MC)

}
regMC_S[[i]] <- regMC
}

# simular respuesta "Y" a partir de probas de emisión del correspondiente estado simulado

RRR <- list()
for (i in 1:4){
  nn <- resultados[[i]][[1]]@ntimes
  RR <- vector()
  for (j in 1:nn){
St   <- regMC_S[[i]][j]
pp <- (resultados[[i]][[2]])[, St]
response   <- rmultinom(n = 1, size = 1, prob = pp)
RR[j] <- rownames(response)[which(response == 1)]
}
RRR[[i]] <- RR
}

resultados.sim <- list(RRR, regMC_S)
save(resultados.sim, file = "resultados.sim.RData")

# comparar con los datos usando  Kullback-Leibler Divergence (KL)
# para dos distribuciones iguales: KL == 0

load("resultados.sim.RData")

RRR <- resultados.sim[[1]]

library(philentropy)

KL <- vector()
for (i in 1:4){
 distr_datos <-  prop.table(table(data.motif[data.motif$ID==i, "STRINGE"]))
 distr_sim   <-  prop.table(table(RRR[[i]]))
 if ( length(distr_sim) < 8) {distr_sim <- c(distr_sim,0)} else {distr_sim <- distr_sim}
 KL[i]   <- KL(rbind(distr_datos, distr_sim), unit = "log2")
}
  
round(KL,4) 
# 0.0357 0.0387 0.0886 0.0430


  # cuadro comparando las distribuciones de respuestas en datos y simulación
  P_datos_sim_1  <- rbind(prop.table(table(data.motif[data.motif$ID==1, "STRINGE"])), 
                          prop.table(table(RRR[[1]])) )
  rownames(P_datos_sim_1) <- c("datos SONIA", "sim SONIA")
  P_datos_sim_2  <- rbind(prop.table(table(data.motif[data.motif$ID==2, "STRINGE"])), 
                          prop.table(table(RRR[[2]])) )
  rownames(P_datos_sim_2) <- c("datos NATI", "sim NATI")
  P_datos_sim_3  <- rbind(prop.table(table(data.motif[data.motif$ID==3, "STRINGE"])), 
                         prop.table(table(RRR[[3]]))  )
  rownames(P_datos_sim_3) <- c("datos JAKO", "sim JAKO")
  P_datos_sim_4  <- rbind(prop.table(table(data.motif[data.motif$ID==4, "STRINGE"])), 
                          prop.table(table(RRR[[4]])) )
  rownames(P_datos_sim_4) <- c("datos MIRKA", "sim MIRKA")
  
  
  P_datos_sim <- rbind(round(P_datos_sim_1,4), round(P_datos_sim_2,4), round(P_datos_sim_3,4), round(P_datos_sim_4,4))
  
  xtable::xtable(P_datos_sim, digits = 3, caption = "Probabilidades empiricas y simuladas de cada estado observado.")
    
  
  # graficar distribuciones de datos y simuladas  
  

suj <- c("SONIA","SONIA", "NATI","NATI", "JAKO","JAKO", "MIRKA", "MIRKA")    

for (i in seq(1,8,2)){
    
df1 <- data.frame(P = c(P_datos_sim[i,], P_datos_sim[i+1,]) , 
                    S = rep(colnames(P_datos_sim),2), 
                    Clase = c(rep("datos",8), rep("simulación",8))   )
  
  
  
Plot1  <- ggplot(df1, aes(x = S, y = P, fill = Clase)) + 
  geom_bar(stat = "identity", position="dodge") +
  labs(title = suj[i])+ ylab("Probabilidad") + xlab("Símbolo") +
  theme(plot.title = element_text(size = rel(2), colour = "black")) + 
  scale_fill_manual(values = alpha(c("black", "gray"), 1)) +
  theme(axis.title.x = element_text(face="bold", size=10))   
  
  ggsave(filename = paste0("datos_sim",i,".png"), plot = Plot1,  width = 5,
         height = 5)  
  
}
  
    


