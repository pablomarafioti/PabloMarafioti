 ## DINAMICA DE SERIES DE ERROR


load("datos_articulo.RData")

datos <- datos_articulo

###PERMUTATION ENTROPY (PE)


#devtools::install_version('statcomp', '0.0.1.1000')
library(statcomp)

# plano entropia-complejidad por alumno


Signal <- round(as.numeric(as.character(datos$RES_BIN)) + runif(nrow(datos),0.01, 0.09),digits = 4)
df.Signal <- data.frame(Signal = Signal, ID = datos$ID, SESION = datos$ID.SESION,
                        B = as.numeric(as.character(datos$RES_CAT)))

DATA.PE <- list()
DATA.C <- list()
SESION.INDEX  <- unique(datos$ID.SESION)
D = 2

for(i in 1:52){
  GLOBAL  <- global_complexity(df.Signal[df.Signal$SESION==SESION.INDEX[i],"B"], ndemb = D)
  
  DATA.PE[[i]] <- GLOBAL[[1]]
  DATA.C[[i]] <- GLOBAL[[2]]
}

freq.abs.sesion <- as.matrix(tapply(as.numeric(as.character(datos$RES_BIN)),
                                    datos$ID.SESION, function(x){sum(x)}))
freq.abs.index <- match(as.character(SESION.INDEX), rownames(freq.abs.sesion))
freq_abs_sesion <- round(freq.abs.sesion[freq.abs.index], digits = 2)


freq.rel.sesion <- as.matrix(tapply(as.numeric(as.character(datos$RES_BIN)),
                                    datos$ID.SESION, function(x){sum(x)/length(x)}))
freq.rel.index <- match(as.character(SESION.INDEX), rownames(freq.rel.sesion))
freq_rel_sesion <- round(freq.rel.sesion[freq.rel.index], digits = 2)

n.sesion <- as.matrix(tapply(as.numeric(as.character(datos$RES_BIN)),
                             datos$ID.SESION, function(x){length(x)}))
n.index <- match(as.character(SESION.INDEX), rownames(n.sesion))
n_sesion <- n.sesion[n.index]

cuts  <- match(c("1.1","1.12","2.1","2.14","3.1","3.14","4.1","4.12"),as.character(SESION.INDEX))
lineas_corte <- cumsum(n_sesion)


knoise.complexity = sapply(X= n_sesion, FUN=function(x){global_complexity(x=powernoise(k=0.5, N=x)[[1]], ndemb=2)})

source("t2p.r")
source("perm_2samples.r")
source("pan_test.r")
source("permsign.r")


e.cor <- vector()
c.cor <- vector()
e.pan <- vector()
c.pan <- vector()
e.perm <- vector()
c.perm <- vector()

g = 0
for (i in c(1,3,5,7)){
  g  = g+1
  df.test <- data.frame(x= (knoise.complexity["PE",])[cuts[i]:cuts[i+1]], 
                        y= (knoise.complexity["MPR_Cjs",])[cuts[i]:cuts[i+1]],
                        e=unlist(DATA.PE)[cuts[i]:cuts[i+1]],
                        c=unlist(DATA.C)[cuts[i]:cuts[i+1]]) 
  e.cor[g] <- cor.test(df.test$x,df.test$e, method="spearman", alternative="two.sided", exact=F)$p.value
  c.cor[g] <- cor.test(df.test$y,df.test$c, method="spearman", alternative="two.sided", exact =F)$p.value
  e.pan[g] <- pan(df.test$x,df.test$e,alt="two.sided",B=10000)$p.value
  c.pan[g] <- pan(df.test$y,df.test$c,alt="two.sided",B=10000)$p.value
  
  e.perm[g]   <- perm.2samples(cbind(rep(1:2,each=length(df.test$e)),c(df.test$x, df.test$e)), alt="two.sided",B=10000)$p.value
  c.perm[g]   <- perm.2samples(cbind(rep(1:2,each=length(df.test$c)),c(df.test$y, df.test$c)),alt="two.sided",B=10000)$p.value
}

data.test <- data.frame(id = 1:4, e.cor = e.cor, e.perm = e.perm, 
                        c.cor = c.cor ,c.perm = c.perm)
data.test <- apply(data.test,2,function(x){round(x, digits=3)})
data.test <- apply(data.test,2,function(x){ifelse(x<0.001,"<0.001",x)})


xtable.test <- xtable::xtable(data.test, 
                              caption = "Test de permutación: e.cor= p valor de S de Spearman para PE; 
e.perm= p valor del test para PE; c.cor= p valor de S de Spearman para C;
e.perm= p valor del test para C",
                              label = "xtable:example")

data.PE.C <- data.frame(PE.Sonia = c(round(unlist(DATA.PE)[cuts[1]:cuts[2]],digits = 3),"",""),
                        C.Sonia = c(round(unlist(DATA.C)[cuts[1]:cuts[2]], digits = 3),"",""),
                        PE.Nati = round(unlist(DATA.PE)[cuts[3]:cuts[4]], digits = 3),
                        C.Nati = round(unlist(DATA.C)[cuts[3]:cuts[4]], digits = 3),
                        PE.Jako = round(unlist(DATA.PE)[cuts[5]:cuts[6]], digits = 3),
                        C.Jako = round(unlist(DATA.C)[cuts[5]:cuts[6]], digits = 3),
                        PE.Mirka = c(round(unlist(DATA.PE)[cuts[7]:cuts[8]], digits = 3),"",""),
                        C.Mirka = c(round(unlist(DATA.C)[cuts[7]:cuts[8]], digits = 3),"",""))



xtable.PE.C <- xtable::xtable(data.PE.C, 
                              caption = "Medidas de Entropía de Permutación (PE) y Complejidad (C) para cada aprendiente.",
                              label = "xtable:example")


df.data.pe <- data.frame(PE = unlist(DATA.PE), 
                         ID = c(rep(1,12), rep(2,14), rep(3,14), rep(4,12)),
                         SESION = c(1:12,1:14,1:14,1:12))

at.no_error <- df.data.pe[df.data.pe$PE <=0.55, c("ID","SESION")]
at.error <- df.data.pe[df.data.pe$PE >=0.80, c("ID","SESION")]
trans <- df.data.pe[-as.numeric(c(rownames(at.no_error), rownames(at.error))),
                    c("ID","SESION")]

at.sonia <-  c(paste(as.character(at.no_error[at.no_error$ID==1,"SESION"]),collapse=","),
               paste(as.character(trans[trans$ID==1,"SESION"]), collapse=","),
               paste(as.character(at.error[at.error$ID==1,"SESION"]), collapse=",")
)

at.nati <-  c(paste(as.character(at.no_error[at.no_error$ID==2,"SESION"]),collapse=","),
              paste(as.character(trans[trans$ID==2,"SESION"]), collapse=","),
              paste(as.character(at.error[at.error$ID==2,"SESION"]), collapse=",")
)

at.jako <-  c(paste(as.character(at.no_error[at.no_error$ID==3,"SESION"]),collapse=","),
              paste(as.character(trans[trans$ID==3,"SESION"]), collapse=","),
              paste(as.character(at.error[at.error$ID==3,"SESION"]), collapse=",")
)

at.mirka <-  c(paste(as.character(at.no_error[at.no_error$ID==4,"SESION"]),collapse=","),
               paste(as.character(trans[trans$ID==4,"SESION"]), collapse=","),
               paste(as.character(at.error[at.error$ID==4,"SESION"]), collapse=",")
)

df.at <- data.frame(rbind(at.sonia, at.nati, at.jako, at.mirka))
row.names(df.at) <- c("SONIA","NATI","JAKO","MIRKA")
colnames(df.at) <- c("NO_ERROR", "TRANSICION", "ERROR")

xtable.at <- xtable::xtable(df.at, 
                            caption = "Sesiones en cada región del plano entropía-complejidad",
                            label = "xtable:example")



# graficas

library(ggplot2)
library(ggrepel)
library(gridExtra)


min_d2 <- limit_curves(2, fun="min")
max_d2 <- limit_curves(2, fun="max")

mind2 <- data.frame(x = min_d2[[1]],
                    y = min_d2[[2]])
maxd2 <- data.frame(x = max_d2[[1]],
                    y = max_d2[[2]])

#"C:/Users/Juan/Desktop/PABLO/doctorado/articulo complejidad/capitulo_introduccion/publicacion_IV"

for (i in c(1,3,5,7)) {
  
  df.plot <- data.frame(x = unlist(DATA.PE[cuts[i]:cuts[i+1]]), 
                        y = unlist(DATA.C[cuts[i]:cuts[i+1]]))
  png(paste0("E_C_plot",i,".png"))
  
  print(
    ggplot(df.plot, aes(x=x, y=y))+ geom_point(color = "red") + 
      geom_line(aes(x=x, y=y), data=mind2) + geom_line(aes(x=x, y=y), data=maxd2) +
      geom_text_repel(aes(label=rownames(df.plot)), size = 3) +
      scale_x_continuous(name = "PE",
                         limits = c(0, 1),
                         breaks = seq(0, 1, by = 0.1)) +
      scale_y_continuous(name = "C",
                         limits = c(0,0.3)) +
      ggtitle("Mapa Entropia-Complejidad") 
  )
  dev.off()
}

df.plot.t <- data.frame(x = unlist(DATA.PE), y = unlist(DATA.C),
                        Frecuencia_Error = ifelse(freq_rel_sesion<0.2,"0-0.2", 
                                                  ifelse(freq_rel_sesion>0.21 & freq_rel_sesion<0.35,"0.21-0.35","> 0.35"))
)


png(paste0("E_C_plot_todos.png"))
print(
  ggplot(df.plot.t, aes(x=x, y=y)) + geom_point(aes(color=Frecuencia_Error)) + 
    geom_line(aes(x=x, y=y), data=mind2) + geom_line(aes(x=x, y=y), data=maxd2) +
    geom_text_repel(aes(label=as.character(freq_rel_sesion), size = 3)) +
    scale_x_continuous(name = "PE",
                       limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.1)) +
    scale_y_continuous(name = "C",
                       limits = c(0, 0.3)) +
    ggtitle("Mapa Entropia-Complejidad")                     
)
dev.off()

## Entropias con ventana movil : deteccion de cambio de dinamica

#install.packages("robts", repos=c("http://R-Forge.R-project.org", "http://CRAN.R-project.org"))

library(robts)

DATA.SIGNAL <- list()
ID.INDEX <- 1:4
v = 80
D = 3
for(i in 1:4){
  NPE  <-   zoo::rollapply(df.Signal[df.Signal$ID==ID.INDEX[i],"B"],v,
                           function(x){permutation_entropy(ordinal_pattern_distribution(x, ndemb = D))},
                           fill=NA, align="right")
  DATA.SIGNAL[[i]] <- NPE[-which(is.na(NPE))]
}


results <- vector()
for (i in 1:4){
  results[i] <- changerob(DATA.SIGNAL[[i]], plot = F)$estimate[[1]]
}


v=80
l = 1 
for (i in 1:4){
  df.plot.2 <-  data.frame(y = c(rep(NA,80), DATA.SIGNAL[[i]]))
  df.plot.2$x <- 1:dim(df.plot.2)[1]
  lineas <- lineas_corte[cuts[l]:cuts[l+1]]
  if (i==1) {vlinea=lineas} else {vlinea=lineas-lineas_corte[cuts[l-1]]}
  
  png(paste0("PE_plot",i,".png"))
  
  print(ggplot(df.plot.2, aes(y=y, x=x)) + geom_line(color="blue") +
          scale_y_continuous(name = "EP",
                             limits = round(range(df.plot.2$y[-c(1:80)]), digits=2),
                             breaks = seq(0, 1, by = 0.1)) +
          scale_x_continuous(name = "Sesiones",
                             breaks =  vlinea) + 
          geom_vline(xintercept=results[i]+v, color = "red") +
          geom_vline(xintercept= vlinea, linetype="dashed", color = "black") +
          ggtitle("EP - Ventana = 80")                     
  )
  dev.off()
  
  l = l+2
} 


### recurrence plots

#library(devtools)
#packageurl <- "http://cran.r-project.org/src/contrib/Archive/fractal/fractal_2.0-4.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")
#install_github("FredHasselman/casnet", build_vignettes = FALSE)


#library(crqa)


rqa.C1 <- as.numeric(as.character(datos$C_1))

rqa.C2 <- ifelse(as.numeric(as.character(datos$C_2))!=0,
                 as.numeric(as.character(datos$C_2))+0.1,0)

rqa.df <- data.frame(id = datos$ID, id.sesion = datos$ID.SESION, 
                     rqa.C1, rqa.C2)


library(casnet)

# ejemplos para sesiones 2 y 12 de SONIA

RP_ejemplo <- list()

j=1

rqa.data <- rqa.df[rqa.df$id==j,]
rownames(rqa.data) <- 1:dim(rqa.data)[1]
SESION.INDEX.RQA  <- unique(rqa.data$id.sesion)

z = 0
for (k in c(2,12)){
  z = z + 1
  array1 = rqa.data[rqa.data$id.sesion==SESION.INDEX.RQA[k],"rqa.C1"]
  array2 = rqa.data[rqa.data$id.sesion==SESION.INDEX.RQA[k],"rqa.C2"]
  
  RP_ejemplo[[z]] <- casnet::rp(y1 = array1, y2 =  array2, emRad = 0.0001, 
                                VLmin = l_min, DLmin = l_min, HLmin = l_min)
} 

g_1 <- rp_plot(RP_ejemplo[[1]], plotDimensions = TRUE, returnOnlyObject = TRUE, 
               title = "SONIA, SESION 2")

g_2 <- rp_plot(RP_ejemplo[[2]], plotDimensions = TRUE, returnOnlyObject = TRUE, 
               title = "SONIA, SESION 12")

cowplot::plot_grid(g_1, g_2)


#
# Medidas para cada sesion por alumno 
#
#
# las siguientes lineas comentadas calculan las medidas de recurrencia,
# por sesión y alumno. El resultado se guarda en RQA.V y RQA.H. 
# Basta cargar esos datos para generar las gráficas del artículo.

# RQA.V <-list(list(), list(), list(), list())
# 
# l_min = 3
# 
# path_out = "C:/Users/Juan/Desktop/PABLO/casnet"
# 
# for (j in 1:4){
#   rqa.data <- rqa.df[rqa.df$id==j,]
#   rownames(rqa.data) <- 1:dim(rqa.data)[1]
#   SESION.INDEX.RQA  <- unique(rqa.data$id.sesion)
#   rqa.res <- list()
#   for (k in 1:length(SESION.INDEX.RQA)){
#     array1 = rqa.data[rqa.data$id.sesion==SESION.INDEX.RQA[k],"rqa.C1"]
#     array2 = rqa.data[rqa.data$id.sesion==SESION.INDEX.RQA[k],"rqa.C2"] 
#     
#     RQA.V[[j]][[k]]  <- rp_cl_main(data = data.frame(y1 = array1, y2 =  array2), emRad = 0.0001, 
#                         VLmin = l_min, DLmin = l_min, returnLineDist = T, path_out = path_out)
#     
#     
#   }
# }
# 
# detach(package:casnet)
# library(casnet)
# 
# RQA.H <-list(list(), list(), list(), list())
# 
# 
# for (j in 1:4){
#   rqa.data <- rqa.df[rqa.df$id==j,]
#   rownames(rqa.data) <- 1:dim(rqa.data)[1]
#   SESION.INDEX.RQA  <- unique(rqa.data$id.sesion)
#   rqa.res <- list()
#   for (k in 1:length(SESION.INDEX.RQA)){
#     array1 = rqa.data[rqa.data$id.sesion==SESION.INDEX.RQA[k],"rqa.C1"]
#     array2 = rqa.data[rqa.data$id.sesion==SESION.INDEX.RQA[k],"rqa.C2"] 
#     
#     RQA.H[[j]][[k]]  <- rp_cl_main(data = data.frame(y1 = array2, y2 =  array1), emRad = 0.0001, 
#                         VLmin = l_min, DLmin = l_min, returnLineDist = T, path_out = path_out)
#     
#   }
# }

#save(RQA.V, file="RQA.V.RData")
#save(RQA.H, file="RQA.H.RData")


# cargar los resultados de las medidas de recurrencia

load("RQA.V.RData")
load("RQA.H.RData")


RR <- list()
DET <- list()
LL <- list()
LAM_v <- list()
LAM_h <- list()
TT_v <- list()
TT_h <- list()
#ANI <- list()
#ENTR_v <- list()
#ENTR_h <- list()


for (j in 1:4){
  
  RR[[j]] <- sapply(1:length(RQA.V[[j]]), function(x){RQA.V[[j]][[x]]$measures$RR})
  DET[[j]] <- sapply(1:length(RQA.V[[j]]), function(x){RQA.V[[j]][[x]]$measures$DET})
  LL[[j]] <- sapply(1:length(RQA.V[[j]]), function(x){RQA.V[[j]][[x]]$measures$L})
  LAM_v[[j]] <- sapply(1:length(RQA.V[[j]]), function(x){RQA.V[[j]][[x]]$measures$LAM})
  LAM_h[[j]] <- sapply(1:length(RQA.V[[j]]), function(x){RQA.H[[j]][[x]]$measures$LAM})
  TT_v[[j]] <- sapply(1:length(RQA.V[[j]]), function(x){RQA.V[[j]][[x]]$measures$TT})
  TT_h[[j]] <- sapply(1:length(RQA.V[[j]]), function(x){RQA.H[[j]][[x]]$measures$TT})
  #ANI[[j]] <- sapply(1:length(RQA.H[[j]]), function(x){RP_m[[j]][[x]]$ANI})
  #ENTR_v[[j]] <- sapply(1:length(RP_m[[j]]), function(x){RP_m[[j]][[x]]$ENT_vl})
  #ENTR_h[[j]] <- sapply(1:length(RP_m[[j]]), function(x){RP_m[[j]][[x]]$ENT_hl})
}

all_dist_d <- list()
all_dist_d_bind <- list()
all_pl_d <- list()

all_dist_h <- list()
all_dist_h_bind <- list()
all_pl_h <- list()

all_dist_v <- list()
all_dist_v_bind <- list()
all_pl_v <- list()

for (j in 1:4){
  
  all_dist_d[[j]] <- data.frame(sapply(1:length(RQA.V[[j]]), function(x){RQA.V[[j]][[x]]$diag_disthist}))
  all_dist_d_bind[[j]]   <- data.table::rbindlist(all_dist_d[[j]])
  all_pl_d[[j]] <- tapply(all_dist_d_bind[[j]]$V2, factor(all_dist_d_bind[[j]]$V1), sum)   
  
  all_dist_h[[j]] <- data.frame(sapply(1:length(RQA.H[[j]]), function(x){RQA.H[[j]][[x]]$hori_disthist}))
  all_dist_h_bind[[j]]   <- data.table::rbindlist(all_dist_h[[j]])
  all_pl_h[[j]] <- tapply(all_dist_h_bind[[j]]$V2, factor(all_dist_h_bind[[j]]$V1), sum)   
  
  all_dist_v[[j]] <- data.frame(sapply(1:length(RQA.V[[j]]), function(x){RQA.V[[j]][[x]]$hori_disthist}))
  all_dist_v_bind[[j]]   <- data.table::rbindlist(all_dist_v[[j]])
  all_pl_v[[j]] <- tapply(all_dist_v_bind[[j]]$V2, factor(all_dist_v_bind[[j]]$V1), sum)   
  
}


n_boot_d <- vector() 
n_boot_h <- vector()
n_boot_v <- vector() 

boot.d <- list(
  list(list(), list(), list()), 
  list(list(), list(), list()),
  list(list(), list(), list()),
  list(list(), list(), list())
)

boot.h <- list(
  list(list(), list(), list()), 
  list(list(),list(), list()), 
  list(list(), list(), list()),
  list(list(), list(), list())
  
)

boot.v <- list(
  list(list(), list(), list()), 
  list(list(),list(), list()), 
  list(list(), list(), list()),
  list(list(), list(), list())
  
)



ci.boot.d <- list(list(), list(), list(), list())
ci.boot.h <- list(list(), list(), list(), list())
ci.boot.v <- list(list(), list(), list(), list())

for (j in 1:4){
  
  n_boot_d[j] <- round((1/length(RQA.V[[j]]))*(sum(all_pl_d[[j]])))
  n_boot_h[j] <- round((1/length(RQA.H[[j]]))*(sum(all_pl_h[[j]])))
  n_boot_v[j] <- round((1/length(RQA.V[[j]]))*(sum(all_pl_v[[j]])))
  
  for (i in 1:2000){
    
    sample.l.d <- sample(x=as.numeric(names(all_pl_d[[j]])), 
                         size=n_boot_d[j], replace = TRUE, prob = prop.table(all_pl_d[[j]]))
    sample.l.h <- sample(x=as.numeric(names(all_pl_h[[j]])), 
                         size=n_boot_h[j], replace = TRUE, prob = prop.table(all_pl_h[[j]]))    
    sample.l.v <- sample(x=as.numeric(names(all_pl_v[[j]])), 
                         size=n_boot_v[j], replace = TRUE, prob = prop.table(all_pl_v[[j]]))      
    
    boot.d[[j]][[1]][[i]] <- sum(sample.l.d[which(sample.l.d>=l_min)])/
      sum(sample.l.d)
    boot.d[[j]][[2]][[i]] <- mean(sample.l.d[which(sample.l.d>=l_min)])
    
    boot.h[[j]][[1]][[i]] <- sum(sample.l.h[which(sample.l.h>=l_min)])/
      sum(sample.l.h)
    boot.h[[j]][[2]][[i]] <- mean(sample.l.h[which(sample.l.h>=l_min)])
    
    boot.v[[j]][[1]][[i]] <- sum(sample.l.v[which(sample.l.v>=l_min)])/
      sum(sample.l.v)
    boot.v[[j]][[2]][[i]] <- mean(sample.l.v[which(sample.l.v>=l_min)])
    
  }
  
  alfa = 0.05
  
  ci.boot.d[[j]][[1]] <- quantile(unlist(boot.d[[j]][[1]]), c(alfa/2, 1-(alfa/2)))
  ci.boot.d[[j]][[2]] <- quantile(unlist(boot.d[[j]][[2]]), c(alfa/2, 1-(alfa/2)))
  ci.boot.h[[j]][[1]] <- quantile(unlist(boot.h[[j]][[1]]), c(alfa/2, 1-(alfa/2)))
  ci.boot.h[[j]][[2]] <- quantile(unlist(boot.h[[j]][[2]]), c(alfa/2, 1-(alfa/2)))
  ci.boot.v[[j]][[1]] <- quantile(unlist(boot.v[[j]][[1]]), c(alfa/2, 1-(alfa/2)))
  ci.boot.v[[j]][[2]] <- quantile(unlist(boot.v[[j]][[2]]), c(alfa/2, 1-(alfa/2)))
  
}


#library(ggplot2)


OBS <- list(DET, LL, LAM_h, TT_h  , LAM_v, TT_v)

est <- c("DET" , "D", "LAM_h", "TT_h" , "LAM_v", "TT_v")

for (k in 1:6){
  
  for (j in 1:4){
    
    if (k == 1) {BOUND = ci.boot.d[[j]][[1]][2]}
    if (k == 2) {BOUND = ci.boot.d[[j]][[2]][2]}
    if (k == 3) {BOUND = ci.boot.h[[j]][[1]][2]}
    if (k == 4) {BOUND = ci.boot.h[[j]][[2]][2]}
    if (k == 5) {BOUND = ci.boot.v[[j]][[1]][2]}
    if (k == 6) {BOUND = ci.boot.v[[j]][[2]][2]}   
    
    #print(paste0(est[k],"_",j,":",BOUND))
    
    png(paste0(est[k],"_",j,".png"))
    
    print(
      
      ggplot(data.frame(x=1:length(OBS[[k]][[j]]), y=OBS[[k]][[j]]), aes(x = x, y = y)) +
        geom_point(color="green") + geom_line(color = "green") +
        #      geom_hline(yintercept=ci.boot[1], color = "red", linetype = "dashed") +
        
        geom_hline(yintercept=BOUND, color = "red", linetype = "dashed") + 
        
        #labs(x="Sesión", y=est[k]) +
        scale_x_continuous(name = "Sesiones", breaks = 1:length(OBS[[k]][[j]]))  +
        scale_y_continuous(name = est[k]) +
        ggtitle(est[k])  + theme_bw() +
        geom_text(aes(1,as.numeric(round(BOUND,digits = 2)),
                      label = as.numeric(round(BOUND,digits = 2)), 
                      vjust = -1))
      
    )
    
    dev.off()
    
  }
  
}


# discretizacion y suma de medidas de recurrencia.

DET_D <- list()
LAM_V <- list()
LAM_H <- list()
LL_D <- list()
TT_V <- list()
TT_H <- list()


for (j in 1:4){
  DET_D[[j]] <- ifelse(DET[[j]]>=ci.boot.d[[j]][[1]][2],1,0)
  LL_D[[j]] <- ifelse(LL[[j]]>=ci.boot.d[[j]][[2]][2],1,0)
  LAM_H[[j]] <- ifelse(LAM_h[[j]]>=ci.boot.h[[j]][[1]][2],1,0)
  TT_H[[j]] <- ifelse(TT_h[[j]]>=ci.boot.h[[j]][[2]][2],1,0)
  LAM_V[[j]] <- ifelse(LAM_v[[j]]>=ci.boot.v[[j]][[1]][2],1,0)
  TT_V[[j]] <- ifelse(TT_v[[j]]>=ci.boot.v[[j]][[2]][2],1,0)
}


medidas.suma <- list()

for (j in 1:4){
  
  medidas.suma[[j]]  <- apply(
    rbind(DET_D[[j]], LAM_V[[j]], LAM_H[[j]], LL_D[[j]], TT_V[[j]], TT_H[[j]]),2,
    sum)
  
}


alumnos <- c("SONIA","NATI", "JAKO", "MIRKA")

medidas.rqa <- list()
medidas.rqa.names <- list()
medidas.rqa.suma <- list()

for (j in 1:4){
  
  medidas.rqa.names[[j]] <- c(paste0("Alumno: ", alumnos[[j]]) ,"DET", "D",
                              "LAM vertical","LAM horizontal ", "TT vertical",
                              "TT horizontal" ) 
  
  medidas.rqa.suma[[j]] <- paste(as.character(medidas.suma[[j]]), collapse = " ")
  
  medidas.rqa[[j]] <- c(" ",
                        paste(as.character(DET_D[[j]]), collapse = " "),
                        paste(as.character(LL_D[[j]]), collapse = " "),
                        paste(as.character(LAM_V[[j]]), collapse = " "),
                        paste(as.character(LAM_H[[j]]), collapse = " "),
                        paste(as.character(TT_V[[j]]), collapse = " "),
                        paste(as.character(TT_H[[j]]), collapse = " ")
  )
}

df.medidas  <- rbind(
  data.frame(medidas = medidas.rqa.names[[1]], binario =  medidas.rqa[[1]]),
  data.frame(medidas = "SUMA: ", binario =  medidas.rqa.suma[[1]]),
  data.frame(medidas = medidas.rqa.names[[2]], binario =  medidas.rqa[[2]]),
  data.frame(medidas = "SUMA: ", binario =  medidas.rqa.suma[[2]]),
  data.frame(medidas = medidas.rqa.names[[3]], binario =  medidas.rqa[[3]]),
  data.frame(medidas = "SUMA: ", binario =  medidas.rqa.suma[[3]]),
  data.frame(medidas = medidas.rqa.names[[4]], binario =  medidas.rqa[[4]]),
  data.frame(medidas = "SUMA: ", binario =  medidas.rqa.suma[[4]])
)


SESION.INDEX  <- unique(datos$ID.SESION)
freq.rel.sesion <- as.matrix(tapply(as.numeric(as.character(datos$RES_BIN)),
                                    datos$ID.SESION, function(x){sum(x)/length(x)}))
freq.rel.index <- match(as.character(SESION.INDEX), rownames(freq.rel.sesion))
freq_rel_sesion <- round(freq.rel.sesion[freq.rel.index], digits = 2)

df.medidas.suma    <- data.frame(id.sesion = SESION.INDEX, 
                                 id = c(rep(1,12), rep(2,14), rep(3,14), rep(4,12)),
                                 freq = freq_rel_sesion, 
                                 medidas = c(medidas.suma[[1]], medidas.suma[[2]], 
                                             medidas.suma[[3]], medidas.suma[[4]])
)


df.medidas

print(xtable::xtable(df.medidas, 
                     caption = "Discretización y suma de medidas de recurrencia",
                     label = "xtable:example"), include.rownames = FALSE)






