
library(igraph)
library(statnet)
#devtools::install_github("DougLuke/UserNetR")
library(UserNetR)
library(Matrix)
library(RSiena)
library(intergraph)
library(ggplot2)
library(dplyr, warn.conflicts=FALSE)
library(tidyr, warn.conflicts=FALSE)
library(ggplot2)

### Evolucion de redes.

load("datos_articulo.RData")
load("token_matrix.RData")
load("redes.RData")


# características de las redes globales

datas.cut.nets <- list(list(), list(), list(), list())

n_sesion <- c(12,14,14,12)


for ( k in 1:4){
  
  datas <- data.frame(token.matrix.list[[k]], datos[datos$ID==k, c("SESION", "C_1", "C_2", "CONTR_WHERE", "RES_BIN")])  
  datas$SESION <- as.numeric(as.character(datas$SESION))
  datas$RES_BIN <- as.numeric(as.character(datas$RES_BIN))
  datas$C1 <- ifelse(datas$C_1!=0,1,0)
  datas$C2 <- ifelse(datas$C_2!=0,1,0)
  
  for (i in 1:n_sesion[k]) {
    
    #datas.cut <- datas[datas$SESION==i,]   
    datas.cut <- datas[datas$SESION<=i,]    
    datas.cut.X1 <- ifelse(datas.cut$CONTR_WHERE==2, datas.cut$C1, datas.cut$C2) 
    datas.cut.X2 <- ifelse(datas.cut$CONTR_WHERE==2, datas.cut$C2, datas.cut$C1) 
    datas.cut$X1.C <- datas.cut.X1
    datas.cut$X2.C <- datas.cut.X2
    
    X1.COUNT <- tapply(datas.cut$X1.C, as.character(datas.cut$X1), sum)
    X2.COUNT <- tapply(datas.cut$X2.C, as.character(datas.cut$X2), sum)
    X1.X2.COUNT <- c(X1.COUNT, X2.COUNT)
    X.COUNT <- table(c(as.character(datas.cut$X1), as.character(datas.cut$X2)))
    
    #datas.cut.net <- network(datas.cut[,1:2], matrix.type="edgelist", directed = T)
    
    datas.cut.net <- graph_from_edgelist(as.matrix(datas.cut[,1:2]), directed = T)
    datas.cut.nodos  <-  V(datas.cut.net)$name
    error_count  <- X1.X2.COUNT[match(datas.cut.nodos,  names(X1.X2.COUNT))] 
    word_total  <- X.COUNT[match(datas.cut.nodos,  names(X.COUNT))] 
    
    error_measure  <- as.vector(error_count)
    V(datas.cut.net)$error_measure <- error_measure
    
    E(datas.cut.net)$weight  <- datas.cut$RES_BIN
    
    datas.cut.nets[[k]][[i]] <- datas.cut.net
    
  }
}


library(poweRlaw)

df.desc.list <- list()

for ( k in 1:4){
  
  deg_mean <- vector()
  cc_mean <- vector()
  v_n <- vector()
  e_n <- vector()
  densidad <- vector()
  islas <- vector()
  deg_max <- vector()
  asort <- vector()
  asort_deg <- vector()
  cc_mean <- vector()
  L_mean <- vector()
  S_in_mean <- vector()
  S_out_mean <- vector()
  gama <- vector()
  
  for (i in 1:n_sesion[k]) {
    
    RED  <- datas.cut.nets[[k]][[i]]
    
    dat_k <- igraph::degree(RED)
    deg_mean[i] <- mean(dat_k)
    deg_max[i] <- max(igraph::degree(RED))
    cc_mean[i] <- igraph::transitivity(RED, type = "average")
    v_n[i] <-  vcount(RED)
    e_n[i] <-  ecount(RED)
    densidad[i]  <- edge_density(RED)
    islas[i] <- sum(igraph::degree(RED) == 0)
    asort[i] <- assortativity(RED, types1 = V(RED)$error_measure, directed = F)
    asort_deg[i] <- assortativity_degree(RED, directed = F)
    L_mean[i]  <- average.path.length(RED)
    S_in_mean[i] <- mean(strength(RED, mode="in"))
    #S_out_mean[i] <- mean(strength(RED, mode="out"))
    #gama[i] <-  attributes(fit_power_law(igraph::degree(RED, mode = "all"), implementation="R.mle"))$coef
    dat_k <- dat_k[dat_k>0]
    if (i == 1 & k ==1) {gama[i] == NULL} else {
      gama[i] <- estimate_xmin(displ$new(dat_k))$pars}
  }
  
  net_names <- sapply(1: n_sesion[k] , function(x){paste0("red_", x)})
  
  df.desc   <-  data.frame(K = deg_mean, 
                           K_max = deg_max, 
                           N = v_n, 
                           E = e_n, 
                           D = densidad, 
                           #islas = islas, 
                           C = cc_mean,
                           L = L_mean,
                           S = S_in_mean,
                           #S_out_mean = S_out_mean,
                           gama = gama,
                           R_deg = asort_deg,
                           R = asort) 
  row.names(df.desc) <- net_names
  
  
  df.desc.list[[k]] <- apply(df.desc,2, function(x){round(x, digits = 3)})
  
}

df.total <- rbind(df.desc.list[[1]][12,], df.desc.list[[2]][14,],
                  df.desc.list[[3]][14,], df.desc.list[[4]][12,])

suj_names <- c("SONIA", "NATI", "JAKO", "MIRKA")

row.names(df.total) <- suj_names 


xtable::xtable(df.total, caption = "Medidas globales de las redes. 
K = grado medio; K.max = grado máximo, N = número de nodos; E = número de aristas;
D = densidad; C = coeficiente de clustering medio; L =diámetro medio;
S = strenght medio; gama = exponente de power law; 
R = coeficiente (cuantitativo) de asortatividad; 
R.deg = correlación de grado")


## distribución de grado

for (i in 1:4){
  
  g <- datas.cut.nets[[i]][[n_sesion[i]]]
  dat <- igraph::degree(g)
  dat <- dat[dat>0]
  
  data.dist <- data.frame(k=0:max(dat),p_k=degree_distribution(g))
  data.dist <- data.dist[data.dist$p_k>0,]
  
  plot.3  <- ggplot(data.dist) + geom_point(aes(x=k, y=p_k), color = "blue") + 
    theme_bw() + geom_line(aes(x=k, y=p_k), color = "blue", linetype = "dashed") +
    ggtitle(paste0(suj_names[i]," : ", "Distribución de grado"))
  
  m_pl <- displ$new(dat)
  est_pl <- estimate_xmin(m_pl)
  #est_pl$pars
  m_pl$setXmin(est_pl)
  plot.dat <- plot(m_pl, draw = F)
  fit.dat <- lines(m_pl, draw = F)
  
  plot.4 <- ggplot(plot.dat) + geom_point(aes(x=log(x), y=log(y))) + 
    labs(x="log(k)", y="log(CDF)") + theme_bw() + 
    geom_line(data=fit.dat, aes(x=log(x), y=log(y)), colour="red")  + 
    ggtitle(paste0(suj_names[i]," : ", "Power Law")) + 
    annotate(geom="text", x=3, y=-1, color="red",
             label=paste0("gamma: ", as.numeric(round(est_pl$pars,digits = 3)))) +
    annotation_logticks() 
  
  png(paste0("power_law_", suj_names[i],".png"))
  
  print(
    
    cowplot::plot_grid(plot.3, plot.4,  nrow = 1, ncol = 2)
  )
  dev.off()
  
}

## graficas de redes

# ejemplo SONIA

net0 <- asNetwork(datas.cut.nets[[1]][[3]])
nodos0 <- network.vertex.names(net0)
edge_cat <- get.edge.attribute(net0,"weight")
edge_cat[c(19,20)] <- 2 
node_size <- get.vertex.attribute(net0,"error_measure")
#node_cat <- factor(ERROR1[match(nodos0,nodos)])
linecol_pal <- c("green","yellow","red")


png("ejemplo.png")

gplot(net0, 
      gmode="digraph",thresh=0, mode='fruchtermanreingold',
      vertex.col="gray",vertex.cex=node_size+1, displaylabels = T,
      label.cex=0.9,label.pos=5, usearrows = T, interactive = F,
      displayisolates=TRUE, pad=0.4,label.col="blue", jitter = T, 
      edge.col=linecol_pal[factor(edge_cat)], edge.lwd=1)
legend("bottomleft",legend=c("0","1","2"),
       col="gray",pch=19,pt.cex=c(1,2,3), bty="n",
       title="error nodo")
legend("bottomright",legend=c("0","1","2"),
       col=linecol_pal,pch=19,pt.cex=1,5, bty="n",
       title="error enlace")

dev.off()

graphics.off()

# Redes globales

op <- par(mar = c(0,0,2,1),mfrow=c(2,2))
gplot(asNetwork(datas.cut.nets[[1]][[12]]), main="SONIA")
gplot(asNetwork(datas.cut.nets[[2]][[14]]), main="NATI")
gplot(asNetwork(datas.cut.nets[[3]][[14]]), main="JAKO")
gplot(asNetwork(datas.cut.nets[[4]][[12]]), main="MIRKA")
par(op)



# graficas de acumulación de error

s.names <- c(1,5,9)

for (i in 1:3){
  
  
  #wave <- as.network(waves.list.list[[3]][[3]])
  
  wave <- as.network(net.adj.alumnos[[1]][[i]])
  
  png(paste("sonia_graph_",i,".png"))
  
  par(mar = c(1, 1, 4, 1), oma = c(0, 0, 0, 0))
  set.seed(0)
  
  plot(wave, vertex.cex = 2.3, 
       vertex.col = "#32648E", vertex.border = "white", 
       edge.col = "red", edge.lwd = 0.01,
       label = get.vertex.attribute(wave, "vertex.names"),
       label.pos = 5, # center the label 
       label.col = "black",
       label.cex = ifelse(nchar(get.vertex.attribute(wave, "vertex.names")) > 4, 0.8, 0.9), 
       pad = -1) 
  title(main = paste0("SONIA ( sesiones ", s.names[i]," a ", s.names[i]+3, " )"), line = 1, cex.main = 0.9)
  
  dev.off()
}



s.names <- c(1,5,9)

for (i in 1:3){
  
  
  wave <- as.network(net.adj.alumnos[[2]][[i]])
  
  
  png(paste("nati_graph_",i,".png"))
  
  par(mar = c(1, 1, 4, 1), oma = c(0, 0, 0, 0))
  set.seed(0)
  
  plot(wave, vertex.cex = 2.3, 
       vertex.col = "#32648E", vertex.border = "white", 
       edge.col = "red", edge.lwd = 0.01,
       label = get.vertex.attribute(wave, "vertex.names"),
       label.pos = 5, # center the label 
       label.col = "black",
       label.cex = ifelse(nchar(get.vertex.attribute(wave, "vertex.names")) > 4, 0.8, 0.9), 
       pad = -1) 
  title(main = paste0("NATI ( sesiones hasta ", i," )"), line = 1, cex.main = 0.9)
  
  dev.off()
}


s.names <- c(1,5,9)

for (i in 1:3){
  
  
  wave <- as.network(net.adj.alumnos[[3]][[i]])
  
  
  png(paste("jako_graph_",i,".png"))
  
  par(mar = c(1, 1, 4, 1), oma = c(0, 0, 0, 0))
  set.seed(0)
  
  plot(wave, vertex.cex = 2.3, 
       vertex.col = "#32648E", vertex.border = "white", 
       edge.col = "red", edge.lwd = 0.01,
       label = get.vertex.attribute(wave, "vertex.names"),
       label.pos = 5, # center the label 
       label.col = "black",
       label.cex = ifelse(nchar(get.vertex.attribute(wave, "vertex.names")) > 4, 0.8, 0.9), 
       pad = -1) 
  title(main = paste0("JAKO ( sesiones hasta ", i," )"), line = 1, cex.main = 0.9)
  
  dev.off()
}

s.names <- c(4,8,12)

for (i in 1:3){
  
  
  wave <- as.network(net.adj.alumnos[[4]][[i]])
  
  png(paste("mirka_graph_",i,".png"))
  
  par(mar = c(1, 1, 4, 1), oma = c(0, 0, 0, 0))
  set.seed(0)
  
  plot(wave, vertex.cex = 0.8, 
       vertex.col = "#32648E", vertex.border = "white", 
       edge.col = "red", edge.lwd = 0.01,
       label = get.vertex.attribute(wave, "vertex.names"),
       label.pos = 5, # center the label 
       label.col = "black",
       label.cex = ifelse(nchar(get.vertex.attribute(wave, "vertex.names")) > 4, 0.8, 0.9), 
       pad = -1) 
  title(main = paste0("MIRKA ( sesiones hasta ", i," )"), line = 1, cex.main = 0.9)
  
  dev.off()
}

############################ Análisis estadístico temporal
### RSiena

# funciones a usar

outTable <- function(x) {
  coef <- abs(x$theta)
  coefPretty <- sprintf("%.3f", round(coef,3))
  se <- diag(x$covtheta)**.5
  sePretty <- sprintf("%.3f", round(se,3))
  stat <- coef/se
  pval <- 2*pnorm(-abs(stat))
  symp <- symnum(pval, corr = FALSE,
                 cutpoints = c(0,  .001,.01,.05, .1, 1),
                 symbols = c("***","**","*","."," "))
  convPretty <- sprintf("%.3f", round(abs(x$tconv),3))
  #out1 <- noquote(cbind(
  out1 <- data.frame(cbind(
    Function = x$effects[[1]], 
    Effect = x$effects[[2]], 
    Coef = coefPretty, 
    StEr = sePretty, 
    Pval = round(pval,4),
    Sig = symp, 
    Conv = convPretty))
  out2 <- paste("Maximum Convergence Ratio:", round(x$tconv.max,3))
  return(list(out1,out2))
}


siena07ToConvergence <- function(alg,dat,eff){
  numr <- 0
  ans <- siena07(alg, data=dat, effects=eff,
                 useCluster=TRUE, initC=TRUE, nbrNodes=8) 
  repeat {
    numr <- numr+1 
    maxt <- max(abs(ans$tconv [!eff$fix[eff$include]]))
    cat (numr , maxt , "\n") 
    if (maxt < 0.10 ) { break } 
    if (maxt > 5) { break }
    if (numr > 10) { break } 
    ans <-  siena07(alg, data=dat, effects=eff , prevAns=ans ,
                    useCluster=TRUE, initC=TRUE, nbrNodes=8)
  }
  ans  
}


testear.efecto  <- function(ans){
  
  for (i in which(ans$test)){
    sct <- score.Test(ans,i)
    cat(ans$requestedEffects$effectName[i], '\n')
    print(sct)}
}

############### preparacion 

waves <- 3

NETdata.list <- list()


scores.list <- list()

for (k in 1:4){
  
  zeros.list  <- apply(chart.list.a[[k]],1,function(x){which(x==0)})
  
  scores <- chart.list.b[[k]]
  
  for (i in 1:length(zeros.list)){scores[i, zeros.list[[i]]] <- NA}
  
  scores.list[[k]] <- scores 
  
}


# grados de nodos


grados.in <- list(list(), list(), list(), list())
grados.out <- list(list(), list(), list(), list())
grados.tot <- list(list(), list(), list(), list())

for (i in 1:4){
  for (j in 1:3){
    grados.in[[i]][[j]] <-  degree(as.network(net.adj.alumnos[[i]][[j]]), cmode="indegree")
    grados.out[[i]][[j]] <-  degree(as.network(net.adj.alumnos[[i]][[j]]), cmode="outdegree")
    grados.tot[[i]][[j]] <-  degree(as.network(net.adj.alumnos[[i]][[j]]))
  }}

# siena preparacion


redes.adj <- net.adj.alumnos


for (k in 1:4){
  
  #behabior <-  t(apply(scores.list[[k]],1, function(x){ifelse(x>0, rank(x, ties.method = "min"),x) }))
  #behabior <-  apply(scores.list[[k]],2, function(x){ifelse(x>6, 6,x) })
  
  behabior <- scores.list[[k]]
  
  nrpalabras <- dim(behabior)[1] 
  
  comp  <- apply(scores.list[[k]],1, function(x){ c(1,2,3)[-which(is.na(x))]})
  
  comp <- lapply(comp, function(x){if (length(x)==0) {x <- c(1,3)} else 
  { if (length(x)==1) {x <- c(2.99, 3)} else {x}  }})
  
  
  red <- sienaDependent(array(unlist(redes.adj[[k]]),
                              dim=c(nrpalabras, nrpalabras,waves)), allowOnly = T)
  error <-  sienaDependent(behabior, type="behavior", allowOnly = T)
  
  changes <- sienaCompositionChange(comp)
  
  f  <- match(row.names(redes.adj[[k]][[3]]), network.vertex.names(redes.list[[k]]))
  
  mod.unique <- sort(unique(atributos.v.list[[k]][f,"MOD"]))[-1]
  
  mod.dummy   <- model.matrix(~factor(atributos.v.list[[k]][f,"MOD"]))[,-1]
  
  for (i in 1:length(mod.unique)){
    
    xx <- paste0("mod_", mod.unique[i])
    
    mod_ <- coCovar(ifelse(as.vector(mod.dummy[,i])==0,0,1), centered=T)
    assign(xx, mod_, envir = .GlobalEnv)
    
  }
  
  grados_in <- varCovar((matrix(unlist(grados.in[[k]]), ncol=3)), centered=T)
  grados_out <- varCovar((matrix(unlist(grados.out[[k]]), ncol=3)), centered=T)
  grados_tot <- varCovar((matrix(unlist(grados.tot[[k]]), ncol=3)), centered=T)
  
  mod <- coCovar(atributos.v.list[[k]][f,"MOD"], centered=T)
  anim <- coCovar(atributos.v.list[[k]][f,"ANIM"], centered=T)
  stem <- coCovar(atributos.v.list[[k]][f,"STEM"], centered=T)
  morf <- coCovar(atributos.v.list[[k]][f,"MORF"], centered=T)
  es <- coCovar(atributos.v.list[[k]][f,"ES"], centered=T)
  esp <- coCovar(atributos.v.list[[k]][f,"ESP"], centered=T)
  
  esp_1 <- coCovar(ifelse(atributos.v.list[[k]][f,"ESP"]==1,1,0), centered=T)
  esp_2 <- coCovar(ifelse(atributos.v.list[[k]][f,"ESP"]==2,1,0), centered=T)
  esp_3 <- coCovar(ifelse(atributos.v.list[[k]][f,"ESP"]==3,1,0), centered=T)
  esp_4 <- coCovar(ifelse(atributos.v.list[[k]][f,"ESP"]==4,1,0), centered=T)
  
  fam.lex <- coCovar(atributos.v.list[[k]][f,"FAM.LEX"], centered=T)
  ima.conc <- coCovar(atributos.v.list[[k]][f,"IMA.CONC"], centered=T)
  
  #cual_mod <- sort(match(c("mod_1","mod_2","mod_3","mod_4") , names(as.list(.GlobalEnv))))
  #env_mod <- as.list(.GlobalEnv)[cual_mod]
  
  dum2 <- matrix(c(0,1,0), nrow=nrpalabras, ncol=3, byrow=TRUE)
  dum3 <- matrix(c(0,0,1), nrow=nrpalabras, ncol=3, byrow=TRUE)
  time2 <- varCovar(dum2)
  time3 <- varCovar(dum3)
  
  
  if (  length(mod.unique) == 3 ){
    
    NETdata.list[[k]] <- sienaDataCreate(red, error,mod, time2, time3,
                                         mod_1, mod_3, mod_4, anim, stem, esp_1, esp_2, esp_3, esp_4,
                                         morf, es, fam.lex, ima.conc,esp, grados_in, grados_out, grados_tot
                                         ,changes)
    
  } else {
    
    NETdata.list[[k]] <- sienaDataCreate(red, error, mod,time2, time3,
                                         mod_1, mod_2, mod_3, mod_4, esp_1, esp_2, esp_3, esp_4, 
                                         anim, stem, morf, es, fam.lex, ima.conc,esp, grados_in, grados_out, grados_tot 
                                         , changes)  
  }
  
}


# cuantos nodos quedan igual en "error"


charts <- scores.list

( totalChange <- table(charts[[1]][,1],charts[[1]][,3],useNA='always') ) 
sum(diag(totalChange)) / sum(totalChange) 

( totalChange <- table(charts[[2]][,1],charts[[2]][,3],useNA='always') ) 
sum(diag(totalChange)) / sum(totalChange) 

( totalChange <- table(charts[[3]][,1],charts[[3]][,3],useNA='always') ) 
sum(diag(totalChange)) / sum(totalChange) 

( totalChange <- table(charts[[4]][,1],charts[[4]][,3],useNA='always') ) 
sum(diag(totalChange)) / sum(totalChange) 


#### ESTIMACION

## resultados finales

#library(RSiena)

myRes <- list()

load("siena_sonia_ans2.RData")
myRes[[1]] <- ans2
load("siena_nati_ans1.RData")
myRes[[2]] <- ans1
load("siena_jako_ans3.RData")
myRes[[3]] <- ans3
load("siena_mirka_ans3.RData")
myRes[[4]] <- ans3

suj <- c("SONIA", "NATI", "JAKO", "MIRKA")

for (k in 1:4){
  
  res_table <- outTable(myRes[[k]])[[1]]
  res_table$Pval <- as.character(res_table$Pval) 
  res_table$Pval <- ifelse(res_table$Pval=="0","p<0.001",res_table$Pval)
  
  print(res_table)
  print(xtable::xtable(res_table, caption = paste0(suj[k],". ", outTable(myRes[[k]])[[2]])))
  
  
  gof.id <- RSiena::sienaGOF(myRes[[k]], verbose=TRUE, varName="red", IndegreeDistribution,
                             join=T, cumulative=F, levls=0:3) 
  
  plot <- plot(gof.id, scale=T, center=T) 
  
  x <- paste0("plot_",k)
  assign(x, plot, envir = .GlobalEnv)
  
}

cowplot::plot_grid(plot_1, plot_2, plot_3, plot_4, nrow = 2, ncol = 2)


################################################################
## pruebas hechas para los modelos

#####SONIA

# poner efectos basicos y de red.

k = 1

effects <- getEffects(NETdata.list[[k]])
effects$include[effects$shortName=="quad"] <- FALSE
#effects <- includeEffects(effects,outRate, name="error", type = "rate" ,interaction1="red",fix=F, test=F)
#effects[effects$shortName=="outRate" & effects$name=="red", "include"] <- TRUE
#effects[effects$shortName=="outRate" & effects$name=="red", "initialValue"] <- 0.8
#effects[effects$shortName=="outRate" & effects$name=="error", "initialValue"] <- 0.6
effects[effects$shortName=="density", "include"] <- TRUE
effects$include[effects$effectName=="error linear shape" & effects$type=="eval"] <- FALSE

#effects <- includeEffects(effects, RateX, name="red", type="rate",interaction1="mod_3",  fix=F, test=F)
#effects <- includeEffects(effects, RateX, name="red", type="rate",interaction1="mod_4",  fix=F, test=F)
#effects <- includeEffects(effects, RateX, name="red", type="rate",interaction1="grados_out",  fix=F, test=F)

#effects <- includeEffects(effects, outPop, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects, inAct, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects, outAct, name="red", type="eval",  fix=F, test=F)
effects <- includeEffects(effects, inPop, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects, antiInIso, name="red", type="eval",  fix=F, test=F) 
effects <- includeEffects(effects,antiIso, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects,outInAss, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects,inOutAss, name="red", type="eval",  fix=T, test=T)
#effects <- includeEffects(effects, outTrunc, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects, inIsDegree, name="red", type="eval",  fix=T, test=T)
#effects <- includeEffects(effects, isolatePop, name="red", type="eval",  fix=T, test=T)
#effects <- includeEffects(effects, inStructEq, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects, isolateNet, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects, gwdspFB, name="red", type="eval",  fix=F, test=F)


#effects <- includeEffects(effects, name="red", altX, interaction1="anim", fix=T, test=T)
effects <- includeEffects(effects, name="red", altX, interaction1="fam.lex", fix=F, test=F)
#effects <- includeEffects(effects, name="red", altX, interaction1="ima.conc", fix=T, test=T)
#effects <- includeEffects(effects, name="red", altX, interaction1="esp_3", fix=T, test=T)
#effects <- includeEffects(effects, name="red", altX, interaction1="es", fix=T, test=T)

#effects <- includeEffects(effects, name="red", simRecipX, interaction1="esp_3", fix=T, test=T)
#effects <- includeEffects(effects, name="red", simRecipX, interaction1="es", fix=T, test=T)

effects <- includeEffects(effects, name="red", sameX,interaction1="es", fix=F, test=F)
effects <- includeEffects(effects, name="red", sameX, interaction1="esp_3", fix=F, test=F)
#effects <- includeEffects(effects, name="red", simX,interaction1="morf", fix=T, test=T)
effects <- includeEffects(effects, name="red", simX,interaction1="stem", fix=F, test=F)


#effects <- includeEffects(effects, name="red", egoX, interaction1="anim", fix=F, test=F)
effects <- includeEffects(effects, name="red", egoX, interaction1="fam.lex", fix=F, test=F)
#effects <- includeEffects(effects, name="red", egoX, interaction1="ima.conc", fix=F, test=F)

#effects <- includeEffects(effects,name="error", avRecAlt, interaction1="red")
#effects <- includeEffects(effects,name="error", avInAlt, interaction1="red")
#effects <- includeEffects(effects,name="error", totAltPop, interaction1="red")

#effects <- includeEffects(effects,name="error", popAlt, interaction1="red")

#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "stem",fix=T, test=T)
#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "es", fix=T, test=T)
#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "anim", fix=F, test=F)
#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "esp_3", fix=T, test=T)
#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "morf", fix=T, test=T)
#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "ima.conc", fix=F, test=F)
#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "fam.lex", fix=F, test=F)

#effects <- includeTimeDummy(effects, antiIso,timeDummy="all")

model <- sienaModelCreate(useStdInits=TRUE, projname='results',
                          firstg = 0.05, diagonalize=.2, seed=786840, cond=FALSE, n3 = 3000,
                          doubleAveraging = 0, nsub = 4, MaxDegree =  c(red=8))


ans1 <- siena07(model,data=NETdata.list[[k]],effects=effects,useCluster=TRUE, batch=FALSE, 
                initC=TRUE,nbrNodes=8, prevAns=NULL, returnDeps=TRUE)

outTable(ans1)



testear.efecto(ans1)

siena07ToConvergence(model,NETdata.list[[k]], effects)


#effects <- includeTimeDummy(effects, effFrom, interaction1 = "stem", timeDummy="all", name="error")


ans2 <- siena07(model,data=NETdata.list[[k]],effects=effects,useCluster=TRUE, batch=FALSE, 
                initC=TRUE,nbrNodes=8, prevAns=ans1, returnDeps=TRUE)

outTable(ans2)


ans3 <- siena07(model,data=NETdata.list[[k]],effects=effects,useCluster=TRUE, batch=FALSE, 
                initC=TRUE,nbrNodes=8, prevAns=ans2, returnDeps=TRUE)

outTable(ans3)

#save(ans2, file = "siena_sonia_ans2.RData")


### NATI

k = 2

effects <- getEffects(NETdata.list[[k]])
effects$include[effects$shortName=="quad"] <- FALSE
#effects$include[effects$shortName=="recip"] <- FALSE
#effects <- includeEffects(effects,outRate, name="error", type = "rate" ,interaction1="red",fix=F, test=F)
#effects[effects$shortName=="outRate" & effects$name=="red", "include"] <- TRUE
#effects[effects$shortName=="outRate" & effects$name=="red", "initialValue"] <- 0.8
#effects[effects$shortName=="outRate" & effects$name=="error", "initialValue"] <- 0.6
#effects[effects$shortName=="density", "include"] <- FALSE
#effects$include[effects$effectName=="error linear shape" & effects$type=="eval"] <- FALSE

#effects <- includeEffects(effects, RateX, name="red", type="rate",interaction1="mod_3",  fix=F, test=F)
#effects <- includeEffects(effects, RateX, name="red", type="rate",interaction1="mod_4",  fix=F, test=F)
#effects <- includeEffects(effects, RateX, name="red", type="rate",interaction1="mod",  fix=F, test=F)
#effects <- includeEffects(effects, RateX, name="red", type="rate",interaction1="grados_out",  fix=F, test=F)

#effects <- includeEffects(effects, cycle4, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects, sharedPop, name="red", type="eval",  fix=F, test=F)
effects <- includeEffects(effects, outPop, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects, inAct, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects, outAct, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects, inPop, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects, antiInIso, name="red", type="eval",  fix=F, test=F) 
#effects <- includeEffects(effects,antiIso, name="red", type="eval",  fix=T, test=T)
#effects <- includeEffects(effects,outInAss, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects,inOutAss, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects,outOutAss, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects,inInAss, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects, outTrunc2, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects, inIsDegree, name="red", type="eval",  fix=F, test=F)
effects <- includeEffects(effects, isolatePop, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects, inStructEq, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects, isolateNet, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects, gwdspFB, name="red", type="eval",  fix=F, test=F)


#effects <- includeEffects(effects, name="red", altX, interaction1="anim", fix=T, test=T)
#effects <- includeEffects(effects, name="red", altX, interaction1="fam.lex", fix=T, test=T)
#effects <- includeEffects(effects, name="red", altX, interaction1="ima.conc", fix=T, test=T)
#effects <- includeEffects(effects, name="red", altX, interaction1="esp_3", fix=F, test=F)
#effects <- includeEffects(effects, name="red", altX, interaction1="es", fix=T, test=T)

#effects <- includeEffects(effects, name="red", simRecipX, interaction1="esp_3", fix=T, test=T)
#effects <- includeEffects(effects, name="red", simRecipX, interaction1="es", fix=T, test=T)
#effects <- includeEffects(effects, name="red", simRecipX, interaction1="morf", fix=T, test=T)
#effects <- includeEffects(effects, name="red", simRecipX, interaction1="stem", fix=T, test=T)

#effects <- includeEffects(effects, name="red", sameX,interaction1="es", fix=T, test=T)
effects <- includeEffects(effects, name="red", sameX, interaction1="esp_3", fix=F, test=F)
#effects <- includeEffects(effects, name="red", simX,interaction1="morf", fix=F, test=F)
#effects <- includeEffects(effects, name="red", simX,interaction1="stem", fix=T, test=T)


#effects <- includeEffects(effects, name="red", egoX, interaction1="anim", fix=T, test=T)
#effects <- includeEffects(effects, name="red", egoX, interaction1="fam.lex", fix=F, test=F)
#effects <- includeEffects(effects, name="red", egoX, interaction1="ima.conc", fix=F, test=F)
effects <- includeEffects(effects, name="red", egoX, interaction1="mod_4", fix=F, test=F)
effects <- includeEffects(effects, name="red", egoX, interaction1="es", fix=F, test=F)

#effects <- includeEffects(effects,name="error", avRecAlt, interaction1="red")
#effects <- includeEffects(effects,name="error", avInAlt, interaction1="red")
#effects <- includeEffects(effects,name="error", totAltPop, interaction1="red")

#effects <- includeEffects(effects,name="error", popAlt, interaction1="red")

#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "stem",fix=T, test=T)
#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "es", fix=T, test=T)
#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "anim", fix=F, test=F)
#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "esp_3", fix=F, test=F)
#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "morf", fix=T, test=T)
#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "ima.conc", fix=F, test=F)
#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "fam.lex", fix=F, test=F)

#effects <- includeTimeDummy(effects, inPop, isolatePop, timeDummy="all")

model <- sienaModelCreate(useStdInits=TRUE, projname='results',
                          firstg = 0.05, diagonalize=.2, seed=786840, cond=FALSE, n3 = 3000,
                          doubleAveraging = 0, nsub = 4, MaxDegree =  c(red=8))


ans1 <- siena07(model,data=NETdata.list[[k]],effects=effects,useCluster=TRUE, batch=FALSE, 
                initC=TRUE,nbrNodes=8, prevAns=NULL, returnDeps=TRUE)

outTable(ans1)



testear.efecto(ans1)

siena07ToConvergence(model,NETdata.list[[k]], effects)


#effects <- includeTimeDummy(effects, effFrom, interaction1 = "stem", timeDummy="all", name="error")


ans2 <- siena07(model,data=NETdata.list[[k]],effects=effects,useCluster=TRUE, batch=FALSE, 
                initC=TRUE,nbrNodes=8, prevAns=ans1, returnDeps=TRUE)

outTable(ans2)


ans3 <- siena07(model,data=NETdata.list[[k]],effects=effects,useCluster=TRUE, batch=FALSE, 
                initC=TRUE,nbrNodes=8, prevAns=ans2, returnDeps=TRUE)

outTable(ans3)

#save(ans1, file = "siena_nati_ans1.RData")

##JAKO


k = 3

effects <- getEffects(NETdata.list[[k]])
effects$include[effects$shortName=="quad"] <- FALSE
#effects$include[effects$shortName=="recip"] <- FALSE
#effects <- includeEffects(effects,outRate, name="error", type = "rate" ,interaction1="red",fix=F, test=F)
#effects[effects$shortName=="outRate" & effects$name=="red", "include"] <- TRUE
#effects[effects$shortName=="inRate" & effects$name=="red", "include"] <- TRUE
#effects[effects$shortName=="outRate" & effects$name=="red", "initialValue"] <- 0.8
#effects[effects$shortName=="outRate" & effects$name=="error", "initialValue"] <- 0.6
#effects[effects$shortName=="density", "include"] <- FALSE
#effects$include[effects$effectName=="error linear shape" & effects$type=="eval"] <- FALSE

#effects <- includeEffects(effects, RateX, name="red", type="rate",interaction1="mod_3",  fix=F, test=F)
#effects <- includeEffects(effects, RateX, name="red", type="rate",interaction1="mod_4",  fix=F, test=F)
#effects <- includeEffects(effects, RateX, name="red", type="rate",interaction1="mod",  fix=F, test=F)
#effects <- includeEffects(effects, RateX, name="red", type="rate",interaction1="grados_out",  fix=F, test=F)

# effects <- includeEffects(effects, cycle4, name="red", type="eval",  fix=T, test=T)
# effects <- includeEffects(effects, sharedPop, name="red", type="eval",  fix=T, test=T)
#effects <- includeEffects(effects, outPop, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects, inAct, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects, outAct, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects, inPop, name="red", type="eval",  fix=F, test=F)
# effects <- includeEffects(effects, antiInIso, name="red", type="eval",  fix=T, test=T) 
#effects <- includeEffects(effects,antiIso, name="red", type="eval",  fix=F, test=F)
effects <- includeEffects(effects,outInAss, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects,inOutAss, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects,outOutAss, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects,inInAss, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects, outTrunc2, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects, inIsDegree, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects, isolatePop, name="red", type="eval",  fix=F, test=F)
effects <- includeEffects(effects, inStructEq, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects, isolateNet, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects, gwdspFB, name="red", type="eval",  fix=T, test=T)


#effects <- includeEffects(effects, name="red", altX, interaction1="anim", fix=T, test=T)
#effects <- includeEffects(effects, name="red", altX, interaction1="fam.lex", fix=T, test=T)
#effects <- includeEffects(effects, name="red", altX, interaction1="ima.conc", fix=T, test=T)
#effects <- includeEffects(effects, name="red", altX, interaction1="esp_3", fix=F, test=F)
#effects <- includeEffects(effects, name="red", altX, interaction1="es", fix=T, test=T)

#effects <- includeEffects(effects, name="red", simRecipX, interaction1="esp_3", fix=T, test=T)
#effects <- includeEffects(effects, name="red", simRecipX, interaction1="es", fix=T, test=T)
#effects <- includeEffects(effects, name="red", simRecipX, interaction1="morf", fix=T, test=T)
#effects <- includeEffects(effects, name="red", simRecipX, interaction1="stem", fix=T, test=T)

#effects <- includeEffects(effects, name="red", sameX,interaction1="es", fix=T, test=T)
#effects <- includeEffects(effects, name="red", sameX, interaction1="esp_3", fix=T, test=T)
effects <- includeEffects(effects, name="red", simX,interaction1="morf", fix=F, test=F)
#effects <- includeEffects(effects, name="red", simX,interaction1="stem", fix=T, test=T)


#effects <- includeEffects(effects, name="red", egoX, interaction1="anim", fix=F, test=F)
#effects <- includeEffects(effects, name="red", egoX, interaction1="fam.lex", fix=F, test=F)
effects <- includeEffects(effects, name="red", egoX, interaction1="ima.conc", fix=F, test=F)
#effects <- includeEffects(effects, name="red", egoX, interaction1="mod_4", fix=T, test=T)
#effects <- includeEffects(effects, name="red", egoX, interaction1="es", fix=T, test=T)
#effects <- includeEffects(effects, name="red", egoX, interaction1="esp_3", fix=T, test=T)

#effects <- includeEffects(effects,name="error", avRecAlt, interaction1="red")
#effects <- includeEffects(effects,name="error", avInAlt, interaction1="red")
#effects <- includeEffects(effects,name="error", totAltPop, interaction1="red")

#effects <- includeEffects(effects,name="error", popAlt, interaction1="red")

#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "stem",fix=F, test=F)
#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "es", fix=T, test=T)
#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "anim", fix=F, test=F)
#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "esp_3", fix=T, test=T)
#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "morf", fix=F, test=F)
#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "ima.conc", fix=F, test=F)
#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "fam.lex", fix=F, test=F)

#effects <- includeTimeDummy(effects, inPop, isolatePop, timeDummy="all")




model <- sienaModelCreate(useStdInits=TRUE, projname='results',
                          firstg = 0.01, diagonalize=.2, seed=786840, cond=FALSE, n3 = 1000,
                          doubleAveraging = 0, nsub = 4, MaxDegree =  c(red=14))

ans1 <- siena07(model,data=NETdata.list[[k]],effects=effects,useCluster=TRUE, batch=FALSE,
                initC=TRUE, nbrNodes=8, returnDeps=TRUE)

outTable(ans1)  


ans2 <- siena07(model,data=NETdata.list[[k]],effects=effects,useCluster=TRUE, batch=FALSE, 
                initC=TRUE,nbrNodes=8, prevAns=ans1, returnDeps=TRUE)

outTable(ans2)

ans3 <- siena07(model,data=NETdata.list[[k]],effects=effects,useCluster=TRUE, batch=FALSE, 
                initC=TRUE,nbrNodes=8, prevAns=ans2, returnDeps=TRUE)

outTable(ans3)

ans4 <- siena07(model,data=NETdata.list[[k]],effects=effects,useCluster=TRUE, batch=FALSE, 
                initC=TRUE,nbrNodes=8, prevAns=ans3, returnDeps=TRUE)

outTable(ans4)




testear.efecto(ans4)



#save(ans3, file = "siena_jako_ans3.RData")
#load("siena_jako_ans3.RData")


##MIRKA



k = 4

#print01Report(NETdata.list[[k]], modelname = 'report_mirka')


effects <- getEffects(NETdata.list[[k]])

effects$include[effects$shortName=="quad"] <- FALSE
#effects <- includeEffects(effects,outRate, name="error", type = "rate" ,interaction1="red",fix=F, test=F)
#effects[effects$shortName=="outRate" & effects$name=="red", "include"] <- TRUE
#effects[effects$shortName=="outRate" & effects$name=="red", "initialValue"] <- 0.8
#effects[effects$shortName=="outRate" & effects$name=="error", "initialValue"] <- 0.6
#effects[effects$shortName=="density", "include"] <- FALSE
effects$include[effects$effectName=="error linear shape" & effects$type=="eval"] <- FALSE


#effects <- includeEffects(effects, RateX, name="red", type="rate",interaction1="grados_out",  fix=F, test=F)

#effects <- includeEffects(effects, RateX, name="red", type="rate",interaction1="mod",  fix=T, test=T)
# effects <- includeEffects(effects, RateX, name="red", type="rate",interaction1="stem",  fix=T, test=T)
#effects <- includeEffects(effects, RateX, name="red", type="rate",interaction1="fam.lex",  fix=F, test=F)
#effects <- includeEffects(effects, RateX, name="red", type="rate",interaction1="ima.conc",  fix=F, test=F)
# effects <- includeEffects(effects, RateX, name="red", type="rate",interaction1="morf",  fix=T, test=T)
#effects <- includeEffects(effects, RateX, name="red", type="rate",interaction1="es",  fix=T, test=T)
#effects <- includeEffects(effects, RateX, name="red", type="rate",interaction1="anim",  fix=T, test=T)
#effects <- includeEffects(effects, RateX, name="red", type="rate",interaction1="esp",  fix=T, test=T)
#effects <- includeEffects(effects, RateX, name="red", type="rate",interaction1="mod_1",  fix=T, test=T)
#effects <- includeEffects(effects, RateX, name="red", type="rate",interaction1="mod_3",  fix=T, test=T)
#effects <- includeEffects(effects, RateX, name="red", type="rate",interaction1="mod_4",  fix=T, test=T)

effects <- includeEffects(effects, balance, name="red", type="eval",  fix=F, test=F)
#effects <- includeEffects(effects, transTies, name="red", type="eval",  fix=T, test=T)
#  effects <- includeEffects(effects, cycle4, name="red", type="eval",  fix=F, test=F)
#  effects <- includeEffects(effects, sharedPop, name="red", type="eval",  fix=F, test=F)
#   effects <- includeEffects(effects, outPop, name="red", type="eval",  fix=F, test=F)
#  effects <- includeEffects(effects, inAct, name="red", type="eval",  fix=T, test=T)
#  effects <- includeEffects(effects, outAct, name="red", type="eval",  fix=T, test=T)
#  effects <- includeEffects(effects, inPop, name="red", type="eval",  fix=F, test=F)
effects <- includeEffects(effects, antiInIso, name="red", type="eval",  fix=F, test=F)
#  effects <- includeEffects(effects,antiIso, name="red", type="eval",  fix=T, test=T)
#effects <- includeEffects(effects,outInAss, name="red", type="eval",  fix=T, test=T)
#  effects <- includeEffects(effects,inOutAss, name="red", type="eval",  fix=T, test=T)
#effects <- includeEffects(effects,outOutAss, name="red", type="eval",  fix=T, test=T)
#effects <- includeEffects(effects,inInAss, name="red", type="eval",  fix=T, test=T)
#effects <- includeEffects(effects, outTrunc2, name="red", type="eval",  fix=T, test=T)
#effects <- includeEffects(effects, inIsDegree, name="red", type="eval",  fix=T, test=T)
#  effects <- includeEffects(effects, isolatePop, name="red", type="eval",  fix=T, test=T)
#  effects <- includeEffects(effects, inStructEq, name="red", type="eval",  fix=T, test=T)
#effects <- includeEffects(effects, isolateNet, name="red", type="eval",  fix=T, test=T)
#  effects <- includeEffects(effects, gwdspFB, name="red", type="eval",  fix=F, test=F)


#effects <- includeEffects(effects, name="red", altX, interaction1="anim", fix=F, test=F)
#effects <- includeEffects(effects, name="red", altX, interaction1="fam.lex", fix=F, test=F)
effects <- includeEffects(effects, name="red", altX, interaction1="ima.conc", fix=F, test=F)
#effects <- includeEffects(effects, name="red", altX, interaction1="esp_3", fix=F, test=F)
#effects <- includeEffects(effects, name="red", altX, interaction1="es", fix=F, test=F)

#effects <- includeEffects(effects, name="red", simRecipX, interaction1="esp_3", fix=T, test=T)
#effects <- includeEffects(effects, name="red", simRecipX, interaction1="es", fix=T, test=T)
#effects <- includeEffects(effects, name="red", simRecipX, interaction1="morf", fix=T, test=T)
#effects <- includeEffects(effects, name="red", simRecipX, interaction1="stem", fix=T, test=T)

effects <- includeEffects(effects, name="red", sameX,interaction1="es", fix=F, test=F)
effects <- includeEffects(effects, name="red", sameX, interaction1="esp_3", fix=F, test=F)
#effects <- includeEffects(effects, name="red", simX,interaction1="morf", fix=T, test=T)
#effects <- includeEffects(effects, name="red", simX,interaction1="stem", fix=T, test=T)


#effects <- includeEffects(effects, name="red", egoX, interaction1="anim", fix=T, test=T)
#effects <- includeEffects(effects, name="red", egoX, interaction1="fam.lex", fix=T, test=T)
#effects <- includeEffects(effects, name="red", egoX, interaction1="ima.conc", fix=T, test=T)
#effects <- includeEffects(effects, name="red", egoX, interaction1="mod_4", fix=T, test=T)
#effects <- includeEffects(effects, name="red", egoX, interaction1="es", fix=T, test=T)
#effects <- includeEffects(effects, name="red", egoX, interaction1="esp_3", fix=T, test=T)

#effects <- includeEffects(effects,name="error", avRecAlt, interaction1="red")
#effects <- includeEffects(effects,name="error", avInAlt, interaction1="red")
#effects <- includeEffects(effects,name="error", totAltPop, interaction1="red")

#effects <- includeEffects(effects,name="error", popAlt, interaction1="red")

#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "stem",fix=T, test=T)
#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "es", fix=T, test=T)
#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "anim", fix=T, test=T)
#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "esp_3", fix=T, test=T)
#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "morf", fix=F, test=F)
effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "ima.conc", fix=F, test=F)
#effects <- includeEffects(effects, name = "error", effFrom, interaction1 = "fam.lex", fix=T, test=T)

#effects <- includeTimeDummy(effects, inPop, timeDummy="all")

#effects <- includeEffects(effects, inPop, interaction1='time2')
#effects <- includeEffects(effects, inPop, interaction1='time3')

model <- sienaModelCreate(useStdInits=TRUE, projname='results',
                          firstg = 0.01, diagonalize=.2, seed=786840, cond=FALSE,
                          doubleAveraging = 0, nsub = 4, MaxDegree =  c(red=14), n3 = 3000)

ans1 <- siena07(model,data=NETdata.list[[k]],effects=effects,useCluster=TRUE, batch=FALSE,
                initC=TRUE, nbrNodes=8, returnDeps=TRUE)

outTable(ans1)  


ans2 <- siena07(model,data=NETdata.list[[k]],effects=effects,useCluster=TRUE, batch=FALSE, 
                initC=TRUE,nbrNodes=8, prevAns=ans1, returnDeps=TRUE)

outTable(ans2)

ans3 <- siena07(model,data=NETdata.list[[k]],effects=effects,useCluster=TRUE, batch=FALSE, 
                initC=TRUE,nbrNodes=8, prevAns=ans2, returnDeps=TRUE)

outTable(ans3)

ans4 <- siena07(model,data=NETdata.list[[k]],effects=effects,useCluster=TRUE, batch=FALSE, 
                initC=TRUE,nbrNodes=8, prevAns=ans3, returnDeps=TRUE)

outTable(ans4)


ans5 <- siena07(model,data=NETdata.list[[k]],effects=effects,useCluster=TRUE, batch=FALSE, 
                initC=TRUE,nbrNodes=8, prevAns=ans4, returnDeps=TRUE)

outTable(ans5)



testear.efecto(ans1)

#save(ans3, file = "siena_mirka_ans3.RData")






## Goodness of fit

myRes <- ans1

# time heterogeinity

timetest <- sienaTimeTest(myRes)
summary(timetest)
plot(timetest)

# indegree


gof.id <- sienaGOF(myRes, verbose=TRUE, varName="red", IndegreeDistribution,
                   join=T, cumulative=F, levls=0:3) 

plot(gof.id, scale=T, center=T)  


# outdegree

gof.od <- sienaGOF(myRes, verbose=TRUE, varName="red", OutdegreeDistribution,
                   join=T, cumulative=F, levls=0:3) 
plot(gof.od, scale=T, center=T)  


gof.behaviourS <- sienaGOF(myRes,BehaviorDistribution, levls=0:6,
                           verbose=TRUE,join=TRUE,varName="error") 
plot(gof.behaviourS, scale=T, center=T) 


gof.gd <- sienaGOF(myRes, verbose=TRUE, varName="red", GeodesicDistribution,
                   join=T, cumulative=F) 
plot(gof.gd)




siena07ToConvergence <- function(alg,dat,eff){
  numr <- 0
  ans <- siena07(alg, data=dat, effects=eff,
                 useCluster=TRUE, initC=TRUE, nbrNodes=8) 
  repeat {
    numr <- numr+1 
    maxt <- max(abs(ans$tconv [!eff$fix[eff$include]]))
    cat (numr , maxt , "\n") 
    if (maxt < 0.10 ) { break } 
    if (maxt > 5) { break }
    if (numr > 10) { break } 
    ans <-  siena07(alg, data=dat, effects=eff , prevAns=ans ,
                    useCluster=TRUE, initC=TRUE, nbrNodes=8)
  }
  ans  
}

GeodesicDistribution <- function (i, data, sims, period, groupName,
                                  varName, levls=c(1:5,Inf), cumulative=TRUE, ...) {
  x <- networkExtraction(i, data, sims, period, groupName, varName)
  require(sna)
  a <- sna::geodist(symmetrize(x))$gdist
  if (cumulative)
  {
    gdi <- sapply(levls, function(i){ sum(a<=i) })
  }
  else
  {
    gdi <- sapply(levls, function(i){ sum(a==i) })
  }
  names(gdi) <- as.character(levls)
  gdi
}


