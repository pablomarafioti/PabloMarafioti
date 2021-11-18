
### Modelación de errores de concordancia en español LE usando sistemas dinámicos
### simulaciones

library(deSolve)
library(ggplot2)
library(phaseR)
library(binom)

# datos y contextos


load("datos_articulo_6.RData")


dif.data <-  datos[, c("RES_CAT", "ID", "SESION" , "INSTANCIA", "CUMRES", "TIME",
"MOD", "ANIM", "ES","MORF.f", "Fabs.SC.f", "CUMRES.f", "FAM.LEX.f", "EST1", "EST5", "EST2")]

# contextos 

error.genero <- data.frame(MOD = "3", Fabs.SC.f = "1", EST1 = "1")
error.epente <- data.frame(MORF.f = "1", FAM.LEX.f = "1", EST5 = "1", EST2 = "1")  
error.plural <- data.frame(ES = "2", MORF.f = "1/2", Fabs.SC.f = "1", FAM.LEX.f = "1")  
error.mixto  <- data.frame(Fabs.SC.f = "1", FAM.LEX.f = "1", EST1 = "1", EST5 = "1" )


########
### modelo de Lokta-Volterra
########


lvcomp <- function(t, n, parms) {
   with(as.list(parms), {
     dn1dt <- r1 * n[1] * (1 - a11 * n[1] - a12 * n[2])
     dn2dt <- r2 * n[2] * (1 - a22 * n[2] - a21 * n[1])
     list(c(dn1dt, dn2dt))
     })
   }

parms <- c(r1 = 0.5, r2 = 0.5, a11 = 1, a21 = 1.5, a22 = 1,  a12 = 1.5)

t = 100



w <- seq(0.01,0.99, by = 0.01)


df.res <- list()
res.loss <- vector()

# cuánto vale w ?

for (h in 1:length(w)) { 

  print(h)
  
res <- list(vector(), vector(), vector(), vector())

for  (l in 1:4){
  
  lv.data <- dif.data[dif.data$ID==l,]
  error_antes <-  (exp(lv.data$CUMRES)-1) / (lv.data$TIME-1)
  error_antes[is.na(error_antes)] <- 0
  
  
  for (i in 1:dim(lv.data)[1]){
    
    conc = data.frame(lv.data[i,-c(1:4)])
    comp = list(conc, error.genero, error.epente, error.plural, error.mixto)
    
    
    m = rep(0,5)  
    
    for (j in 1:5){ 
      
      A = conc[colnames(comp[[j]])]
      
      if (j == 4 ){
        B =  comp[[j]]
        
        index.morf  = which(unlist(strsplit(as.character(comp[[j]]$MORF.f), "[/]")) ==  as.character(comp[[1]]$MORF.f))
        if (length(index.morf) == 0) {B$MORF.f = "00"} else {B$MORF.f =  ifelse(index.morf == 1, "1", "2") }
      } else {B = comp[[j]]}
      
      A <- apply(A, 2, function(x){as.numeric(as.character(x))})
      B <- apply(B, 2, function(x){as.numeric(as.character(x))})
      
      
       
      m[j] <- length(A) - sum(A==B) 
      
    }
    
  
    
    E <- (w[h])*mean(crch::ptlogis(m, location = 2, left = -1, right = 4)) + (1-w[h])*error_antes[i]  
     
    
    if (E==0.5){E <- ifelse(runif(1)<0.5,0.25,0.75)}
     
      
    initialN <- c(1-E,  E)
    
    out <- ode(y = initialN, times = 1:t, func = lvcomp, parms = parms)
   
    
    res[[l]][i] <- which.max(out[t,c("1","2")])
    
  }
   
}

# resultados ("1" = correcto, "2" = error)

sim  = rbind(prop.table(table(res[[1]])), prop.table(table(res[[2]])),
      prop.table(table(res[[3]])), prop.table(table(res[[4]])))

real  = rbind(prop.table(table(datos[datos$ID==1, "RES_BIN"])),
          prop.table(table(datos[datos$ID==2, "RES_BIN"])),
          prop.table(table(datos[datos$ID==3, "RES_BIN"])),
          prop.table(table(datos[datos$ID==4, "RES_BIN"]))
)

res_df <- data.frame(error.sim = sim[,2], error.real = real[,2],
          correcto.sim = sim[,1], correcto.real = real[,1])

row.names(res_df) <- c("SONIA","NATI","JAKO","MIRKA")


df.res[[h]]  <- res_df

res.loss[h] <-  MLmetrics::RMSE(df.res[[h]][,1], df.res[[h]][,2])

}


df.res[[which.min(res.loss)]]

# h = 71 --> w = 0.71


# Simulación con w = 0.71

h = 71
res <- list(vector(), vector(), vector(), vector())

for  (l in 1:4){
  
  lv.data <- dif.data[dif.data$ID==l,]
  error_antes <-  (exp(lv.data$CUMRES)-1) / (lv.data$TIME-1)
  error_antes[is.na(error_antes)] <- 0
  
  
  for (i in 1:dim(lv.data)[1]){
    
    conc = data.frame(lv.data[i,-c(1:4)])
    comp = list(conc, error.genero, error.epente, error.plural, error.mixto)
    
    
    m = rep(0,5)  
    
    for (j in 1:5){ 
      
      A = conc[colnames(comp[[j]])]
      
      if (j == 4 ){
        B =  comp[[j]]
        
        index.morf  = which(unlist(strsplit(as.character(comp[[j]]$MORF.f), "[/]")) ==  as.character(comp[[1]]$MORF.f))
        if (length(index.morf) == 0) {B$MORF.f = "00"} else {B$MORF.f =  ifelse(index.morf == 1, "1", "2") }
      } else {B = comp[[j]]}
      
      A <- apply(A, 2, function(x){as.numeric(as.character(x))})
      B <- apply(B, 2, function(x){as.numeric(as.character(x))})
      
      
      
      m[j] <- length(A) - sum(A==B) 
      
    }
    
    
    
    E <- (w[h])*mean(crch::ptlogis(m, location = 2, left = -1, right = 4)) + (1-w[h])*error_antes[i]  
    
    
    if (E==0.5){E <- ifelse(runif(1)<0.5,0.25,0.75)}
    
    
    initialN <- c(1-E,  E)
    
    out <- ode(y = initialN, times = 1:t, func = lvcomp, parms = parms)
    
    
    res[[l]][i] <- which.max(out[t,c("1","2")])
    
  }
  
}

# resultados ("1" = correcto, "2" = error)

sim  = rbind(prop.table(table(res[[1]])), prop.table(table(res[[2]])),
             prop.table(table(res[[3]])), prop.table(table(res[[4]])))

real  = rbind(prop.table(table(datos[datos$ID==1, "RES_BIN"])),
              prop.table(table(datos[datos$ID==2, "RES_BIN"])),
              prop.table(table(datos[datos$ID==3, "RES_BIN"])),
              prop.table(table(datos[datos$ID==4, "RES_BIN"]))
)

res_df <- data.frame(error.sim = sim[,2], error.real = real[,2],
                     correcto.sim = sim[,1], correcto.real = real[,1])

row.names(res_df) <- c("SONIA","NATI","JAKO","MIRKA")




dist.sonia <- vector()
dist.nati <- vector()
dist.jako <- vector()
dist.mirka <- vector()


  
  error_por_sesiones_sim <- list()
  error_por_sesiones_real <- list()
  
  
  for (i in 1:4){
    
    error_por_sesiones_sim[[i]] <-
      
      tapply(ifelse(as.numeric(as.character(res[[i]]))==2,1,0), datos[datos$ID==i,"SESION"], sum)
    
    error_por_sesiones_real[[i]] <-  
      
      tapply(as.numeric(as.character(datos[datos$ID==i, "RES_BIN"])), datos[datos$ID==i,"SESION"], sum)
    
  }
  
  dist_error_sesion <- c(
  TSdist::MinkowskiDistance(as.vector(error_por_sesiones_sim[[1]][-c(13,14)]),as.vector(error_por_sesiones_real[[1]][-c(13,14)]), p = 3),
  TSdist::MinkowskiDistance(as.vector(error_por_sesiones_sim[[2]]),as.vector(error_por_sesiones_real[[2]]), p = 3),
  TSdist::MinkowskiDistance(as.vector(error_por_sesiones_sim[[3]]),as.vector(error_por_sesiones_real[[3]]), p = 3),
  TSdist::MinkowskiDistance(as.vector(error_por_sesiones_sim[[4]][-c(13,14)]),as.vector(error_por_sesiones_real[[4]][-c(13,14)]), p = 3)
  )



NN <- tapply(datos$RES_BIN, datos$ID, length)


prop_test_error   <- c(
prop.test(c(res_df[1,1]*NN[1], res_df[1,2]*NN[1]), rep(NN[1],2))$p.value,
prop.test(c(res_df[2,1]*NN[2], res_df[2,2]*NN[2]), rep(NN[2],2))$p.value,
prop.test(c(res_df[3,1]*NN[3], res_df[3,2]*NN[3]), rep(NN[3],2))$p.value,
prop.test(c(res_df[4,1]*NN[4], res_df[4,2]*NN[4]), rep(NN[4],2))$p.value
)

resultados_df <- data.frame(res_df[,-c(3,4)], p.value = prop_test_error, dist = dist_error_sesion)

# xtable::xtable(resultados_df, digits = 3,
# caption = "Resultados de la simulación: 
# (i) proporciones entre errores reales y simulados; 
# (ii)  p-valores de un test para diferencia entre dichas proporciones;
# (iii) distancia de Minkowski entre los errores simulados y reales por sesiones.")


# plot por sesiones

sesion.tot <- list()

for (i in 1:4){
  sesion.tot[[i]]  <- tapply(datos[datos$ID==i,"RES_BIN"], datos[datos$ID==i,"SESION"], length)
}

df.plot.1 <- data.frame(
  SESIONES = c(1:12, 1:12), 
  ERROR.PROP = c(error_por_sesiones_sim[[1]][1:12]/sesion.tot[[1]][1:12]
                 , error_por_sesiones_real[[1]][1:12]/sesion.tot[[1]][1:12]),
  CLASE = c(rep("SONIA.SIM", 12), rep("SONIA.REAL", 12))
)

p.1 <- ggplot(data = df.plot.1, aes(x=SESIONES, y=ERROR.PROP, group=CLASE, color=CLASE, shape=CLASE)) +
  geom_point() +   theme_bw() + ggtitle("Errores reales y simulados (SONIA)") +
  scale_color_manual(values=c("black", "black")) +
  scale_x_continuous(breaks = 1:12)


df.plot.2 <- data.frame(
  SESIONES = c(1:14, 1:14), 
  ERROR.PROP = c(error_por_sesiones_sim[[2]]/sesion.tot[[2]]
                 , error_por_sesiones_real[[2]]/sesion.tot[[2]]),
  CLASE = c(rep("NATI.SIM", 14), rep("NATI.REAL", 14))
  
)

p.2 <- ggplot(data = df.plot.2, aes(x=SESIONES, y=ERROR.PROP, group=CLASE, color=CLASE, shape=CLASE)) +
  geom_point() +   theme_bw() + ggtitle("Errores reales y simulados (NATI)") +
  scale_color_manual(values=c("black", "black")) +
  scale_x_continuous(breaks = 1:14)


df.plot.3 <- data.frame(
  SESIONES = c(1:14, 1:14), 
  ERROR.PROP = c(error_por_sesiones_sim[[3]]/sesion.tot[[3]]
                 , error_por_sesiones_real[[3]]/sesion.tot[[3]]),
  CLASE = c(rep("JAKO.SIM", 14), rep("JAKO.REAL", 14))
  
)

p.3 <- ggplot(data = df.plot.3, aes(x=SESIONES, y=ERROR.PROP, group=CLASE, color=CLASE, shape=CLASE)) +
  geom_point() +   theme_bw() + ggtitle("Errores reales y simulados (JAKO)") +
  scale_color_manual(values=c("black", "black")) +
  scale_x_continuous(breaks = 1:14)


df.plot.4 <- data.frame(
  SESIONES = c(1:12, 1:12), 
  ERROR.PROP = c(error_por_sesiones_sim[[4]][1:12]/sesion.tot[[4]][1:12]
                 , error_por_sesiones_real[[4]][1:12]/sesion.tot[[4]][1:12]),
  CLASE = c(rep("MIRKA.SIM", 12), rep("MIRKA.REAL", 12))
)

p.4 <- ggplot(data = df.plot.4, aes(x=SESIONES, y=ERROR.PROP, group=CLASE, color=CLASE, shape=CLASE)) +
  geom_point() +   theme_bw() + ggtitle("Errores reales y simulados (MIRKA)") +
  scale_color_manual(values=c("black", "black")) +
  scale_x_continuous(breaks = 1:12)

png("simulacion_lv_sesiones.png")
gridExtra::grid.arrange(p.1, p.2, p.3, p.4)
dev.off()


# plot ejemplos (SONIA: "muchos pueblos", "les idiomas")

# generar out con l = 1,  i = 2 ("muchos pueblos"), i = 10 ("les idiomas")  

h = 71

lv.data <- dif.data[dif.data$ID==l,]
error_antes <-  (exp(lv.data$CUMRES)-1) / (lv.data$TIME-1)
error_antes[is.na(error_antes)] <- 0

  
  conc = data.frame(lv.data[i,-c(1:4)])
  comp = list(conc, error.genero, error.epente, error.plural, error.mixto)
  
  
  m = rep(0,5)  
  
  for (j in 1:5){ 
    
    A = conc[colnames(comp[[j]])]
    
    if (j == 4 ){
      B =  comp[[j]]
      
      index.morf  = which(unlist(strsplit(as.character(comp[[j]]$MORF.f), "[/]")) ==  as.character(comp[[1]]$MORF.f))
      if (length(index.morf) == 0) {B$MORF.f = "00"} else {B$MORF.f =  ifelse(index.morf == 1, "1", "2") }
    } else {B = comp[[j]]}
    
    A <- apply(A, 2, function(x){as.numeric(as.character(x))})
    B <- apply(B, 2, function(x){as.numeric(as.character(x))})
    
    
    
    m[j] <- length(A) - sum(A==B) 
    
  }
  
  
  
  E <- (w[h])*mean(crch::ptlogis(m, location = 2, left = -1, right = 4)) + (1-w[h])*error_antes[i]  
  
  
  if (E==0.5){E <- ifelse(runif(1)<0.5,0.25,0.75)}
  
  
  initialN <- c(1-E,  E)
  
  out <- ode(y = initialN, times = 1:t, func = lvcomp, parms = parms)
  


df.plot <- data.frame(TIME = rep(1:(t),2), 
                      ACT = c(out[,"1"], out[,"2"]),
                      CLASE = c(rep("correcto", t), rep("error",t))
)

ggplot(data = df.plot, aes(x=TIME, y=ACT, group=CLASE, color=CLASE)) +
  geom_line(aes(linetype=CLASE)) +   theme_bw() + ggtitle("Flujo en el tiempo") +
  xlab("Tiempo") + ylab("Activación") +
  scale_color_manual(values=c("black", "black")) +
  scale_linetype_manual(values=c("dashed", "dotted"))


#### diagramas de flujo de ejemplos

parms <- c(r1 = 0.5, r2 = 0.5, a11 = 1, a21 = 1.5, a22 = 1,  a12 = 1.5)

 lotkaVolterra.flowField <-
   flowField(lvcomp, xlim = c(0, 1), ylim = c(0, 1), state.names = c("corrrecto","error"),
               parameters = parms, points = 19, add = FALSE)
 grid() 
 
  lotkaVolterra.nullclines <-
   nullclines(lvcomp, xlim = c(0, 1), ylim = c(0, 1), state.names = c("corrrecto","error"),
                parameters = parms, points = 500, add.legend = F, col = c("black", "black"))

# "los pueblos":   
 lotkaVolterra.trajectory.1 <-
   trajectory(lvcomp, y0 = c(1-0.492579, 0.492579), tlim = c(1,100), state.names = c("corrrecto","error"),
                parameters = parms, col = "black", lty = "dashed")
# "les idiomas":
 lotkaVolterra.trajectory.2 <-
   trajectory(lvcomp, y0 = c(1-0.5616874, 0.5616874), tlim = c(1,100), state.names = c("corrrecto","error"),
              parameters = parms, col = "black", lty = "dotted")
  
 

 ##########################
 ## evolutionary games
 

 ## modelo 
 
 rd <- function(t, n, parms) {
   with(as.list(parms), {
     dn1dt <- n[1] * (1 - n[1]) * (n[1] * (a11 + a22 - a12 - a21) + a12 - a22)
     dn2dt <- -dn1dt
     list(c(dn1dt, dn2dt))
   })
 }
 
 
 parms <- c(a11 = 3, a21 = 0, a22 = 1,  a12 = 0)
 
 
 error_por_sesiones_real <- list()
 
 for (z in 1:4){
  
   error_por_sesiones_real[[z]] <-  
     
  tapply(as.numeric(as.character(datos[datos$ID==z, "RES_BIN"])), datos[datos$ID==z,"SESION"], sum)
   
 }
 
# cuánto vale theta ?
 
 par.grid <- expand.grid(mu1 = seq(0.2, 0.3, by = 0.01) , mu2 = seq(0.8,0.9, by = 0.01))
 
 df.res <- list()
 res.loss <- vector()
 
 
 t = 50
 
 for (h in 1:dim(par.grid)[1]) {
 
 print(h)
         
   res <- list(vector(), vector(), vector(), vector())
   
   for  (l in 1:4){
     
     lv.data <- dif.data[dif.data$ID==l,]
     error_antes <-  (exp(lv.data$CUMRES)-1) / (lv.data$TIME-1)
     error_antes[is.na(error_antes)] <- 0
     
     
     for (i in 1:dim(lv.data)[1]){
       
       conc = data.frame(lv.data[i,-c(1:4)])
       comp = list(conc, error.genero, error.epente, error.plural, error.mixto)
       
       
       m = rep(0,5)  
       
       for (j in 1:5){ 
         
         A = conc[colnames(comp[[j]])]
         
         if (j == 4 ){
           B =  comp[[j]]
           
           index.morf  = which(unlist(strsplit(as.character(comp[[j]]$MORF.f), "[/]")) ==  as.character(comp[[1]]$MORF.f))
           if (length(index.morf) == 0) {B$MORF.f = "00"} else {B$MORF.f =  ifelse(index.morf == 1, "1", "2") }
         } else {B = comp[[j]]}
         
         A <- apply(A, 2, function(x){as.numeric(as.character(x))})
         B <- apply(B, 2, function(x){as.numeric(as.character(x))})
         
         
         m[j] <- length(A) - sum(A==B) 
         
       }
       
       
       
       #mu1 <- 0.26
       #mu2 <- 0.88
       
       
       mu1 <- par.grid[h ,"mu1"]
       mu2 <- par.grid[h ,"mu2"] 
       
       theta <- c(mu1, mu2)
       
       theta_est <- vector()
       theta_est2 <- vector()
       
       prior1 <- round(error_antes[i], 2)
       prior2 <- 1 - prior1
        
       
       for (u in 2:5){
        
       Likxprior1 = dbinom(x=m[u],size=length(comp[[u]]), p=theta[1])*prior1
       Likxprior2 = dbinom(x=m[u],size=length(comp[[u]]), p=theta[2])*prior2
       marginal = Likxprior1 + Likxprior2
      
       post_1 = Likxprior1 / marginal
       post_2 = Likxprior2 / marginal
       
      posterior <- c(post_1, post_2)
      theta_est[u] = which.max(posterior)
      
      theta_est2[u] = max(posterior)
       }
        
       tabla <- data.frame(error = 0, correcto = 0)
       tabla[, as.numeric(names(table(theta_est[-1])))] <- table(theta_est[-1])
       
       
        ifelse(tabla[,1]==tabla[,2], E <- mean(theta_est2[2:5]), E <- theta[which.max(tabla)]) 
       
      
         
       initialN <- c(1 - E, E)
       
       out <- ode(y = initialN, times = 1:t, func = rd, parms = parms)
       
       # 1 = error, 2 = correcto 
       res[[l]][i] <- which.max(out[t,c("1","2")])
       
     }
     
   }
    
   
   
   # resultados
   
   sim  = rbind(prop.table(table(res[[1]])), prop.table(table(res[[2]])),
                prop.table(table(res[[3]])), prop.table(table(res[[4]])))
   real  = rbind(prop.table(table(datos[datos$ID==1, "RES_BIN"])),
                 prop.table(table(datos[datos$ID==2, "RES_BIN"])),
                 prop.table(table(datos[datos$ID==3, "RES_BIN"])),
                 prop.table(table(datos[datos$ID==4, "RES_BIN"]))
   )
   
   
   res_df <- data.frame(error.sim = sim[,1], error.real = real[,2],
                        correcto.sim = sim[,2], correcto.real = real[,1])
   
   row.names(res_df) <- c("SONIA","NATI","JAKO","MIRKA")
   
   
   df.res[[h]]  <- res_df
   
     res.loss[h] <-  MLmetrics::RMSE(df.res[[h]][,1], df.res[[h]][,2])
     
   
 }
 
 
 #res_ev_games <- list(df.res, res.loss)
 
 #save(res_ev_games, file = "res_ev_games.RData")
 
 df.res[[which.min(res.loss)]]
 
 par.grid[95,]
 #m1      m2    a1   a2
 # 0.26  0.88    3    1
 
# simulación con theta = (0.26, 0.88)
 

 res <- list(vector(), vector(), vector(), vector())
 
 for  (l in 1:4){
   
   lv.data <- dif.data[dif.data$ID==l,]
   error_antes <-  (exp(lv.data$CUMRES)-1) / (lv.data$TIME-1)
   error_antes[is.na(error_antes)] <- 0
   
   
   for (i in 1:dim(lv.data)[1]){
     
     conc = data.frame(lv.data[i,-c(1:4)])
     comp = list(conc, error.genero, error.epente, error.plural, error.mixto)
     
     
     m = rep(0,5)  
     
     for (j in 1:5){ 
       
       A = conc[colnames(comp[[j]])]
       
       if (j == 4 ){
         B =  comp[[j]]
         
         index.morf  = which(unlist(strsplit(as.character(comp[[j]]$MORF.f), "[/]")) ==  as.character(comp[[1]]$MORF.f))
         if (length(index.morf) == 0) {B$MORF.f = "00"} else {B$MORF.f =  ifelse(index.morf == 1, "1", "2") }
       } else {B = comp[[j]]}
       
       A <- apply(A, 2, function(x){as.numeric(as.character(x))})
       B <- apply(B, 2, function(x){as.numeric(as.character(x))})
       
       
       m[j] <- length(A) - sum(A==B) 
       
     }
     
     
     mu1 <- 0.26
     mu2 <- 0.88
     
     
     #mu1 <- par.grid[h ,"mu1"]
     #mu2 <- par.grid[h ,"mu2"] 
     
     
     theta <- c(mu1, mu2)
     
     theta_est <- vector()
     theta_est2 <- vector()
     
     prior1 <- round(error_antes[i], 2)
     prior2 <- 1 - prior1
     
     
     for (u in 2:5){
       
       Likxprior1 = dbinom(x=m[u],size=length(comp[[u]]), p=theta[1])*prior1
       Likxprior2 = dbinom(x=m[u],size=length(comp[[u]]), p=theta[2])*prior2
       marginal = Likxprior1 + Likxprior2
       
       post_1 = Likxprior1 / marginal
       post_2 = Likxprior2 / marginal
       
       posterior <- c(post_1, post_2)
       theta_est[u] = which.max(posterior)
       
       theta_est2[u] = max(posterior)
     }
     
     tabla <- data.frame(error = 0, correcto = 0)
     tabla[, as.numeric(names(table(theta_est[-1])))] <- table(theta_est[-1])
     
     
     ifelse(tabla[,1]==tabla[,2], E <- mean(theta_est2[2:5]), E <- theta[which.max(tabla)]) 
     
     
     
     initialN <- c(1 - E, E)
     
     out <- ode(y = initialN, times = 1:t, func = rd, parms = parms)
     
     # 1 = error, 2 = correcto 
     res[[l]][i] <- which.max(out[t,c("1","2")])
     
   }
   
 }
 
 
 
 # resultados
 
 sim  = rbind(prop.table(table(res[[1]])), prop.table(table(res[[2]])),
              prop.table(table(res[[3]])), prop.table(table(res[[4]])))
 real  = rbind(prop.table(table(datos[datos$ID==1, "RES_BIN"])),
               prop.table(table(datos[datos$ID==2, "RES_BIN"])),
               prop.table(table(datos[datos$ID==3, "RES_BIN"])),
               prop.table(table(datos[datos$ID==4, "RES_BIN"]))
 )
 
 
 res_df <- data.frame(error.sim = sim[,1], error.real = real[,2],
                      correcto.sim = sim[,2], correcto.real = real[,1])
 
 row.names(res_df) <- c("SONIA","NATI","JAKO","MIRKA")
 
 
 
  
 
# resultados 
 
 error_por_sesiones_sim <- list()
 error_por_sesiones_real <- list()
 
 for (i in 1:4){
 
 error_por_sesiones_sim[[i]] <-
 
  tapply(ifelse(as.numeric(as.character(res[[i]]))==1,1,0), datos[datos$ID==i,"SESION"], sum)
 
 error_por_sesiones_real[[i]] <-  
 
    tapply(as.numeric(as.character(datos[datos$ID==i, "RES_BIN"])), datos[datos$ID==i,"SESION"], sum)

 }

 
 dist_error_sesion <- c(
   TSdist::MinkowskiDistance(as.vector(error_por_sesiones_sim[[1]][-c(13,14)]),as.vector(error_por_sesiones_real[[1]][-c(13,14)]), p = 3),
   TSdist::MinkowskiDistance(as.vector(error_por_sesiones_sim[[2]]),as.vector(error_por_sesiones_real[[2]]), p = 3),
   TSdist::MinkowskiDistance(as.vector(error_por_sesiones_sim[[3]]),as.vector(error_por_sesiones_real[[3]]), p = 3),
   TSdist::MinkowskiDistance(as.vector(error_por_sesiones_sim[[4]][-c(13,14)]),as.vector(error_por_sesiones_real[[4]][-c(13,14)]), p = 3)
 )
 
  
 NN <- tapply(datos$RES_BIN, datos$ID, length)
 
 prop_test_error <- c(
 prop.test(c(res_df[1,1]*NN[1], res_df[1,2]*NN[1]), rep(NN[1],2))$p.value,
 prop.test(c(res_df[2,1]*NN[2], res_df[2,2]*NN[2]), rep(NN[2],2))$p.value,
 prop.test(c(res_df[3,1]*NN[3], res_df[3,2]*NN[3]), rep(NN[3],2))$p.value,
 prop.test(c(res_df[4,1]*NN[4], res_df[4,2]*NN[4]), rep(NN[4],2))$p.value
 )
 
 resultados_df <- data.frame(res_df[,-c(3,4)], p.value = prop_test_error, dist = dist_error_sesion)
 
#  xtable::xtable(resultados_df, digits = 3,
#                 caption = "Resultados de la simulación: 
# (i) proporciones entre errores reales y simulados; 
# (ii)  p-valores de un test para diferencia entre dichas proporciones;
# (iii) distancia de Minkowski entre los errores simulados y reales por sesiones.")
#  
 
 
 # plot por sesiones
 
 sesion.tot <- list()
 
 for (i in 1:4){
sesion.tot[[i]]  <- tapply(datos[datos$ID==i,"RES_BIN"], datos[datos$ID==i,"SESION"], length)
 }
 
 df.plot.1 <- data.frame(
 SESIONES = c(1:12, 1:12), 
 ERROR.PROP = c(error_por_sesiones_sim[[1]][1:12]/sesion.tot[[1]][1:12]
           , error_por_sesiones_real[[1]][1:12]/sesion.tot[[1]][1:12]),
 CLASE = c(rep("SONIA.SIM", 12), rep("SONIA.REAL", 12))
 )
 
p.1 <- ggplot(data = df.plot.1, aes(x=SESIONES, y=ERROR.PROP, group=CLASE, color=CLASE, shape=CLASE)) +
   geom_point() +   theme_bw() + ggtitle("Errores reales y simulados (SONIA)") +
   scale_color_manual(values=c("black", "black")) +
   scale_x_continuous(breaks = 1:12)
 
 
df.plot.2 <- data.frame(
  SESIONES = c(1:14, 1:14), 
  ERROR.PROP = c(error_por_sesiones_sim[[2]]/sesion.tot[[2]]
                 , error_por_sesiones_real[[2]]/sesion.tot[[2]]),
  CLASE = c(rep("NATI.SIM", 14), rep("NATI.REAL", 14))
  
)

p.2 <- ggplot(data = df.plot.2, aes(x=SESIONES, y=ERROR.PROP, group=CLASE, color=CLASE, shape=CLASE)) +
  geom_point() +   theme_bw() + ggtitle("Errores reales y simulados (NATI)") +
  scale_color_manual(values=c("black", "black")) +
  scale_x_continuous(breaks = 1:14)


df.plot.3 <- data.frame(
  SESIONES = c(1:14, 1:14), 
  ERROR.PROP = c(error_por_sesiones_sim[[3]]/sesion.tot[[3]]
                 , error_por_sesiones_real[[3]]/sesion.tot[[3]]),
  CLASE = c(rep("JAKO.SIM", 14), rep("JAKO.REAL", 14))
  
)

p.3 <- ggplot(data = df.plot.3, aes(x=SESIONES, y=ERROR.PROP, group=CLASE, color=CLASE, shape=CLASE)) +
  geom_point() +   theme_bw() + ggtitle("Errores reales y simulados (JAKO)") +
  scale_color_manual(values=c("black", "black")) +
  scale_x_continuous(breaks = 1:14)


df.plot.4 <- data.frame(
  SESIONES = c(1:12, 1:12), 
  ERROR.PROP = c(error_por_sesiones_sim[[4]][1:12]/sesion.tot[[4]][1:12]
                 , error_por_sesiones_real[[4]][1:12]/sesion.tot[[4]][1:12]),
  CLASE = c(rep("MIRKA.SIM", 12), rep("MIRKA.REAL", 12))
)

p.4 <- ggplot(data = df.plot.4, aes(x=SESIONES, y=ERROR.PROP, group=CLASE, color=CLASE, shape=CLASE)) +
  geom_point() +   theme_bw() + ggtitle("Errores reales y simulados (MIRKA)") +
  scale_color_manual(values=c("black", "black")) +
  scale_x_continuous(breaks = 1:12)

png("simulacion_ev_games_sesiones.png")
gridExtra::grid.arrange(p.1, p.2, p.3, p.4)
dev.off()


 
# plot ejemplo ("los estudio")
# generar out con l = 1 , i = 80 

lv.data <- dif.data[dif.data$ID==l,]
error_antes <-  (exp(lv.data$CUMRES)-1) / (lv.data$TIME-1)
error_antes[is.na(error_antes)] <- 0


  
conc = data.frame(lv.data[i,-c(1:4)])
comp = list(conc, error.genero, error.epente, error.plural, error.mixto)
  
  
  m = rep(0,5)  
  
  for (j in 1:5){ 
    
    A = conc[colnames(comp[[j]])]
    
    if (j == 4 ){
      B =  comp[[j]]
      
      index.morf  = which(unlist(strsplit(as.character(comp[[j]]$MORF.f), "[/]")) ==  as.character(comp[[1]]$MORF.f))
      if (length(index.morf) == 0) {B$MORF.f = "00"} else {B$MORF.f =  ifelse(index.morf == 1, "1", "2") }
    } else {B = comp[[j]]}
    
    A <- apply(A, 2, function(x){as.numeric(as.character(x))})
    B <- apply(B, 2, function(x){as.numeric(as.character(x))})
    
    
    m[j] <- length(A) - sum(A==B) 
    
  }
  
  
  mu1 <- 0.26
  mu2 <- 0.88
  
  
  #mu1 <- par.grid[h ,"mu1"]
  #mu2 <- par.grid[h ,"mu2"] 
  
  
  theta <- c(mu1, mu2)
  
  theta_est <- vector()
  theta_est2 <- vector()
  
  prior1 <- round(error_antes[i], 2)
  prior2 <- 1 - prior1
  
  
  for (u in 2:5){
    
    Likxprior1 = dbinom(x=m[u],size=length(comp[[u]]), p=theta[1])*prior1
    Likxprior2 = dbinom(x=m[u],size=length(comp[[u]]), p=theta[2])*prior2
    marginal = Likxprior1 + Likxprior2
    
    post_1 = Likxprior1 / marginal
    post_2 = Likxprior2 / marginal
    
    posterior <- c(post_1, post_2)
    theta_est[u] = which.max(posterior)
    
    theta_est2[u] = max(posterior)
  }
  
  tabla <- data.frame(error = 0, correcto = 0)
  tabla[, as.numeric(names(table(theta_est[-1])))] <- table(theta_est[-1])
  
  
  ifelse(tabla[,1]==tabla[,2], E <- mean(theta_est2[2:5]), E <- theta[which.max(tabla)]) 
  
  
  
  initialN <- c(1 - E, E)
  
  out <- ode(y = initialN, times = 1:t, func = rd, parms = parms)
  

 df.plot <- data.frame(TIME = rep(1:(t),2), 
                       ACT = c(out[,"1"], out[,"2"]),
                       CLASE = c(rep("error", t), rep("correcto",t))
 )
 
 
 
 ggplot(data = df.plot, aes(x=TIME, y=ACT, group=CLASE, color=CLASE)) +
   geom_line(aes(linetype=CLASE)) +   theme_bw() + ggtitle("Flujo en el tiempo") +
   xlab("Tiempo") + ylab("Activación") +
   scale_color_manual(values=c("black", "black")) +
   scale_linetype_manual(values=c("dashed", "dotted"))
 
 
 
 # espacio de fase:

 library(EvolutionaryGames) 
 
 A = matrix(c(0, 0, 0, 0),2,2)
 A[1,1] =  3
 A[2,2] =  1
 ESS(A, c("X1", "X2"))
 phaseDiagram2S(A, Replicator, strategies = c("1", "2"))
 
 
  
#####################
##GSC: Gradient Symbolic Computation

##funcion adaptada para R a partir de:  
##https://cocalc.com/projects/18c77389-a5c3-49de-946a-7593b53d3fb2/
 
 #install.packages("pracma")
 #install.packages("rlist")

 library("pracma")
 library("rlist")

params = c(      # parameters in conceptual coordinates
  'W' = 0,     # (self-)weight
  'b' =  0,     # bias
  'z' = 0,     # bowl center
  'beta' = 4   # bowl strength
)


harmony1D <- function(x, ex, q, params){
  # retornar valor de función de armonía evaluado en el estado x
  
H_0 = 0.5 * x * params['W'] * x + params['b'] * x + ex * x
H_1 = -(params['beta']/2) * (x-params['z'])^2
Q_0 = -(x-1)^2 * (x+1)^2
return( (H_0 + H_1) + q * Q_0 )

}

harmonyGrad1D <- function(x, ex, q, params){
  # Retornar el valor del gradiente de la función de armonía en x
  
dH_0 = params['W'] * x + params['b'] + ex
dH_1 = -params['beta'] * (x-params['z'])
dQ_0 = -4 * x * (x-1) * (x+1)
return ( (dH_0 + dH_1) + q * dQ_0 )

}



 run1D <- function(tspan, x_init, params ,opts=NULL){
 # tspan: list [t_init, t_end]

   
  opts0 = c(
    'T_decay_rate' = 0,
    'T_init' = 0.01,
    'q_init' = 0,
    'dq' = 1,
    'max_dt' = 0.05,
    'ex' = 0
  )

  
if  (is.null(opts)==FALSE) {
  
 for (i in 1:length(opts)){
   
  index <- which(names(opts)[i]==names(opts0))
  opts0[index] <- opts[i]
 }  
  
  } 
  
  opts = opts0
  
if (is.null(x_init)){
  x = 0.1 * np.random.randn()} else {x = x_init}
  
  
t = tspan[1]
q = opts['q_init']
T = opts['T_init']
ex = opts['ex']

t_trace = vector()
x_trace = vector()
h_trace = vector()
T_trace = vector()
q_trace = vector()
surf_trace = list()
 
xvals = pracma::linspace(-1.2, 1.2, 100)

# Log estado inicial

x_trace  <- append(x_trace,x) 
h_trace <- append(h_trace, harmony1D(x, ex, q, params))
T_trace <- append(T_trace, T)
t_trace  <- append(t_trace,t) 
q_trace  <- append(q_trace,q) 
surf_trace  <- rlist::list.append(surf_trace, harmony1D(xvals, ex, q, params)) 



while (t <= tspan[2]){
  
  grad = harmonyGrad1D(x, ex, q, params)

# Update step size to prevent overflow
  
if (grad != 0){
  dt = min(opts['max_dt'], opts['max_dt']/abs(grad))}else{
  dt = opts['max_dt']}

# Update state
  
x = x + dt * grad + sqrt(2*T*dt) * rnorm(1)

# update temperatura (exponential decay)

T = exp(-opts['T_decay_rate'] * dt) * T

# update time

t = t + dt

# update q

q = q + opts['dq'] * dt

# Log

x_trace  <- append(x_trace,x) 
h_trace <- append(h_trace, harmony1D(x, ex, q, params))
T_trace <- append(T_trace, T)
t_trace  <- append(t_trace,t) 
q_trace  <- append(q_trace,q) 
surf_trace  <- rlist::list.append(surf_trace, harmony1D(xvals, ex, q, params))

}

len <- length(x_trace)

return( list(surf.trace = surf_trace[-len], x.trace = x_trace[-len], 
             h.trace = h_trace[-len],  T.trace = T_trace[-len], 
             t.trace = t_trace[-len],  q.trace = q_trace[-len]) )

 }
 
 
mm <- seq(0,1,by=0.01)
escala <- scales::rescale(mm, c(-1,1))
 

 error_por_sesiones_real <- list()
 
 for (z in 1:4){
   
   error_por_sesiones_real[[z]] <-  
     
     tapply(as.numeric(as.character(datos[datos$ID==z, "RES_BIN"])), datos[datos$ID==z,"SESION"], sum)
   
 }
 
 
 
# Establecer el valor de  "tresh" a partir del cual se declara "error"
 
 
 tresh = seq(0.85, 0.95, by = 0.01)
 
 df.res <- list()
 res.loss <- vector()
 
 vector.dist <- list()
 res.out <- list()
 
 

 for (h in 1:length(tresh)) {

   print(h) 
   res <- list(vector(), vector(), vector(), vector())
   
   for  (l in 1:4){
     
     lv.data <- dif.data[dif.data$ID==l,]
     error_antes <-  (exp(lv.data$CUMRES)-1) / (lv.data$TIME-1)
     error_antes[is.na(error_antes)] <- 0.01
     error_antes[error_antes==0] <- 0.01
     
     for (i in 1:dim(lv.data)[1]){
       
       conc = data.frame(lv.data[i,-c(1:4)])
       comp = list(conc, error.genero, error.epente, error.plural, error.mixto)
       
       
       m = rep(0,5)  
       
       for (j in 1:5){ 
         
         A = conc[colnames(comp[[j]])]
         
         if (j == 4 ){
           B =  comp[[j]]
           
           index.morf  = which(unlist(strsplit(as.character(comp[[j]]$MORF.f), "[/]")) ==  as.character(comp[[1]]$MORF.f))
           if (length(index.morf) == 0) {B$MORF.f = "00"} else {B$MORF.f =  ifelse(index.morf == 1, "1", "2") }
         } else {B = comp[[j]]}
         
         A <- apply(A, 2, function(x){as.numeric(as.character(x))})
         B <- apply(B, 2, function(x){as.numeric(as.character(x))})
         
         m[j] <- philentropy::distance(rbind(A,B), method = "jaccard") 
         
         
         
       }
       
        
        
        
  
 
 
  KK    <- crch::qtlogis(round(exp(-((median(m[2:5]))^2/(2*(error_antes[i])))),2), left = -1, right = 1)
 
              
   out <- run1D(tspan = c(0,5), x_init = -1, params = params ,
   opts= c(ex = KK))
       
    x.trace.last   <- out$x.trace[length(out$x.trace)]
       
       #res[[l]][i] <- ifelse(x.trace.last>0.89,1,0) 
       res[[l]][i] <- ifelse(x.trace.last>tresh[h],1,0) 
       
    
     }
     
   }
    
   # resultados
   
   res.out[[h]] <- res
   
   sim  = rbind(prop.table(table(res[[1]])), prop.table(table(res[[2]])),
                prop.table(table(res[[3]])), prop.table(table(res[[4]])))
   real  = rbind(prop.table(table(datos[datos$ID==1, "RES_BIN"])),
                 prop.table(table(datos[datos$ID==2, "RES_BIN"])),
                 prop.table(table(datos[datos$ID==3, "RES_BIN"])),
                 prop.table(table(datos[datos$ID==4, "RES_BIN"]))
   )
   
   
res_df <- data.frame(error.sim = sim[,2], error.real = real[,2],
                        correcto.sim = sim[,1], correcto.real = real[,1])
   
   row.names(res_df) <- c("SONIA","NATI","JAKO","MIRKA")
   
   
   df.res[[h]]  <- res_df

   
     
   res.loss[h] <-  sum((df.res[[h]][,1] - df.res[[h]][,2])^2)*1/4
   
   
  
   
 }
  
# h ("tresh") minimo es h = 5 --> s = 0.89
 
 df.res[[which.min(res.loss)]]
 
 tresh[5]
 
 # 50 simulaciones con s = 0.89
 
 
 df.res.sim <- list()
 res.out.sim <- list()
 

 
 
 
 n_sim = 50
 
 for (v in 1:n_sim) {
   
   print(v) 
   res.sim <- list(vector(), vector(), vector(), vector())
   
   for  (l in 1:4){
     
     lv.data <- dif.data[dif.data$ID==l,]
     error_antes <-  (exp(lv.data$CUMRES)-1) / (lv.data$TIME-1)
     error_antes[is.na(error_antes)] <- 0.01
     error_antes[error_antes==0] <- 0.01
     
     for (i in 1:dim(lv.data)[1]){
       
       conc = data.frame(lv.data[i,-c(1:4)])
       comp = list(conc, error.genero, error.epente, error.plural, error.mixto)
       
       
       m = rep(0,5)  
       
       for (j in 1:5){ 
         
         A = conc[colnames(comp[[j]])]
         
         if (j == 4 ){
           B =  comp[[j]]
  
           index.morf  = which(unlist(strsplit(as.character(comp[[j]]$MORF.f), "[/]")) ==  as.character(comp[[1]]$MORF.f))
           if (length(index.morf) == 0) {B$MORF.f = "00"} else {B$MORF.f =  ifelse(index.morf == 1, "1", "2") }
         } else {B = comp[[j]]}
         
         A <- apply(A, 2, function(x){as.numeric(as.character(x))})
         B <- apply(B, 2, function(x){as.numeric(as.character(x))})
         
         m[j] <- philentropy::distance(rbind(A,B), method = "jaccard") 
         
        
         
       }
       
      
       
       
       KK    <- crch::qtlogis(round(exp(-((median(m[2:5]))^2/(2*(error_antes[i])))),2), left = -1, right = 1)
       
       
       
       out <- run1D(tspan = c(0,5), x_init = -1, params = params ,
                    opts= c(ex = KK))
       
       x.trace.last   <- out$x.trace[length(out$x.trace)]
       
       res.sim[[l]][i] <- ifelse(x.trace.last>0.89,1,0) 
       
     }
     
   }
   
   # resultados
   
   res.out.sim[[v]] <- res.sim
   
   sim  = rbind(prop.table(table(res.sim[[1]])), prop.table(table(res.sim[[2]])),
                prop.table(table(res.sim[[3]])), prop.table(table(res.sim[[4]])))
   real  = rbind(prop.table(table(datos[datos$ID==1, "RES_BIN"])),
                 prop.table(table(datos[datos$ID==2, "RES_BIN"])),
                 prop.table(table(datos[datos$ID==3, "RES_BIN"])),
                 prop.table(table(datos[datos$ID==4, "RES_BIN"]))
   )
   
   
   res_df <- data.frame(error.sim = sim[,2], error.real = real[,2],
                        correcto.sim = sim[,1], correcto.real = real[,1])
   
   row.names(res_df) <- c("SONIA","NATI","JAKO","MIRKA")
   
   
   df.res.sim[[v]]  <- res_df
   
 }
 
   
 # grabar los resultados de las 50 simulaciones
 #res_gsc <- list(df.res.sim, res.out.sim)
 
 #save(res_gsc, file = "res_gsc.RData")
 
 load("res_gsc.RData")
 
 
 # Patrón por sesiones (media de las 50 simulaciones)
 
 dist.sonia <- vector()
 dist.nati <- vector()
 dist.jako <- vector()
 dist.mirka <- vector()

 res.out <- res_gsc[[2]]
 
  
 for (k in 1:length(res.out)){ 
 
 error_por_sesiones_sim <- list()
 error_por_sesiones_real <- list()
 
 res <- res.out[[k]]
 
 for (i in 1:4){
   
   error_por_sesiones_sim[[i]] <-
     
     tapply(ifelse(as.numeric(as.character(res[[i]]))==1,1,0), datos[datos$ID==i,"SESION"], sum)
   
   error_por_sesiones_real[[i]] <-  
     
     tapply(as.numeric(as.character(datos[datos$ID==i, "RES_BIN"])), datos[datos$ID==i,"SESION"], sum)
   
 }
 

dist.sonia[k]  <- TSdist::MinkowskiDistance(as.vector(error_por_sesiones_sim[[1]][-c(13,14)]),as.vector(error_por_sesiones_real[[1]][-c(13,14)]), p = 3)
dist.nati[k] <-  TSdist::MinkowskiDistance(as.vector(error_por_sesiones_sim[[2]]),as.vector(error_por_sesiones_real[[2]]), p = 3)
dist.jako[k] <-  TSdist::MinkowskiDistance(as.vector(error_por_sesiones_sim[[3]]),as.vector(error_por_sesiones_real[[3]]), p = 3) 
dist.mirka[k]  <- TSdist::MinkowskiDistance(as.vector(error_por_sesiones_sim[[4]][-c(13,14)]),as.vector(error_por_sesiones_real[[4]][-c(13,14)]), p = 3)
    
 }

 dist_error_sesion_mean <- c(
 mean(dist.sonia),
 mean(dist.nati),
 mean(dist.jako),
 mean(dist.mirka)
 )
 
 
 # Patrón de error global 
 
 
 NN <- tapply(datos$RES_BIN, datos$ID, length)
 

p.valor.sonia <- vector()
p.valor.nati <- vector()
p.valor.jako <- vector()
p.valor.mirka <- vector()

df.res <- res_gsc[[1]]

for (k in 1:length(df.res)){ 

 resp_df <- df.res[[k]]  
  
p.valor.sonia[k] <- prop.test(c(resp_df[1,1]*NN[1], resp_df[1,2]*NN[1]), rep(NN[1],2))$p.value
p.valor.nati[k] <- prop.test(c(resp_df[2,1]*NN[2], resp_df[2,2]*NN[2]), rep(NN[2],2))$p.value
p.valor.jako[k] <- prop.test(c(resp_df[3,1]*NN[3], resp_df[3,2]*NN[3]), rep(NN[3],2))$p.value
p.valor.mirka[k] <- prop.test(c(resp_df[4,1]*NN[4], resp_df[4,2]*NN[4]), rep(NN[4],2))$p.value
 
}
 
prop_test_error_mean <- c(
sum(p.valor.sonia<0.05)/length(df.res), 
sum(p.valor.nati<0.05)/length(df.res), 
sum(p.valor.jako<0.05)/length(df.res), 
sum(p.valor.mirka<0.05)/length(df.res) 
)

mean_error_sim <- rowMeans(sapply(1:50, function(k) {df.res[[k]][,1]}))
se_error_sim <- apply(sapply(1:50, function(k) {df.res[[k]][,1]}),1, function(x){sd(x)/sqrt(50)})
error_real <- df.res[[1]][,2]


# Resultados 

resultados_df <- data.frame(error.sim.mean = round(mean_error_sim,4) , 
error.sim.se = round(se_error_sim,4),
error.real = round(error_real,4),
p.value.freq = prop_test_error_mean, 
dist.mean = round(dist_error_sesion_mean,4))

# xtable::xtable(resultados_df, digits = 4,
#                caption = "Resultados de la simulación: 
# (i) proporciones de errores reales y proporciones medias de simulados; 
# (ii) test para diferencia entre dichas proporciones: frecuencias de p.valores < 0.05;
# (iii) distancia media de Minkowski entre los errores simulados y reales por sesiones.")


 # plot por sesiones
 
 sesion.tot <- list()
 
 for (i in 1:4){
   sesion.tot[[i]]  <- tapply(datos[datos$ID==i,"RES_BIN"], datos[datos$ID==i,"SESION"], length)
 }
 

 res.out <- res_gsc[[2]]
 
 error_por_sesiones_real <- list()
 
 for (i in 1:4){
 error_por_sesiones_real[[i]] <-  
   
   tapply(as.numeric(as.character(datos[datos$ID==i, "RES_BIN"])), datos[datos$ID==i,"SESION"], sum)
 }

     error_por_sesiones_sim_list <- list()   
     
   for (k in 1:length(res.out)){ 
       
       error_por_sesiones_sim <- list()
       res <- res.out[[k]]
       
   for (i in 1:4){
     
     error_por_sesiones_sim[[i]] <-
       
       tapply(ifelse(as.numeric(as.character(res[[i]]))==1,1,0), datos[datos$ID==i,"SESION"], sum)/sesion.tot[[i]]
     
     
     }
   
   
  error_por_sesiones_sim_list[[k]]   <-  error_por_sesiones_sim
   
 }
 
 error_por_sesiones_simulacion_sonia_mean <-
   
  rowMeans(sapply(1:50, function(x){error_por_sesiones_sim_list[[x]][[1]]})[-(13:14),])
 
 error_por_sesiones_simulacion_nati_mean <-
   
   rowMeans(sapply(1:50, function(x){error_por_sesiones_sim_list[[x]][[2]]}))
 
 error_por_sesiones_simulacion_jako_mean <-
   
   rowMeans(sapply(1:50, function(x){error_por_sesiones_sim_list[[x]][[3]]}))
 
 error_por_sesiones_simulacion_mirka_mean <-
   
   rowMeans(sapply(1:50, function(x){error_por_sesiones_sim_list[[x]][[4]]})[-(13:14),])
 
 
 
 error_por_sesiones_simulacion_sonia_se <-
   
   matrixStats::rowSds(sapply(1:50, function(x){error_por_sesiones_sim_list[[x]][[1]]})[-(13:14),])/sqrt(50)
 
 error_por_sesiones_simulacion_nati_se <-
   
   matrixStats::rowSds(sapply(1:50, function(x){error_por_sesiones_sim_list[[x]][[2]]}))/sqrt(50)
 
 error_por_sesiones_simulacion_jako_se <-
   
   matrixStats::rowSds(sapply(1:50, function(x){error_por_sesiones_sim_list[[x]][[3]]}))/sqrt(50)
 
 error_por_sesiones_simulacion_mirka_se <-
   
   matrixStats::rowSds(sapply(1:50, function(x){error_por_sesiones_sim_list[[x]][[4]]})[-(13:14),])/sqrt(50)
 
 
 
 
 df.plot.1 <- data.frame(
   SESIONES = c(1:12, 1:12), 
   ERROR.PROP.MEAN = c(error_por_sesiones_simulacion_sonia_mean
                  , error_por_sesiones_real[[1]][1:12]/sesion.tot[[1]][1:12]),
   ERROR.PROP.SE = c(error_por_sesiones_simulacion_sonia_se
                       , rep(NA,12)),
   CLASE = c(rep("SONIA.SIM", 12), rep("SONIA.REAL", 12))
 )
 pd <- position_dodge(0.1)
 
 
 p.1 <- ggplot(data = df.plot.1, aes(x=SESIONES, y=ERROR.PROP.MEAN, group=CLASE, color=CLASE)) +
   theme_bw() + ggtitle("Errores reales y simulados (SONIA)") +
   geom_errorbar(aes(ymin=ERROR.PROP.MEAN-ERROR.PROP.SE, ymax=ERROR.PROP.MEAN+ERROR.PROP.SE), colour="black", width=1, position=pd) +
   geom_line(position=pd, aes(linetype=CLASE)) +
   geom_point(position=pd, size=1)+
   scale_color_manual(values=c("black", "black")) +
   scale_x_continuous(breaks = 1:12) +
   scale_linetype_manual(values=c("dashed", "dotted"))
 
 df.plot.2 <- data.frame(
   SESIONES = c(1:14, 1:14), 
   ERROR.PROP.MEAN = c(error_por_sesiones_simulacion_nati_mean
                  , error_por_sesiones_real[[2]]/sesion.tot[[2]]),
   ERROR.PROP.SE = c(error_por_sesiones_simulacion_nati_se
                     , rep(NA,14)),
   CLASE = c(rep("NATI.SIM", 14), rep("NATI.REAL", 14))
   
 )
 
 p.2 <- ggplot(data = df.plot.2, aes(x=SESIONES, y=ERROR.PROP.MEAN, group=CLASE, color=CLASE)) +
   theme_bw() + ggtitle("Errores reales y simulados (NATI)") +
   geom_errorbar(aes(ymin=ERROR.PROP.MEAN-ERROR.PROP.SE, ymax=ERROR.PROP.MEAN+ERROR.PROP.SE), colour="black", width=1, position=pd) +
   geom_line(position=pd, aes(linetype=CLASE)) +
   geom_point(position=pd, size=1)+
   scale_color_manual(values=c("black", "black")) +
   scale_x_continuous(breaks = 1:14) +
   scale_linetype_manual(values=c("dashed", "dotted"))
 
 df.plot.3 <- data.frame(
   SESIONES = c(1:14, 1:14), 
   ERROR.PROP.MEAN = c(error_por_sesiones_simulacion_jako_mean,
                       error_por_sesiones_real[[3]]/sesion.tot[[3]]),
   ERROR.PROP.SE = c(error_por_sesiones_simulacion_jako_se
                     , rep(NA,14)),
   CLASE = c(rep("JAKO.SIM", 14), rep("JAKO.REAL", 14))
   
 )
 
 p.3 <- ggplot(data = df.plot.3, aes(x=SESIONES, y=ERROR.PROP.MEAN, group=CLASE, color=CLASE)) +
      theme_bw() + ggtitle("Errores reales y simulados (JAKO)") +
   geom_errorbar(aes(ymin=ERROR.PROP.MEAN-ERROR.PROP.SE, ymax=ERROR.PROP.MEAN+ERROR.PROP.SE), colour="black", width=1, position=pd) +
   geom_line(position=pd, aes(linetype=CLASE)) +
   geom_point(position=pd, size=1)+
   scale_color_manual(values=c("black", "black")) +
   scale_x_continuous(breaks = 1:14) +
   scale_linetype_manual(values=c("dashed", "dotted"))
 
 df.plot.4 <- data.frame(
   SESIONES = c(1:12, 1:12), 
   ERROR.PROP.MEAN = c(error_por_sesiones_simulacion_mirka_mean
                       , error_por_sesiones_real[[1]][1:12]/sesion.tot[[1]][1:12]),
   ERROR.PROP.SE = c(error_por_sesiones_simulacion_mirka_se
                     , rep(NA,12)),
   CLASE = c(rep("MIRKA.SIM", 12), rep("MIRKA.REAL", 12))
 )
 
 
 
 p.4 <- ggplot(data = df.plot.4, aes(x=SESIONES, y=ERROR.PROP.MEAN, group=CLASE, color=CLASE)) +
   theme_bw() + ggtitle("Errores reales y simulados (MIRKA)") +
   geom_errorbar(aes(ymin=ERROR.PROP.MEAN-ERROR.PROP.SE, ymax=ERROR.PROP.MEAN+ERROR.PROP.SE), colour="black", width=1, position=pd) +
   geom_line(position=pd, aes(linetype=CLASE)) +
   geom_point(position=pd, size=1)+
   scale_color_manual(values=c("black", "black")) +
   scale_x_continuous(breaks = 1:12) +
   scale_linetype_manual(values=c("dashed", "dotted"))
 
 
 png("gsc_sesiones.png")
 gridExtra::grid.arrange(p.1, p.2, p.3, p.4)
 dev.off()
 

 # ejemeplo con e = 1
 
 valores_de_x <-  pracma::linspace(-1.2, 1.2, 100)
 
 
 out <- run1D(tspan = c(0,5), x_init = -1, params = params ,
              opts= c(ex = 1))
 
 c = c(0, 2.76, 3.97, 4.96)
 
 for (k in 1:4){
   
   i = which(round(out$q.trace,2)==c[k])
   
   out_df <-      data.frame(xvals = valores_de_x,
                             yvals = out$surf.trace[[i]], 
                             x.trace = round(rep(out$x.trace[i],100),2) ,
                             h.trace = round(rep(out$h.trace[i],100),2),
                             q.trace = round(rep(out$q.trace[i],100),2),
                             t.trace = round(rep(out$t.trace[i],100),2),
                             T.trace = round(rep(out$T.trace[i],100),2)
   ) 
   
   
   p <- ggplot(out_df, aes(x=xvals, y=yvals)) +
     geom_line(color = "black") +
     geom_point(aes(x= x.trace, y=h.trace), colour="black") +
     geom_text(aes(label = paste0("t: ", t.trace), x = -1, y = 0)) +
     geom_text(aes(label = paste0("T: ", T.trace), x = -1, y = 0.5)) +
     geom_text(aes(label = paste0("q: ", q.trace), x = -1, y = 1)) +
     labs(title = 'Ascenso de gradiente',x = 'x', y = 'H(x)')
   
   assign(paste0("p",k, sep=""), p)
   
 }
 
 png("gsc_ejemplo_con_e_1.png")
 gridExtra::grid.arrange(p1, p2, p3, p4)
 dev.off() 
 
 
 
 ## ejemplo : "muchas lugares" [SONIA, sesión 4, línea 66]
 
 
 valores_de_x <-  pracma::linspace(-1.2, 1.2, 100)
 
 
 out <- run1D(tspan = c(0,5), x_init = -1, params = params ,
              opts= c(ex = 0.32))
 
 c = c(0, 1.32, 2.22, 4.93)
 
 for (k in 1:4){
   
   i = which(round(out$q.trace,2)==c[k])
   
   out_df <-      data.frame(xvals = valores_de_x,
                               yvals = out$surf.trace[[i]], 
                               x.trace = round(rep(out$x.trace[i],100),2) ,
                               h.trace = round(rep(out$h.trace[i],100),2),
                               q.trace = round(rep(out$q.trace[i],100),2),
                               t.trace = round(rep(out$t.trace[i],100),2),
                               T.trace = round(rep(out$T.trace[i],100),2)
                    ) 
   
   
 p <- ggplot(out_df, aes(x=xvals, y=yvals)) +
   geom_line(color = "blue") +
   geom_point(aes(x= x.trace, y=h.trace), colour="red") +
   geom_text(aes(label = paste0("t: ", t.trace), x = -1, y = 0)) +
   geom_text(aes(label = paste0("T: ", T.trace), x = -1, y = 0.5)) +
   geom_text(aes(label = paste0("q: ", q.trace), x = -1, y = 1)) +
   labs(title = 'Ascenso de gradiente',x = 'x', y = 'H(x)')
 
 assign(paste0("p",k, sep=""), p)
 
 }
 
 png("gsc_ejemplo_muchas_lugares_SONIA_i_56.png")
 gridExtra::grid.arrange(p1, p2, p3, p4)
 dev.off()
 
 
 ### animación de "muchas lugares"
 # intalar "Rtools": https://cran.r-project.org/bin/windows/Rtools/
 #devtools::install_github("thomasp85/transformr")
 library(gganimate)
 library(gifski)
 
 valores_de_x <-  pracma::linspace(-1.2, 1.2, 100)
 
 
 out <- run1D(tspan = c(0,5), x_init = -1, params = params ,
              opts= c(ex = 0.32))
 
 out_df <- data.frame()
 
 for (i in 1:length(out$x.trace)){
   
   out_df <-  rbind(out_df,  
                    data.frame(xvals = valores_de_x,
                               yvals = out$surf.trace[[i]], 
                               x.trace = round(rep(out$x.trace[i],100),2) ,
                               h.trace = round(rep(out$h.trace[i],100),2),
                               q.trace = round(rep(out$q.trace[i],100),2),
                               t.trace = round(rep(out$t.trace[i],100),2),
                               T.trace = round(rep(out$T.trace[i],100),2)
                    ) 
   )
   
 }
 
 
 
 
 p <- ggplot(out_df, aes(x=xvals, y=yvals)) +
   geom_line(color = "black") +
   geom_point(aes(x= x.trace, y=h.trace), colour="black") +
   geom_text(aes(label = paste0("t: ", t.trace), x = -1, y = 0)) +
   geom_text(aes(label = paste0("T: ", T.trace), x = -1, y = 0.5)) +
   geom_text(aes(label = paste0("q: ", q.trace), x = -1, y = 1)) +
   labs(title = 'Ascenso de gradiente',x = 'x', y = 'H(x)') +
   transition_time(q.trace) +
   ease_aes('linear')
 
 
 animate(p, duration = 10, fps = 20, width = 400, height = 400, renderer = gifski_renderer())
 
 
 #anim_save("output.gif")
 
 
 
 
 
########## gráfico de bifurcación tridente super-crítica imperfecta

library(latex2exp) 
 
# caso r < 0  
 
e <- 1
q <- 0.5 
a <- 1/(2*q^(1/2))
r <- -4+4*q # r = -2
h <- (1/a)*e # 1.414214
#hc <- ((2/3)*r)*sqrt(r/3) # INF: no hay bifuración 

x <- seq(-1,1,0.01)
u <- x/a
y1 <-r*u - u^3


g1 <- data.frame(u, y1)
#f <- function(x) r*x - x^3

p <- 
  ggplot(g1, aes(u,y1)) +  geom_line(color="red") +
  #geom_point()+
  #stat_function(fun=f, colour="red") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_hline(yintercept = -h, color = "blue") +
  geom_hline(yintercept = -3, color = "green") +
  geom_hline(yintercept = 1.5, color = "green") +
  geom_hline(yintercept = 4, color = "green") +
  geom_point(data = data.frame(x= -1.18, y = 4), aes(x,y),shape = 21, colour = "black", fill = "black", size = 3) +
  geom_point(data = data.frame(x= -0.64, y = 1.5), aes(x,y),shape = 21, colour = "black", fill = "black", size = 3) +
  geom_point(data = data.frame(x= 0.6, y = -h), aes(x,y),shape = 21, colour = "black", fill = "black", size = 3) +
  geom_point(data = data.frame(x= 1, y = -3), aes(x,y),shape = 21, colour = "black", fill = "black", size = 3) +
  geom_text(aes(label = TeX("$f\\left(u\\right)=ru-u^{3}$", output = "character"), x = -1, y = 6), color = "red", size = 4, parse = TRUE)  +  
  geom_text(aes(label = TeX("$g\\left(u\\right)=\\;-h$", output = "character"), x = 0.8, y = -h+0.5), color = "blue", size = 4, parse = TRUE)  +
  theme_bw() + ggtitle(TeX("$r<0\\;(q<1)$")) + xlab(TeX("$u$")) + ylab(TeX("$\\frac{dx}{dt}=ru-u^{3}$"))
 

model <- function(u) c(F1 = r*u - u^3 + h)

ss <- rootSolve::uniroot.all(f = model, c(min(u),max(u)))
ss # 0.5993471

r - 3*(0.5993471)^2 # -3.077651 --> estable

# caso: r > 0

#library("gghalves")
#library(extrafont)
#font_import()
#loadfonts(device = "win")

library(showtext)
# descargar tipo de fuente para símbolos
download.file("https://github.com/gearit/ttf-symbola",
              "Symbola.ttf", mode="wb")
# intalarla 
# verificar path de carpeta de fuentes
font_paths() # "C:\\Windows\\Fonts"
# Agregar fuente
font_add(family = "symbola", regular = "C:\\Windows\\Fonts\\Symbola.ttf")



e <- 1
q <- 3
a <- 1/(2*q^(1/2))
r <- -4+4*q # r = 8
h <- (1/a)*e # 3.464102
hc <- ((2/3)*r)*sqrt(r/3) # 8.709297
u_max <- sqrt(r/3) # 1.632993


x <- seq(-1,1,0.01)
u <- x/a

y1 <-r*u - u^3
g1 <- data.frame(u, y1)
#f <- function(x) r*x - x^3



p1 <- ggplot(g1, aes(u,y1)) + geom_line(color = "red") +
  #stat_function(fun=f, colour="red") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_hline(yintercept = 11, color = "green") +
  geom_hline(yintercept = -h, color = "blue") +
  geom_hline(yintercept = hc, color = "green") +
  geom_hline(yintercept = 5, color = "green") +
  geom_hline(yintercept = -hc, color = "green") +
  geom_hline(yintercept = -11, color = "green") +
  
  geom_point(data = data.frame(x= -3.35, y = 11), aes(x,y),shape = 21, colour = "black", fill = "black", size = 3) +
  
  geom_text(data = data.frame(x= u_max, y = hc), aes(x,y), label = "\u25D7", size=3) +
  geom_point(data = data.frame(x= -3.25, y = hc), aes(x,y),shape = 21, colour = "black", fill = "black", size = 3) +
  
  geom_point(data = data.frame(x= 0.7, y = 5), aes(x,y),shape = 21, colour = "black", fill = "white", size = 3) +
  geom_point(data = data.frame(x= 2.45, y = 5), aes(x,y),shape = 21, colour = "black", fill = "black", size = 3) +
  geom_point(data = data.frame(x= -3.1, y = 5), aes(x,y),shape = 21, colour = "black", fill = "black", size = 3) +
  
  geom_point(data = data.frame(x= -2.6, y = -h), aes(x,y),shape = 21, colour = "black", fill = "black", size = 3) +
  geom_point(data = data.frame(x= -0.45, y = -h), aes(x,y),shape = 21, colour = "black", fill = "white", size = 3) +
  geom_point(data = data.frame(x= 3, y = -h), aes(x,y),shape = 21, colour = "black", fill = "black", size = 3) +
  
  geom_text(data = data.frame(x= -u_max, y = -hc), aes(x,y), label = "\u25D6", size=3) +
  geom_point(data = data.frame(x= 3.25, y = -hc), aes(x,y),shape = 21, colour = "black", fill = "black", size = 3) +
  
  geom_point(data = data.frame(x= 3.35, y = -11), aes(x,y),shape = 21, colour = "black", fill = "black", size = 3) +
  
  
  geom_text(aes(label = TeX("$|h|>h_{c}$", output = "character"), x = 3, y = 11-0.8), color = "blue", size = 4, parse = TRUE)  +
  geom_text(aes(label = TeX("$|h|=h_{c}$", output = "character"), x = 3, y = hc-0.8), color = "blue", size = 4, parse = TRUE)  +
  geom_text(aes(label = TeX("$|h|<h_{c}$", output = "character"), x = 3, y = 5-0.8), color = "blue", size = 4, parse = TRUE)  +
  
  geom_text(aes(label = TeX("$f\\left(u\\right)=ru-u^{3}$", output = "character"), x = -2.5, y = 13), color = "red", size = 4, parse = TRUE)  +  
  geom_text(aes(label = TeX("$g\\left(u\\right)=\\;-h$", output = "character"), x = 1, y = -h+0.8), color = "blue", size = 4, parse = TRUE)  +
  theme_bw() + ggtitle(TeX("$r>0\\;(q>1)$")) + xlab(TeX("$u$")) + ylab(TeX("$\\frac{dx}{dt}=ru-u^{3}$"))
  
  

cairo_pdf("p1.pdf", family="symbola")
  p1
  dev.off()

  # raices (puntos fijos)
  
  model <- function(u) c(F1 = r*u - u^3 + h)
  
  ss <- rootSolve::uniroot.all(f = model, c(min(u),max(u)))
  ss
# -2.5801997 (estable) -0.4437909 (inestable)  3.0241436 (estable)

 