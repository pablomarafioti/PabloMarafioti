## ratings
## Estadística descriptiva

# recordar que hay que intslar los paquetes !!
# install.packages("languageR")

library(languageR)

# datos de tipos de ratings para 81 palabras sobre plantas y animales
# Word: token (palabra)
# Class: tipo semántico (planta o animal)
# Length: largo de la palabra (en cantidad de letras)
# Frequency: log(frecuencia) de la palabra
# meanWeightRating: rating subjetivo medio del peso del referente
# meanSizeRating: rating subjetivo medio del tamaño del referente
# meanFamiliarity: rating subjetivo medio de la familiaridad de la 
#                  palabra (frecuencia subjetiva)
# SynsetCount: logaritmo del número de Synsets (conjuntos de sinonimia) 
#              en la cual está listada en WordNet 
# FamilySize: logaritmo de número de palabras complejas en las que ocurre
# DerivEntropy: entropía derivacional de la palabra (mide complejidad morfológica) 
# Complex: si la palabra es compleja o simple
# FreqSingular: frecuencia de ocurrecias de la palabra en singular
# FreqPlural: frecuencia de ocurrecias de la palabra en plural
# ringfl:  log(FreqSingular / FreqPlural)
 
colnames(ratings) # nombres de las variables
head(ratings)

str(ratings)  # tipos de variables

#library(modeest) # paquete con función mfv() para calcular la moda
#library(moments) # paquete para calcular asimetría y curtosis

# calcular estadísticos descriptivos para la variable "Length"

table(ratings$Length) # distribución (freq. absolutas)
prop.table(table(ratings$Length)) # distribución (freq. relativas)

mean(ratings$Length) # media
median(ratings$Length) # mediana
modeest::mfv(ratings$Length) # moda
mean(ratings$Length, trim = 0.1) # media podada al 10 %
min(ratings$Length) # el mínimo
max(ratings$Length) # el máximo
range(ratings$Length) # el rango
IQR(ratings$Length) # rango intercuartil
summary(ratings$Length) # min, 1Q, mediana (2Q), media, 3Q, max
quantile(ratings$Length) # percentiles del 0%, 25%, 50%, 75%, 100% 
quantile(ratings$Length, 0.25) # percentil del 25 % = 1 cuartil
quantile(ratings$Length, 0.5) # percentil del 50 % = 2 cuartil o mediana
quantile(ratings$Length, 0.75) # percentil del 75 % = 3 cuartil
quantile(ratings$Length, 0.1) # percentil del 10 % (primer decil)
quantile(ratings$Length, 0.9) # percentil del 90 % (noveno decil)
quantile(ratings$Length, seq(0, 1, 0.1)) # percentiles de 0 % a 100 %
var(ratings$Length) # varianza
sd(ratings$Length) # desvío
mad(ratings$Length) # median absolute deviation
mad(ratings$Length)/0.6745 # MADN
moments::kurtosis(ratings$Length) # curtosis
moments::skewness(ratings$Length) # asimetría

# para todas las variables cuantitativas


psych::describe(ratings[,-c(1,5,6,10)]) # sacar las cuali


library(Rcmdr)

numSummary(ratings[,-c(1,5,6,10)], statistics=c("mean", "sd", "IQR", "quantiles", "cv", "skewness", "kurtosis"), 
                   quantiles=c(0,.25,.5,.75,1), type="2")


# calcular un estadístico segmentado por grupo (Class)


numSummary(ratings[,-c(1,5,6,10)], statistics=c("mean", "sd"), 
                 groups = ratings$Class)

psych::describeBy(ratings[,-c(1,5,6,10)], group = ratings$Class, mat=T, digits = 2)

# cruzando "class" con "complex"
aggregate(ratings[,-c(1,5,6,10)], list(ratings$Class, ratings$Complex), c("mean"))


### Usando funciones apply, tapply


apply(ratings[,-c(1,5,6,10)],2, mean)

tapply(ratings$Length, ratings$Class, mean)

# creando una función

medidas <- function(x){M  <- c(media=mean(x), 
              mediana=median(x), 
              desvio=sd(x), 
              mad = mad(x),
              q1=quantile(x,.25), 
              q3=quantile(x,.75), 
              min=min(x), 
              max=max(x),
              rango=range(x),
              iqr = IQR(x),
              curtosis  = moments::kurtosis(x),
              asimetria = moments::skewness(x)
              )
return(M)
}

medidas(ratings$Length)

moda <- function(x){
  
  return(as.numeric(names(which.max(table(x)))))
}

moda(ratings$Length)

apply(ratings[,-c(1,5,6,10)],2, medidas)

tapply(ratings$Length, ratings$Class, medidas)

## Gráficos

library(ggplot2)

# las variables "Class" y "Complex" son cualitativas 
# --> gráfico de barras 

table(ratings$Class)

ggplot(data=ratings, aes(x=Class, fill=Class)) +
       geom_bar(stat="count") +  ggtitle("Class")

table(ratings$Complex)

datos <- data.frame(Tipo = ratings$Complex)
levels(datos$Tipo) <- c("complejo", "simple")
# o bien
#datos$Tipo <- recode_factor(datos$Tipo, complex = "complejo", simplex = "simple")

ggplot(data=datos, aes(x=Tipo, fill=Tipo)) +
  geom_bar(stat="count") + ggtitle("Morfología")


# La variable "Length" es cuantitativa discreta  --> gráfico de bastones

# con frecuencia absoluta

ggplot(ratings, aes(x = factor(Length))) + 
geom_bar(color="blue", fill="green", width = 0.1) +
  ylab("Frecuencia") + xlab("Largo")

# con frecuencia relativa

conteo <- prop.table(table(ratings$Length))
df <- data.frame(conteo)

ggplot(df, aes(x=Var1, y=Freq)) +
  geom_segment( aes(xend=Var1, yend=0)) +
  geom_point( size=2, color="orange") +
  theme_bw() + 
  ylab("Freq. relativa") + xlab("Largo")

# Las variables "Frequency", "SynsetCount", "FamilySize", "DerivEntropy"
# son cuantitativas continuas --> histograma
# automáticamente divide en 30 bins, si no funciona, usar, por ejemplo:

# v = ratings$Frequency 
# n = dim(ratings)[1]
# b = round(10*log(n))
# b = round(2*sqrt(n))
# b = round(1 + log2(n)) 
# h = 3.49*sd(v)*(n^(-1/3))
# h = 2*IQR(v)*(n^(-1/3))

n = dim(ratings)[1]

b = round(1 + log2(n))

ggplot(ratings, aes( x = Frequency) ) +
  geom_histogram( aes(y = ..density..), color = "black", fill = "beige",
                  binwidth = NULL, bins = b) + xlab("Log Frecuencia") + 
  ylab("Densidad") + geom_density(kernel = "epanechnikov", color="red")


b = round(1 + log2(n))

ggplot(ratings, aes( x = SynsetCount) ) +
  geom_histogram( aes(y = ..density..), color = "black", fill = "beige",
                  binwidth = NULL, bins = b) + xlab("Log Synset") + 
  ylab("Densidad") + geom_density( kernel = "epanechnikov", color="red")


v = ratings$FamilySize
h = 2*IQR(v)*(n^(-1/3))

ggplot(ratings, aes( x = FamilySize ) ) +
  geom_histogram( aes(y = ..density..), color = "black", fill = "beige",
                  binwidth = h, bins = NULL) + xlab("Log FamilySize") + 
  ylab("Densidad") + geom_density( kernel = "epanechnikov", color="red")


b = round(1 + log2(n))

ggplot(ratings, aes( x = DerivEntropy ) ) +
  geom_histogram( aes(y = ..density..), color = "black", fill = "beige",
                  binwidth = NULL, bins = b) + xlab("DerivEntropy") + 
  ylab("Densidad") + geom_density( kernel = "epanechnikov", color="red")


##  Box plots

# Length
ggplot(ratings, aes(x = "" , y= Length)) + 
geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) + 
xlab("") 

# Log Frequency
ggplot(ratings, aes(x = "" , y= Frequency)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) + 
  xlab("") 

# Log SynsetCount
ggplot(ratings, aes(x = "" , y= SynsetCount)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) + 
  xlab("") 

# Log FamilySize
ggplot(ratings, aes(x = "" , y= FamilySize)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) + 
  xlab("") 

# DerivEntropy
ggplot(ratings, aes(x = "" , y= DerivEntropy)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) + 
  xlab("") 


## Dispersogramas


library("PerformanceAnalytics")

datos <- ratings[,c("Frequency", "FamilySize","SynsetCount", "Length",
                  "DerivEntropy")]

chart.Correlation(datos, histogram=TRUE, pch=19, method="spearman")


# scatter plot por sujeto

xylowess.fnc(Rating ~  Frequency | Subject, data = weightRatings,
              xlab = "log Frequency", ylab = "Weight Rating")



## mosaic plot (dative)


library(vcd)


 xtabs( ~ AnimacyOfRec + AccessOfRec + RealizationOfRecipient,
        data = dative)


mosaic(~ AnimacyOfRec + AccessOfRec + RealizationOfRecipient, 
       data = dative, main = "Animicidad vs Realización del Rec. vs Accesib. del Rec.", 
               shade = FALSE, legend = TRUE)       
        
##  Pym

#install.packages("Rling_1.0.tar.gz", repos = NULL, type = "source")

library(Rling)
library(ggplot2)
        
data(pym_high)
data(pym_low)

# 51 palabras
# syl: número de sílabas
# let: número de letras
# imag: índice subjetivo (1 a 7) de imaginabilidad (score promedio)
# conc: índice subjetivo (1 a 7) de concretud (score promedio)
# assoc: número promedio de asociaciones para cada palabra

# ¿la "asociatividad" de una palabra depende de su nivel de frecuencia?


psych::describe(pym_high$assoc)
psych::describe(pym_low$assoc)


# parece haber diferencias entre las medias / medianas de ambos grupos
# no parece haber mucha diferencia entre los desvíos
# "high" tiene una leve simetría a derecha y 
# "low", una leve simetría a izquierda

pym_assoc <- data.frame(assoc = c(pym_high$assoc, pym_low$assoc), 
                    freq = c(rep("high", 50), rep("low", 51)))


ggplot(pym_assoc, aes(x = freq, y = assoc)) + geom_boxplot() +
  xlab("Frecuencia") + ylab("Número promedio de asociaciones")

IQR(pym_high$assoc)
IQR(pym_low$assoc)

# simetría a derecha en "high" y a izquierda en "low"
# poca diferencia en rango inytercuartil (IQR)
# el número de asociaciones es más bajo el "low" respecto a "high"


# low tiene un otlier 

pym_low$assoc[pym_low$assoc <  4] # 3


# ¿la "concretud" de una palabra depende de su nivel de frecuencia?

psych::describe(pym_high$conc)
psych::describe(pym_low$conc)

# la diferencia entre medianas de ambos grupos es más alta que en las medias
# poca diferencia entre desvíos
# ambos grupos con distribuciones asimétricas a izquierda 
# ("high" más que "low")

pym_conc <- data.frame(conc = c(pym_high$conc, pym_low$conc), 
              freq = c(rep("high", 50), rep("low", 51)))

ggplot(pym_conc, aes(x = freq, y = conc)) + geom_boxplot() +
       xlab("Frecuencia") + ylab("concretud promedio")

IQR(pym_high$conc)
IQR(pym_low$conc)  
  
# no hay outliers
# poca diferencia en la dispersión
# la mediana de la concretud es más alta en "high" que en "low" 

# GJT: 
#datos <- read.spss("C:/Users/Juan/Desktop/PABLO/Notas_de_clase/DATA/SSPS/DeKeyser2000.sav", to.data.frame=TRUE)
#save(datos, file = "DeKeyser_2000.RData")

load("DeKeyser_2000.RData")

library(ggplot2)

 ggplot(datos, aes(x = Status, y = GJTScore)) + geom_boxplot() +
  xlab("Status") + ylab("GJT Score") 

 psych::describeBy(datos$GJTScore, datos$Status)
 
 # Ambos grupos son asimétricos a izquierda. El primero tiene un outliaer
 # la dispersion (IQR) parece diferente en cada grupo
 # la mediana parece difrente




### distribuciones de probabilidad

# default es N(mean= 0, sd = 1)
#1) dnorm(x,mean,sd): la funcion de densidad / masa de probabilidad
#2) pnorm(q, mean, sd, lower.tail = TRUE): 
#la función de distribución acumulada. 
#Con lower.tail = TRUE [default] calcula  P(X <= x)
#Con lower.tail = FALSE calcula P(X >= x) o sea: 1 - P(X <= x)
#3) qnorm(p, mean, sd): el percentil p-ésimo de la distribución acumulada 
#4) rnorm(n, mean, sd): da números aleatorios de una Normal (muestra aleatoria)

##Ej 1:

# distribución binomial Bi(n,p), n = 8, con p = 0.1,0.5,0.8


library(ggplot2)


df1 <- data.frame(anim = 0:8, prob = dbinom(x = 0:8, size = 8, prob = 0.1))

plot1 <- ggplot(df1, aes(x=anim, y=prob)) +
  geom_segment(aes(xend=anim, yend=0)) +
  geom_point(size=2, color="darkred") +
  theme_bw() +
  #xlab("éxitos") + ylab("proba") +
  geom_text(
    aes(label = round(prob, 2), y = prob + 0.03), # coordenada y
    position = position_dodge(0.9), # ajusta posición horizontal
    size = 3,  # tamaño del texto
    vjust = 0  # justificacion vertical = 0 = bottom
  ) +
  labs(title = "Probabilidad de X = x animados.",
       subtitle = "bi(8, 1/10)",
       x = "animados (x)",
       y = "probabilidad") 

df2 <- data.frame(anim = 0:8, prob = dbinom(x = 0:8, size = 8, prob = 0.5)) 
   
plot2 <- ggplot(df2, aes(x=anim, y=prob)) +
  geom_segment(aes(xend=anim, yend=0)) +
  geom_point(size=2, color="darkred") +
  theme_bw() +
#  xlab("éxitos") + ylab("proba") +
  geom_text(
    aes(label = round(prob,2), y = prob + 0.03),
    position = position_dodge(0.9),
    size = 3,
    vjust = 0
  ) +
  labs(title = "Probabilidad de X = x animados.",
       subtitle = "bi(8, 1/2)",
       x = "animados (x)",
       y = "probabilidad") 


 df3 <- data.frame(anim = 0:8, prob = dbinom(x = 0:8, size = 8, prob = 0.8))
   
 plot3 <- ggplot(df3, aes(x=anim, y=prob)) +
   geom_segment(aes(xend=anim, yend=0)) +
   geom_point(size=2, color="darkred") +
   theme_bw() +
  # xlab("éxitos") + ylab("proba") +
  geom_text(
    aes(label = round(prob, 2), y = prob + 0.03), # coordenada y
    position = position_dodge(0.9), # ajusta posición horizontal
    size = 3,  # tamaño del texto
    vjust = 0 # justificacion vertical = 0 = bottom
  ) +
  labs(title = "Probabilidad de X = x animados.",
       subtitle = "bi(8, 8/10)",
       x = "animados (x)",
       y = "probabilidad") 

cowplot::plot_grid(plot1, plot2, plot3, nrow = 2, ncol = 2)


# distribución Poisson(5/3)

df4 <- data.frame(conteo = 0:8, prob = dpois(x = 0:8, lambda = 5/3))

plot4 <- ggplot(df4, aes(x=conteo, y=prob)) +
  geom_segment(aes(xend=conteo, yend=0)) +
  geom_point(size=2, color="darkred") +
  theme_bw() +
  # xlab("éxitos") + ylab("proba") +
  labs(title = "Probabilidad de X = x",
       subtitle = "Poisson(5/3)",
       x = "conteo de errores (x)",
       y = "probabilidad") 



# distribución Normal: N(0,0.81); mean = 0, sd = sqrt(0.81) = 0.9

library(ggfortify)

ggdistribution(dnorm, seq(-3, 3, 0.1), mean = 0, sd = 0.9, fill = "blue")

# distribución Gamma(alfa = 3, beta = 2), alfa (shape); beta (scale)

ggdistribution(dgamma, seq(0.1, 18, 0.1), shape = 3, scale = 2, fill = "blue")

# distribución chi-cuadrado: chisq(v = 3)

ggdistribution(dchisq, seq(0.1, 18, 0.1), df = 3, fill = "blue")


##Ej 2: Z es normal estándar: N(0,1) 

ggdistribution(dnorm, seq(-4, 4, 0.1), mean = 0, sd = 1, fill = "blue")


#P(Z >= 1.23) = 0.1093486
1 - pnorm(1.23)     
pnorm(1.23, lower.tail = FALSE)
#P(Z <= 1.23) = 0.8906514
pnorm(1.23)
#P(Z >= -0.45) = 0.6736448
1- pnorm(-0.45)
pnorm(-0.45, lower.tail = FALSE)
# P(0.5 <= Z <= 1.45) = 0.2350083
pnorm(1.45) - pnorm(0.5)
# P(Z = 1.53) = 0 (¡la densidad/area bajo la curva en un punto es cero!)
# P(-1 <= Z <= 1) = 0.6826895
pnorm(1) - pnorm(-1)
# P(-2 <= Z <= 2) = 0.9544997
pnorm(2) - pnorm(-2)
# P(-3 <= Z <= 3) = 0.9973002
pnorm(3) - pnorm(-3)

##Ej. 3: X es binomial:  bi(5,0.7)

proba <- dbinom(x = 0:5, size = 5, prob = 0.7)
df.5 <- data.frame(x = 0:5, y = proba)

ggplot(df.5, aes(x=x, y=y)) +
  geom_segment(aes(xend=x, yend=0)) +
  geom_point(size=4, color="darkred") +
  theme_bw() +
  xlab("éxitos") + ylab("proba") 
   


# P(X = 3) =  0.3087
dbinom(x = 3, size = 5, prob = 0.7)
pbinom(q = 3, size = 5, prob = 0.7) - pbinom(q = 2, size = 5, prob = 0.7)
# P(X > 3) = 0.52822
1 - pbinom(3,5,0.7)
1 - sum(dbinom(0:3,5,0.7))
pbinom(3,5,0.7, lower.tail = F)
sum(dbinom(4:5,5,0.7))
# P(X < 3) = 0.16308
pbinom(2,5,0.7)
sum(dbinom(0:2,5,0.7))
# P(X <= 3) = 0.47178
pbinom(3,5,0.7)
sum(dbinom(0:3,5,0.7))
# P(X >= 3) = 0.83692
1 - pbinom(2,5,0.7) 
1 - sum(dbinom(0:2,5,0.7))
pbinom(2,5,0.7, lower.tail = F) 
# P(2 <= X <= 4) = 0.80115
sum(dbinom(2:4,5,0.7)) 
pbinom(4,5,0.7) -  pbinom(2,5,0.7) + dbinom(2,5,0.7)

## Ej 4. X es Gamma(3,2)

ggdistribution(dgamma, seq(0.1, 18, 0.1), shape = 3, scale = 2, fill = "blue")


# P(X >= 6) = 0.4231901
1 - pgamma(6, shape = 3, scale = 2)
pgamma(6, shape = 3, scale = 2, lower.tail = F)
# P(6 <= X <= 8) = 0.1850868
pgamma(8,shape = 3, scale = 2) - pgamma(6, shape = 3, scale = 2)
# P(X <= 12) = 0.9380312
pgamma(12, shape = 3, scale = 2)

## Ej 5. X es Poisson: pois(10)

proba <- dpois(x = 0:25, lambda = 10)
df.6 <- data.frame(x = 0:25, y = proba)

ggplot(df.6, aes(x=x, y=y)) +
  geom_segment(aes(xend=x, yend=0)) +
  geom_point(size=3, color="darkred") +
  theme_bw() +
  xlab("conteo") + ylab("proba") 


# P(X = 4) = 0.01891664
dpois(4, lambda = 10)
ppois(4, 10) - ppois(3,10)
# P(X < 7) = 0.1301414
ppois(6 ,10)
sum(dpois(0:6, lambda = 10))
# P(X <= 9) = 0.4579297
ppois(9, 10)
sum(dpois(0:9, lambda = 10))
# P(X >= 3) = 0.9972306
1 - ppois(2, 10)
1 - sum(dpois(0:2, lambda = 10))
ppois(2, 10, lower.tail = F)
# P(X > 6) = 0.8698586
1 - ppois(6, 10)
1 - sum(dpois(0:6, lambda = 10))
ppois(6, 10, lower.tail = F)

## Ej. 6. percentiles bajo la distribución normal estándar: N(0,1)

ggdistribution(dnorm, seq(-4, 4, 0.1), mean = 0, sd = 1, fill = "blue")

# P(Z >= z) = 0.005 --> z = 2.575829
qnorm(0.005, lower.tail = F)
qnorm(1 - 0.005)
# P(Z <= z) = 0.025 --> z = -1.959964
qnorm(0.025)
# P(Z >= z) = 0.05 --> z =  1.644854
qnorm(0.05, lower.tail = F)
qnorm(1 - 0.05)
# P(Z <= z) = 0.975 --> z =  1.959964
qnorm(0.975)
qnorm(1-0.025)
qnorm(0.025, lower.tail = F)
# P(Z >= z) = 0.5 --> z = 0 (mediana = media = 0)
qnorm(0.5)
qnorm(0.5, lower.tail = F)
qnorm(1 - 0.5)

##Ej. 7: sacar muestras de n=300 a  partir de funciones de probabilidad.

n = 300
N <- rnorm(n,0,2) # N(0,4)
E <- rexp(n, rate = 0.5) # Exp(0.5)
F <- rf(n, df1 = 5, df2 = 4) # F(5,4)
T <- rt(n, df = 4) # t(4)

b = round(2*sqrt(n))

H1 <- ggplot(data.frame(N), aes( x = N) ) +
  geom_histogram( aes(y = ..density..), color = "black", fill = "beige",
                  binwidth = NULL, bins = b) + xlab("x") + ggtitle("N(0,4)") +
  ylab("Densidad") + geom_density(kernel = "epanechnikov", color="red")

H2 <- ggplot(data.frame(E), aes( x = E) ) +
  geom_histogram( aes(y = ..density..), color = "black", fill = "beige",
                  binwidth = NULL, bins = b) + xlab("x") + ggtitle("Exp(0.5)")+
  ylab("Densidad") + geom_density(kernel = "epanechnikov", color="red")

H3 <- ggplot(data.frame(F), aes( x = F) ) +
  geom_histogram( aes(y = ..density..), color = "black", fill = "beige",
                  binwidth = NULL, bins = b) + xlab("x") + ggtitle("F(5,4)") +
  ylab("Densidad") + geom_density(kernel = "epanechnikov", color="red")

H4 <- ggplot(data.frame(T), aes( x = T) ) +
  geom_histogram( aes(y = ..density..), color = "black", fill = "beige",
                  binwidth = NULL, bins = b) + xlab("x") +  ggtitle("t(4)") +
  ylab("Densidad") + geom_density(kernel = "epanechnikov", color="red")

  gridExtra::grid.arrange(H1, H2, H3, H4, ncol = 2, nrow = 2)

### Distribuciones muestrales, intervalos de confianza, muestreo al azar

#Ej. 1. X es normal N(mu=72,sigma^2=16), sigma = 4.
# 
#(i) P(X >= 74) = P( (x-72)/4 >= (74-72)/4) = P(Z >= 1/2) = 0.3085375
1 - pnorm(0.5); pnorm(0.5, lower.tail = F)
# (ii) m.a. de n = 64; 
# La media muestral X_barra se distribuye N(72, (4^2)/64), 
# desvio = sqrt((sigma^2)/n) = sigma/sqrt(n)= 4/(sqrt(64)) = 4/8 = 1/2
#
# P(X_barra >= 74) = P( (X_barra-72)/(1/2) >= (74-72)/(1/2) )     
#                  = P(Z >= 4) = 3.17*10^(-5) --> P(Z >= 4) < 0.0001
(74-72)/(1/2)
1 - pnorm(4); pnorm(4, lower.tail = F)

#(iii) si la distribución de X no es normal:
# En (i) hay que usar la distribución que tenga X para sacar la probabilidad.
# (¡Ojo! el TCL habla sobre la media muestral de X, no de X)
#En (ii) como n >= 30 , por el TCL, se calcula igual. 

# Ej. 2
# (i)
# X tiene cualquier distribución con esperanza mu = 13.3, sigma = 1.12
# m.a. de n = 36 
# Por el TCL la media muestral se distribuye como 
# N(13.3, 1.12^2/36) ; desvio_media = 1.12/sqrt(36) = 0.186
1.12/sqrt(36)

# P(13 <= X_barra <= 13.6) 
# = P( (13 - 13.3)/0.186 <= (X_barra - 13.3)/0.186 <= (13.6 - 13.3)/0.186)
(13 - 13.3)/0.186
(13.6 - 13.3)/0.186
#P(-1.61 <= Z <= 1.61) = 0.8926
pnorm(1.61) -  pnorm(-1.61)

# (ii) m.a es de n = 36 (grande), X tiene cualquier distribución, 
# sigma conocido; aplicar el T.C.L.
# confianza 1 - alfa = 0.95 --> alfa = 0.05 --> alfa/2 = 0.025
# X_barra = 13
# ES = 1.12/sqrt(36) = 0.186 
# z(alfa/2) = qnorm(0.025, lower.tail = F) = 1.96
# IC(mu) lower = x_barra - 1.96*ES
# IC(mu) upper = x_barra + 1.96*ES

c(13 - 1.96*0.186, 13 + 1.96*0.186)  
# [12.63544, 13.36456] contiene a mu = 13.3

# (iii) el intervalo es IC 95 % = [13.1, 13.5]
# su longitud es L = 13.5 - 13.1 = 0.4
# sigma = 1.12
# para sacar el n muestral necesario aplicamos la fórmula:
# n = ((2*z(alfa/2)*sigma)/L)^2

n = ((2*1.96*1.12)/0.4)^2
# n = 120

# Ej. 3.
# X es Normal: N(16,4), sigma = 2.

# (i) P(3 <=  X  <= 8) = 3.17*10^(-5)
pnorm(8,16,2) - pnorm(3,16,2)

# (ii) Como X es Normal --> X_barra es N(16, sigma^2/n) 
#n = 5; 
#desvio_media = sigma/sqrt(n) = 2/sqrt(5) = 0.89
# P(14 <=  X_barra  <= 18) =  
# = P( (14-16)/0.89 <=  (X_barra-16)/0.89  <= (18-16)/0.89)
# = P(-2.25 <= Z <= 2.25) = 0.9755511
pnorm(2.25)-pnorm(-2.25)

# (iii) muestra es chica pero sigma es conocido (X es Normal)
# confianza 1 - alfa = 0.95 --> alfa = 0.05 --> alfa/2 = 0.025
# x_barra = 16
# ES = 2/sqrt(5) = 0.89 
# z(alfa/2) = qnorm(0.025, lower.tail = F) = 1.96
# IC(mu) lower = x_barra - 1.96*ES
# IC(mu) upper = x_barra + 1.96*ES

c(16 - 1.96*0.89, 16 + 1.96*0.89) 
# IC = [14.2556, 17.7444]

#Ej 4. 
# 5 m.a. de n = 10, X es Normal con mu = 112 , sigma desconocido
# --> IC usando distribución t(n-1)

m1 <- c(97, 117, 140, 78, 99, 148, 108, 135, 126, 121)
m2 <- c(177, 198, 107, 99, 104, 121, 148, 133, 126, 115)
m3 <- c(97, 125, 62, 120, 132, 135, 118, 137, 126, 118, 117)
m4 <- c(101, 114, 79, 120, 115, 117, 106, 86, 110, 119)
m5 <- c(137, 118, 78, 129, 87,110, 106, 116, 140, 98)

# confianza 1 - alfa = 0.95 --> alfa = 0.05 --> alfa/2 = 0.025
# distribución t(n-1) = t(10-1) = t(9)
# t(n-1,alfa/2) = qt(0.025, df = 9 ,lower.tail = F) = 2.262
# ES = sd(muestra) / sqrt(n) = sd(muestra) / sqrt(10)
# IC lower = mean(muestra) - 2.262*ES
# IC upper = mean(muestra) + 2.262*ES

# Entonces los IC de nivel 0.95 para mu = 112 son

# muestra 1
c(mean(m1) - 2.262*sd(m1)/sqrt(10), mean(m1) + 2.262*sd(m1)/sqrt(10))
# IC = [101.3799, 132.4201]

# muestra 2
c(mean(m2) - 2.262*sd(m2)/sqrt(10), mean(m2) + 2.262*sd(m2)/sqrt(10))
# IC = [109.4679, 156.1321]

# muestra 3

c(mean(m3) - 2.262*sd(m3)/sqrt(10), mean(m3) + 2.262*sd(m3)/sqrt(10))
# IC = [101.7755, 132.2245]

# muestra 4

c(mean(m4) - 2.262*sd(m4)/sqrt(10), mean(m4) + 2.262*sd(m4)/sqrt(10))
# IC = [96.59499, 116.80501]

# muestra 5

c(mean(m5) - 2.262*sd(m5)/sqrt(10), mean(m5) + 2.262*sd(m5)/sqrt(10))
# IC = [97.26599, 126.53401]

# (ii) todos los intervalos contienen a mu = 112

#Ej. 5
#(i)
#m.a. de n= 100; sigma desconocido --> aplicar T.C.L.
# confianza 1 - alfa = 0.99 --> alfa = 0.01 --> alfa/2 = 0.005
# x_barra = 36.22
# s = 1.105
# ES = s/sqrt(n) = 1.105/sqrt(100) = 0.1105
# z(alfa/2) = qnorm(0.005, lower.tail = F) = 2.576
# IC(mu) lower = x_barra - 2.576*ES
# IC(mu) upper = x_barra + 2.576*ES

c(36.22 - 2.576*0.1105, 36.22 + 2.576*0.1105)
#IC = [35.93535, 36.50465] 

#(ii) longitud es L = 0.4. Aplicar fórmula para el n muestral

# n >= ((2*z(alfa/2)*s)/L)^2

n =  ((2*2.576*1.105)/0.4)^2
n
# n = 202 

# Ej. 6

# X = cantidad de errores en 50 blancos a llenar; 
# X es binomial: bi(n= 50, p = 0.3); 
# E[X] = n*p = 50*(0.3) = 15
# V[X] = n*p*(1-p) = 50*(0.3)*(0.7) = 10.5

# siguiendo la distribución binomial:
# P(X >= 10) = 0.959
1 - pbinom(q = 9, size = 50, prob = 0.3)
pbinom(q = 9, size = 50, prob = 0.3, lower.tail = F)

# siguiendo la aproximación a la normal por n grande:
 
# Y es N(mu = np, sigma^2 = np(1-p)) --> sigma = sqrt(np(1-p))

# P(Y >= 10) =  P(Z >= (10 - 15)/sqrt(10.5) ) = P(Z >= -1.54) = 0.938
pnorm(-1.54, lower.tail = F)

# siguiendo la aproximación a la normal con corrección de continuidad:
# P(X >= x) = P(Y >= x-0.5)

# P(Y >= 9.5) =  P(Z >= (9.5 - 15)/sqrt(10.5) ) = P(Z >= -1.69) = 0.954
pnorm(-1.69, lower.tail = F)

# Ej. 7.
# n = 200
# confianza 1 - alfa = 0.95 --> alfa = 0.05 --> alfa/2 = 0.025
# x_barra = 180 / 200 = 0.9
# ES = sqrt((p*(1-p))/n) = sqrt((0.9*0.1)/200) = 0.021
# z(alfa/2) = qnorm(0.025, lower.tail = F) = 1.96
# IC(mu) lower = x_barra - 1.96*ES
# IC(mu) upper = x_barra + 1.96*ES

c(0.9 - 1.96*0.021, 0.9 + 1.96*0.021) 
# IC = [0.85884, 0.94116]; p = 96.3 no está en el intervalo

#Ej. 8:

U <- c(12,10,13,11,14,12,14,8,15,13,13,11,6,8,9,10,12,7,11,12,7,14,14,13,
  13,15,7,15,6,10,8,6,8,15,10,10,6,12,10,14,11,13,5,11,14,7,9,13,16,15)

# (i) seleccionar 20 valores por muestreo aleatorio simple

N = 50
n = 20

indices <- sample(1:20, replace = FALSE)
muestra <- U[indices]

#(ii) calcular la media y la varianza

y_barra <- mean(muestra)
f = n/N
y_var <- (var(muestra)/n)*(1-f)


#(iii) calcular el n muestral con un error del 5 %

n_o = ((1.96^2)*(y_var))/(0.05^2)

n = n_o / (1 + (n_o/N))
n
# n = 42





















