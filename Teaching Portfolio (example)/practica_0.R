###### PRACTICA 0: INTRODUCCIÓN A R.

# R es sensible a la diferencia entre mayúsculas y minúsculas.

#  Los directorios se expresan con "/"

#  La asignación de contenido a un objeto se realiza mediante: <-
  
# instalación de paquetes:               install.packages("paquete")      

# para usar paquetes instalados:     library("paquete") o library(paquete)

# La forma sintáctica de escribir en R es: comando()

#  para fijar directorio de trabajo: setwd("C:/Users/Pablo/APUNTES")

# para ver cuál es el directorio de trabajo: getwd()

# para conocer el contenido dentro de un directorio de trabajo: dir()

# Para pedir ayuda sobre algún comando o función de R: ?setwd

# para borrar todos los objetos de la memoria: rm(list=ls(all=T))

# Para actualizar R:
#install.packages("installr")
#library(installr)
#updateR()


##########################
###########TIPOS DE ESTRUCTURAS DE DATOS.

# (1) VECTORES: colección ordenada de valores, todos del MISMO tipo.

Freq <- c(10, 2, 2,  0, 5) # numérico
S <- c("yo" ,"tengo", "una", "vaca", "lechera") # caracteres
L <- c(TRUE, FALSE, TRUE, TRUE, FALSE)   # lógico
alfa <- 2 # constante (un vector de longitud uno)

Freq; S; L; alfa

# (2) FACTORES: para representar variables cualitativas. 
# Tiene una cantidad limitada de niveles observados.

Anim_1 <- factor(c("animado", "animado", "no animado", "animado", "no animado"),
                 levels = c("animado", "no animado"))

Anim_2 <- factor(c(1, 1, 0, 1), levels = c(1, 0))

Anim_1; Anim_2

levels(Anim_1) # para ver los niveles del factor 

# si no establecemos los niveles, R elige como referencia el
# que empieza con la letra menor del alfabeto

# para cambiar el nivel de referencia de un factor

Anim_1 <- relevel(Anim_1, ref = "no animado")
levels(Anim_1)

# (3) MATRICES: arreglo de números en filas y columnas.

matriz <- matrix(c(517,300, 33, 47),2,2, byrow=T, 
          dimnames = list(c("Animado","No animado"), 
                          c("NP","PP")))
matriz

View(matriz)


# (4) DATA.FRAME:  una base de datos plana, que permite combinar 
#              vectores con diferente tipo de datos.

df <- data.frame(Oración = S, Frecuencia = Freq, Lógico =  L)

df

View(df)

# (5) LISTAS: Es un vector genérico que contiene otros objetos creados antes.
# Elemenos pueden ser de distinto tipo.

Lista <- list(S, Freq, L, Anim_1)

# (6) SERIES TEMPORALES: es un vector de datos que contiene información 
# sobre la distribución  temporal de esos datos.

serie_temporal <-ts(1:20, start = 1980) 

## mode: para ver el tipo de los vectores.
## str ("structure"): características de los vectores. 

mode(Freq); mode(S); mode(L); mode(serie_temporal)

# info sobre cada vector: str
str(Freq); str(S); str(L); str(serie_temporal)

# PROCEDIMIENTOS PARA CREAR VECTORES O PARTES DE ELLOS

a<-rep(8,4)     # Crea un vector de números 8, con 4 elementos idénticos.
a


b<-seq(1,10,2)   # Crea una secuencia que va de 1 a 10, salteando de a dos.
b


c<-c(1:20)      # Crea un vector con números consecutivos del 1 al 20.
c

f<-c(seq(2,20,3),rep(8,3),a,b,c, Freq) # se pueden concatenar elementos
f                                      # del mismo tipo.  


########################

### OPERADORES ARITMÉTICOS

#   +   SUMA
#   -   RESTA
#   *   MULTIPLICACIÓN
#   /   DIVISIÓN
#   ^   POTENCIA 

### OPERADORES COMPARATIVOS

#   < MENOR
#   > MAYOR
#   <=  MENOR O IGUAL
#   >=  MAYOR O IGUAL
#   ==  IGUAL (OJO: no es " = ")
#   != DISTINTO          

###  OPERADORES LÓGICOS
#   !     NO
#   &     Y
#   |     O      


#################################################
# OPERACIONES CON VECTORES ##

Freq

length(Freq)  # Número de elementos
max(Freq) # valor máximo
min(Freq) # vaor mínimo 
sum(Freq) # sumatoria
diff(Freq) # diferencia entre elementos consecutivos
cumsum(Freq) # suma acumulada
sort(Freq) # ordena los valores de menor a mayor
sort(Freq, decreasing = T) # ordena los valores de mayor a menor
order(Freq) # los índices de los valores en manera ascendente; o sea,
            # cómo obtener los valores en orden ascendente 
unique(Freq) # valores no repetidos

# selección de observaciones (índices). El índice menor es 1 !!

Freq

Freq[3] # selecciona el tercer elemento
Freq[1:4] # selecciona los primeros cuatro elementos
Freq[c(1,3,5)] # selecciona el elemento primero, tercero y quinto
which.max(Freq) # seleciona el ÍNDICE del máximo valor
which.min(Freq) # seleciona el ÍNDICE del mínimo valor
Freq[Freq==2]  # selecciona el/los elementos de valor igual a 2
Freq[Freq<5]   # selecciona el/los elementos de valor menor a 5
Freq[Freq<=2]  # selecciona el/los elementos de valor menor o igual a 2
Freq[Freq>0]   # selecciona el/los elementos de valor mayor a 0
Freq[Freq>=5]  # selecciona el/los elementos de valor menor o igual a 5
Freq[Freq!=5]  # selecciona el/los elementos de valor distinto a 5

## operaciones

# calculamos la media muestral

n <- length(Freq)
media_muestral <- sum(Freq)/n

mean(Freq)

# calculamos la varianza muestral

varianza_muestral <- sum((Freq - media_muestral)^2)/(n-1)

var(Freq)

# calculamos el desvío 

desvío <- sqrt(varianza_muestral)

sd(Freq)

# calculamos el CV

CV <- desvío / media_muestral

sd(Freq) / mean(Freq)

# calculamos la mediana (para impar)

mediana_muestral <- sort(Freq)[(n+1)/2]

median(Freq)

## Reemplazar valores

Freq_2 <- Freq

Freq_2[Freq_2==2] <- 3
Freq_2[Freq_2==5] <- NA # reemplazar por un valor faltante
Freq; Freq_2


## borramos los elementos en la memoria

rm(list=ls(all = TRUE))

## CARGA DE DATOS.

V <- c("Trans","Trans","Ditr","Trans","Intrans","Intrans","Trans",
       "Trans","Ditr","Intrans", "Intrans","Trans","Intrans","Intrans",
       "Trans","Trans","Intrans","Intrans","Intrans", "Intrans")

S <- c("Hum", "Abstr","Abstr","Hum","Abstr","Hum","MatObj", "Abstr",
       "Animal", "Abstr", "Hum", "Hum", "Hum", "MatObj", "Hum",
        "Hum", "Hum", "Hum", "MatObj", "Hum")

C <- rpois(20, 4) # 20 valores muestreados de una distribución Poisson con lambda = 4

V <- factor(V)
S <- factor(S)


# concatenar columnas: cbind(V,S)
# concatenar filas: rbind(V,S)

df <- data.frame(verbo = V, sujeto = S, conteo = C)

# exportar y cargar un archivo excel.

library("readxl")
library("xlsx")


write.xlsx(df, file = "ejemplo.xlsx",
           sheetName = "ejemplo", append = FALSE)

ejemplo_df <- read_excel("ejemplo.xlsx")

ejemplo_df <- ejemplo_df[, c(2:4)]  # tomar las columnas 2 y 3

# Operaciones con data frames

dim(ejemplo_df) # cantidad de filas y columnas 
colnames(ejemplo_df) # nombres de las columnas
row.names(ejemplo_df) # nombres de las filas
head(ejemplo_df) # las primeras observaciones
tail(ejemplo_df) # las últimas observaciones
ejemplo_df$verbo; ejemplo_df[, "verbo"]; ejemplo_df[,1]  # la columna "verbo"
ejemplo_df[, c("verbo", "sujeto")] # las columnas "verbo" y "sujeto"
ejemplo_df[, -c(1)] # todas las columnas MENOS la columna "verbo"
ejemplo_df[1:8,] # las primeras 8 observaciones
ejemplo_df[-c(1:8),] # todas las observaciones MENOS las primeras 8 observaciones

# las frecuencias de las oraciones con verbo transitivo Y sujetos "humano".

ejemplo_df[ejemplo_df$verbo=="Trans" &  ejemplo_df$sujeto=="Hum","conteo"]

# las observaciones de las oraciones con verbo transitivo O distransitivo Y sujetos "humano" o "abstractos".

ejemplo_df[(ejemplo_df$verbo=="Trans"|ejemplo_df$verbo=="Ditr") &  (ejemplo_df$sujeto=="Hum"|ejemplo_df$sujeto=="Abstr"),]

# las observaciones con frecuencia menores o iguales a 5 y mayor a 2 (2<x<5).

ejemplo_df[ejemplo_df$conteo <= 5 & ejemplo_df$conteo > 2,]


# la distribución de las frecuencias

freq_dist <- table(ejemplo_df$conteo) # absolutas
prop.table(freq_dist) # relativas 


# Agregar una columna con conteos en escala logarítmica.

ejemplo_df$log_conteo <- log(ejemplo_df$conteo + 1)


# Agregar una columna con conteos estandarizados.

ejemplo_df$z_conteo <- scale(ejemplo_df$conteo)

# Agregar una columna con conteos recodificados como categorías en 
# frecuencias: baja=[0,3), media=[3,6), alta=[6,8]

#library(car)
ejemplo_df$cat_conteo <- recode(ejemplo_df$conteo, "0:2=1; 3:5=2; 6:9=3")
ejemplo_df$cat_conteo <- factor(ejemplo_df$cat_conteo,  labels = c("bajo", "medio", "alto"))

# o bien usando la función "cut":
# cut(x, breaks, labels = NULL, include.lowest = FALSE, right = TRUE)

ejemplo_df$cat_conteo <- factor(cut(ejemplo_df$conteo, breaks = c(0,3,6,9),  labels = c("bajo", "medio", "alto"), include.lowest = TRUE, right = FALSE)) 
ejemplo_df$cat_conteo <- factor(cut(ejemplo_df$conteo, breaks = c(-1,2,5,9),  labels = c("bajo", "medio", "alto")))   
  

table(ejemplo_df$cat_conteo)

# otras formas de discretizar.
#library(arules)

ejemplo_df$cat_conteo <- arules::discretize(ejemplo_df$conteo, method="fixed", breaks = c(0,3,6,9), labels = c("bajo","medio", "alto"))
# [0,3), [3,6), [6,9]

ejemplo_df$cat_conteo <- arules::discretize(ejemplo_df$conteo, method="frequency", breaks = 3, labels = c("bajo","medio", "alto"))


ejemplo_df$cat_conteo <- arules::discretize(ejemplo_df$conteo, method="cluster", breaks = 3, labels = c("bajo","medio", "alto"))

table(ejemplo_df$cat_conteo)
unique(ejemplo_df[ejemplo_df$cat_conteo=="alto", "conteo"])
unique(ejemplo_df[ejemplo_df$cat_conteo=="medio", "conteo"])
unique(ejemplo_df[ejemplo_df$cat_conteo=="bajo", "conteo"])

# Unir categorías de una variable categórica: 
# Verbo en transitivo y no transitivo

table(ejemplo_df$verbo)

ejemplo_df$verbo_binario  <- car::recode(ejemplo_df$verbo, "'Trans'='Trans'; c('Ditr', 'Intrans')='notTrans'")

table(ejemplo_df$verbo_binario)

# funciones "apply" y "tapply"

# apply: aplicar una funcion a las filas o a las columnas de un data.frame
# apply(data.frame, fila = 1 / columna = 2, función)

data <- data.frame(a=c(88, 85, 82, 97, 67, 77, 74, 86, 81, 95),
                   b=c(77, 88, 85, 76, 81, 82, 88, 91, 92, 99),
                   c=c(67, 68, 68, 74, 74, 76, 76, 77, 78, 84))

#Ej: calcular el coeficiente de variación para cada columna

apply(data,2, function(x){sd(x) / mean(x)})

#Ej: calcular la media para cada columna

apply(data,2, mean)

#Ej: calcular el máximo en cada fila

apply(data,1, max)


## tapply: calcular una función para datos agrupados según un criterio
## tapply(vector de datos, criterio de agrupación, función)

#Ej: calcular la media para cada nivel de un factor

tapply(ejemplo_df$conteo, ejemplo_df$verbo, mean)

tapply(ejemplo_df$conteo, ejemplo_df$sujeto, mean)


# libreria dplyr

# install.packages("dplyr")

library(dplyr)

nettle <- read.csv('nettle_1999_climate.csv')

head(nettle)

# seleccionar filas según algún criterio
filter(nettle, Langs > 500) # los paises con más de 500 lenguas
# nettle[nettle$Langs > 500,]
filter(nettle, Country == 'Nepal') # los datos de Nepal
# nettle[nettle$Country == 'Nepal',]

# seleccionar columnas según algún criterio

select(nettle, Langs, Country) # seleccionar las columnas 'Langs' y 'Country'
# nettle[, c("Langs", "Country")]; nettle[, c(5,1)]

 select(nettle, -Country) # selecionar todas las columnas menos 'Country'

#nettle[,-1]

select(nettle, Area:Langs) # seleccionar consecutivamente de la columna 'Area' a la columna "Langs"
# nettle[, c(3:5)]

# crear nuevas columnas

mutate(nettle, Lang100 = Langs / 100)
#nettle$Lang100 <- nettle$Langs / 100

mutate(nettle, Lang_z = scale(Langs)) # estandarizar una variable
#nettle$Lang_z <- scale(nettle$Langs)

# renombrar columnas

rename(nettle, Pop = Population)
# colnames(nettle)[2] <- "Pop"

# re-organizar filas de modo ascendente o descendente según los valores de alguna columna

arrange(nettle, Langs) # ascendente
# nettle[sort(nettle$Langs),]

arrange(nettle, desc(Langs)) # descendente
# nettle[sort(nettle$Langs, decreasing = TRUE),]

# scatter plot de Lang versus MGS

library(ggplot2)

ggplot(nettle, aes(x = MGS, y = Langs)) + geom_point()

# poner nombre a los puntos
ggplot(nettle, aes(x = MGS, y = Langs, label = Country)) +
  geom_text(check_overlap = TRUE) # check_overlap: para que no se superpongan


# juntar dos gráficos

plot1 <- ggplot(nettle, aes(x = MGS, y = Langs)) + geom_point()


plot2 <- ggplot(nettle, aes(x = MGS, y = Langs, label = Country)) +
  geom_text(check_overlap = TRUE)

plot_both <- gridExtra::grid.arrange(plot1, plot2, ncol = 2)

# grabar un gráfico

ggsave(file = 'nettle.png', plot = plot_both, width = 8, height = 6)

# otro ejemplo
library(dplyr)

icon <- read.csv('perry_winter_2017_iconicity.csv')
mod <-  read.csv('lynott_connell_2009_modality.csv')

head(icon)
head(mod)

icon <- select(icon, Word, POS, Iconicity) # seleccionar columnas

# histograma para variable iconicidad (distribución asimétrica a derecha)

ggplot(icon, aes(x = Iconicity)) +
  geom_histogram(fill = 'peachpuff3') +
  geom_vline(aes(xintercept = 0), linetype = 2) +
  theme_minimal()

head(mod)  


mod <- select(mod, Word, DominantModality:Smell) # seleccionar columnas
# mod[, c("Word","DominantModality", "Sight", "Touch", "Sound","Taste",   
#"Smell")]; mod[,c(2, 3:8)]

mod <- rename(mod, Modality = DominantModality) # renombrar columna
# colnames(mod)[2] <- "Modality"

both <- left_join(icon, mod) # unir data.frames según Word

head(both)

# seleccionar las clases Adj, Verb, Noun
both <- filter(both, POS %in% c('Adjective', 'Verb', 'Noun'))
# both[both$POS %in% c('Adjective', 'Verb', 'Noun'),]

# filtrar los NA

both <- filter(both, !is.na(Modality)) # sacar los NA de Modality
# both[!is.na(both$Modality),]

# boxplot de iconcidad según tipo de modalidad

count(both, Modality) # cuántas obs por modalidad
# table(both$Modality)

# hacer bar plot

ggplot(both, aes(x = Modality, fill = Modality)) +
  geom_bar(stat = 'count') + theme_minimal()

# barplot con porcentajes:

ggplot(data = both, aes(x = Modality, y = ..prop.., group = 1)) + 
  geom_bar(stat = "count") + 
  scale_y_continuous(labels = scales::percent_format())


