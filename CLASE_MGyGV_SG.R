#clase MGyGV
#GWAS

# 1. INSTALAR Y CARGAR LOS PAQUETES NECESARIOS
# 2. SETEAR UNA SEMILLA ASI TODOS TIENEN LOS MISMOS DATOS
# 3. SIMULAR MATRIZ DE DATOS GENOMICOS Y FENOTIPICOS (QTLs y EFECTOS)
# 4. PONERLE NOMBRE A LOS INDIVIDUOS Y MARCADORES
# 5. RANKEAR MARCADORES ASOCIADOS Y EFECTOS: EXPORTAR COMO TABLA
# 6. REALIZAR ANALISIS GWAS CON DISTINTAS METODOLOGIAS Y COMPARAR RESULTADOS CON NUESTROS "QTLs REALES"

# instalacion de paquetes
#PAQUETE PARA REALIZAR SELECCION GENOMICA MEDIANTE RRBLUP
install.packages("rrBLUP")


#CARGA DE PAQUETES
library(simulMGF)
library(rrBLUP)
source("https://raw.githubusercontent.com/mngar/forest/main/rrblup.R")

#setear la semilla (VAMOS A SIMULAR LOS MISMOS DATOS QUE UTILIZAMOS EN LA PRACTICA DE GWAS)
set.seed(1234)


#SIMULACION DE DATOS  (MISMOS PARAMETROS UTILIZADOS EN PRACTICA DE GWAS)
Nind <- 1000
Nmarkers <- 10000
Nqtl <- 50
Esigma <- .5
Pmean <- 25
Perror <- .25

simulN(Nind, Nmarkers, Nqtl, Esigma, Pmean, Perror)
 str(nsimout)

#los QTLs
  QTL <- cbind(nsimout$QTN, nsimout$Meffects)
  QTL <- cbind(QTL,abs(nsimout$Meffects))
  colnames(QTL) <- c("marker", "effect", "effabs")
  QTL <- as.data.frame(QTL)
  QTL <- QTL[order(-QTL$effabs),]

#matriz de datos genomicos
population <- nsimout$geno
colnames(population) <- c(paste("M", 1:Nmarkers,sep = ""))
rownames(population) <- c(paste("IND", 1:Nind,sep = ""))

#matriz de datos fenotipicos
popfeno <- nsimout$pheno
colnames(popfeno) <- "PHENO"
rownames(popfeno) <- c(paste("IND", 1:Nind,sep = ""))
feno <- data.frame ( IND = rownames(popfeno),
                     Pheno = popfeno)



#mapa
map <- data.frame (SNP = colnames(population),
                  Chromosome = c(rep(1,(Nmarkers/5)),rep(2,(Nmarkers/5)),rep(3,(Nmarkers/5)),rep(4,(Nmarkers/5)),rep(5,(Nmarkers/5))),
                  Position = c(1:(Nmarkers/5),1:(Nmarkers/5),1:(Nmarkers/5),1:(Nmarkers/5),1:(Nmarkers/5))
                  )






##########
##  SG  ##
##########
#setear el directorio de trabajo
setwd("D:/MGyGV/R/SG")


rrblup(population, feno$PHENO, 10, 90, "MGyGV")


#los archivos de salida los podemos visualizar en Excel o en R

presicion <- read.csv("MGyGV_salida_PRECISION_SG.csv", header = FALSE)
efectos <- read.csv("MGyGV_SALIDA_EFECTOS_RR.csv", header = FALSE)
obspred <- read.csv("MGyGV_observado_vs_predicho.csv", header = FALSE)
index <- read.csv("MGyGV_Index.csv", header = FALSE)

# presicion es la correlacion entre los valores predichos y observados
presicion
boxplot(presicion)

#aqui tenemos los efectos para cada uno de los marcadores
head(efectos)
#tomamos los de la primer validacion cruzada
datas = efectos[1:10000,]
colnames(datas) <- c("fold", "marker", "efecto", "efectoabs")
datas <- datas[order(-datas$efectoabs),]
#al graficar podemos comprobar que en RR-BLUP los efectos son muestreados desde una distribucion normal
plot(1:10000,datas$efecto)
#también podemos corroborar que asigne un mayor efecto absoluto a los marcadores que estan asociados a un QTL
head(datas)
head(QTL)


# valores observados vs. predichos
head(obspred)
plot(obspred$V2,obspred$V3)                                        h
plot(obspred$V2[1:100],obspred$V3[1:100])

#coloreo los individuos con valores mayores a la media poblacional para visualizar como se compone lo que estoy seleccionando     
data<- obspred[1:100,2:3]
colnames(data) = c("obs","pred")
#color por defecto
data$Colour="black"
# colores segun valor
data$Colour[data$obs>mean(feno$PHENO)]="green"
data$Colour[data$obs<mean(feno$PHENO)]="red"
# graficar
plot(data$obs,data$pred, xlim=c(15,35), col=data$Colour, ylim=c(22,28))


#en azul los que seleccionaria si hiciera un corte en el 25% superior
data<- obspred[1:100,2:3]
colnames(data) = c("obs","pred")
#color por defecto
data$Colour="black"
# colores segun valor
data$Colour[data$obs>mean(feno$PHENO)]="green"
data$Colour[data$obs<mean(feno$PHENO)]="red"
data$Colour[data$pred>quantile(data$pred, c(0.75), type = 6)]="blue"
# graficar
plot(data$obs,data$pred, xlim=c(15,35), col=data$Colour, ylim=c(22,28))



#en azul los que seleccionaria si hiciera un corte en el 25% superior, ahora en verde el 25% superior real
data<- obspred[1:100,2:3]
colnames(data) = c("obs","pred")
#color por defecto
data$Colour="black"
# colores segun valor
data$Colour[data$obs>quantile(feno$PHENO, c(0.75), type = 6)]="green"
data$Colour[data$obs<mean(feno$PHENO)]="red"
data$Colour[data$pred>quantile(data$pred, c(0.75), type = 6)]="blue"
# graficar
plot(data$obs,data$pred, xlim=c(15,35), col=data$Colour, ylim=c(22,28))


sum(data$Colour == "blue")

data[data$Colour == "blue",]

sum(data[data$Colour == "blue",]$obs>quantile(feno$PHENO, c(0.75), type = 6))
sum(data[data$Colour == "blue",]$obs<mean(feno$PHENO))

sum(data$obs>quantile(feno$PHENO, c(0.75), type = 6))

eficiencia = sum(data[data$Colour == "blue",]$obs>quantile(feno$PHENO, c(0.75), type = 6))/sum(data$obs>quantile(feno$PHENO, c(0.75), type = 6))*100


