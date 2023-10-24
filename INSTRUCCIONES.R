
# 1. En caso de que no los tengan instalados previamente: Instalar R y Rtools desde https://www.r-project.org/
# 2. INSTALAR Y CARGAR LOS SIGUIENTES PAQUETES copiando y pegando el siguiente código en la consola de R:                                                          

#Para la clase de mapeo por asociación:
#PAQUETE DE SIMULACION DE DATOS
install.packages("simulMGF")

#PAQUETES PARA CHECKEAR NORMALIDAD DEL FENOTIPO (GRAFICA Y ANALITICAMENTE)
install.packages("ggplot2")
install.packages("nortest")

#PAQUETE PARA MAPEO POR ASOCIACION
install.packages("devtools")
devtools::install_github("jiabowang/GAPIT", force=TRUE)

#SOLO SI LO ANTERIOR NO FUNCIONA:
# copiar el codigo presente en la siguiente página al block de notas y guardar como .txt
# http://zzlab.net/GAPIT/gapit_functions.txt

#comprobar que los paquetes se hayan cargado correctamente:
library(simulMGF)
library(ggplot2)
library(nortest)
library(GAPIT)


#Para la clase de SG
install.packages("rrBLUP")
install.packages("BLR")


# 3. DESCARGAR LOS ARCHIVOS CON EXTENSION ".R" PRESENTES EN  https://github.com/mngar/CURSO
# AQUI SE ENCUENTRAN LOS CODIGOS A UTILIZAR EN LA PARTE PRACTICA
