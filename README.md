# MGyGV
## Mejoramiento Genético y Genómico Vegetal 
- Asignatura de grado de la carrera de Ingeniería en Agrobiotecnología, de la Universidad Nacional de San Martín.
- Curso de posgrado de la Facultad de Exactas y Naturales de la Universidad de Buenos Aires.

### Material adicional de las clases de Mejoramiento Genético Forestal, GWAS y Selección Genómica en el segundo cuatrimestre del ciclo lectivo 2023.
Dictadas por el Dr. Martín García y la Lic. Catalina Molina.


Clase práctica de mapeo de asociación de genoma completo (GWAS)

1.	INSTALAR Y CARGAR LOS PAQUETES NECESARIOS
•	SIMULACION DE DATOS: simulMGF
•	CHECKEO DE NORMALIDAD DEL FENOTIPO (GRAFICA Y ANALITICAMENTE): ggplot2, nortest
•	ANALISIS DE ASOCIACION: GAPIT3
2.	SETEAR UNA SEMILLA ASI TODOS OBTIENEN LOS MISMOS DATOS
3.	SIMULAR MATRIZ DE DATOS GENOMICOS Y FENOTIPICOS (QTLs y EFECTOS):
•	Población de estudio
•	Población de validación
4.	PONERLE NOMBRE A LOS INDIVIDUOS Y MARCADORES
5.	RANKEAR MARCADORES ASOCIADOS Y EFECTOS: EXPORTAR COMO TABLA
6.	REALIZAR ANALISIS GWAS CON DISTINTAS METODOLOGIAS Y COMPARAR RESULTADOS CON NUESTROS "QTLs REALES"
7.	MEDIANTE LOS MARCADORES ASOCIADOS DE UNA DE LAS METODOLOGIAS AJUSTAR UN MODELO LINEAL Y PREDECIR FENOTIPOS EN POBLACION DE VALIDACION

Clase práctica de Selección Genómica (SG)

1.	INSTALAR Y CARGAR LOS PAQUETES NECESARIOS
•	SIMULACION DE DATOS: simulMGF
•	SELECCIÓN GENOMICA: rrBLUP
2.	SETEAR UNA SEMILLA ASI TODOS OBTIENEN LOS MISMOS DATOS
3.	SIMULAR MATRIZ DE DATOS GENOMICOS Y FENOTIPICOS (QTLs y EFECTOS):
•	Población de estudio
4.	CORRER 1 VALIDACION CRUZADA 90:10 MEDIANTE LA FUNCION “AD HOC” rrblup (originalmente eran 10 validaciones cruzadas pero los tiempos varían mucho dependiendo de las capacidades de cada PC, el año próximo se puede implementar GBLUP que es computacionalmente más rápida y equivalente en precisión).
5.	ESTUDIAR LOS ARCHIVOS DE SALIDA (precisión, valores predichos y observados de la población de validación, efectos asignados a cada marcador).
6.	GRAFICAR VALORES OBSERVADOS Y PREDICHOS
7.	ESTIMAR LA EFICIENCIA DE LA SELECCIÓN BAJO UN SUPUESTO VALOR DE CORTE DEL 25% SUPERIOR.
