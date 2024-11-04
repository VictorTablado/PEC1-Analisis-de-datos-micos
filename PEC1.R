if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
#BiocManager::install ("SummarizedExperiment")
library(SummarizedExperiment)
library (ggplot2)
print ("ahora si?")
#install.packages("readxl")
library (readxl)
datos <- read_excel("C:/Users/Victor/Desktop/R/datos/TIO2+PTYR-human-MSS+MSIvsPD.xlsx")
head (datos)
matriz_datos <- as.matrix(datos[,c("M1_1_MSS","M1_2_MSS","M5_1_MSS","M5_2_MSS", "T49_1_MSS", "T49_2_MSS", "M42_1_PD", "M42_2_PD", "M43_1_PD", "M43_2_PD", "M64_1_PD", "M64_2_PD")])

Metadata <- data.frame (sample_id = colnames(matriz_datos), condition = c(rep("MSS", 6), rep ("PD", 6)), strigsAsFactors = FALSE)

CaracteristicasMetadata <- data.frame(
  SequenceModifications = datos$SequenceModifications, 
  stringsAsFactors = FALSE
)
CaracteristicasMetadata
se <- SummarizedExperiment(
  assays = list(counts = matriz_datos),
  rowData = CaracteristicasMetadata,
  colData = Metadata,
)

#Comprobamos el Summarized Experiment
assay(se, "counts")
colData(se)
rowData(se)

#Primero analizamos la estructura del dataset
str(se)

#Obtenemos las dimensiones dle dataset
dim(se)

#Representamos la distribución de los datos
boxplot (assay(se, "counts"), main = "Distribución de los datos", xlab = "muestras", ylab = "niveles de expresión", las = 2, col = "lightblue")


#Representamos el histograma de las 5 primeras muestras para ver su distribución más detalladamente 
par (mfrow = c(1, 5))
for (i in 1:5) {
  hist(assay(se, "counts")[, i], main = paste("Histograma de", colnames(se)[i]),
       xlab = "Nivel de expresión", col = "lightgreen", breaks = 30)
}


#Normalizamos los datos en base de logaritmo base 2, sumamos uno ya que algunos valores son 0 y de esta forma evitamos errores en la normalización
log_se <- log2(assay(se, "counts") + 1)

#Representamos los datos para ver la nueva distribución de los datos
colores <- c(rep("red", 2), rep("blue", 2), rep("green", 2), rep("orange", 2), rep("grey", 2), rep("yellow", 2))
boxplot (log_se, main = "Distribución logaritmica", xlab = "Muestras", ylab = "log2(Expresión", las = 2, col = colores)


#Calculamos la matriz de correlaciones
Cor_se <- cor (log_se)

#Y la visualizamos mediante un mapa de distancias 
heatmap (Cor_se, sym = TRUE, col = colorRampPalette(c("blue", "white", "red"))(100), main = "Matriz de correlación")
legend ("topright",
        legend = c("Correlación negativa", "Sin correlación", "Correlación positiva"),
        fill = c("blue", "white", "red"),
        title = "correlación",
        cex = 0.5)


PCA <- prcomp(t(log_se), scale = FALSE)
loads <- round (PCA$sdev^2/sum (PCA$sdev^2)*100, 1)
#Visualizamos el PCA
xlab <- c(paste("PCA1", loads[1], "%"))
ylab <- c(paste ("PCA2", loads [2], "%"))
plot(PCA$x[, 1:2], col = colores, pch = 19, xlab = xlab, ylab = ylab, main ="PCA")

legend ("topright", 
        legend = unique(Metadata$sample_id),
        col = colores,
        pch = 19,
        title = "Condiciones",
        cex = 0.6)

#T test
t_test <- apply (log_se, 1, function(x){
  grupo_MSS <- x[Metadata$condition == "MSS"]
  grupo_PD <- x[Metadata$condition == "PD"]
  t.test(grupo_MSS, grupo_PD)$p.value
})

#calculamos cuantos de los valores obtenidos son significativos
suma <- sum ((t_test) < 0.05, na.rm = TRUE)
suma

#Creamos un data frame con los resultados significativos
p_value_Significativo <- data.frame (
  Fosfopeptido = CaracteristicasMetadata$SequenceModifications,
  p_value = t_test < 0.05,
  valores = t_test
)

Significativos <- p_value_Significativo [p_value_Significativo$p_value == TRUE,]
print (Significativos)
