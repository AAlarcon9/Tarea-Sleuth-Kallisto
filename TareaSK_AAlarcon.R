####Instalar librerias necesarias
devtools::install_github("pachterlab/sleuth")
install.packages("devtools")

library("sleuth")

BiocManager::install("biomaRt")
library(biomaRt)
library(dplyr)

#Esta funci√≥n permite mapear, a partir de la base de datos de
tx2gene <- function(){
  
#Seleccionar la base de datos que se quiere utilizar y sus atributos
  
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = "ensembl_transcript_id",
                       ens_gene = "ensembl_gene_id", ext_gene = "external_gene_name")
  return(t2g)
}

t2g <- tx2gene()

t2g

#El directorio del cual voy a recuperaar los datos
base_dir<-"~/6to semestre/Genomica/"


#Generar un objeto con los nombres de las carpetas de abundancia que yo quiera analizar

samples <- paste0("sample", c("1","2","3",
                               "10","11", "12"))
                               

#Indicar de donde y que archivos recuperar

kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))

#Establecer los grupos, es decir, las condiciones experimentales. 

s2c <- data.frame(path=kal_dirs, sample=samples, muestras = c("control","control","control","mutante",
                        
                                                                                                    "mutante", "mutante"), stringsAsFactors=FALSE)
#Crear un objeto con los resultados de kallisto

so <- sleuth_prep(s2c, ~muestras, target_mapping = t2g,extra_bootstrap_summary = TRUE)


#Ver el ajuste y estimar varianza

so <- sleuth_fit(so)


#Indicar nuestra conficiÛn de referencia contra la que vamos a comparar 

so <- sleuth_wt(so, which_beta="muestrasmutante") 


#Visualizar de manera interactiva los resultados de los an·lisis, gr·ficas y tablas

library(shiny)
sleuth_live(so)


#Te posicionas en el directorio de tu preferencia

setwd("~/6to semestre/Genomica/")


#Despues de guardar tu tabla del live, la cargas en R

resultados<-read.table("test_table.csv",sep=",",
                       header=TRUE)


#Creas un objeto que contenga los genes significativos 

significativos<-which(resultados$qval<0.1)
significativos<-resultados[significativos,]


#De los genes con valores significativos seleccionas los upregulated y los metes en un objeto

upregulated<-which(significativos$b>0)
upregulated<-significativos[upregulated,]


#Seleccionas y metes en objeto los downregulated

downregulated<-which(significativos$b<0)
downregulated<-significativos[downregulated,]


#Los gruadas en una tabla en un archivo de texto

write.table(upregulated,file="~/6to semestre/Genomica/UpregulatedCvsM.txt",sep="\t")
write.table(downregulated,file="~/6to semestre/Genomica/DownregulatedCvsM.txt",sep="\t")


#Visualizar los genes upregulated

upregulated$ext_gene


