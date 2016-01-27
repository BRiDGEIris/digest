# DiGeST :Distributed Gene/Variant Scoring Tool

DiGeST is a gene/variant scoring tool, that takes as input two populations (case/control) of SNPs, and outputs a ranking of gene/variant according to their suspected pathogeneticity. It is composed of a front-end which allows a user to browse,  filter and create groups of case and control SNPs, and a back-end that perform scoring and ranking in a distributed way, using the Map/Reduce programming model. Scalable storage and computing rely on Parquet file and Spark distributed computing framework, respectively. 


## Install

DiGeST can run on a desktop, or in a distributed environment with a Hadoop cluster back-end. 

### Requirements
 
* R version > 3.2.2
* httr (1.0.0),  ggplot2 (1.0.1), shinyjs (0.2.3), shinyBS (0.61), DT (0.1.39), queryBuildR (0.1), RCurl (1.95-4.7),  bitops (1.0-6),
RJDBC (0.2-5), rJava (0.9-7), plyr (1.8.3), rpivotTable (0.1.5.2), jsonlite (0.9.17), RMySQL (0.10.3), DBI (0.3.1), shiny (0.12.2.9000)  

### Stand alone

* Clone GitHub repository
```
git clone https://github.com/Yannael/digest
```
* In RStudio, open 'digest/code/digest.Rproj'
* Open 'global.R'
* Change 'USE_CLUSTER<-TRUE' to 'USE_CLUSTER<-FALSE'
* Set environment variables to right locations:
```
SPARK_HOME<-"/path/to/spark/folder" #Where Spark is installed
VARIANTS<-"/path/to/parquet/files" #Where SNPs data are, in parquet format
```
* Run
```
library(shiny)
runApp(launch.browser = T)
```

### Hadoop cluster with R Studio

 
 



 