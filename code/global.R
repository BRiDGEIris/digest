
require('RMySQL')
require("jsonlite")
require('rpivotTable')
require('plyr')
require('shiny')
require('RJDBC')
require('RCurl')
require('queryBuildR')

CliniPhenomeAPI<-"http://bridgeiris.ulb.ac.be:81/bridgeirisportal/search.php"

pathVariants<-"/Users/yalb/Projects/Github/digest/variantsulb"

USE_CLUSTER<-TRUE

if (USE_CLUSTER) {
  SPARK_HOME<-"/home/shiny/spark"
  SPARK_HOME<-"/Users/yalb/spark"
  Sys.setenv(SPARK_HOME=SPARK_HOME)
  Sys.setenv(PATH=paste0(SPARK_HOME,"/bin:",SPARK_HOME,"/sbin:",Sys.getenv("PATH")))
  
  IMPALA_CLASSPATH<-"impala-jdbc-0.5-2"
  IMPALA_SERVER<-"jdbc:hive2://127.0.0.1:21050/;auth=noSasl"
  
  drv <- JDBC(driverClass = "org.apache.hive.jdbc.HiveDriver",
              classPath = list.files(IMPALA_CLASSPATH,pattern="jar$",full.names=T),
              identifier.quote="`")

  VARIANTS_TABLE<-'highlander.exomes_hc'
#  VARIANTS_TABLE<-'gdegols.exomes_1000g'
  
} else {
  #SPARK_HOME<-"/home/docker/spark"
  #VARIANTS<-"/home/docker/variantsulb"
  SPARK_HOME<-"/Users/yalb/spark"
  VARIANTS<-"/Users/yalb/Projects/Github/digest/variantsulb"
  
  Sys.setenv(SPARK_HOME=SPARK_HOME)
  Sys.setenv(PATH=paste0(SPARK_HOME,"/bin:",SPARK_HOME,"/sbin:",Sys.getenv("PATH")))
  .libPaths(c(file.path(Sys.getenv("SPARK_HOME"), "R", "lib"), .libPaths()))
  library(SparkR)
  require('shiny') #Necessary for column function
  
  sparkEnvir <- list('spark.sql.parquet.binaryAsString'='true') #Needed to read strings from Parquet
  sc<-sparkR.init(master="local[4]",sparkEnvir=sparkEnvir)
  sqlContext <- sparkRSQL.init(sc) #sparkRHive.init(sc)#
  
  VARIANTS_TABLE<-"variants"
  
  df <- read.df(sqlContext, VARIANTS, "parquet")
  registerTempTable(df, VARIANTS_TABLE);
  
}

if (!file.exists("/tmp/spark-events")) {
  dir.create("/tmp/spark-events")
}

######################
#CliniPhenome
#####################

#This function retrieve the phenotype data from CliniPhenome API
updatePhenotypeTable<-function() {
  phenotypes<-fromJSON(getURL("http://bridgeiris.ulb.ac.be/cliniphenome/searchcp?ont_keyword=all"))
  phenotypes<-phenotypes[which(phenotypes[,1]!="UZ Brussel"),]
  phenotypes<-phenotypes[which(phenotypes[,1]!=""),]
  
  phenotypesdb <- dbConnect(RSQLite::SQLite(), "phenotypes.db")
  colnames(phenotypes)<-c("Data_Source","Sample_ID","Pathology","Gender","Super_Population")
  dbWriteTable(phenotypesdb,"phenotypes",phenotypes,overwrite=T,row.names=F)
  dbDisconnect(phenotypesdb)
  
  data<-phenotypes
  data[,1]<-as.factor(data[,1])
  data[,2]<-as.character(data[,2])
  data[,3]<-as.factor(data[,3])
  data[,4]<-as.factor(data[,4])
  data[,5]<-as.factor(data[,5])
  filters<-getFiltersFromTable(data)
  save(file="filterPhenotypeSpec.Rdata",filters)
}

#updatePhenotypeTable()

#Retrieve phenotype table from local DB
loadPhenotypes<-function(sql) {
  
  if (sql!="") {
    sql<-paste0(" where ",sql)
    sql<-gsub(',',"','",sql)
  }
  condb<-dbConnect(RSQLite::SQLite(), paste0("phenotypes.db"))
  data<-dbGetQuery(condb,paste0("select * from phenotypes ",sql))
  dbDisconnect(condb)
  data
}

######################
#Highlander
#####################


#Only keep a subset of Highlander fields
fieldsToKeep<-c("patient","chr","pos","reference","alternative","zygosity","read_depth","genotype_quality","filters",
                "gene_symbol", "gene_ensembl", 
                "consensus_MAF","consensus_MAC", "snpeff_effect","snpeff_impact", 
                "allelic_depth_proportion_ref","allelic_depth_proportion_alt",
                "downsampled","mapping_quality_zero_reads","allele_num",
                "change_type", "transcript_ensembl","num_genes","clinvar_rs", "dbsnp_id", 
                "cadd_phred","cadd_raw","vest_score","pph2_hdiv_score", "pph2_hdiv_pred", 
                "pph2_hvar_score", "pph2_hvar_pred", "sift_score", "sift_pred", "short_tandem_repeat",
                "variant_confidence_by_depth"
                #,"change_num_blood_filt"
                )
fields_select<-paste(unique(c(fieldsToKeep)),collapse=",")

######################
#Generic functions
#####################

getWidgetTable<-function(data,session,selection='none') {
  action <- dataTableAjax(session, data,rownames=F)
  widget<-datatable(data, 
                    rownames=F,
                    escape=F,
                    selection = selection,
                    options = list(
                      ajax = list(url = action),
                      dom= 'lipt',
                      lengthMenu = list(c(20, 100, 1000), c('10', '100','1000')),pageLength = 20,
                      columnDefs = list(
                        list(
                          targets = c(1),
                          render = JS(
                            "function(data, type, row, meta) {",
                            "return type === 'display' && data.length > 15 ?",
                            "'<span title=\"' + data + '\">' + data.substr(0, 11) + '...</span>' : data;",
                            "}")
                        ),
                        list(className="dt-right",targets="_all")
                      )
                    )
  )
  widget
}

preprocSQL<-function(sql) {
  if (sql!="") {
    sql<-paste0(" where ",sql)
    sql<-gsub(' *, *',"','",sql)
  }
  sql
}

loadData<-function(sql,noLimit=F,maxRows=1000) {
  sql<-preprocSQL(sql)
  if (USE_CLUSTER) {
    condb <- dbConnect(drv,IMPALA_SERVER )
    nbrows<-dbGetQuery(condb,paste0("select count(*) from ",VARIANTS_TABLE," ",sql))
  }
  else {
    nbrows<-collect(sql(sqlContext,paste0("select count(*) from ",VARIANTS_TABLE," ",sql)))
  }
  
  if (noLimit) limit<-""
  else limit<-paste0(" limit ",maxRows)
  nbRowsExceededWarningMessage<-""
  if (nbrows>maxRows) {
    nbRowsExceededWarningMessage<-paste0("Warning: Query returns <b>",nbrows," records</b>. First ",maxRows," retrieved.")
  }
  
  if (USE_CLUSTER) {
    data<-dbGetQuery(condb,paste0("select ",fields_select," from ",VARIANTS_TABLE," ",sql,limit))
    dbDisconnect(condb)
  }
  else {
    data<-collect(sql(sqlContext,paste0("select ",tolower(fields_select)," from ",VARIANTS_TABLE," ",sql,limit)))
  }
  
  results<-list(data=data,nbRowsExceededWarningMessage=nbRowsExceededWarningMessage)
  results
}

######################
#Process results
#####################


get4<-function(l) {l[[4]]}
get1<-function(l) {l[[1]]}
get2<-function(l) {l[[2]]}

procRes<-function(results) {
  res<-list()
  res$name<-results[[1]]
  res$scale<-results[[2]]
  res$scope<-results[[3]]
  res$start_time<-as.POSIXct(results[[4]],origin = "1970-01-01",tz="Europe/Brussels")
  res$end_time<-as.POSIXct(results[[5]],origin = "1970-01-01",tz="Europe/Brussels")
  res$run_time<-round(results[[6]],digits=2)
  res$scores<-sapply(results[[7]],get2)
  res$locus<-sapply(results[[7]],get1)
  res$caseSampleID=results[[8]]
  res$controlSampleID=results[[9]]
  res$group1name=results[[10]]
  res$group2name=results[[11]]
  if (res$scale=="variant") {
    if (res$scope=="monogenic") {
      variantsID<-res$locus[1,]
      scoreSummary<-cbind(t(data.frame(sapply(variantsID,strsplit,':'))),res$locus[2,])
      colnames(scoreSummary)<-c("Chr","Position","Reference","Alternative",'Gene_Symbol')
      rownames(scoreSummary)<-NULL
      scoreSummary[,'Gene_Symbol']<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",scoreSummary[,'Gene_Symbol'],"' target='_blank'>",scoreSummary[,'Gene_Symbol'],"</a>")
      #if (results[[2]]=="pairVariantsMonogenic") {
      #  uniqueid2<-apply(res$locus[5:7,],2,paste,collapse=":")
      #  to.keep<-match(uniqueid2,dbvariants$uniqueid)
      #  res$infovariants2<-dbvariants[to.keep,]
      #  
      #  scoreSummary2<-cbind(res$locus[5,],res$infovariants2[,'reference'],res$infovariants2[,'alternative'],res$infovariants1[,'gene_symbol'])
      #  colnames(scoreSummary2)<-c("Locus2","Reference2", "Alternative2","Gene symbol")
      #}
    }
    
    if (results[[2]]=="pairVariantsDigenic") {
      uniqueid2<-apply(res$locus[6:8,],2,paste,collapse=":")
      to.keep<-match(uniqueid2,dbvariants$uniqueid)
      res$infovariants2<-dbvariants[to.keep,]
      scoreSummary2<-cbind(res$infovariants1[,'gene_symbol'],res$locus[6,],res$infovariants2[,'reference'],res$infovariants2[,'alternative'],res$infovariants2[,'gene_symbol'])
      colnames(scoreSummary2)<-c("Gene symbol1","Locus2","Reference2", "Alternative2","Gene symbol2")
      
    }
    
    res$scoreSummary<-cbind(Score=t(res$scores),scoreSummary)
    colnames(res$scoreSummary)[1:3]<-c("Score","Score_Case","Score_Control")
    
  }
  
  if (res$scale=="gene") {
    if (res$scope=="monogenic") {
      #browser()
      geneID<-res$locus
      geneID_Link<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",geneID,"' target='_blank'>",geneID,"</a>")
      res$scores[1,]<-format(as.numeric(res$scores[1,]),scientific=FALSE,digits=3)
      res$scores[2,]<-format(as.numeric(res$scores[2,]),scientific=TRUE,digits=3)
      res$scores[3,]<-format(as.numeric(res$scores[3,]),scientific=FALSE,digits=3)
      res$scores[4,]<-format(as.numeric(res$scores[4,]),scientific=FALSE,digits=3)
      res$scoreSummaryRaw<-cbind(t(res$scores),geneID)
      colnames(res$scoreSummaryRaw)<-c("Ratio_Difference","P_Value","Ratio_Case","Ratio_Control","Score_Case","Score_Control","Gene_Symbol")
      res$scoreSummary<-cbind(t(res$scores),geneID_Link)
      colnames(res$scoreSummary)<-c("Ratio_Difference","P_Value","Ratio_Case","Ratio_Control","Score_Case","Score_Control","Gene_Symbol")
    }
    
    if (res$scope=="digenic") {
      genes<-ldply(res$scores[1,])
      res$scores<-data.matrix(data.frame(res$scores[2:7,]))
      res$scores[1,]<-format(as.numeric(res$scores[1,]),scientific=FALSE,digits=3)
      res$scores[2,]<-format(as.numeric(res$scores[2,]),scientific=TRUE,digits=3)
      res$scores[3,]<-format(as.numeric(res$scores[3,]),scientific=FALSE,digits=3)
      res$scores[4,]<-format(as.numeric(res$scores[4,]),scientific=FALSE,digits=3)
      
      geneID1<-genes[,1]
      geneID1_Link<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",geneID1,"' target='_blank'>",geneID1,"</a>")
      geneID2<-genes[,2]
      geneID2_Link<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",geneID2,"' target='_blank'>",geneID2,"</a>")
      res$scoreSummaryRaw<-cbind(t(res$scores),geneID1,geneID2)
      colnames(res$scoreSummaryRaw)<-c("Ratio_Difference","P_Value","Ratio_Case","Ratio_Control","Score_Case","Score_Control","Gene_Symbol1","Gene_Symbol2")
      res$scoreSummary<-cbind(t(res$scores),geneID1_Link,geneID2_Link)
      colnames(res$scoreSummary)<-c("Ratio_Difference","P_Value","Ratio_Case","Ratio_Control","Score_Case","Score_Control","Gene_Symbol1","Gene_Symbol2")
    }
    #browser()
  }
  res
}

analysesFiles<-dir("users/analyses/","*.*")
analysesNames<-as.vector(unlist(sapply(analysesFiles,strsplit,'.txt')))

dummy<-function() {
  
  condb<-dbConnect(RSQLite::SQLite(), paste0("digest/phenotypes.db"))
  data<-dbGetQuery(condb,paste0("select * from phenotypes "))
  dbDisconnect(condb)
  
  data<-data[which(data[,1]=="1000 Genomes"),]
  dataEUR<-data[(which(data[,5]=="EUR"))[1:25],]
  dataEAS<-data[(which(data[,5]=="EAS"))[1:25],]
  
  phenotypes<-rbind(dataEUR,dataEAS)
  
  condb<-dbConnect(RSQLite::SQLite(), paste0("digest/phenotypesdemo.db"))
  dbWriteTable(condb,"phenotypes",phenotypes,overwrite=T,row.names=F)
  dbDisconnect(condb)
  
  data[which(data[,1]=="ULB"),5]<-""
  write.table(file="phenotypes.csv",data,col.names=T,row.names=F,sep=',')

  gene_list<-read.table("Guillaume_genes_new.txt",stringsAsFactor=F)[,1]
  gene_list_str<-paste0(gene_list,collapse=',')
  
  gene_list<-read.table("annick.txt",stringsAsFactor=F)[,1]
  gene_list_str<-paste0(gene_list,collapse=',')
  
  #ARHGAP11B,ASPM,ATR,ATRIP,BLM,BRAT1,C7orf27,BUB1B,CASC5,CASK,CCDC7,CDC6,CDK5RAP2,CDT1,CENPF,CENPJ,CEP135,CEP152,CEP250,CEP63,CIT,COX7B,DYRK1A,EFTUD2,EIF2AK3,ERCC3,ERCC4,ERCC5,ERCC6,ERCC8,IER3IP1,KIF11,KMT2B,MLL4,MLL2,LIG4,MCPH1,MYCN,NBN,NDE1,NIN,NIPBL,ORC1,ORC1L,ORC4,ORC4L,ORC6,ORC6L,PCNT,PHC1,PLK4,PNKP,PPM1D,RAD50,RBBP8,RNU4ATAC,SASS6,SLC25A19,SLC9A6,SMC1A,SMC3,STAMBP,STIL,TRMT10A,RG9MTD2,TUBA1A,TUBB,TUBB2B,TUBB3,TUBG1,TUBGCP4,TUBGCP6,UBE3A,WDR62,ZEB2
  
  
}



