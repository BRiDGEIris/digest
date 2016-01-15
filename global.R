SPARK_HOME<-"/Users/yalb/spark"
#SPARK_HOME<-"/home/docker/spark"

#Sys.setenv(SPARK_HOME=SPARK_HOME)
#Sys.setenv(PATH=paste0(SPARK_HOME,"/bin:",SPARK_HOME,"/sbin:",Sys.getenv("PATH")))
#.libPaths(c(file.path(Sys.getenv("SPARK_HOME"), "R", "lib"), .libPaths()))
#library(SparkR)

require('RMySQL')
require("jsonlite")
require('rpivotTable')
require('plyr')
require('shiny')
require('RJDBC')
require('RCurl')
require('queryBuildR')

#VARIANTS<-"/home/docker/variantsulb"
VARIANTS<-"/Users/yalb/Projects/Github/variant-ranking/variantsulb2"

#VARIANTS_TABLE<-"gvr4.variantsulb"
VARIANTS_TABLE<-'highlander.exomes_hc'
##VARIANTS_TABLE<-"variants"
#VARIANTS_TABLE<-"highlander"

CliniPhenomeAPI<-"http://bridgeiris.ulb.ac.be:81/bridgeirisportal/search.php"

IMPALA_CLASSPATH<-"impala-jdbc-0.5-2"
IMPALA_SERVER<-"jdbc:hive2://127.0.0.1:21050/;auth=noSasl"


if (!file.exists("/tmp/spark-events")) {
  dir.create("/tmp/spark-events")
}

drv <- JDBC(driverClass = "org.apache.hive.jdbc.HiveDriver",
                        classPath = list.files(IMPALA_CLASSPATH,pattern="jar$",full.names=T),
            identifier.quote="`")


#sparkEnvir <- list('spark.sql.parquet.binaryAsString'='true') #Needed to read strings from Parquet
#sc<-sparkR.init(master="local[12]",sparkEnvir=sparkEnvir)
#sqlContext <- sparkRHive.init(sc)#sparkRSQL.init(sc)

#df <- read.df(sqlContext, VARIANTS, "parquet")
#registerTempTable(df, VARIANTS_TABLE);

######################
#CliniPhenome
#####################

#This function retrieve the phenotype data from CliniPhenome API
updatePhenotypeTable<-function() {
  phenotypes<-fromJSON(getURL("http://bridgeiris.ulb.ac.be/cliniphenome/searchcp?ont_keyword=all"))
  phenotypes<-phenotypes[which(phenotypes[,2]>7),]
  
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

updatePhenotypeTable()

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
                "pph2_hvar_score", "pph2_hvar_pred", "sift_score", "sift_pred", "short_tandem_repeat","variant_confidence_by_depth")
fields_select<-paste(unique(c(fieldsToKeep)),collapse=",")

updateVariantFilter<-function() {
  
  condb <- dbConnect(drv,IMPALA_SERVER )
  rs<-dbSendQuery(condb,paste0("select ",fields_select," from ",VARIANTS_TABLE," limit 1"))
  column_info<-dbColumnInfo(rs)
  column_info[,2]<-as.character(column_info[,2])
  column_info[,1]<-as.character(column_info[,1])
  for (i in 1:nrow(column_info)) {
    if (column_info[i,2]=="string") {
      print(i)
      nbFields<-dbGetQuery(condb,paste0("select count(distinct ",column_info[i,1],") from ",VARIANTS_TABLE))
      if (nbFields<30) column_info[i,2]<-"factor"
    }
  }
  
  filters<-list()
  for (i in 1:nrow(column_info)) {
    filterCol<-
      switch(column_info[i,2],
             string=list(
               id= column_info[i,1],
               label= column_info[i,1],
               type= 'string',
               default_value="",
               operators=list('equal','not_equal','contains', 'in', 'not_in','begins_with', 'ends_with','is_null', 'is_not_null')),
             factor={
               values<-sort(dbGetQuery(condb,paste0("select distinct ",column_info[i,1]," from ",VARIANTS_TABLE))[,1])
               list(
                 id= column_info[i,1],
                 label= column_info[i,1],
                 type= 'string',
                 input='select',
                 values=values,
                 default_value=values[1],
                 operators=list('equal','not_equal','contains', 'in', 'not_in','is_null', 'is_not_null'))
             },
             tinyint=list(
               id= column_info[i,1],
               label= column_info[i,1],
               type= 'integer',
               default_value="",
               operators=list('equal','not_equal','less', 'less_or_equal', 'greater','greater_or_equal','between','in', 'not_in','is_null', 'is_not_null')),
             int=list(
               id= column_info[i,1],
               label= column_info[i,1],
               type= 'integer',
               default_value="",
               operators=list('equal','not_equal','less', 'less_or_equal', 'greater','greater_or_equal','between','in', 'not_in','is_null', 'is_not_null')),
             double=list(
               id= column_info[i,1],
               label= column_info[i,1],
               type= 'double',
               default_value="",
               operators=list('equal','not_equal',  'less', 'less_or_equal', 'greater','greater_or_equal','between','is_null', 'is_not_null'))
      )
    filters<-c(filters,list(filterCol))
  }
  save(file="filterVariantSpec.Rdata",filters)
}

saveGeneListHighlander<-function() {
  connectFile<-"../../connectHighlander.R"
  source(connectFile)
  
  dbGetQuery(highlanderdb,paste0("show tables"))
  data<-dbGetQuery(highlanderdb,paste0("select * from users_data where username='yleborgn' and analysis='exomes_hc'"))
  
  data<-dbGetQuery(highlanderdb,paste0("select * from users_data"))
                   
  
  userdatadb <- dbConnect(RSQLite::SQLite(), "userdata.db")
  dbWriteTable(userdatadb,"users_data",data,overwrite=T,row.names=F)
  dbDisconnect(userdatadb)
  
  userdatadb <- dbConnect(RSQLite::SQLite(), "userdata.db")
  data<-dbGetQuery(userdatadb,paste0("select * from users_data where username='yleborgn' and analysis='exomes_hc'"))
  
  dbGetQuery(userdatadb,"INSERT INTO users_data ('username', 'type','key','value') VALUES ('yleborgn', 'VALUES','SHARE|yleborgn|(digest)','CUSTOM$gene_symbol!EQUAL!A1CF?A2M-AS1?AAAS!!0!0');")
  
  
  dbGetQuery(highlanderdb,"INSERT INTO users_data (`username`,`type`,`analysis`,`key`,`value`) VALUES ('yleborgn', 'FILTER','exomes_hc','SHARE|yleborgn|(digest)','CUSTOM$gene_symbol!EQUAL!A1CF?A2M-AS1?AAAS!!0!0');")
  
  dbGetQuery(highlanderdb,"DELETE FROM users_data WHERE id=102697;")
  
  data<-dbGetQuery(userdatadb,paste0("select * from users_data where username='yleborgn' and analysis='exomes_hc' and id=200000"))

  max_id<-dbGetQuery(highlanderdb,paste0("select max(id) from users_data"))[1,1]
  dbGetQuery(highlanderdb,paste0("INSERT INTO users_data ('id', 'username', 'type','analysis','key','value') VALUES ('",max_id+1,"', 'yleborgn', 'FILTER','test3','CUSTOM$gene_symbol!EQUAL!A1CF?A2M-AS1?AAAS!!0!0');"))

}

#########
#LDAP
#########

connectLDAP<-function() {
  #curl -u cn=yann,dc=example,dc=org 'ldap://192.168.99.100/cn=yann,dc=example,dc=org'
  library("RCurl")
  getURL('ldap://192.168.99.100/cn=yann,dc=example,dc=org',.opts=list(userpwd = "cn=yann,dc=example,dc=org:test"))
  
}



getWidgetTable<-function(data,session,selection='none',targetsShort="_all") {
  action <- dataTableAjax(session, data,rownames=F)
  widget<-datatable(data, 
                    rownames=F,
                    escape=F,
                    selection = selection,
                    options = list(
                      ajax = list(url = action),
                      dom= 'lipt',
                      lengthMenu = list(c(10, 100, 1000), c('10', '100','1000')),pageLength = 10,
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
  
  condb <- dbConnect(drv,IMPALA_SERVER )
  nbrows<-dbGetQuery(condb,paste0("select count(*) from ",VARIANTS_TABLE," ",sql))
  #nbrows<-collect(sql(sqlContext,paste0("select count(*) from ",VARIANTS_TABLE," ",sql)))
  
  if (noLimit) limit<-""
  else limit<-paste0(" limit ",maxRows)
  nbRowsExceededWarningMessage<-""
  if (nbrows>maxRows) {
    nbRowsExceededWarningMessage<-paste0("Warning: Query returns ",nbrows," records. First ",maxRows," retrieved.")
  }
  
  data<-dbGetQuery(condb,paste0("select ",fields_select," from ",VARIANTS_TABLE," ",sql,limit))
  #data<-collect(sql(sqlContext,paste0("select * from ",VARIANTS_TABLE," ",sql,limit)))
    
  dbDisconnect(condb)
  
  results<-list(data=data,nbRowsExceededWarningMessage=nbRowsExceededWarningMessage)
  results
}

getNbRows<-function(sql) {
  sql<-preprocSQL(sql)
  condb <- dbConnect(drv,IMPALA_SERVER )
  nbrows<-dbGetQuery(condb,paste0("select count(*) from ",VARIANTS_TABLE," ",sql))
  dbDisconnect(condb)
  nbrows
}


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
      geneID<-res$locus
      geneID<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",geneID,"' target='_blank'>",geneID,"</a>")
      res$scoreSummary<-cbind(t(res$scores),geneID)
      colnames(res$scoreSummary)<-c("Score","Score_Case","Score_Control","Gene_Symbol")
    }
    
    if (res$scope=="digenic") {
      genes<-ldply(res$scores[1,])
      res$scores<-res$scores[2:4,]
      geneID1<-genes[,1]
      geneID1<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",geneID1,"' target='_blank'>",geneID1,"</a>")
      geneID2<-genes[,2]
      geneID2<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",geneID2,"' target='_blank'>",geneID2,"</a>")
      res$scoreSummary<-cbind(t(res$scores),geneID1,geneID2)
      colnames(res$scoreSummary)<-c("Score","Score_Case","Score_Control","Gene_Symbol1","Gene_Symbol2")
    }
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
  
  
}


