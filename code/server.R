library(shiny)
library(DT)
library(plyr)
library(ggplot2)
library(queryBuildR)
library(shinyBS)
library(rpivotTable)
library(httr)
library(shinyjs)
require(jsonlite)
library(RCurl)

source("filterPhenotypes.R")
source("filterVariants.R")

shinyServer(function(input, output,session) {
  sessionvalues <- reactiveValues()
  data<-loadData("")
  sessionvalues$variants<-data$data
  sessionvalues$nbRowsExceededWarningMessage<-data$nbRowsExceededWarningMessage
  sessionvalues$phenotypes<-loadPhenotypes("")
  
  sessionvalues$analysesNames<-analysesNames
  sessionvalues$variantDataGene<-NULL
  
  sessionvalues$logged_user<-""
  
  #Add UI components for filtering phenotypes/variant tabs
  output<-createFilterPhenotype(input,output,session,sessionvalues)
  output<-createFilterVariant(input,output,session,sessionvalues)
  
  observe({
    query<-input$sampleid
    if (length(query)>0) {
      query<-substr(query,2,nchar(query))
      if (query!="") {
        updateTabsetPanel(session, "tabset", selected = "Gene & variant filtering manager")
      }
    }
  })
  
  
  ####################################################
  #User login and UI controls
  ####################################################
  
  disableButtons<-function() {
    shinyjs::disable("filterPhenotypeSave")
    shinyjs::disable("filterPhenotypeDelete")
    shinyjs::disable("filterVariantSave")
    shinyjs::disable("filterVariantDelete")
    shinyjs::disable("startAnalysisButton")
  }
  
  enableButtons<-function() {
    shinyjs::enable("filterPhenotypeSave")
    shinyjs::enable("filterPhenotypeDelete")
    shinyjs::enable("filterVariantSave")
    shinyjs::enable("filterVariantDelete")
    shinyjs::enable("startAnalysisButton")
  }
  
  output$loginUI<-renderUI({
    if (length(sessionvalues$logged_user)>0) {
      if (sessionvalues$logged_user=="") {
        disableButtons()
        div(div(
          strong("Not logged in. "),
          actionButton("loginButton", icon("log-in","fa-2x",lib="glyphicon"),tooltip="Log in"),
          bsAlert("alertLogin"),
          align="right"),
          bsModal("modalLogin", "Login", "loginButton", 
                  size = "small",
                  textInput("userIDLogin", "User ID:", value = ""),
                  passwordInput("passwordLogin", "Password:", value = ""),
                  actionButton("confirmLogin", label = "Log in")
          ))
      }
      else {
        enableButtons()
        div(div(
          strong(paste0("Logged in as ",sessionvalues$logged_user)),
          actionButton("logoutButton", icon("log-out","fa-2x",lib="glyphicon")),  
          align="right"))
      }
    }
  })
  
  observe({
    input$filterPhenotypeSave
    input$filterPhenotypeDelete
    input$filterVariantSave
    input$filterVariantDelete
    if (sessionvalues$logged_user=="") {
      disableButtons()
    }
    else {
      enableButtons()
    }
  })
  
  observe({
    if (length(input$confirmLogin)>0) {
      if (input$confirmLogin!=0) {
        isolate({
          users<-read.table("users.csv",header=T,stringsAsFactors=F,colClasses=c("character","character"))
          i.users<-which(users$login==input$userIDLogin)
          error_msg<-""
          if (length(i.users)>0) {
            if (input$passwordLogin==users$password[i.users]) {
              sessionvalues$logged_user<-input$userIDLogin
              runjs("$('.modal-backdrop').remove();")
              runjs('$("body").removeClass("modal-open")')
            }
            else {
              error_msg<-"Incorrect password"
            }
          }
          else {
            error_msg<-"Unknown login"
          }
          if (error_msg!="") {
            toggleModal(session, "modalLogin", toggle = "close")
            createAlert(session, "alertLogin", title = "Oops",
                        content = error_msg, append = FALSE)
          }
        })
      }
    }
  }) 
  
  observe({
    if (length(input$logoutButton)>0) {
      if (input$logoutButton) 
        sessionvalues$logged_user<-""
    }
  }) 
  
  ####################################################
  #Phenotypes group manager
  ####################################################
  
  
  #Get sample IDs button
  output$listSamplesIDs<-renderText({
    listIDs<-sessionvalues$phenotypes$Sample_ID
    result<-paste(listIDs,sep="",collapse=" , ")
    result
  })
  
  #If apply filter, update phenotype table
  observe({
    if (length(input$filterPhenotypeQueryBuilderSQL)>0)
      sessionvalues$phenotypes<-loadPhenotypes(input$filterPhenotypeQueryBuilderSQL)
  })
  
  output$showVarPhenotypesUI<-renderUI({
    isolate({
      niceNames<-as.vector(sapply(colnames(sessionvalues$phenotypes),idToName))
      selectInput('showVarPhenotypes', 'Select variables to display', niceNames, 
                  selected=niceNames,multiple=TRUE, selectize=TRUE,width='1050px')
    })
  })
  
  output$phenotypesTable<-DT::renderDataTable({
    if (length(input$showVarPhenotypes)>0) {
      data<-sessionvalues$phenotypes[,sapply(input$showVarPhenotypes,nameToId)]
      data[is.na(data)]<-''
      colnames(data)<-input$showVarPhenotypes
      getWidgetTable(data,session)
    }
  },server=T)
  
  output$downloadPhenotypesSelection <- downloadHandler(
    filename = function() {
      paste('phenotypes.zip', sep='')
    },
    content = function(con) {
      write.csv(sessionvalues$phenotypes, file="phenotypes.csv", row.names=F,quote=T)
      zip(con,c('phenotypes.csv'))
    }
  )
  
  output$pivotTablePhenotypes<-renderRpivotTable({
    if (length(input$showVarPhenotypes)>0) {
      data<-sessionvalues$phenotypes[,sapply(input$showVarPhenotypes,nameToId)]
      colnames(data)<-input$showVarPhenotypes
      rpivotTable(data,width='1050px')
    }
  })
  
  ####################################################
  #Variant group manager
  ####################################################
  
  output$nbRowsExceededWarningMessage<-renderText({
    sessionvalues$nbRowsExceededWarningMessage
  })
  
  #Apply filters
  observe({
    if (length(input$filterVariantQueryBuilderSQL)) {
      withProgress(min=1, max=3, expr={
        setProgress(message = 'Retrieving data, please wait...',
                    value=2)
        data<-loadData(input$filterVariantQueryBuilderSQL)
        sessionvalues$variants<-data$data
        sessionvalues$nbRowsExceededWarningMessage<-data$nbRowsExceededWarningMessage
      })
    }
  })
  
  output$showVarVariantsUI<-renderUI({
    isolate({
      niceNames<-as.vector(sapply(colnames(sessionvalues$variants),idToName))
      niceNames<-colnames(sessionvalues$variants)
      selectInput('showVarVariants', 'Select variables to display', niceNames, 
                  selected=niceNames[c(1,2:7)],multiple=TRUE, selectize=TRUE,width='1050px')
    })
  })
  
  output$variantsTable<-DT::renderDataTable({
    if (length(input$showVarVariants)>0) {
      data<-sessionvalues$variants[,input$showVarVariants]
      data[is.na(data)]<-''
      colnames(data)<-input$showVarVariants
      getWidgetTable(data,session)
    }
  },server=T)
  
  output$downloadVariantsSelection <- downloadHandler(
    filename = function() {
      paste('variantSelection.zip', sep='')
    },
    content = function(con) {
      write.csv(sessionvalues$variants, file="variantSelection.csv", row.names=F,quote=T)
      zip(con,c('variantSelection.csv'))
    }
  )
  output$pivotTableVariants<-renderRpivotTable({
    if (length(input$showVarVariants)>0) {
      data<-sessionvalues$variants[,sapply(input$showVarVariants,nameToId)]
      colnames(data)<-input$showVarVariants
      rpivotTable(data,width='1050px')
    }
  })
  
  ####################################################
  #Scoring tool
  ####################################################
  
  output$selectSampleGroup1UI<-renderUI({
    input$filterVariantConfirmSave
    variantsGroup<-read.table(sessionvalues$variantGroupFile,header=T,stringsAsFactors=F,colClasses=c("character","character"))
    selectInput('selectSampleGroup1', 'Control group', 
                choices = list("Groups"=variantsGroup$Name), 
                selected=variantsGroup$Name[1],
                selectize = FALSE)
  })
  
  output$selectSampleGroup2UI<-renderUI({
    input$filterVariantConfirmSave
    variantsGroup<-read.table(sessionvalues$variantGroupFile,header=T,stringsAsFactors=F,colClasses=c("character","character"))
    selectInput('selectSampleGroup2', 'Case group', 
                choices = list("Groups"=variantsGroup$Name), 
                selected=variantsGroup$Name[1],
                selectize = FALSE)
  })
  
  observe({
    if (input$startAnalysisButton>0) {
      isolate({
        analysisName<-input$analysisName
        scope<-input$rankingScope
        scale<-input$rankingScale
      })
      
      if (analysisName!='') {
        analysis<-list()
        variantsGroup<-read.table(sessionvalues$variantGroupFile,header=T,stringsAsFactors=F,colClasses=c("character","character"))
        
        isolate({
          sampleGroup1name<-input$selectSampleGroup1
          sampleGroup2name<-input$selectSampleGroup2
          
          selectSampleGroupIndex2<-which(variantsGroup$Name==sampleGroup2name)
          group2sql<-variantsGroup$SQL[selectSampleGroupIndex2]
          group2sql<-preprocSQL(group2sql)
          
          if (!is.null(sampleGroup1name)) { 
            selectSampleGroupIndex1<-which(variantsGroup$Name==sampleGroup1name)
            group1sql<-variantsGroup$SQL[selectSampleGroupIndex1]
            group1sql<-preprocSQL(group1sql)
          }
          else {
            sampleGroup1name<-"NULL"
            group1sql<-"NULL"
          }
          
          controlMAF<-input$controlGroupMAF
          caseMAF<-input$caseGroupMAF
          
          jobArguments<-rbind(analysisName,scope,scale,group1sql,group2sql,sampleGroup1name,sampleGroup2name,controlMAF,caseMAF,VARIANTS)
          setwd("spark")
          write.table(file="jobsArguments.conf",jobArguments,quote=F,col.names=F,row.names=F)
          #startCommand<-paste('spark-submit --name ",analysisName," --master local --conf spark.eventLog.enabled=true --conf spark.eventLog.dir=hdfs://node001:8020/user/yleborgn/logs ../../spark/GVR.py &')
          #startCommand<-paste('spark-submit --name ",analysisName," --master local --conf spark.eventLog.enabled=true ../../spark/GVR.py &')
          #system("../../connect.sh")
          system(paste0("./run_local.sh ",analysisName," ../users/",sessionvalues$logged_user,"/analyses"))
          setwd("..")
        })
      }  
      else {
        createAlert(session, "alertStartAnalysis", title = "Oops",
                    content = "Specify a name for your analysis", append = FALSE)
      }
    }
  })
  
  
  ####################################################
  #Results explorer
  ####################################################
  
  observe({
    input$refreshResultsButton
    if (sessionvalues$logged_user=="") folderAnalyses<-"users/analyses/"
    else folderAnalyses<-paste0("users/",sessionvalues$logged_user,"/analyses/")
    analysesFiles<-dir(folderAnalyses,"*.txt")
    sessionvalues$analysesNames<-as.vector(unlist(sapply(analysesFiles,strsplit,'.txt')))
  })
  
  output$selectAnalysisUI<-renderUI({
    selectInput('selectAnalysis', 'Select analysis', choices = sessionvalues$analysesNames
                , selected=sessionvalues$analysesNames[1],selectize = FALSE)
  })
  
  output$showVarResultsUI<-renderUI({
    nameAnalysis<-input$selectAnalysis
    
    if (length(sessionvalues$results)>0) {
      
      if (sessionvalues$results$scale=="variant") {
        if (sessionvalues$results$scope=="monogenic") {
          niceNames<-as.vector(sapply(colnames(sessionvalues$results$scoreSummary),idToName))
          initialSelect<-niceNames[c(1:5,8)]
        }
      }
      
      if (sessionvalues$results$scale=="gene") {
        if (sessionvalues$results$scope=="monogenic") {
          niceNames<-as.vector(sapply(colnames(sessionvalues$results$scoreSummary),idToName))
          initialSelect<-niceNames
        }
        if (sessionvalues$results$scope=="digenic") {
          niceNames<-as.vector(sapply(colnames(sessionvalues$results$scoreSummary),idToName))
          initialSelect<-niceNames
        }
      }
      selectInput('showVarResults', 'Select variables to display', niceNames, 
                  selected=initialSelect,multiple=TRUE, selectize=TRUE,width='1050px')
    }
  })
  
  
  output$resultsTable<-DT::renderDataTable({
    if ((length(input$showVarResults)>0) & (length(sessionvalues$results)>0)) {
      nameAnalysis<-input$selectAnalysis
      if (length(setdiff(sapply(input$showVarResults,nameToId),colnames(sessionvalues$results$scoreSummary)))==0) {
        data<-sessionvalues$results$scoreSummary[,sapply(input$showVarResults,nameToId)]
        #targetsShort<-which(colnames(data)!="Gene_Symbol")
        colnames(data)<-input$showVarResults
        getWidgetTable(data,session,selection='single')#,targetsShort=targetsShort)
      }
    }
  },server=TRUE)
  
  observe({
    nameAnalysis<-input$selectAnalysis
    if (length(nameAnalysis)>0) {
      if (sessionvalues$logged_user=="") folderAnalyses<-"users/analyses/"
      else folderAnalyses<-paste0("users/",sessionvalues$logged_user,"/analyses/")
      if (file.exists(paste0(folderAnalyses,nameAnalysis,".txt"))) {
        sessionvalues$results<-procRes(fromJSON(txt=paste0(folderAnalyses,nameAnalysis,".txt")))
        if (length(input$resultsTable_rows_selected)) {
          withProgress(min=1, max=4, expr={
            setProgress(message = 'Retrieving control data, please wait...',
                        value=2)
            
#             if (sessionvalues$results$scale=="variant") {
#               if (sessionvalues$results$scope=="monogenic") {
#                 variantsData<-sessionvalues$results$scoreSummary[input$resultsTable_rows_selected,]
#                 data<-read.table(("filterVariant.csv"),header=T,stringsAsFactors=F,colClasses=c("character","character"))
#                 sqlControl<-data$SQL[which(data$Name==sessionvalues$results$group1name)]
#                 if (length(sqlControl)>0) {
#                   sqlControl<-paste0(sqlControl," and chr='",variantsData['Chr'],"' and position=",variantsData['Position']," and reference='",variantsData['Reference'],"' and alternative='",variantsData['Alternative'],"'")
#                   variantsControl<-loadData(sqlControl)$data
#                   nControl<-nrow(variantsControl)
#                 }
#                 else {
#                   variantsControl<-NULL
#                   nControl<-0
#                 }
#                 setProgress(message = 'Retrieving case data, please wait...',
#                             value=3)
#                 sqlCase<-data$SQL[which(data$Name==sessionvalues$results$group1name)]
#                 sqlCase<-paste0(sqlCase," and chr='",variantsData['Chr'],"' and position=",variantsData['Position']," and reference='",variantsData['Reference'],"' and alternative='",variantsData['Alternative'],"'")
#                 variantsCase<-loadData(sqlCase)$data
#                 variants<-rbind(variantsControl,variantsCase)
#                 variants<-cbind("Group"=c(rep(sessionvalues$results$group1name,nControl),rep(sessionvalues$results$group2name,nrow(variantsCase))),variants)
#                 sessionvalues$variantDataGene<-variants
#               }
#             }
#             
            if (sessionvalues$results$scale=="gene") {
              if (sessionvalues$results$scope=="monogenic") {
                #browser()
                geneID<-sessionvalues$results$scoreSummaryRaw[input$resultsTable_rows_selected,'Gene_Symbol']
                data<-read.table(paste0("users/",sessionvalues$logged_user,"/filterVariant.csv"),header=T,stringsAsFactors=F,colClasses=c("character","character"))
                sqlControl<-data$SQL[which(data$Name==sessionvalues$results$group2name)]
                if (length(sqlControl)>0) {
                  sqlControl<-paste0(sqlControl," and gene_symbol='",geneID,"'")
                  variantsControl<-loadData(sqlControl,noLimit=T)$data
                  nControl<-nrow(variantsControl)
                }
                else {
                  variantsControl<-NULL
                  nControl<-0
                }
                setProgress(message = 'Retrieving case data, please wait...',
                            value=3)
                data<-read.table(paste0("users/",sessionvalues$logged_user,"/filterVariant.csv"),header=T,stringsAsFactors=F,colClasses=c("character","character"))
                sqlCase<-data$SQL[which(data$Name==sessionvalues$results$group1name)]
                sqlCase<-paste0(sqlCase," and gene_symbol='",geneID,"'")
                variantsCase<-loadData(sqlCase,noLimit=T)$data
                variants<-rbind(variantsControl,variantsCase)
                variants<-cbind("Group"=c(rep(sessionvalues$results$group2name,nControl),rep(sessionvalues$results$group1name,nrow(variantsCase))),variants)
                sessionvalues$variantDataGene<-variants
              }
              if (sessionvalues$results$scope=="digenic") {
                geneID1<-sessionvalues$results$scoreSummaryRaw[input$resultsTable_rows_selected,'Gene_Symbol1']
                geneID2<-sessionvalues$results$scoreSummaryRaw[input$resultsTable_rows_selected,'Gene_Symbol2']
                data<-read.table(paste0("users/",sessionvalues$logged_user,"/filterVariant.csv"),header=T,stringsAsFactors=F,colClasses=c("character","character"))
                sqlControl<-data$SQL[which(data$Name==sessionvalues$results$group2name)]
                if (length(sqlControl)>0) {
                  sqlControl<-paste0(sqlControl," and (gene_symbol='",geneID1,"' or gene_symbol='",geneID2,"')")
                  variantsControl<-loadData(sqlControl,noLimit=T)$data
                  nControl<-nrow(variantsControl)
                }
                else {
                  variantsControl<-NULL
                  nControl<-0
                }
                setProgress(message = 'Retrieving case data, please wait...',
                            value=3)
                data<-read.table(sessionvalues$variantGroupFile,header=T,stringsAsFactors=F,colClasses=c("character","character"))
                sqlCase<-data$SQL[which(data$Name==sessionvalues$results$group1name)]
                sqlCase<-paste0(sqlCase," and (gene_symbol='",geneID1,"' or gene_symbol='",geneID2,"')")
                variantsCase<-loadData(sqlCase,noLimit=T)$data
                variants<-rbind(variantsControl,variantsCase)
                variants<-cbind("Group"=c(rep(sessionvalues$results$group2name,nControl),rep(sessionvalues$results$group1name,nrow(variantsCase))),variants)
                sessionvalues$variantDataGene<-variants
              }
            }
          }
          )
        }
      }
    }
  })
  
  output$showVarMetadataUI<-renderUI({
    nameAnalysis<-input$selectAnalysis
    if (length(sessionvalues$variantDataGene)>0) {
      niceNames<-as.vector(sapply(colnames(sessionvalues$variantDataGene),idToName))
      selectInput('showVarMetadata', 'Select variables to display', niceNames, 
                  selected=niceNames[c(1:8,24:26,17)],multiple=TRUE, selectize=TRUE,width='1050px')
    }
  })
  
  output$variantsMetadataTable<-DT::renderDataTable({
    nameAnalysis<-input$selectAnalysis
    if (length(sessionvalues$variantDataGene)>0 & length(input$showVarMetadata)>0) {
      data<-sessionvalues$variantDataGene[,sapply(input$showVarMetadata,nameToId)]
      getWidgetTable(data,session)
    }
  },server=TRUE)
  
  output$variantsMetadataPivotTable<-renderRpivotTable({
    #input$checkboxAddMetadata
    nameAnalysis<-input$selectAnalysis
    if (length(sessionvalues$variantDataGene)>0 & length(input$showVarMetadata)>0) {
      data<-sessionvalues$variantDataGene[,sapply(input$showVarMetadata,nameToId)]
      rpivotTable(data)
    }
  })
  
  
  output$resultsMetadata<-renderUI({
    fluidRow(
      column(3,
             strong("Control group: "),br(),
             strong("Pathological group: "),br(),
             strong("Start time: "),br(),
             strong("End time: "),br(),
             strong("Total run time: ")
      ),
      column(4,
             paste0(sessionvalues$results$group1name," (n=",length(sessionvalues$results$controlSampleID),")"),br(),
             paste0(sessionvalues$results$group2name," (n=",length(sessionvalues$results$caseSampleID),")"),br(),
             sessionvalues$results$start_time,br(),
             sessionvalues$results$end_time,br(),
             paste(sessionvalues$results$run_time, "seconds")
      )
      
    )
  })
  
  #Download score list as CSV
  output$downloadScoreTable <- downloadHandler(
    filename = function() {
      paste0(input$selectAnalysis,'.zip')
    },
    content = function(con) {
      filename=paste0(input$selectAnalysis,'.csv')
      write.csv(sessionvalues$results$scoreSummaryRaw, file=filename, row.names=F,quote=F)
      zip(con,c(filename))
    }
  )
  
  observe({
    input$shareHighlanderConfirm
    isolate({
      if ((length(input$shareHighlanderPatientSelect)>0) & length(input$shareHighlanderGeneSelect)>0) {
        patients<-input$shareHighlanderPatientSelect
        genes<-input$shareHighlanderGeneSelect
        username<-input$shareHighlanderUserSelect
        
        result<-paste0("CUSTOM$&[patient!EQUAL!",paste0(patients,collapse="?"),"!!0!0#gene_symbol!EQUAL!",
                        paste0(genes,collapse="?"),"!!0!0#]")
        #browser()
        connectFile<-"../connectHighlander.R"
        source(connectFile)
        dbGetQuery(highlanderdb,paste0("INSERT INTO users_data (`username`,`type`,`analysis`,`key`,`value`) VALUES ('yleborgn', 'FILTER','exomes_hc','SHARE|",username,"|DiGeST','",result,"');"))
        dbDisconnect(highlanderdb)
        
        toggleModal(session, "shareHighlanderModal", toggle = "close")
      }
    })
  })
  
  
  output$shareHighlanderUI<-renderUI({
    connectFile<-"../connectHighlander.R"
    browser()
    source(connectFile)
    users<-sort(dbGetQuery(highlanderdb,paste0("select * from users"))$username[-1])
    dbDisconnect(highlanderdb)
    genes<-sessionvalues$results$scoreSummaryRaw[,'Gene_Symbol']
    patients<-c(sessionvalues$results$caseSampleID)
    fluidRow(
      column(12,
      selectInput('shareHighlanderUserSelect', 'Share with user', 
                  choices = users, 
                  selected=users[1],
                  selectize = FALSE),
      selectInput('shareHighlanderPatientSelect', 'Sample list', 
                  choices = patients, 
                  selected=patients,
                  selectize = TRUE,multiple=TRUE,width='900px'),
      selectInput('shareHighlanderGeneSelect', 'Gene list', 
                  choices = genes, 
                  selected=genes[1],
                  selectize = TRUE,multiple=TRUE,width='900px')
      )
    )
  })
  
  output$resultsPanel<-renderUI({
    fluidRow(
      column(12,
             uiOutput("resultsMetadata"),
             hr(),
             h3("Scoring results"),
             div(
               downloadButton('downloadScoreTable', label = "Download score table (CSV)",class = NULL),
               actionButton("shareHighlanderButton", label = "Share in Highlander"),
               align="right"),
             bsModal("shareHighlanderModal", "Share in Highlander", "shareHighlanderButton", 
                     size = "large",
                     uiOutput("shareHighlanderUI"),
                     actionButton("shareHighlanderConfirm", label = "Share")
             ),
             uiOutput("showVarResultsUI"),
             DT::dataTableOutput('resultsTable'),
             hr(),
             fluidRow(
               column(12,
                      fluidRow(
                        uiOutput("showVarMetadataUI"),
                        dataTableOutput('variantsMetadataTable')
                      ),
                      fluidRow(
                        h4("Pivot Table"),
                        rpivotTableOutput('variantsMetadataPivotTable')
                      )
               )
             )
      )
    )
    
  })
  
})

