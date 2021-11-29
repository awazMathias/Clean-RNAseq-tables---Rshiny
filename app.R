## Author : Peries Mathias
## Date : Jully 2021


library(shiny)
library(shinydashboard)
library(data.table)
library(readxl)
library(tools)
library(ggplot2)
library(reshape2)
library(plotly)

ui <- dashboardPage(
  dashboardHeader(title = "Cleaning RNAseq Tables"),
  dashboardSidebar(disable = TRUE,
                   sidebarMenu(
                     menuItem("Cleaning fir RNAseq table", tabName = "Cleaning")
                   )),
  dashboardBody(
    tabItems(
      uiOutput("tab"),
      tabItem(tabName = "Cleaning",
            fluidRow(
              column(width = 12,
                     tabItem(tabName = "space1",
                             h2("                                                                                          ")
                     )
              ),
              sidebarLayout(
                sidebarPanel(
                  tabItem(tabName = "Uploader",

                          ## file input for upload only csv, xlx, xlsx or text file.
                          fileInput("datafile",label = "1 - Upload your .txt, .csv and/or .xls table of RNAseq (one or several files can be uploaded but, if multiple files, upload them one by one)",
                                    buttonLabel = "Browse...",
                                    placeholder = "No file selected",
                                    accept = c(".csv",".txt",".xls",".xlsx"),
                                    multiple = TRUE),

                          ## numeric input for define  minimal value.
                          numericInput("values", label = "2 - Give a specific value for delete rows with one values in a row less than this value (clic on the button 'Remove lines values')",
                                       value = 0),
                          ## remove all the rows with at least, one value lower than the value give before.
                          actionButton("supr0", "Remove lines < value"),

                          ## For show the number of rows deleted
                          textOutput("supr0return"),

                          ## Column for improve display
                          column(width = 12,
                                 tabItem(tabName = "space1",
                                         h2("                                                                                          ")
                                 )
                          ),
                          column(width = 12,
                                 tabItem(tabName = "space1",
                                         h2("                                                                                          ")
                                 )
                          ),

                          # Determine the specific threshold
                          numericInput("threshold", label = "3 - Give a specific threshold for delete rows with an average below this threshold (clic on the button 'Clean by specific treshold')",
                                       value = 10),

                          # Delete all rows with a mean lower than the threshold
                          actionButton("suprtrshld", "Clean by specific threshold"),

                          ## For show the number of rows deleted
                          textOutput("suprtrshldreturn")
                  ),


                  tabItem(tabName = "downloader",

                          column(width = 12,
                                 tabItem(tabName = "space1",
                                         h2("                                                                                          ")
                                 )
                          ),
                          column(width = 12,
                                 tabItem(tabName = "space1",
                                         h2("                                                                                          ")
                                 )
                          ),

                          ## For choose the ouput format of the table(s)
                          radioButtons("format_choices", label = "4 -Choose the format you want for the download table(s) (extension, separator and quoted or not",
                                       choices = list("txt" = ".txt",
                                                      "csv" = ".csv"
                                                      ),selected = ".txt"),

                          h3("Parameters"),

                          ## For choose the separator of the output file
                          radioButtons(inputId = "sep",
                                       label = "5 - Separator",
                                       choices = c(Virgule = ",",
                                                   Point_virgule = ";",
                                                   Tab = "\t",
                                                   Space = " "),
                                       selected = ",", inline=T),

                          ## For choose the quote option of the output file
                          radioButtons(inputId = "quote",
                                       label= "6 - Quote",
                                       choices = c(None = "",
                                                   "Double Quote" = '"double"',
                                                   "Single Quote" = '"escape"'),
                                       selected = "", inline=T),

                          ## Download button of an unique table
                          downloadButton("dlcomplete", label = "7 - Download the complete table"),

                          ## Download button of zip archive containing one file by column of the table create.
                          downloadButton("dlcondition", label = "7 - Download condition by condition"),

                          column(width = 12,
                                 tabItem(tabName = "space1",
                                         h2("                                                                                          ")
                                 )
                          )
                      )
                ),
                mainPanel(
                  ## Output of the table create by fusion of all files give by the user
                  DT::dataTableOutput("bigtable"),
                  ## Output of the plot that give, by condition, the repartition of RNAseq count values.
                  plotlyOutput("renderplott")
                )
              )
              )))))

server <- function(input, output, session) {
  options(shiny.maxRequestSize = 1000000*1024^2) ## Increase the possible max weight of the uploaded file (here 1T)

  ## Create reactive values for the boxplot and the table output (live updated)
  rendertable <- reactiveValues(site = NULL)
  df <- reactiveValues(site = NULL)

  observeEvent(input$datafile, {
    req(input$datafile)

    ## For determine the file extension and how the file going to be read. Save in the rendertable reactive variable.
    extension <- file_ext(input$datafile$name)
    if(is.null(rendertable$site)){ # if its the first file upload
      if(extension == 'txt'){
        rendertable$site <- read.table(input$datafile$datapath)
      }else if(extension == 'csv'){
        rendertable$site <- read.csv(input$datafile$datapath)
      }else{
        rendertable$site <- read_excel(input$datafile$datapath)
      }
      names(rendertable$site)[1] <- "genes" ## For named the first column "genes"
      names(rendertable$site)[2] <- paste("'",input$datafile$name,"'") # For give at the second column the name of the uploaded file

    }else{ ## if its not the fist file uploaded
      if(extension == 'txt'){
        saved <- read.table(input$datafile$datapath)
      }else if(extension == 'csv'){
        saved <- read.csv(input$datafile$datapath)
      }else{
        saved <- read_excel(input$datafile$datapath)
      }

      ## For show if the different file dont have the same number of genes and show how many rows are deleted.
      if(nrow(saved) > nrow(rendertable$site)){
        showNotification(paste("The number of genes between the actual and precedents files is different,", nrow(saved) - nrow(rendertable$site),"missing genes were deleted"), duration =NULL, type = "warning")
      }
      if(nrow(saved) < nrow(rendertable$site)){
        showNotification(paste("The number of genes between the actual and precedents files is different,", nrow(rendertable$site) - nrow(saved),"missing genes were deleted"), duration =NULL, type = "warning")
      }

      names(saved)[1] <- "genes"
      names(saved)[2] <- paste("'",input$datafile$name,"'")

      # Merge the file saved data table and the rendertable reactive data.
      rendertable$site <- merge(rendertable$site, saved, by = "genes")

      ## Transform the table in order that the gene column become the row names.
      saved <- data.frame(rendertable$site, row.names =  1)
      ## make a log10 transformation
      saved <- log10(saved)
      ## Transform the table by condition
      df$site = reshape2::melt(saved, variable.name = "Samples")
      ## Another transformation for create condition by samples.
      df$site <- data.frame(df$site, Condition = substr(df$site$Samples, 1, as.numeric(ncol(saved))))
    }
  })

  output$bigtable <- DT::renderDataTable({
    ## render the rendertable in the ui
    rendertable$site
  })

  ## remove all the rows with at least, one value lower than the value give by the user.
  observeEvent(input$supr0, {
      numberrows <- as.numeric(nrow(rendertable$site)) ## determine the actual number of rows of the rendertable
      rendertable$site <- rendertable$site[apply(rendertable$site[,2:ncol(rendertable$site)] > input$values, 1, all),] ## delete rows with value lower than the value give by the user.
      numberrows <- numberrows - as.numeric(nrow(rendertable$site)) ## determine the number of rows deleted
      ###
      saved <- data.frame(rendertable$site, row.names =  1)
      saved <- log10(saved)
      df$site = reshape2::melt(saved, variable.name = "Samples")
      df$site <- data.frame(df$site, Condition = substr(df$site$Samples, 1, as.numeric(ncol(saved))))
      ###

      ## Show in the ui the number of rows deleted.
      output$supr0return <- renderText({
    paste(numberrows,"genes have been deleted") ## gi
  })
    }
  )

  ## Deletion all rows with a mean lower than the threshold
  observeEvent(input$suprtrshld, {
    numberrows <- as.numeric(nrow(rendertable$site))
    rendertable$site <- rendertable$site[rowMeans(rendertable$site[2:ncol(rendertable$site)]) >= input$threshold,] ## delete rows with mean lower than the threshold.
    numberrows <- numberrows - as.numeric(nrow(rendertable$site))
    ###
    saved <- data.frame(rendertable$site, row.names =  1)
    saved <- log10(saved)
    df$site = reshape2::melt(saved, variable.name = "Samples")
    df$site <- data.frame(df$site, Condition = substr(df$site$Samples, 1, as.numeric(ncol(saved))))
    ###

    ## Show in the ui the number of rows deleted.
    output$suprtrshldreturn <- renderText({
      paste(numberrows,"genes have been deleted")
    })
    }
  )

  output$renderplott <- renderPlotly({
    ## render the boxlpot.
    if (!is.null(df$site)){
      ggplotly(ggplot(df$site, aes(x = Samples, y = value, fill = Samples)) +
                 ## Add boxplot
                 geom_boxplot() +
                 ## Title of the boxplot
                 ggtitle("Boxplot of the count repartition by genes and by conditions") +
                 ## Axes Name
                 xlab("Conditions") + ylab("log10") +
                 #For showing the legend in the bottom of the graph
                 theme(legend.position = 'bottom')
               )%>%
                  layout(legend=list(orientation = "h", y = -0.3))
    }})

  ## Save the rendertable entirely for download the complete table.
  observeEvent(input$format_choices, {
  output$dlcomplete <- downloadHandler(
    filename = function() { ## determine the filename and the extension give by the user.
      if (input$format_choices == ".txt"){
        paste("Complete_table","txt", sep=".")
      }else if (input$format_choices == ".csv"){
        paste("Complete_table","csv", sep=".")
      }
    },
    content = function(file){ ## Download the table compared to he's extension
      if (input$format_choices == ".txt"){
        if(input$quote == "" || is.null(input$quote)){
        write.table(rendertable$site, file = file, row.names = FALSE, col.names = FALSE, sep = input$sep)
        }else{
          write.table(rendertable$site, file = file, row.names = FALSE, col.names = FALSE, sep = input$sep, qmethod = input$quote)
        }
      }else if (input$format_choices == ".csv"){
        if(input$quote == "" || is.null(input$quote)){
        write.csv(rendertable$site, file = file, row.names = FALSE, col.names = FALSE, sep = input$sep)
        }else{
          write.csv(rendertable$site, file = file, row.names = FALSE, col.names = FALSE, sep = input$sep, qmethod = input$quote)
        }
      }}
    )

  ## Save the rendertable column by column and create an zip archive with one file by column (and so by condition)
  output$dlcondition <- downloadHandler(
    filename = function(){ ## Determine the name of the zip archive
      paste0("Tables_RNAseq_conditions_by_conditions",".zip")
    },
    content = function(con){
      ## Create a tempory directory
      tmpdir <- setwd(tempdir())
     on.exit(setwd(tmpdir))
      filesToSave <-NULL;

      # For all column of the rendertable create one file compared to the extension choose by the user.
      for (i in 2:ncol(rendertable$site)){
        if (input$format_choices == ".txt"){
          fileName <- paste(names(rendertable$site[i]),"txt", sep=".")
          temp_tab = data.frame("gene" = rendertable$site[,1], "Count" = rendertable$site[,i])
          if(input$quote == "" || is.null(input$quote)){
            write.table(temp_tab, file = fileName, row.names = FALSE, col.names = FALSE, sep = input$sep)
          }else{
            write.table(temp_tab, file = fileName, row.names = FALSE, col.names = FALSE, sep = input$sep, qmethod = input$quote)
          }
        }else if (input$format_choices == ".csv" || input$format_choices == ".xls" || input$format_choices == ".xlsx"){
          fileName <- paste(names(rendertable$site[i]),"csv", sep=".")
          temp_tab = data.frame("gene" = rendertable$site[,1], "Count" = rendertable$site[,i])
          if(input$quote == "" || is.null(input$quote)){
            write.csv(temp_tab, file =  fileName,  row.names = FALSE, col.names = FALSE, sep = input$sep)
          }else{
            write.csv(temp_tab, file =  fileName,  row.names = FALSE, col.names = FALSE, sep = input$sep, qmethod = input$quote)
          }
        }
        filesToSave <- c(fileName,filesToSave) ## Save the file and the filename in a string
      }
      # Create the zip archive containing all files create before
      zip(con, files = filesToSave)
    },
    contentType = "application/zip"
    )
  })
}


shinyApp(ui, server)
