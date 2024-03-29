##### CyTOF App
### PACKAGE INSTALLATION AND MANAGEMENT 
#list of packages required
list.of.packages <- c('shiny', 'gateR', 'dplyr', 'shinyFiles')

#checking missing packages from list
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

#install missing ones
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)


## Bioconductor packages
if(!('flowCore' %in% installed.packages())){
    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("flowCore")
}



### ATTACHING NECESSARY PACKAGES

library(flowCore)
library(tools)
library(shiny)
library(gateR)
library(ggplot2)
library(dplyr)
library(shinyFiles)


### TIMEOUT

 options(shiny.maxRequestSize=500*1024^2)
# timeoutSeconds <- 60*15 # 15'
# inactivity <- sprintf("function idleTimer() {
#   var t = setTimeout(logout, %s);
#   window.onmousemove = resetTimer; // catches mouse movements
#   window.onmousedown = resetTimer; // catches mouse movements
#   window.onclick = resetTimer;     // catches mouse clicks
#   window.onscroll = resetTimer;    // catches scrolling
#   window.onkeypress = resetTimer;  //catches keyboard actions
#   function logout() {
#   Shiny.setInputValue('timeOut', '%ss')
#   }
#   function resetTimer() {
#   clearTimeout(t);
#   t = setTimeout(logout, %s);  // time is in milliseconds (1000 is 1 second)
#   }
#   }
#   idleTimer();", timeoutSeconds*1000, timeoutSeconds, timeoutSeconds*1000)



### GLOBAL VARS

#Establishing a blacklist for flowcell data we will not need
blacklist <- c("Time", "Cell_length", "beadDist", 'Event_length', 'Center', 'Residual', 'Width', 'Offset')
manual_colnames <- c("ID", "C1", "C2","Dataset")

# A queue of notification IDs
ids <- character(0)
# A counter
n <- 0

### SERVER VS LOCAL

if(strsplit(getwd(), split = '/')[[1]][2] == 'srv'){
  # This is the path if being run on BioServer
  maindir <- '/var/www/download/cytof/' # Needs a path if on the BioServer
}else{
  # This includes being run locally (current default only works for Andrew’s path)
  maindir <- '' # Empty for AG or a path to local version recording
}



### DATAFRAME RESTRUCTURING FUNCTION

move_to_start <- function(x, to_move) {
  x[, c(to_move, setdiff(colnames(x), to_move))]
} 


### UI

ui <- fluidPage(
  #tags$script(inactivity),  
  tabsetPanel(  
    tabPanel("Processor", fluid = TRUE,
             titlePanel("CyTOF Processor"),
             # Copy the line below to make a text input box
             sidebarPanel(
               textInput("name", label = h4("Name of Analysis"), placeholder = "Please label your analysis"),
               br(),
               shinyDirButton('dir_btn', 'Select directory', 'Please select a directory to save the analysis files to', multiple = FALSE),
               textOutput('sel_dir'),
               hr(),
               
               fileInput("fcs", "Choose .fcs files",
                         multiple = TRUE,
                         accept = c('.fcs')),
               
               fileInput("csv", "Choose metadata .csv file",
                         multiple = FALSE,
                         accept=c('text/csv', 'text/comma-separated- values,text/plain', '.csv')),
               
               selectInput("corr_val", label = "Correction:",
                            choices = c("None" = "none",
                              "False Discovery Rate" = "FDR",
                              "Correlated Sidak" = "correlated Sidak",
                              "Uncorrelated Sidak" = "uncorrelated Sidak",
                              "Correlated Bonferroni" = "correlated Bonferroni",
                              "Uncorrelated Bonferroni" = "uncorrelated Bonferroni",
                              "Adler and Hasofer" = "Adler and Hasofer")),
               
               radioButtons("alpha", "Alpha:",
                            c("0.001" = 0.001,
                              "0.01" = 0.01,
                              "0.05" = 0.05)),
               checkboxInput("numerator", "Extract cells from case clusters", TRUE),
               checkboxInput("arcsinh", "Arcsinh Transform", FALSE),
               selectizeInput('plotted_markers',
                              choices = NULL,
                              multiple = T,
                              label = HTML('<p style=“color:black;“>Markers for Gating
                                           </p>')),
               
               # textInput('bl',
               #                placeholder = "Enter non-analysis columns"
               #                ),
               hr(),
               actionButton("script", "Send")
             ),
             
             mainPanel(
               h3("App Overview"),
               h5("This application is designed to accept a user's .fcs files and transform the data into a format that the gateR package will use to perform a gating analysis. The user will need to provide a name for the analysis and a directory location for the final results to be saved to. The application will create a directory at this location and will have the name specified earlier. The final directory will contain the original metadata file, the gated .RDS object, and any generated RRS plots. Use the next tab to view density plots of the markers from a completed object. Continue reading for information on how to use this application."),
               tags$a(href="https://cran.r-project.org/web/packages/gateR/gateR.pdf", target="_blank",
                      "Documentation for the gateR package is available here."),
               
               h3("Metadata"),
               p("The .csv metadata file must contain the filenames, conditions, and indicate which of those conditons should be treated as the control:"),
               tableOutput('mtable'),
               p("The first column contains the filename without the extension. The C1 and C2 columns contain the conditions in the experiment.The C1-Control and C2-Control columns contain a boolean value specifying which of the contion values was the control in the experiment. Only a single TRUE and FALSE value is needed for each of these columns. In this example, \'untreated\' is the control, so there is a TRUE in the cell adjacent to it, and a FALSE adjacent to the \'treated\' value.", span(strong("The application currently only supports 2 conditions max.")), "If there is only a single condition, the C2 and C2. Control columns may be ommitted from the metadata file. A metadata sample is available for download below:"),
               downloadButton("downloadData", "Download Metadata Example"),
               
               h3("Selecting Markers for Gating"),
               p("An even number of markers must be selected for the analysis.", span(strong(" The order of the marker selection matters. From the gateR documentation: the odd-numbered elements are the markers used on the x-axis of a gate and the even-numbered elements are the markers used on the y-axis of a gate."))," Markers can be repeated in successive gates.")
               
             ) # End tab one main
             
    ), # end tab one
    
    tabPanel("Plotting", fluid = TRUE,
             titlePanel("CyTOF Plotter"),
             sidebarPanel(
            fileInput("rds", "Choose .RDS object",
                         multiple = FALSE,
                         accept = c('.RDS')),
               
               
             ), #end of sidebar
             
             mainPanel(
               
               # conditionalPanel(condition = "length(conditions()) == 2",
               #                   radioButtons("cview", "Condition:",
               #                                c("0.001" = 0.001,
               #                                  "0.01" = 0.01,
               #                                  "0.05" = 0.05))),
               
               #uiOutput("markers_to_plot_ui"),
               selectizeInput('markers_to_plot',
                              choices = NULL,
                              multiple = F,
                              label = HTML('<p style=“color:black;“>Marker for Plotting:</p>')),
               plotOutput("plot1"),
               downloadButton("downloadPlot", "Download Plot")
               
             )      
             
    ) ### End tab two
    
    # tabPanel("FAQ", fluid=TRUE,
    #          #https://github.com/Waller-SUSAN/gateR
    #          mainPanel(
    #            h3("FAQ")
    #           
    #          )
      
    #)
    
    
  ) # End Tabset Panel
)  ## End of UI



### PATH DETECTION 

if(strsplit(getwd(), split = "/")[[1]][2] == "srv"){
  # This is the path if being run on BioServer
  #sessionDF <- versionDF(version_df_dir = "/home/goodspeed/App_package_versions") # Needs a path if on the BioServer
  download_path <- "/var/www/download/cytof/"
}else{
  # This includes being run locally (current default only works for Andrew's path)
  #sessionDF <- versionDF() # Empty for AG or a path to local version recording
  download_path <- "./"
}


####### SERVER

server <- function(input, output, session) {
  
  # observeEvent(input$timeOut, {
  #   # Modified from: https://stackoverflow.com/questions/33839543/shiny-server-session-time-out-doesnt-work
  #   showModal(modalDialog(
  #     title = "Inactivity Timeout",
  #     paste("Session timeout due to", input$timeOut, "inactivity -", Sys.time()),
  #     footer = NULL
  #   ))
  #   stopApp() # I made this change so that the app closes instead of the window only
  # })
  
  
  ## This is all for the metadata example and download
  data_example <- read.csv("metadata-template.csv", stringsAsFactors = FALSE, check.names=FALSE)
  
  output$mtable <- renderTable(data_example)
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("metadata-template-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(data_example, file)
    }
  )
  
  
  ### First tab analysis 
  shinyDirChoose(input, 'dir_btn', roots = c(home = '~'), session=session)
  dir <- reactive({parseDirPath(c(home = '~'), input$dir_btn)})
  
  observe({
    if(!is.null(dir())){
      output$sel_dir <- renderText(dir())
    }
  })
  
  
  #Get the name of the run and create a directory with that name
  fcsdir <- reactive({
    
    #Create the run name here
    dirname <- paste(dir(), '/', paste(input$name, format(Sys.time(), "%H%M%S_%m-%d-%y"), sep = "_"), sep='')
    print(dirname)
    
    # if statement, dirname == strsplit(dirname, /), then lapply

    lapply(dirname, function(x) if(!dir.exists(x)) dir.create(x))

    
    #create a new directory
    
    
    return(dirname)
    
  })

  
  set <- reactive({
    fs <- NULL
    req(input$fcs$datapath)
    #print(input$fcs$datapath)
    if(length(grep("\\.fcs$", input$fcs$datapath)) == 0){
        id <- showNotification("Incorrect file type",type = "error", duration = 10)
        ids <<- c(ids, id)
        n <<- n + 1
      }
    
   else { 
      fs <- read.flowSet(files = input$fcs$datapath, pattern = ".fcs")
      
      #Insert catch error for mismatched files here
      sampleNames(fs) <- input$fcs$name
      return(fs)
      
    }
  })
  
  
  channel_obj <- reactive({
    channels <- colnames(set())
    channels <- channels[!(channels %in% blacklist)]
  })
  
  
  #Build out a dataframe for the markers -- same as the one in the script
  marker_object <- reactive({
    req(set())
    m_dim <- length(channel_obj())
    
    markers <- as.data.frame(matrix(nrow = 2, ncol = m_dim))
    for (ch in 1:m_dim) {
      
      if(!is.na(match(getChannelMarker(set()[[1]],channel_obj()[ch])$desc, markers[1,]))) {
        
        firstindex <- match(getChannelMarker(set()[[1]],channel_obj()[ch])$desc, markers[1,])
        suffix <- gsub("[()]", "", channel_obj()[ch])
        orig_suffix <- gsub("[()]", "", channel_obj()[firstindex])
        
        markers[1, ch] <- paste(getChannelMarker(set()[[1]],channel_obj()[ch])$desc, suffix, sep = "-")
        markers[1, firstindex] <- paste(getChannelMarker(set()[[1]],channel_obj()[ch])$desc, orig_suffix, sep = "-")
        markers[2, ch] <- channel_obj()[ch]
        
      }
      else {
        markers[1, ch] <- getChannelMarker(set()[[1]],channel_obj()[ch])$desc
        markers[2, ch] <- channel_obj()[ch]
      }
      #gsub("[()]", "", x)
    }
    
    return(markers)
    
  })
  
  
  #Parse the markers inputted by the user and make sure that there are an even number selected
  markerParse <- reactive({
    collapsedmarks <- paste(input$plotted_markers, collapse="xyz")
    return(collapsedmarks)
  })
  
  #output$marktest <- renderText(markerParse())
  
  #This populates the marker dropdown
  observeEvent(set(), {
    
    updateSelectizeInput(session, 'plotted_markers',
                         choices = as.character(marker_object()[1,]),
                         server = TRUE,
                         selected = NULL)
  }, once = FALSE)
  
  observeEvent(input$script, {
    
    if(is.null(input$fcs$datapath)) {
      id <- showNotification("Please upload .fcs files.",type = "error", duration = 10)
      ids <<- c(ids, id)
      n <<- n + 1
    }
    
    else if(is.null(input$dir_btn)) {
      id <- showNotification("Please select a directory.",type = "error", duration = 10)
      ids <<- c(ids, id)
      n <<- n + 1
    }
    
    else if(is.null(input$name)) {
      id <- showNotification("Please name your analysis.",type = "error", duration = 10)
      ids <<- c(ids, id)
      n <<- n + 1
    }
    
    else if(is.null(input$csv$datapath)) {
      id <- showNotification("Please upload .csv file.",type = "error", duration = 10)
      ids <<- c(ids, id)
      n <<- n + 1
    }
    
    else if(length(input$plotted_markers) %% 2 != 0) {
      id <- showNotification("Please enter an even number of markers for gating.",type = "error", duration = 10)
      ids <<- c(ids, id)
      n <<- n + 1
    }
    
    if(length(input$plotted_markers) < 1) {
      id <- showNotification("Not enough markers for gating.",type = "error", duration = 10)
      ids <<- c(ids, id)
      n <<- n + 1
    }
    
    else{
     # with progress -> no fcs files present in directory, data has been preprocessed
   
      if (length(ids) > 0) { 
        removeNotification(ids[1])
        ids <<- ids[-1]
      }
      
      path <- paste(fcsdir(), "/", sep = "")
      if (is.null(input$fcs)) return()
      
      #Copy all the files to the newly created directory
      file.copy(from = input$fcs$datapath, paste0(path, input$fcs$name))
      file.copy(from = input$csv$datapath, paste0(path, input$csv$name))
      
      withProgress(message = 'Running... (this may take a while)',
                   value = 0, {
                     
                     ## Replace with path to the bash script
                     cmd <- paste("./parameters.sh ", 
                                  fcsdir(), #newly created directory for files
                                  input$corr_val, #correlation
                                  input$alpha, #alpha
                                  input$csv$name, #name of the metadata csv file
                                  markerParse(), #selected markers with string 'xyz' to separate them
                                  input$numerator, #with respect to numerator / not
                                  input$arcsinh, #true or false for arcsinh
                                  input$name) #run name
                     incProgress(0.5)
                     system(cmd)
                         
                     invalidateLater(1000, session)
                     fcsfiles <- list.files(path = path, pattern = "\\.fcs$")
                     rdsfiles <- list.files(pattern = "\\.RDS$")
                     
                     if(length(fcsfiles) == 0) {
                       incProgress(0.3, detail = "Beginning Gating")
                     }
                     
                     if(length(rdsfiles) == 1) {
                       incProgress(0.2, detail = "Gating complete")
                       showNotification("Success, you may now close the window", type = "message",  duration = 10)
                     }
                        
                   }) ##withProgress
      
    }

      
  }) ### tab 1 end
  
  
  ###### SECOND TAB
  
  gated_obj <- reactive({
    req(input$rds)
    rds_file <- readRDS(input$rds$datapath)
    full_obj <- rds_file
    

    return(full_obj)
  }) 
  
  conditions <- reactive({
    req(gated_obj())
    if("C2" %in% colnames(gated_obj()$obs)) {
      return(c("C1", "C2"))
    }
    else{return(c("C1"))}
  })
  
  
  
  channel_RDS <- reactive({
    req(gated_obj())
    channel_names <- gated_obj()$markerTable
    
    print(channel_names[1,])
    return(channel_names)
  })
  
  plotter <- reactive({
    req(input$markers_to_plot)
    m_to_plot <- which(channel_RDS()[1,] %in% input$markers_to_plot)
    #print(m_to_plot)
    observations <- gated_obj()$obs#[, names(gated_obj()$obs) %in% channel_RDS()[2,]]
    #channel_to_plot <- channel_RDS()[2, match(m_to_plot, channel_RDS()[1,])]
    channel_to_plot <- channel_RDS()[2, m_to_plot]
    print(channel_to_plot)
    #print(to_plot)
    
    cdn <- conditions()
    
    
    if(length(conditions()) < 2) {
      ggplot(observations, aes(x=(observations[,channel_to_plot]), color=gated_obj()$obs[,cdn])) + geom_density() + theme_classic() + xlab(input$markers_to_plot) + ylab("Density") + scale_colour_discrete("Conditions") + ggtitle(paste("Density of", input$markers_to_plot)) + theme(plot.title = element_text(hjust = 0.5))
    }
    
    else{
      
      cdn_1_obs <- subset(observations, select= -C2) #omit condition 2 from cdn1
      cdn_2_obs <- subset(observations, select = -C1)
      
      print(head(cdn_1_obs))
      print(head(cdn_2_obs))
      
      #ggplot(observations, aes(x=(observations[,channel_to_plot]), color=gated_obj()$obs[,cdn])) + 
      g <- ggplot()
      g <- g + geom_density(data = observations, aes(observations[,channel_to_plot], colour= gated_obj()$obs[,'C1']))          
      g <- g + geom_density(data = observations, aes(observations[,channel_to_plot], colour= gated_obj()$obs[,'C2'])) + theme_classic() + xlab(input$markers_to_plot) + ylab("Density") + scale_colour_discrete("Conditions") +   ggtitle(paste("Density of", input$markers_to_plot)) + theme(plot.title = element_text(hjust = 0.5))}
    
    g
    
  })  
  


  

  observeEvent(input$markers_to_plot, {

    output$plot1 <- renderPlot({

        plotter()

      })
   })

  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("CytofPlot-", Sys.Date(), ".pdf", sep="")
    },
    content = function(file) {
      ggsave(file, plotter(), device = "pdf", width=11, height=8.5)    
    }
  )
  
  
  observeEvent(input$rds, {

    updateSelectizeInput(session, 'markers_to_plot',
                         choices = as.character(channel_RDS()[1,]),
                         server = TRUE,
                         selected = NULL)
  }, once = FALSE)

  
  # output$markers_to_plot_ui <- renderUI({
  #   selectInput('markers_to_plot', NULL, choices = as.character(channel_RDS()[1,]),
  #                selected = as.vector(channel_RDS()[1,])[1], multiple = F ) })
  
  
  # output$image<-renderUI({
  #   
  # })
  
  
} ##### End of Server

shinyApp(ui, server)