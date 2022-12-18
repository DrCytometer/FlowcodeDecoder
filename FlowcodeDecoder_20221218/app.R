# assume group are defined in sample name as xxx group_xxxxxx.csv
# 0 or at least 2 PrmOIs should be choosen to present PrmOI table to be converted as vector => should correct the code to avoid data.frame subset to vector

# launch from command line: R -e "shiny::runApp(path,  launch.browser = TRUE)"

required_packages = c( "shiny", "shinyFiles", "shinyWidgets", "shinybusy", "htmlwidgets", "ggplot2", "gridExtra", "scales", "dplyr", "tidyr", "gtools", "ggrepel", "DT", "rhandsontable")

for (req.package in required_packages){
    if(!requireNamespace(req.package, quietly=TRUE)){
      install.packages(req.package, repos='http://cran.us.r-project.org')
    }
  }


library(shiny)
library(shinyFiles)
library(shinyWidgets)
library(shinybusy)
library(htmlwidgets)
library(ggplot2)
library(gridExtra)
library(scales)
library(dplyr)
library(tidyr)
options(dplyr.summarise.inform = FALSE)
options(shiny.maxRequestSize=1000*1024^2)
library(gtools)
library(ggrepel)
library(DT)
library(rhandsontable)
library(pdftools)



## Shiny app


# UI --------------------
ui <- navbarPage("Flowcode decoder",
                 
 # Files selection =========================

    tabPanel("Files selection",
             
             h3("Flow files to process"),
             hr(),
             fluidRow(column(1),
                      column(4, h4("Use FlowJo population export with csv channel values option")),
                      column(1),
                      column(5, fileInput("flow_files",  label = tags$span(
                        "csv files", 
                        tags$i(
                          class = "fa fa-info-circle", 
                          style = "color:#0072B2;font-size: 12px",
                          title = "0-1023 scale\nexpected format from FlowJo population export with csv channel values option"
                        )), multiple = TRUE, accept = c(".csv")))),
             br(),
             br(),
             br(),
             h3("Procode Combination file"),
             hr(),
             fluidRow(column(1),
                      column(4, h4("Provide a file that specify Id associated to each Procode combination")),
                      column(1),
                      column(5, fileInput("comb_file",  label = tags$span(
                        "csv file", 
                        tags$i(
                          class = "fa fa-info-circle", 
                          style = "color:#0072B2;font-size: 12px",
                          title = "columns should be Id,	Procode_tag1,	Procode_tag2,	Procode_tag3"
                        )), multiple = FALSE, accept = c(".csv")))),
             hr(),

             
             
             ),
 
 # Parameters of interest panel =========================
 
 tabPanel("Parameters of interest",
          fluidRow(column(10,  h3("Select parameters of interest",tags$i(
            class = "fa fa-info-circle", 
            style = "color:#0072B2;font-size: 14px",
            title = "generally exclude parameters that have been used to define the total populations (e.g. viability, CD45)"
          )))),
          hr(),
          fluidRow(column(2, uiOutput("PrmOIs_selector"),offset = 1), 
                   column(1), 
                   column(4, h4(textOutput("PrmOIs_count")),  style="padding:50px;"))
 ),
 
 
 # Threshold setting panel =========================

    tabPanel("Threshold setting",
             
             fluidRow(column(10, h3("Adjust positive population threshold for each parameter",  tags$i(
               class = "fa fa-info-circle", 
               style = "color:#0072B2;font-size: 12px",
               title = "click on a graph to setup threshold value for the x-axis parameter"
             )))),
             hr(),
             br(),
             fluidRow(
               column(6),
               column(2, h4("Import thresholds from file")), column(3, fileInput("imported_threshold_file",  label =  tags$span(
                 "thresholds.csv file", 
                 tags$i(
                   class = "fa fa-info-circle", 
                   style = "color:#0072B2;font-size: 12px",
                   title = "should contains parameter,value\nautomatically generated when a FlowCode decoder analysis has been performed"
                 )), multiple = FALSE, accept = c("thresholds.csv")))),
             
             
             fluidRow(column(1), column(2, uiOutput("threshold_plot_yvar_selector")) ),
             fluidRow(column(1),
                      column(2, uiOutput("threshold_selector")),
                      column(4),
                      column(6, plotOutput("threshold_plot", click="threshold_plot_click",  height="250px")))
                      
             ),
 
 # Analysis panel =========================
    tabPanel("Analysis",
             
    
             
             
             hr(),
             fluidRow(column(1),
                      column(2, h3("Analysis name:"))),
             fluidRow(column(1),
                      column(6, textInput("analysis_name", label = "", value = paste0("FlowcodeDecoder_", format(Sys.time(), "%Y%m%d%H%M"))))),
             hr(),
             fluidRow(column(1),
                      column(2, h3("Save analysis in:"))),
             fluidRow(column(1),
                      column(6, verbatimTextOutput("output_dir_txt")),
                      column(2, shinyDirButton("output_dir_selector", "Browse...", "Select location for ouputs") )
             ),
             
             hr(),
             
             br(),br(),br(),
             fluidRow(column(3, offset = 5, actionButton("run", h4("Run analysis")))) 
             
             
      ),


)
 # ====================


# Server --------------------
server <- function(input, output, session) {

  
 reactV<- reactiveValues()
 
 # Files selection =========================
 
 observeEvent(input$flow_files, {
     reactV$channel_values_matrix = do.call(rbind, lapply(input$flow_files$datapath, read.csv)) # 0-1023 scale ; expected format from FlowJo export population with csv channel values option
     reactV$perSample_cellcount = sapply(lapply(input$flow_files$datapath, read.csv),nrow)
     reactV$rowData = data.frame(sample = rep(sub(".csv$", "",input$flow_files$name),  reactV$perSample_cellcount))
     reactV$colData = data.frame(original_matrix_colnames = colnames(reactV$channel_values_matrix))
     
     if(sum(grepl( "\\.\\.\\.\\.", reactV$colData$original_matrix_colnames ))){ # signal_origin :: Id ; expected format from FlowJo export population with csv include header "Both" stain and parameter
       reactV$colData$signal_origin =  gsub("\\.\\.\\.\\..*", "",  reactV$colData$original_matrix_colnames)
       reactV$colData$parameterId =  gsub(".*\\.\\.\\.\\.", "",  reactV$colData$original_matrix_colnames)
       reactV$colData$is_labeled_parameter = !(reactV$colData$parameterId == reactV$colData$original_matrix_colnames)
     }else{ # for compatibility for non formatted header
       reactV$colData$parameterId = reactV$colData$original_matrix_colnames
     }
     colnames(reactV$channel_values_matrix) = reactV$colData$parameterId
     
     if(!is.null(input$comb_file)){
       if(all(reactV$procodeIds %in% reactV$colData$parameterId) == FALSE){
         report_warning("all Procode Ids do not match flow data channel names", 
                        paste( paste(reactV$procodeIds[reactV$procodeIds %in% reactV$colData$parameterId == FALSE], collapse = ","),
                               "not found among channel names:",
                               paste(reactV$colData$parameterId, collapse = ",")))
       }else{
         reactV$labeled_parameters = sort(reactV$colData$parameterId)
         reactV$preselect = reactV$procodeIds
         report_success("Procode channels recognized", "Set thresholds in the next tab")}
     }
     
     
   })
  
  observeEvent(input$comb_file, {
    reactV$combdf = read.csv(input$comb_file$datapath)
    reactV$procodeIds = unique(unlist(reactV$combdf[,-1]))
    
    if(!is.null(input$flow_files)){
        if(all(reactV$procodeIds %in% reactV$colData$parameterId) == FALSE){
          report_warning("all Procode Ids do not match flow data channel names", 
                         paste( paste(reactV$procodeIds[reactV$procodeIds %in% reactV$colData$parameterId == FALSE], collapse = ","),
                                "not found among channel names:",
                                paste(reactV$colData$parameterId, collapse = ",")))
        }else{
          reactV$labeled_parameters = sort(reactV$colData$parameterId)
          reactV$preselect = reactV$procodeIds
          report_success("Procode channels recognized", "Set thresholds in the next tab")}
    }
  
    })
  
  # Parameters of interest selection ============
  


        # Output dynamic group of buttons to select parameters of interest (PrmOIs) ####################
        
        reactV$preselect = c()

        output$PrmOIs_selector <- renderUI({
          if(is.null(reactV$labeled_parameters))
            return(NULL)
          checkboxGroupButtons("PrmOIs", label = "", 
                               choices = reactV$labeled_parameters ,
                               justified = TRUE,
                               selected = reactV$preselect,
                               direction = "vertical",
                               checkIcon = list(yes = icon("ok",lib = "glyphicon"))
          )
        })
        
        
        # Get PrmOIs from input #################
        observeEvent(input$PrmOIs,{
          reactV$PrmOIs = input$PrmOIs
          reactV$selectedPrmOI = reactV$PrmOIs[1] # initialization
        })
        
        # Output PrmOIs count ###################
        output$PrmOIs_count <- renderText({
          if(is.null(reactV$PrmOIs))
            return(NULL)
          paste(length(reactV$PrmOIs), "parameters selected")
        })
        
        
        
  
  
  
  # Threshold setting panel =========================
    
      reactV$thresholds = data.frame(parameter = character(),  value = numeric())
      
      # Output selector for threshold ###############
      
      output$threshold_selector <- renderUI({
        if(is.null(reactV$procodeIds))
          return(NULL)
        
        PrmOIs = input$PrmOIs
        
        
        choiceNames_list = list()
        for(i in 1:length(PrmOIs)){
          if(PrmOIs[i] %in% reactV$thresholds$parameter){
            choiceNames_list[[i]] =  tags$span(style = "color:#28b78d", PrmOIs[i])
          }else { choiceNames_list[[i]] =  tags$span(style = "color:black", PrmOIs[i])}
        }
        
        
        radioGroupButtons("Selected_PrmOI_for_threshold", label = "x-axis: threshold parameter", 
                          choiceNames = choiceNames_list,
                          choiceValues = PrmOIs,
                          justified = TRUE,
                          direction = "vertical", 
                          selected = reactV$selectedPrmOI
    
        )
      })
      
     observeEvent(input$Selected_PrmOI_for_threshold, {
         reactV$selectedPrmOI = input$Selected_PrmOI_for_threshold
     })
      
      # Output threshold_plot_yvar ###############
      
      output$threshold_plot_yvar_selector<- renderUI({
        if(is.null(reactV$procodeIds))
          return(NULL)
        selectInput("threshold_plot_yvar", label = tags$span("y-axis: contrast parameter", tags$i(
          class = "fa fa-info-circle", 
          style = "color:#0072B2;font-size: 12px",
          title = "change this parameter to better discriminate the x-axis populations"
        )),  choices = sort(colnames(reactV$channel_values_matrix) ), selected = "SSC.A")
      })
      

      # Output threshold_plot ############
      
      plot_threshold_plot <- function(df, title){
        xvar = colnames(df)[1]
        yvar = colnames(df)[2]
        xscale = c(0,1024)
        yscale = c(0,1024)
        
        g = ggplot(df, aes_string(x = xvar, y = yvar) )+
          ggtitle(title)+
          scale_x_continuous(limits = xscale)+  #labels =  math_format(10^.x) , breaks = seq(xscale[1],xscale[2]) 
          scale_y_continuous(limits = yscale)+
          geom_bin2d( aes_string(x = xvar, y = yvar), binwidth = c(2, 5)) 
        
        if(max(ggplot_build(g)$data[[1]]$count) == 1){
          g = g +scale_fill_gradientn(trans = "log10", colours="#0000FF")
        }else{
          g = g +scale_fill_gradientn(trans = "log10", colours=rev(rainbow(10, end = 4/6))) 
        } 
        
        g = g +
          geom_vline(xintercept = reactV$thresholds[reactV$thresholds$parameter == xvar,"value"], size = 1)+
          geom_hline(yintercept = reactV$thresholds[reactV$thresholds$parameter == yvar,"value"], size = 0.5)+
          theme_bw() +
          theme(plot.title = element_text(face="bold"), legend.position = "none")
        
        g
      }   
      
      output$threshold_plot <- renderPlot({
        if(is.null(input$Selected_PrmOI_for_threshold))  
          return(NULL)
        xvar = input$Selected_PrmOI_for_threshold
        yvar = input$threshold_plot_yvar
        df = as.data.frame(reactV$channel_values_matrix[,c(xvar,yvar)])
        plot_threshold_plot(df,"merged samples")
      })
      
      observeEvent(input$threshold_plot_click, {
          if(input$Selected_PrmOI_for_threshold %in% reactV$thresholds$parameter){
            reactV$thresholds[reactV$thresholds$parameter == input$Selected_PrmOI_for_threshold,"value"] = input$threshold_plot_click$x 
          } else{
            new_row = data.frame(parameter = input$Selected_PrmOI_for_threshold,  value = input$threshold_plot_click$x )
            reactV$thresholds = rbind(reactV$thresholds,new_row)
          }        
        })
      
      
      # Import thresholds #####################
      
      observeEvent(input$imported_threshold_file, {
        reactV$thresholds = read.csv(input$imported_threshold_file$datapath)
      })  
      
      # Success when all defined ############
      
      observeEvent(reactV$thresholds, {
        if(!is.null(input$PrmOIs)){
          if(all(length(reactV$thresholds$parameter) == length(input$PrmOIs))){ 
            if(all(sort(reactV$thresholds$parameter) == sort(input$PrmOIs))){
          report_success("All thresholds set", "Start analysis in the next tab")
            }
          }
        }
      })
      

    # Analysis panel =====================
      

     reactV$output_dir = getwd() 

    
      volumes = c(getVolumes()(), "Current working directory" = getwd()  )
      
      shinyDirChoose(input, id = "output_dir_selector", roots = volumes)
      
      observeEvent(input$output_dir_selector,{
        
        reactV$output_dir = parseDirPath(volumes, input$output_dir_selector)
      })
      
      output$output_dir_txt <- renderText({ reactV$output_dir})
      

      
      observeEvent(input$run, {
        
        dir.create(paste0(reactV$output_dir,"/",input$analysis_name))
        reactV$output_dir = paste0(reactV$output_dir,"/",input$analysis_name)
        

        show_modal_spinner(spin = "fading-circle")  
        
        thresholds = setNames(reactV$thresholds$value, reactV$thresholds$parameter) 

        
       #convert combination string data frame to binary #######
        update_modal_spinner(text = "convert combination string data frame to binary")
        
        reactV$bcombdf = as.data.frame(t(apply(reactV$combdf, 1, function(x) reactV$procodeIds %in% x[-1]))*1)
        colnames(reactV$bcombdf) = reactV$procodeIds
        rownames(reactV$bcombdf) = reactV$combdf$Id
        

        #get abovethreshold_matrix ######
        update_modal_spinner(text = "abovethreshold_matrix")
        
        thresholds = setNames(reactV$thresholds$value, reactV$thresholds$parameter)
        aboveThreshold_matrix  = t(apply(reactV$channel_values_matrix[,names(thresholds)],1,function(x) x>thresholds ))*1
        
        reactV$rowData$abovethreshold_ProcodeCount = rowSums(aboveThreshold_matrix[, reactV$procodeIds]  )
        
      
        #discard signal >3rd if signaloverthreshold of 3rd is at least 2times higher than 4th #######
        update_modal_spinner(text = "discard signal >3rd if signaloverthreshold of 3rd is at least 2times higher than 4th")
        signaloverthreshold_matrix = t(apply(reactV$channel_values_matrix[,reactV$procodeIds], 1, function(x) x-thresholds[reactV$procodeIds]))
        
        for(i in which(reactV$rowData$abovethreshold_ProcodeCount>3)){
          x = sort(unlist(signaloverthreshold_matrix[i,]), decreasing = TRUE)
          if(x[3]>2*x[4]){
            aboveThreshold_matrix[i, names(x)[4:length(x)]] = 0
          }
        }
        
        # link Procode_combination to Id#########
        update_modal_spinner(text = "link Procode_combination to Id")
        reactV$rowData$Procode_combination = apply(aboveThreshold_matrix[,reactV$procodeIds],1,function(x) paste(colnames(aboveThreshold_matrix[,reactV$procodeIds])[which(as.logical(x))], collapse =  "_"))
        reactV$combdf$Procode_combination = apply(reactV$combdf, 1, function(x) paste(sort(unlist(x[-1])),collapse = "_"))
        
        reactV$rowData = left_join(reactV$rowData, reactV$combdf[, c(1,ncol(reactV$combdf))], by = c("Procode_combination"))
        
        
        reactV$rowData$Id[is.na(reactV$rowData$Id)&reactV$rowData$abovethreshold_ProcodeCount>=3] = "3 or more unexpected signals"
        reactV$rowData$Id[is.na(reactV$rowData$Id)&reactV$rowData$abovethreshold_ProcodeCount %in% 1:2] = "1 or 2 unexpected signals"
        reactV$rowData$Id[is.na(reactV$rowData$Id)&reactV$rowData$abovethreshold_ProcodeCount == 0] = "no procode signal"
        
        reactV$rowData$Id = factor(reactV$rowData$Id, levels = sort( unique( c(reactV$rowData$Id, reactV$combdf$Id) ) ))
        
        # stats per sample #########
        update_modal_spinner(text = "get stats per sample")
        
        
        reactV$perSample_stats = reactV$rowData %>%
                            group_by(sample, Id) %>%
                            summarize(count = n()) %>%
                            mutate(pct = round(count/sum(count)*100, 2)) %>%
                            complete(Id, fill = list(count = 0, pct = 0)) %>%
                            as.data.frame()
        
        reactV$count_matrix = reactV$perSample_stats %>%
                                  select(sample, Id, count) %>%
                                  pivot_wider(names_from = Id, values_from = count)%>%
                                  as.data.frame()
        
        reactV$pct_matrix = reactV$perSample_stats %>%
                                  select(sample, Id, pct) %>%
                                  pivot_wider(names_from = Id, values_from = pct)%>%
                                  as.data.frame()
        
        
        
        # decoder QC #######
        update_modal_spinner(text = "get decoder QC")
        
        reactV$decoder_QC = reactV$perSample_stats %>%
                                filter(Id %in% c("3 or more unexpected signals", "1 or 2 unexpected signals",  "no procode signal"))  %>%
                                select(sample,Id,pct) %>%
                                pivot_wider( names_from = Id, values_from = pct, values_fill = 0) %>%
                                mutate(expected_signals = 100-`no procode signal`- `1 or 2 unexpected signals` - `3 or more unexpected signals`,
                                       expected_signals_over_all_signals = expected_signals/(expected_signals+`1 or 2 unexpected signals`+`3 or more unexpected signals`)*100, 
                                       expected_signals_over_3plus_signals = expected_signals/(expected_signals+`3 or more unexpected signals`)*100)
        
        colnames(reactV$decoder_QC)[-1] = paste0(colnames(reactV$decoder_QC)[-1], "_pct")
        
        
        expected_signals_count = reactV$perSample_stats %>%
                                filter(!Id %in% c("3 or more unexpected signals_pct", "1 or 2 unexpected signals_pct",  "no procode signal_pct")) %>%
                                    group_by(sample)%>%
                                    summarize(expected_signals_count = sum(count))
        
        total_count = reactV$perSample_stats %>%
          group_by(sample)%>%
          summarize(total_count = sum(count))
        
        
        reactV$decoder_QC = left_join(reactV$decoder_QC,  expected_signals_count, by = c("sample"))
        reactV$decoder_QC = left_join(reactV$decoder_QC,  total_count, by = c("sample"))
        
        reactV$decoder_QC = as.data.frame(reactV$decoder_QC[,c("sample",
                                                               "total_count", 
                                                               "expected_signals_count",
                                                               "expected_signals_pct",
                                                               "3 or more unexpected signals_pct",
                                                               "1 or 2 unexpected signals_pct",
                                                               "no procode signal_pct",
                                                               "expected_signals_over_all_signals_pct",
                                                               "expected_signals_over_3plus_signals_pct")])
        
        
        # troubleshoot tables #####
        update_modal_spinner(text = "get troubleshoot tables")
        
        reactV$perSample_unexpected_3plus = reactV$rowData %>%
                                    filter(Id == "3 or more unexpected signals")%>%
                                    group_by(sample, Procode_combination) %>%
                                    summarize(count = n()) %>%
                                    mutate(pct = round(count/sum(count)*100, 2)) %>%
                                    complete(Procode_combination, fill = list(count = 0, pct = 0)) %>%
                                    as.data.frame()
        
        reactV$unexpected_3plus_count_matrix = reactV$perSample_unexpected_3plus %>%
          select(sample, Procode_combination, count) %>%
          pivot_wider(names_from = Procode_combination, values_from = count, values_fill = 0)%>%
          as.data.frame()
        
        reactV$unexpected_3plus_pct_matrix = reactV$perSample_unexpected_3plus %>%
          select(sample, Procode_combination, pct) %>%
          pivot_wider(names_from = Procode_combination, values_from = pct, values_fill = 0)%>%
          as.data.frame()
        
        
        
        
        reactV$perSample_unexpected_1n2 = reactV$rowData %>%
          filter(Id == "1 or 2 unexpected signals")%>%
          group_by(sample, Procode_combination) %>%
          summarize(count = n()) %>%
          mutate(pct = round(count/sum(count)*100, 2)) %>%
          complete(Procode_combination, fill = list(count = 0, pct = 0)) %>%
          as.data.frame()
        
        reactV$unexpected_1n2_count_matrix = reactV$perSample_unexpected_1n2 %>%
          select(sample, Procode_combination, count) %>%
          pivot_wider(names_from = Procode_combination, values_from = count, values_fill = 0)%>%
          as.data.frame()
        
        reactV$unexpected_1n2_pct_matrix = reactV$perSample_unexpected_1n2 %>%
          select(sample, Procode_combination, pct) %>%
          pivot_wider(names_from = Procode_combination, values_from = pct, values_fill = 0)%>%
          as.data.frame()
        
        
        
        # Other PrmOIs analysis ########
        update_modal_spinner(text = "Other PrmOIs analysis")
        
      
        if(ncol(aboveThreshold_matrix) > length(reactV$procodeIds)  ){
          
          otherPrmOIs = colnames(aboveThreshold_matrix)[!(colnames(aboveThreshold_matrix) %in% reactV$procodeIds)]
          
          
          
          perCell_otherPrmOIs_positiveness = cbind(reactV$rowData, aboveThreshold_matrix[,otherPrmOIs])
          perCell_otherPrmOIs_positiveness$cellId = 1:nrow(perCell_otherPrmOIs_positiveness)

          
          perSample_otherPrmOIs_stats = perCell_otherPrmOIs_positiveness %>% 
                                  pivot_longer(cols = otherPrmOIs,  names_to = "PrmOI", values_to = "positiveness" ) %>% 
                                  group_by(sample, Id, PrmOI ) %>%
                                  summarise(n = n()) %>%
                                  ungroup() %>%
                                  complete(sample,Id,PrmOI  ,  fill = list(n = 0))
          pstat = perCell_otherPrmOIs_positiveness %>% 
                        pivot_longer(cols = otherPrmOIs,  names_to = "PrmOI", values_to = "positiveness" ) %>%
                        filter(positiveness == 1) %>%
                        group_by(sample, Id, PrmOI ) %>%
                        summarise(npos = n()) %>%
                        ungroup() %>%
                        complete(sample,Id,PrmOI  ,  fill = list(npos = 0))  %>%
                        select(npos)
          
          perSample_otherPrmOIs_stats = perSample_otherPrmOIs_stats %>%
                                  mutate(npos = pstat$npos, 
                                         pctpos = npos / n*100) %>%
                                  arrange(PrmOI)
          
          perSample_otherPrmOIs_stats$pctpos[is.na(perSample_otherPrmOIs_stats$pctpos)] = 0
          

          ### define groups ###
          update_modal_spinner(text = "define groups")
          
          
          perSample_otherPrmOIs_stats$group = unname(t(as.data.frame(strsplit(unname(t(as.data.frame((strsplit(perSample_otherPrmOIs_stats$sample, split = " "))))[,2]), split = "_")))[,1])
          
          perGroup_otherPrmOIs_stats = perSample_otherPrmOIs_stats %>%
                                                group_by(group, Id, PrmOI) %>%
                                                summarise(N = n(),
                                                          mean_pctpos = mean(pctpos),
                                                          sd_pctpos = sd(pctpos)) %>%
                                                arrange(PrmOI)
          
          
            }

        

        


        #Graphical output #######
        update_modal_spinner(text = "create graphical output")
        
        dir.create(paste0(reactV$output_dir, "/tmp"))
          
        
        perCell_data = cbind(reactV$rowData, reactV$channel_values_matrix)
        perCell_data$Id = factor(perCell_data$Id, levels = sort(unique(perCell_data$Id)))
        perCell_data$Id_int = as.numeric(perCell_data$Id)
        perCell_data$Id_int_jit = perCell_data$Id_int + runif(length(perCell_data$Id_int), min = -0.2, max = +0.2)
        
        POIs = colnames(reactV$channel_values_matrix)[!(colnames(reactV$channel_values_matrix) %in% c("FSC.A", "FSC.H", "FSC.W", "SSC.A", "SSC.H", "SSC.W", "SSC.B.A", "SSC.B.H", "SSC.B.W", "Time", "AF.A") )]
        
        for( poi in POIs){
          g = ggplot(perCell_data)+
            scale_x_continuous(limits = c(0,1024))+
            scale_y_continuous(breaks = 1:length(levels(perCell_data$Id)), labels = levels(perCell_data$Id) )+
            geom_bin2d( aes_string(x = poi, y = "Id_int_jit"), binwidth = c(5, 0.05), na.rm = TRUE) +
            scale_fill_gradientn(trans = "log10", colours=rev(rainbow(10, end = 4/6))) +
            theme_bw()+
            facet_grid(.~sample)+
            theme(strip.text.x = element_text(size=8, angle=90),
                  panel.grid.minor = element_blank(), 
                  axis.text.x=element_text(angle=45, hjust=1))
          ggsave(paste0(reactV$output_dir, "/tmp/", poi, ".pdf"), g, units = "mm", width = 29*length(unique(perCell_data$sample)), height = length(levels(perCell_data$Id))*8 , limitsize = FALSE ) 
        }
        pdf_combine(paste0(reactV$output_dir, "/tmp/", POIs, ".pdf"), output = paste0(reactV$output_dir,"/markers_graphs.pdf"))
        
        
        #Threshold graphs ####
        update_modal_spinner(text = "create threholds graphs")
        
        plot_threshold_plot <- function(df, xvar){
          yvar = "SSC.A"
          ggplot(df, aes_string(x = xvar, y = yvar) )+
            scale_x_continuous(limits = c(0,1024) )+
            scale_y_continuous(limits = c(0,1024) )+  
            geom_bin2d( aes_string(x = xvar, y = yvar), binwidth = c(3, 30))+ 
            geom_vline(xintercept = thresholds[xvar], size = 1)+
            scale_fill_gradientn(trans = "log10", colours=rev(rainbow(10, end = 4/6)))+
            theme_bw()+
            guides(fill="none")
        }
        

        for( sample in unique(perCell_data$sample)){
          
          threshold_graphs =  vector(mode = "list", length = length(thresholds))
          names(threshold_graphs) = names(thresholds)
          
          for(procode in names(thresholds)){
            threshold_graphs[[procode]] = plot_threshold_plot(perCell_data[perCell_data$sample == sample,], procode)
          }
          
          #grid.arrange(grobs = threshold_graphs, ncol = 1, top = paste0(sample," threshold graphs") ) ## display plot
          ggsave(file = paste0(reactV$output_dir, "/tmp/",sample,"_threshold_graphs.pdf"), arrangeGrob(grobs = threshold_graphs, ncol = 1, top = paste0(sample," threshold graphs")), width = 297, height = 35*length(thresholds), units = "mm")  ## save plot
          
        }
  
        pdf_combine( paste0(reactV$output_dir, "/tmp/",unique(perCell_data$sample),"_threshold_graphs.pdf"), output = paste0(reactV$output_dir,"/thresholds_graphs.pdf"))
        

        #### PrmOI #######
        update_modal_spinner(text = "create PrmOIs graphs")
        
        for(PrmOI in otherPrmOIs ){
          ggplot()+
            ggtitle(PrmOI)+
            geom_point(data = perSample_otherPrmOIs_stats[perSample_otherPrmOIs_stats$PrmOI == PrmOI,], aes(x = pctpos, y = Id)) +
            theme_bw()
          ggsave(file = paste0(reactV$output_dir, "/tmp/",PrmOI,"_perSample.pdf"),   width = 297, height = 210 , units = "mm")  ## save plot
        }
        
        pdf_combine( paste0(reactV$output_dir, "/tmp/",otherPrmOIs,"_perSample.pdf"), output = paste0(reactV$output_dir,"/PrmOI_perSample.pdf"))
        
        
        
        for(PrmOI in otherPrmOIs ){
          ggplot()+
            ggtitle(PrmOI)+
            geom_bar( stat = "identity", data = perGroup_otherPrmOIs_stats[perGroup_otherPrmOIs_stats$PrmOI == PrmOI,], aes(y = mean_pctpos, x = Id ))+
            geom_errorbar(data = perGroup_otherPrmOIs_stats[perGroup_otherPrmOIs_stats$PrmOI == PrmOI,], aes(x = Id, ymin=mean_pctpos-1.96*sd_pctpos/sqrt(N), ymax=mean_pctpos+1.96*sd_pctpos/sqrt(N)), width=.2) +
            geom_jitter(data = perSample_otherPrmOIs_stats[perSample_otherPrmOIs_stats$PrmOI == PrmOI,], aes(y = pctpos, x = Id), height = 0, width = 0.1, shape = 21) +
            theme_bw()+
            coord_flip()+
            facet_grid(. ~ group)
          ggsave(file = paste0(reactV$output_dir, "/tmp/",PrmOI,"_perGroup.pdf"),   width = 297, height = 210 , units = "mm")  ## save plot
        }
        
        pdf_combine( paste0(reactV$output_dir, "/tmp/",otherPrmOIs,"_perGroup.pdf"), output = paste0(reactV$output_dir,"/PrmOI_perGroup.pdf"))
        
        
        reactV$perSample_stats$group = sapply(strsplit(reactV$perSample_stats$sample, "_"), function(x) strsplit(x[1]," ")[[1]][2] )
        reactV$perGroup_stats = reactV$perSample_stats %>%    
          group_by(group, Id) %>%
          summarise(N = n(),
                    mean_count = mean(count),
                    sd_count = sd(count),
                    mean_pct = mean(pct),
                    sd_pct = sd(pct)) %>%
          arrange(group, Id)
        
        
        
        
        
        #save #######
        update_modal_spinner(text = "save tables")
        write.csv(perCell_data , paste0(reactV$output_dir, "/", input$analysis_name, "_perCell_data.csv"), row.names = FALSE)
        write.csv(perSample_otherPrmOIs_stats , paste0(reactV$output_dir, "/", input$analysis_name, "_perSample_otherPrmOIs_stats.csv"), row.names = FALSE)
        write.csv(perGroup_otherPrmOIs_stats , paste0(reactV$output_dir, "/", input$analysis_name, "_perGroup_otherPrmOIs_stats.csv"), row.names = FALSE)
        write.csv(reactV$perSample_stats , paste0(reactV$output_dir, "/", input$analysis_name, "_perSample_stats.csv"), row.names = FALSE)
        write.csv(reactV$perGroup_stats , paste0(reactV$output_dir, "/", input$analysis_name, "_perGroup_stats.csv"), row.names = FALSE)
        write.csv(reactV$count_matrix, paste0(reactV$output_dir, "/", input$analysis_name, "_count_matrix.csv"), row.names = FALSE)
        write.csv(reactV$pct_matrix, paste0(reactV$output_dir, "/", input$analysis_name, "_pct_matrix.csv"), row.names = FALSE)
        write.csv(reactV$decoder_QC, paste0(reactV$output_dir, "/", input$analysis_name, "_decoder_QC.csv"), row.names = FALSE)
        write.csv(reactV$unexpected_3plus_count_matrix, paste0(reactV$output_dir, "/", input$analysis_name, "_unexpected_3plus_count_matrix.csv"), row.names = FALSE)
        write.csv(reactV$unexpected_3plus_pct_matrix, paste0(reactV$output_dir, "/", input$analysis_name, "_unexpected_3plus_pct_matrix.csv"), row.names = FALSE)
        write.csv(reactV$unexpected_1n2_count_matrix, paste0(reactV$output_dir, "/", input$analysis_name, "_unexpected_1n2_count_matrix.csv"), row.names = FALSE)
        write.csv(reactV$unexpected_1n2_pct_matrix, paste0(reactV$output_dir, "/", input$analysis_name, "_unexpected_1n2_pct_matrix.csv"), row.names = FALSE)
        
        
        write.csv(reactV$thresholds, paste0(reactV$output_dir, "/", input$analysis_name, "_thresholds.csv"), row.names = FALSE)
        
        
        
        
        #end ####
        remove_modal_spinner()
        
        
        report_success("Analysis completed", paste("files saved in:", reactV$output_dir ))
        

        
        
        
        
        

      })
      
      
    
    

  }
 
 

    


# shinyApp------------

shinyApp(ui = ui, server = server)