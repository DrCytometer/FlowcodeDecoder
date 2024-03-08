start.time <- Sys.time()
# linearized flowcode decoder for remote computing

required_packages = c( "ggplot2", "gridExtra", "scales", "dplyr", "tidyr", 
                       "gtools", "ggrepel", "DT", "rhandsontable", "pdftools" )

for (req.package in required_packages){
  if(!requireNamespace(req.package, quietly=TRUE)){
    install.packages(req.package, repos='http://cran.us.r-project.org')
  }
}


library(ggplot2)
library(gridExtra)
library(scales)
library(dplyr)
library(tidyr)
options(dplyr.summarise.inform = FALSE)
library(gtools)
library(ggrepel)
library(DT)
library(rhandsontable)
library(pdftools)
#library(data.table)

# required inputs
input_dir <- "./Input"
data_dir <- "./Data"
output_dir <- "./Output"
dir.create(output_dir)
threshold_file <- list.files(input_dir)[grep("thresholds.csv", list.files(input_dir))]
comb_file <- list.files(input_dir)[grep("combinations.csv", list.files(input_dir))]
analysis_name <- paste0("FlowCodeDebarcoder_nonApp_", format(Sys.time(), "%Y%m%d_%H%M%S"))
flow_code_data_files <- list.files(data_dir)

# read threshold file
thresholds = read.csv(file = file.path(input_dir, threshold_file))

# procode combination file
combdf = read.csv(file = file.path(input_dir, comb_file))
procodeIds = unique(unlist(combdf[,-1]))


channel_values_matrix = do.call(rbind, lapply(file.path(data_dir, flow_code_data_files), read.csv))
perSample_cellcount = sapply(lapply(file.path(data_dir, flow_code_data_files), read.csv),nrow)
rowData = data.frame(sample = rep(sub(".csv$", "",flow_code_data_files),  perSample_cellcount))
colData = data.frame(original_matrix_colnames = colnames(channel_values_matrix))
if(sum(grepl( "\\.\\.\\.\\.", colData$original_matrix_colnames ))){ # signal_origin :: Id ; expected format from FlowJo export population with csv include header "Both" stain and parameter
  colData$signal_origin =  gsub("\\.\\.\\.\\..*", "",  colData$original_matrix_colnames)
  colData$parameterId =  gsub(".*\\.\\.\\.\\.", "",  colData$original_matrix_colnames)
  colData$is_labeled_parameter = !(colData$parameterId == colData$original_matrix_colnames)
}else{ # for compatibility for non formatted header
  colData$parameterId = colData$original_matrix_colnames
}
colnames(channel_values_matrix) = colData$parameterId
labeled_parameters = sort(colData$parameterId)

# analysis
thresholds = setNames(thresholds$value, thresholds$parameter) 

bcombdf = as.data.frame(t(apply(combdf, 1, function(x) procodeIds %in% x[-1]))*1)
colnames(bcombdf) = procodeIds
rownames(bcombdf) = combdf$Id
aboveThreshold_matrix  = t(apply(channel_values_matrix[,names(thresholds)],1,function(x) x>thresholds ))*1

rowData$abovethreshold_ProcodeCount = rowSums(aboveThreshold_matrix[, procodeIds]  )

signaloverthreshold_matrix = t(apply(channel_values_matrix[,procodeIds], 1, function(x) x-thresholds[procodeIds]))

for(i in which(rowData$abovethreshold_ProcodeCount>3)){
  x = sort(unlist(signaloverthreshold_matrix[i,]), decreasing = TRUE)
  if(x[3]>2*x[4]){
    aboveThreshold_matrix[i, names(x)[4:length(x)]] = 0
  }
}

#rowData$Procode_combination = apply(aboveThreshold_matrix[,sort(procodeIds)],1,
#                                    function(x) paste(colnames(aboveThreshold_matrix[,sort(procodeIds)])[which(as.logical(x))], collapse =  "_"))

# Vectorized operation
aboveThreshold_ProCode_matrix <- aboveThreshold_matrix[, sort(procodeIds)]
aboveThreshold_matrix_logical <- as.logical(aboveThreshold_ProCode_matrix)
dim(aboveThreshold_matrix_logical) <- dim(aboveThreshold_ProCode_matrix)
colnames(aboveThreshold_matrix_logical) <- colnames(aboveThreshold_ProCode_matrix)
rowData$Procode_combination <- apply(aboveThreshold_matrix_logical, 1, function(x) paste(names(which(x)), collapse = "_"))



combdf$Procode_combination = apply(combdf, 1, function(x) paste(sort(unlist(x[-1])),collapse = "_"))

rowData = left_join(rowData, combdf[, c(1,ncol(combdf))], by = c("Procode_combination"))


rowData$Id[is.na(rowData$Id)&rowData$abovethreshold_ProcodeCount>=3] = "3 or more unexpected signals"
rowData$Id[is.na(rowData$Id)&rowData$abovethreshold_ProcodeCount %in% 1:2] = "1 or 2 unexpected signals"
rowData$Id[is.na(rowData$Id)&rowData$abovethreshold_ProcodeCount == 0] = "no procode signal"

rowData$Id = factor(rowData$Id, levels = sort( unique( c(rowData$Id, combdf$Id) ) ))


perSample_stats = rowData %>%
  group_by(sample, Id) %>%
  summarize(count = n()) %>%
  mutate(pct = round(count/sum(count)*100, 2)) %>%
  complete(Id, fill = list(count = 0, pct = 0)) %>%
  as.data.frame()

count_matrix = perSample_stats %>%
  select(sample, Id, count) %>%
  pivot_wider(names_from = Id, values_from = count)%>%
  as.data.frame()

pct_matrix = perSample_stats %>%
  select(sample, Id, pct) %>%
  pivot_wider(names_from = Id, values_from = pct)%>%
  as.data.frame()

#QC
decoder_QC = perSample_stats %>%
  filter(Id %in% c("3 or more unexpected signals", "1 or 2 unexpected signals",  "no procode signal"))  %>%
  select(sample,Id,pct) %>%
  pivot_wider( names_from = Id, values_from = pct, values_fill = 0) %>%
  mutate(expected_signals = 100-`no procode signal`- `1 or 2 unexpected signals` - `3 or more unexpected signals`,
         expected_signals_over_all_signals = expected_signals/(expected_signals+`1 or 2 unexpected signals`+`3 or more unexpected signals`)*100, 
         expected_signals_over_3plus_signals = expected_signals/(expected_signals+`3 or more unexpected signals`)*100)

colnames(decoder_QC)[-1] = paste0(colnames(decoder_QC)[-1], "_pct")

expected_signals_count = perSample_stats %>%
  filter(!Id %in% c("3 or more unexpected signals", "1 or 2 unexpected signals",  "no procode signal")) %>%
  group_by(sample)%>%
  summarize(expected_signals_count = sum(count))

total_count = perSample_stats %>%
  group_by(sample)%>%
  summarize(total_count = sum(count))

decoder_QC = left_join(decoder_QC,  expected_signals_count, by = c("sample"))
decoder_QC = left_join(decoder_QC,  total_count, by = c("sample"))

decoder_QC = as.data.frame(decoder_QC[,c("sample",
                                                       "total_count", 
                                                       "expected_signals_count",
                                                       "expected_signals_pct",
                                                       "3 or more unexpected signals_pct",
                                                       "1 or 2 unexpected signals_pct",
                                                       "no procode signal_pct",
                                                       "expected_signals_over_all_signals_pct",
                                                       "expected_signals_over_3plus_signals_pct")])

perSample_unexpected_3plus = rowData %>%
  filter(Id == "3 or more unexpected signals")%>%
  group_by(sample, Procode_combination) %>%
  summarize(count = n()) %>%
  mutate(pct = round(count/sum(count)*100, 2)) %>%
  complete(Procode_combination, fill = list(count = 0, pct = 0)) %>%
  as.data.frame()

unexpected_3plus_count_matrix = perSample_unexpected_3plus %>%
  select(sample, Procode_combination, count) %>%
  pivot_wider(names_from = Procode_combination, values_from = count, values_fill = 0)%>%
  as.data.frame()

unexpected_3plus_pct_matrix = perSample_unexpected_3plus %>%
  select(sample, Procode_combination, pct) %>%
  pivot_wider(names_from = Procode_combination, values_from = pct, values_fill = 0)%>%
  as.data.frame()

perSample_unexpected_1n2 = rowData %>%
  filter(Id == "1 or 2 unexpected signals")%>%
  group_by(sample, Procode_combination) %>%
  summarize(count = n()) %>%
  mutate(pct = round(count/sum(count)*100, 2)) %>%
  complete(Procode_combination, fill = list(count = 0, pct = 0)) %>%
  as.data.frame()

unexpected_1n2_count_matrix = perSample_unexpected_1n2 %>%
  select(sample, Procode_combination, count) %>%
  pivot_wider(names_from = Procode_combination, values_from = count, values_fill = 0)%>%
  as.data.frame()

unexpected_1n2_pct_matrix = perSample_unexpected_1n2 %>%
  select(sample, Procode_combination, pct) %>%
  pivot_wider(names_from = Procode_combination, values_from = pct, values_fill = 0)%>%
  as.data.frame()

perCell_data = cbind(rowData, channel_values_matrix)
perCell_data$Id = factor(perCell_data$Id, levels = sort(unique(perCell_data$Id)))
perCell_data$Id_int = as.numeric(perCell_data$Id)
perCell_data$Id_int_jit = perCell_data$Id_int + runif(length(perCell_data$Id_int), min = -0.2, max = +0.2)

perSample_stats$group = sapply(strsplit(perSample_stats$sample, "_"), function(x) strsplit(x[1]," ")[[1]][2] )
perGroup_stats = perSample_stats %>%    
  group_by(group, Id) %>%
  summarise(N = n(),
            mean_count = mean(count),
            sd_count = sd(count),
            mean_pct = mean(pct),
            sd_pct = sd(pct)) %>%
  arrange(group, Id)

# Other PrmOIs analysis ########
if(ncol(aboveThreshold_matrix) > length(procodeIds)  ){
  
  otherPrmOIs = colnames(aboveThreshold_matrix)[!(colnames(aboveThreshold_matrix) %in% procodeIds)]
  
  perCell_otherPrmOIs_positiveness = cbind(rowData, aboveThreshold_matrix[,otherPrmOIs])
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
  perSample_otherPrmOIs_stats$group = unname(t(as.data.frame(strsplit(unname(t(as.data.frame((strsplit(perSample_otherPrmOIs_stats$sample, split = " "))))[,2]), split = "_")))[,1])
  
  perGroup_otherPrmOIs_stats = perSample_otherPrmOIs_stats %>%
    group_by(group, Id, PrmOI) %>%
    summarise(N = n(),
              mean_pctpos = mean(pctpos),
              sd_pctpos = sd(pctpos)) %>%
    arrange(PrmOI)
}
  
#save ######
write.csv(perCell_data , paste0(output_dir, "/", analysis_name, "_perCell_data.csv"), row.names = FALSE)
if(ncol(aboveThreshold_matrix) > length(procodeIds)  ){
  write.csv(perSample_otherPrmOIs_stats , paste0(output_dir, "/", analysis_name, "_perSample_otherPrmOIs_stats.csv"), row.names = FALSE)
  write.csv(perGroup_otherPrmOIs_stats , paste0(output_dir, "/", analysis_name, "_perGroup_otherPrmOIs_stats.csv"), row.names = FALSE)
}
write.csv(perSample_stats , paste0(output_dir, "/", analysis_name, "_perSample_stats.csv"), row.names = FALSE)
write.csv(perGroup_stats , paste0(output_dir, "/", analysis_name, "_perGroup_stats.csv"), row.names = FALSE)
write.csv(count_matrix, paste0(output_dir, "/", analysis_name, "_count_matrix.csv"), row.names = FALSE)
write.csv(pct_matrix, paste0(output_dir, "/", analysis_name, "_pct_matrix.csv"), row.names = FALSE)
write.csv(decoder_QC, paste0(output_dir, "/", analysis_name, "_decoder_QC.csv"), row.names = FALSE)
write.csv(unexpected_3plus_count_matrix, paste0(output_dir, "/", analysis_name, "_unexpected_3plus_count_matrix.csv"), row.names = FALSE)
write.csv(unexpected_3plus_pct_matrix, paste0(output_dir, "/", analysis_name, "_unexpected_3plus_pct_matrix.csv"), row.names = FALSE)
write.csv(unexpected_1n2_count_matrix, paste0(output_dir, "/", analysis_name, "_unexpected_1n2_count_matrix.csv"), row.names = FALSE)
write.csv(unexpected_1n2_pct_matrix, paste0(output_dir, "/", analysis_name, "_unexpected_1n2_pct_matrix.csv"), row.names = FALSE)
write.csv(thresholds, paste0(output_dir, "/", analysis_name, "_thresholds.csv"), row.names = FALSE)

# plotting will be slow with large datasets
#save plots######
dir.create(paste0(output_dir, "/tmp"))

perCell_data = cbind(rowData, channel_values_matrix)
perCell_data$Id = factor(perCell_data$Id, levels = sort(unique(perCell_data$Id)))
perCell_data$Id_int = as.numeric(perCell_data$Id)
perCell_data$Id_int_jit = perCell_data$Id_int + runif(length(perCell_data$Id_int), min = -0.2, max = +0.2)

POIs = colnames(channel_values_matrix)[!(colnames(channel_values_matrix) %in% c("FSC.A", "FSC.H", "FSC.W", "SSC.A", "SSC.H", "SSC.W", "SSC.B.A", "SSC.B.H", "SSC.B.W", "Time", "AF.A") )]

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
  ggsave(paste0(output_dir, "/tmp/", poi, ".pdf"), g, units = "mm", width = 29*length(unique(perCell_data$sample))*2, height = length(levels(perCell_data$Id))*8 , limitsize = FALSE ) 
}
pdf_combine(paste0(output_dir, "/tmp/", POIs, ".pdf"), output = paste0(output_dir,"/markers_graphs.pdf"))


#Threshold graphs ####

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
  
  ggsave(file = paste0(output_dir, "/tmp/",sample,"_threshold_graphs.pdf"), arrangeGrob(grobs = threshold_graphs, ncol = 1, top = paste0(sample," threshold graphs")), width = 297, height = 35*length(thresholds), units = "mm")  ## save plot
    
}

pdf_combine( paste0(output_dir, "/tmp/",unique(perCell_data$sample),"_threshold_graphs.pdf"), output = paste0(output_dir,"/thresholds_graphs.pdf"))


#### PrmOI plots#######
if(ncol(aboveThreshold_matrix) > length(procodeIds)  ){
  for(PrmOI in otherPrmOIs ){
      ggplot()+
        ggtitle(PrmOI)+
        geom_point(data = perSample_otherPrmOIs_stats[perSample_otherPrmOIs_stats$PrmOI == PrmOI,], aes(x = pctpos, y = Id)) +
        theme_bw()
      ggsave(file = paste0(output_dir, "/tmp/",PrmOI,"_perSample.pdf"),   width = 297, height = 210 , units = "mm")  ## save plot
  }
    
  pdf_combine( paste0(output_dir, "/tmp/",otherPrmOIs,"_perSample.pdf"), output = paste0(output_dir,"/PrmOI_perSample.pdf"))


  for(PrmOI in otherPrmOIs ){
      ggplot()+
        ggtitle(PrmOI)+
        geom_bar( stat = "identity", data = perGroup_otherPrmOIs_stats[perGroup_otherPrmOIs_stats$PrmOI == PrmOI,], aes(y = mean_pctpos, x = Id ))+
        geom_errorbar(data = perGroup_otherPrmOIs_stats[perGroup_otherPrmOIs_stats$PrmOI == PrmOI,], aes(x = Id, ymin=mean_pctpos-1.96*sd_pctpos/sqrt(N), ymax=mean_pctpos+1.96*sd_pctpos/sqrt(N)), width=.2) +
        geom_jitter(data = perSample_otherPrmOIs_stats[perSample_otherPrmOIs_stats$PrmOI == PrmOI,], aes(y = pctpos, x = Id), height = 0, width = 0.1, shape = 21) +
        theme_bw()+
        coord_flip()+
        facet_grid(. ~ group)
    ggsave(file = paste0(output_dir, "/tmp/",PrmOI,"_perGroup.pdf"),   width = 297, height = 210 , units = "mm")  ## save plot
  }
    
  pdf_combine( paste0(output_dir, "/tmp/",otherPrmOIs,"_perGroup.pdf"), output = paste0(output_dir,"/PrmOI_perGroup.pdf"))
}

calc.time <- Sys.time() - start.time
calc.time
