#!/usr/bin/env Rscript
#########
rm(list= ls())
options(warn = -1)
library(optparse)
option_list <- list(
  make_option(c("-s","--script_path"),type = 'character',default = NULL,
              help = 'Path to the script', metavar = "/SCRIPTS/PATH"),
  make_option(c("-i","--input_kreport_path"),type = 'character',default = NULL,
              help = 'Path to the krakenuniq report',metavar = "/INPUT/KRAKENUNIQ/REPORT/PATH"),
  make_option(c("-S","--study_name"),type = 'character',default = NULL,
              help = 'Study name',metavar = "PRJNA00000"),
  make_option(c("-u","--unique_kmer_threshold"),type = "integer",default = 1000,
              help = 'Minimal unique kmer number',metavar = '1000'),
  make_option(c("-f","--fold_threshold"),type = "integer",default = 5,
              help = 'Minimal fold of total kmer : total reads',metavar = '5'),
  make_option(c("-n","--min_sample"),type = "integer",default = 3,
              help = 'At least appear in n samples',metavar = '3'),
  make_option(c("-C","--correct"),type = "logical",default = TRUE,
              help = 'Whether to correct for the p value when compute correlations',metavar = 'T'),
  make_option(c("-o","--output_path"),type = 'character',default = NULL,
              help = 'Path to output files',metavar = "/OUTPUT/PATH"),
  make_option(c("-P","--p_threshold"),type = "double",default = 0.05,
              help = "p/FDR value threshold",metavar = "0.05")
)
opt_parser <- OptionParser(usage = "decotaminate.R -s /SCRIPTS/PATH -i /INPUT/KRAKENUNIQ/REPORT/PATH -S PRJNA00000 -u 1000 -f 5 -n 3 -C T -o /OUTPUT/PATH -P 0.05",
                           option_list = option_list)
opt <- parse_args(opt_parser)

script_path <- opt$script_path
#script_path <- "/Users/yinwen/myfiles/Projects/2NAFLD/0.2Decontamination"
kreport_path <- opt$input_kreport_path
#kreport_path <- "/Users/yinwen/myfiles/Projects/2NAFLD/0.1KrakenUniqOutput/117PRJNA682622"
study <- opt$study_name
#study <- 'PRJNA682622'
unique_kmer_threshold <- opt$unique_kmer_threshold
#unique_kmer_threshold <- 1000
fold_threshold <- opt$fold_threshold
#fold_threshold <- 5
min_sample <- opt$min_sample
#min_sample <- 3
is_correct <- opt$correct
#is_correct <- T
output_path <- opt$output_path
p_threshold <- opt$p_threshold
#p_threshold <- 0.05
cell_line <- script_path


#output_path <- "/Users/yinwen/myfiles/Projects/2NAFLD/0.2Decontamination/Output"
output_raw_kreport <- paste0(output_path,"/",study,"/01raw_kreport")
output_unique_kmer <- paste0(output_path,"/",study,"/02unique_kmer_filted")
output_total_kmer <- paste0(output_path,"/",study,"/03total_kmer_filted")
output_study_level <- paste0(output_path,"/",study,"/04study_level_filted")
output_correlation <- paste0(output_path,"/",study,"/05correlation_filted")
output_correlation_plot <- paste0(output_path,"/",study,"/06correlation_plot")
output_cell_line <- paste0(output_path,"/",study,"/07cell_line_plot")
output_cell_line_quantative <- paste0(output_path,"/",study,"/08cell_line_quantative")
output_sample_result <- paste0(output_path,"/",study,"/09sample_decontaminated_results")


if (!dir.exists(output_raw_kreport)) {
  dir.create(output_raw_kreport,recursive = T)
}
if (!dir.exists(output_unique_kmer)) {
  dir.create(output_unique_kmer,recursive = T)
}
if (!dir.exists(output_total_kmer)) {
  dir.create(output_total_kmer,recursive = T)
}
if (!dir.exists(output_study_level)) {
  dir.create(output_study_level,recursive = T)
}
if (!dir.exists(output_correlation)) {
  dir.create(output_correlation,recursive = T)
}
if (!dir.exists(output_correlation_plot)) {
  dir.create(output_correlation_plot,recursive = T)
}
if (!dir.exists(output_cell_line)) {
  dir.create(output_cell_line,recursive = T)
}
if (!dir.exists(output_cell_line_quantative)) {
  dir.create(output_cell_line_quantative,recursive = T)
}
if (!dir.exists(output_sample_result)) {
  dir.create(output_sample_result,recursive = T)
}


suppressPackageStartupMessages({
  library(data.table)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(dplyr)
  library(tools)
})
############
source(paste0(script_path,"/","decotaminate_function.R"))

##step0:读入#######
single_kreport <- read_single_reports(kreport_path,study_name = study)
sample_name_list <- names(single_kreport)
for (i in 1:length(single_kreport)){
 # i = 1
  sample_name <- sample_name_list[[i]]
  write.table(single_kreport[[i]],paste0(output_raw_kreport,"/",sample_name,".txt"),quote = F,row.names = F,sep = '\t')
}

##step1:unique kmer过滤########
unique_kmer_filted_result <- lapply(single_kreport,function(x) unique_kmer_filter(x,unique_kmer_threshold))
###统计结果，合并所有样本
unique_kmer_filted_result_combined <- rbindlist(unique_kmer_filted_result)
write.table(unique_kmer_filted_result_combined,paste0(output_unique_kmer,'/unique_kmer_results.txt'),quote = F,row.names = F,sep = '\t')
###筛选出非污染
unique_kmer_filted_result_combined <- unique_kmer_filted_result_combined[unique_kmer_filted_result_combined$unique_contamination == 'FALSE',]
write.table(unique_kmer_filted_result_combined,paste0(output_unique_kmer,'/unique_kmer_filted.txt'),quote = F,row.names = F,sep = '\t')
###输出提示
num1 <- length(unique(unique_kmer_filted_result_combined$taxid)) #211
print(paste0("Step1: ",num1," taxonomy passed the unique kmer filtering!"))
#unique kmer的数量大于1000的物种还有211个

##step2: total kmer过滤######
total_kmer_filted_result <- total_kmer_filter(unique_kmer_filted_result_combined,fold_threshold)
write.table(total_kmer_filted_result,paste0(output_total_kmer,'/total_kmer_results.txt'),quote = F,row.names = F,sep = '\t')

##筛选并统计非污染结果
total_kmer_filted_result  <- total_kmer_filted_result[total_kmer_filted_result$total_contamination ==F,]
write.table(total_kmer_filted_result,paste0(output_total_kmer,'/total_kmer_filted.txt'),quote = F,row.names = F,sep = '\t')
num2 <- length(unique(total_kmer_filted_result$taxid)) 
print(paste0("Step2: ",num2," taxonomy passed the total kmer filtering!"))

##step3:study level过滤########
study_level_filted_result <- sample_filter(total_kmer_filted_result,min_sample)
write.table(study_level_filted_result,paste0(output_study_level,'/sample_level_results.txt'),quote = F,row.names = F,sep = '\t')

##筛选并统计非污染结果
study_level_filted_result <- study_level_filted_result[study_level_filted_result$sample_contamination == F,]
write.table(study_level_filted_result,paste0(output_study_level,'/sample_level_filted.txt'),quote = F,row.names = F,sep = '\t')

num3 <- length(unique(study_level_filted_result$taxid)) #132
print(paste0("Step3: ",num3," taxonomy passed the sample level filtering!"))


##step4.1: correlation 过滤########
correlation_filted_results <- correlation_filter(study_level_filted_result,correct = is_correct,threshold = p_threshold)
write.table(correlation_filted_results[[1]],paste0(output_correlation,"/correlation_summary.txt"),quote = F,row.names = F,sep = '\t')
write.table(correlation_filted_results[[2]],paste0(output_correlation,"/correlation_filted.txt"),quote = F,row.names = F,sep = '\t')
write.table(correlation_filted_results[[3]],paste0(output_correlation,"/correlation_contamination.txt"),quote = F,row.names = F,sep = '\t')

num4 <- nrow(correlation_filted_results[[2]])
print(paste0("Step4: ",num4," taxonomy passed the correlation spearman filtering!"))
##step4.3:correlation画图########
output_correlation_plot <- paste0(output_correlation_plot,"/correlation plot.pdf")
correlation_plot <- correlation_plot(correlation_filted_results$summary,study_name =study)
##保存
if (!is.null(output_correlation_plot)) {
  ggsave(correlation_plot, filename = output_correlation_plot,width = 8,height = 6)
}

###step5.1: cell line过滤 定性画图#########
cell_line_path <- paste0(cell_line,"/cell.lines.txt")
cell.line_file = read.delim(cell_line_path, header = T,sep = ' ') %>% tibble()
cell.line_file = cell.line_file[,] %>% mutate(study = 'cell lines')
colnames(cell.line_file)[6]  <- "totalkmer"
species_plot_df <- single_kreport %>% rbindlist() %>% rbind(cell.line_file) %>% 
  filter(rank %in% c('S',"species")) %>% subset(name %in% correlation_filted_results$filtered$name)
species_number <- length(unique(species_plot_df$name))
if (species_number < 50){
  cell_line_species_plot <- cell_line_filter(species_plot_df,study)
  output_cell_line_plot <- paste0(output_cell_line,"/cell line plot_species.pdf")
  if (!is.null(output_cell_line_plot)) {
    ggsave(cell_line_species_plot , filename = output_cell_line_plot ,width = 10,height = (species_number %/% 5 + 1)*2)
  }
} else {
  print("too many taxonomies for cell line density plots drawing!")
}

genus_plot_df <- single_kreport %>% rbindlist() %>% rbind(cell.line_file) %>%
  filter(rank %in% c('G',"genus")) %>% subset(name %in% correlation_filted_results$filtered$name)
genus_number <- length(unique(species_plot_df$name))
if (genus_number < 50 ){
  cell_line_genus_plot <- cell_line_filter(genus_plot_df ,study)
  output_cell_line_plot <- paste0(output_cell_line,"/cell line plot_genus.pdf")
 
  if (!is.null(output_cell_line_plot)) {
    ggsave(cell_line_genus_plot , filename = output_cell_line_plot ,width = 10,height = (genus_number %/% 5 + 1)*2)
  }
} else{
  print("too many taxonomies for cell line density plots drawing!")
}

###step5.2: cell line过滤 定量分析#########
correlation_filed_df <- single_kreport %>% 
  rbindlist() %>% 
  filter(rank %in% c('G',"genus",'S',"species")) %>%
  subset(name %in% correlation_filted_results$filtered$name)

cell_line_quantitive_result <- cell_line_quantitive(correlation_filed_df,cell.line_file,qtile = 0.99)
cell_line_quantitive_result.rmna <- na.omit(cell_line_quantitive_result)
cell_line_contaminations <- cell_line_quantitive_result.rmna[cell_line_quantitive_result.rmna$cell_line_contaminate == "TRUE",]$name

num5 <-  length(cell_line_contaminations)
print(paste0("Step5: ",num5," taxonomy failed the cell line filtering!"))


retain_results_list <- correlation_filted_results$filtered[!correlation_filted_results$filtered$name %in% cell_line_contaminations,]
write.table(retain_results_list,paste0(output_cell_line_quantative,"/cell_line_quantative.txt"),quote = F,row.names = F,sep = '\t')

###step6:summary, remove and output########
Nsample <- length(sample_name_list)
pass_decontaminate_statistics  <- data.frame("raw_taxonomies" = rep(0,Nsample),
                                        "unique_kmer_passed" = rep(0,Nsample),
                                        "total_kmer_passed" = rep(0,Nsample),
                                        "study_level_passed" = rep(0,Nsample),
                                        "correlation_spearman_passed" = rep(0,Nsample),
                                        "cell_line_passed" = rep(0,Nsample))
rownames(pass_decontaminate_statistics) <- sample_name_list
unique_kmer_result_list <- list()
total_kmer_result_list <- list()
study_level_result_list <- list()
correlation_result_list <- list()
cell_line_result_list <- list()

raw_kreport_list <- read_raw_kreport(kreport_path)



for (i in 1:Nsample){
  #i = 1
  sample_name <- sample_name_list[[i]]
  unique_kmer_result_list[[sample_name]] <- raw_kreport_list[[sample_name]][raw_kreport_list[[sample_name]]$taxID %in% unique_kmer_filted_result_combined[which(unique_kmer_filted_result_combined$sample == sample_name),]$taxid,]
  total_kmer_result_list[[sample_name]] <- raw_kreport_list[[sample_name]][raw_kreport_list[[sample_name]]$taxID %in% total_kmer_filted_result[which(total_kmer_filted_result$sample == sample_name),]$taxid,]
  study_level_result_list[[sample_name]] <- raw_kreport_list[[sample_name]][raw_kreport_list[[sample_name]]$taxID %in% study_level_filted_result[which(study_level_filted_result$sample == sample_name),]$taxid,]
  correlation_result_list[[sample_name]] <- raw_kreport_list[[sample_name]][raw_kreport_list[[sample_name]]$taxID %in% study_level_filted_result[which(study_level_filted_result$sample == sample_name &
                                                    study_level_filted_result$name %in% correlation_filted_results$filtered$name),]$taxid,]
  cell_line_result_list[[sample_name]] <- raw_kreport_list[[sample_name]][raw_kreport_list[[sample_name]]$taxID %in% study_level_filted_result[which(study_level_filted_result$sample == sample_name &
                                                                            study_level_filted_result$name %in% retain_results_list$name),]$taxid,]
  #write.table(unique_kmer_result_list[[sample_name]],paste0(output_unique_kmer,"/",sample_name,'_unique_kmer_filted.txt'),quote = F,row.names = F,sep = '\t')
  #write.table(total_kmer_result_list[[sample_name]],paste0(output_total_kmer,"/",sample_name,'_total_kmer_filted.txt'),quote = F,row.names = F,sep = '\t')
  #write.table(study_level_result_list[[sample_name]],paste0(output_study_level,"/",sample_name,'_study_level_filted.txt'),quote = F,row.names = F,sep = '\t')
  #write.table(correlation_result_list[[sample_name]],paste0(output_correlation,"/",sample_name,'_correlation_spearman_filted.txt'),quote = F,row.names = F,sep = '\t')
  #write.table(cell_line_result_list[[sample_name]],paste0(output_cell_line_quantative,"/",sample_name,'_cell_line_filted.kreport'),quote = F,row.names = F,sep = '\t')
  
  pass_decontaminate_statistics[i,"raw_taxonomies"] <- nrow(single_kreport[[sample_name]])
  pass_decontaminate_statistics[i,"unique_kmer_passed"] <- length(which(unique_kmer_filted_result_combined$sample == sample_name))
  pass_decontaminate_statistics[i,"total_kmer_passed"] <- length(which(total_kmer_filted_result$sample == sample_name))
  pass_decontaminate_statistics[i,"study_level_passed"] <- length(which(study_level_filted_result$sample == sample_name))
  pass_decontaminate_statistics[i,"correlation_spearman_passed"] <- length(which(study_level_filted_result$sample == sample_name &
                                                                                   study_level_filted_result$name %in% correlation_filted_results$filtered$name))
  pass_decontaminate_statistics[i,"cell_line_passed"] <- length(which(study_level_filted_result$sample == sample_name &
                                                                        study_level_filted_result$name %in% retain_results_list$name))
}


###step7: check######
taxDB <- read.table(paste0(script_path,"/taxDB"),sep = "\t",fill = T,header = F,quote = "")
#remove the line which child and parent line are both 1 to avoid endless loop
hierarchy_df <- taxDB[taxDB$V3 != 'root',] 
#map nodes
colnames(hierarchy_df)[1:2] <- c("Child","Parent")
parent_mapping <- setNames(hierarchy_df$Parent, hierarchy_df$Child)


checked_file_list <- list()
for (i in 1:Nsample){
  #i = 12
  sample_name <- sample_name_list[[i]]
  check_file <- cell_line_result_list[[sample_name]]
  checked_file_list[[sample_name]] <- check_parent_node(check_file)
  pass_decontaminate_statistics$final_result[i] <- length(which(checked_file_list[[sample_name]]$pass_check == T))
  checked_file_list[[sample_name]] <- checked_file_list[[sample_name]][checked_file_list[[sample_name]]$pass_check == T,1:9]
  write.table(checked_file_list[[sample_name]],
              paste0(output_sample_result,"/",sample_name,'_filted.kreport'),quote = F,row.names = F,sep = '\t')
}

write.table(pass_decontaminate_statistics,paste0(output_path,"/",study,"/decontamination_statistic.txt"),sep = '\t',quote = F)
  

print("finish analysing!")


