#函数描述：
##读入单个krakenuniq report
#参数：
##filepath - 储存krakenuniq report的文件路径
##study_name - study名字
#返回值：
##输出成一个列表，每个元素以样本名为名字
##列表中的表格列组成：study_name - study名称; sample - sample名称;rank - 物种层级; taxid - taxid；name - 去除空格后的物种名称; 
##reads - 注释到该物种的reads总数; total kmer - 总kmer数，计算：unique kmer * dup; uniq - unique kmer; rpm - reads/total reads(root+uncultrued); 
##rpmm - reads / total microbiome reads(bacteria + fungi + Viruses)  
#注意事项：
##需要输入为标准的krakennuniq report的格式，即前两行为注释，列名为标准命名规范
##cell line过滤用到的微生物相对丰度计算：某一行的reads数除以（细菌+真菌+病毒）的总reads数
read_single_reports = function(filepath = inputKreport,
                               study_name){
  require(data.table)
  require(tibble)
  require(tools)
  files <- list.files(filepath,full.names = T)
  kreportlist <- list()
  sample_names = file_path_sans_ext(basename(files))
  for (i in 1:length(files)){
    # i = 1
    study_name_list = rep(study_name, length(files))
    sample_name <- file_path_sans_ext(basename(files[[i]]))
    kreport = read.delim(files[[i]], header = T,skip = 2)
    kreport$taxName = trimws(kreport$taxName)
    total_reads = kreport$reads[1] + kreport$reads[2]
    n_microbiome_reads = sum(kreport$reads[kreport$taxName %in% c('Bacteria', 'Fungi', 'Viruses')])
    kreportlist[[sample_name]] = data.frame(study = study_name_list[i], sample = sample_names[i],
                                            rank = kreport$rank, taxid = kreport$taxID, 
                                            name = kreport$taxName, reads = kreport$reads,
                                            totalkmer = kreport$kmers * kreport$dup, uniq = kreport$kmers,
                                            rpm = kreport$reads/total_reads*10^6,
                                            rpmm = kreport$reads/n_microbiome_reads*10^6)
  }
  return(kreportlist)
}

#函数描述：
##第一步过滤：unique kmer至少大于1000
#参数：
##kraken_report - 原始krakenuniq的列表
##threshold - unique kmer的阈值设置，默认为1000
#返回值：
##在原列表的基础上在表格中加上一列代表是否为这一步筛选得到的污染，即如果unique kmer数小于阈值，则认为是污染
#注意事项：
##没有过滤，仅标记
unique_kmer_filter = function(kraken_report,threshold = 1000){
  unique_kmer_filter_result <-  kraken_report %>%
    mutate(unique_contamination = ifelse(uniq < threshold, T, F))
  return(unique_kmer_filter_result)
}

#函数描述：
##第二步过滤:total kmer至少大于n倍的总reads数
#参数：
##unique_kmer_filted_result - 输入合并后的长文件
##threshold -阈值
#返回值：
##在原列表的基础上在表格中加上一列代表是否为这一步筛选得到的污染，即如果total kmer数小于阈值倍的总reads数，则认为是污染
total_kmer_filter = function(unique_kmer_filted_result,threshold){
  total_kmer_filted_result <- unique_kmer_filted_result %>%
    mutate(total_contamination = ifelse(totalkmer < threshold * reads, T, F))
  return(total_kmer_filted_result)
}

#函数描述：
##第三步过滤：物种至少出现在3个样本里
#参数：
##total_kmer_filted_result - 输入合并后的长文件
##min_sample-阈值
#返回值：
##在原列表的基础上在表格中加上一列代表是否为这一步筛选得到的污染，即如果只在1个或2个样本中出现，则认为是污染
sample_filter = function(total_kmer_filted_result,min_sample =3){
  study_level_filted_result <- total_kmer_filted_result  %>%
    group_by(taxid) %>%
    mutate(n()) %>%
    mutate(sample_contamination = ifelse(n() < min_sample, T, F))
  return(study_level_filted_result)
}

#函数描述：
##第四步过滤：相关性显著
#参数：
##study_level_filted_result - 输入合并后的长文件
##correct -是否矫正p值
##threshold -
#返回值：
##列表，第一个为总结，第二个为筛选后得到非污染物种，第三个为污染物种
correlation_filter = function(study_level_filted_result, # input file
                              col1 = 'totalkmer',col2 = 'uniq', col3 = 'reads', 
                              correct, # whether correct for p value
                              threshold){
  correlation_results <- list()
  #calculate correlation
  correlation_summary <- study_level_filted_result %>%
    group_by(rank) %>%
    group_by(name) %>%
    summarize(r1 = cor(get(col1), get(col2), method = 'spearman'),
              r2 = cor(get(col1), get(col3), method = 'spearman'),
              r3 = cor(get(col2), get(col3), method = 'spearman'),
              p1 = cor.test(get(col1), get(col2), method = 'spearman')$p.value,
              p2 = cor.test(get(col1), get(col3), method = 'spearman')$p.value,
              p3 = cor.test(get(col2), get(col3), method = 'spearman')$p.value
    )
  # correction for p value
  if (correct) {
    correlation_summary <- correlation_summary %>%
      mutate(FDR1 = p.adjust(p1, "BH"),
             FDR2 = p.adjust(p2, "BH"),
             FDR3 = p.adjust(p3, "BH"))
  }
  # filt
  filted_results <- correlation_summary %>%
    filter(if (correct) FDR1 < threshold & FDR2 < threshold & FDR3 < threshold 
           else p1 < threshold & p2 < threshold & p3 < threshold)
  correlation_results$summary <- correlation_summary
  correlation_results$filtered <- filted_results
  correlation_results$contamination <- correlation_results$summary[!correlation_results$summary[[1]] %in% correlation_results$filtered[[1]],]
  return(correlation_results)
}


#函数描述：
##相关性绘图
#参数：
##correlation_filter_result - 相关性列表结果中的总结数据框
#返回值：
##图
#correlation_filted_results$summary
#study_name =study
#output_path = output_correlation_plot


correlation_plot <- function(correlation_filter_result =correlation_filted_results$summary,
                             col1 = 'FDR1',col2 = 'FDR2',col3 = 'FDR3',
                             study_name =""){
  #loadfonts(device = "pdf")
  correlation_filter_result = correlation_filted_results$summary
  correlation_filter_result$color_group <- ifelse(correlation_filter_result[[col3]] < 0.05,
                                                  paste0(col3,"<0.05"), paste0(col3,">0.05"))
  #correlation_filter_result <- correlation_filter_result[1:50,]
  plot <- ggplot(correlation_filter_result, aes_string(-log2(as.numeric((correlation_filter_result[[col1]]))+0.00001), 
                                                       -log2(as.numeric((correlation_filter_result[[col2]]))+0.00001),
                                                       color = "color_group")) +
    geom_density_2d(size = 0.75, color = "blue") + 
    geom_point(size = 1.5,alpha = 0.7) +
    geom_vline(xintercept = 2.995732, color = "red", linetype = "dashed") +
    geom_hline(yintercept = 2.995732, color = "red", linetype = "dashed") +
    annotate("text", x = 10, y = 3, label = "FDR2 = 0.05", vjust = -1.5, color = "red") +
    annotate("text", x = 2, y = 12, label = "FDR1 = 0.05", hjust = -0.5, color = "red") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
    scale_color_manual(values = setNames(c("red", "grey"), 
                                         c(paste0(col3, "<0.05"), paste0(col3, ">0.05"))))+
    xlab(paste0("Spearman correlation's ", col1, "\n between total kmer vs unique kmer")) +
    ylab(paste0("Spearman correlation's ", col2, "\n between total kmer vs reads")) +
    labs(color = paste0("Spearman correlation's ",col3, "\n between unique kmer vs reads"))+
    ggtitle(study_name)
  return(plot)
}


#函数描述：
##cell line过滤
#参数：
##correlation_filed_cell_line_combined_df -将已过滤的和cell line的表格结合在一起的表格
cell_line_filter = function(correlation_filed_cell_line_combined_df,study_name){
  #combined_df <- single_kreport %>% rbindlist() %>% rbind(cell.line_file) %>% filter(rank %in% c('S',"species")) 
  #pre_filted_df <- 
  tax_number <- length(unique(correlation_filed_cell_line_combined_df$name))
  if (tax_number < 6){
    n <- 1
  } else{
    n <- tax_number %/% 5 + 1
  }
  color_values <- setNames(c('grey50', 'navyblue'), c('cell lines', study_name))
  plot <- ggplot(correlation_filed_cell_line_combined_df, aes(rpmm, fill = study, ..scaled..)) + 
    geom_density(alpha = 0.5, color = NA) +
    facet_wrap(~name, scales='free',nrow = n) + 
    scale_fill_manual(values = color_values) +
    theme_minimal() + 
    xlab('Microbiome reads per million') +
    ylab('Density') +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x)),
                  breaks = trans_breaks("log10", function(x) 10^x, n=4),
                  oob = scales::squish, expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(fill = NA, color = 'black'),
          strip.background = element_blank(),
          axis.ticks.x = element_line(size=0),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank(), 
          strip.text = element_text(color = 'black', size = 10),
          axis.text.x = element_text(color = 'black', size = 10),
          axis.text.y = element_text(color = 'black', size = 10),
          axis.title.y = element_text(color = 'black', size = 10),
          axis.title.x = element_text(color = 'black', size = 10),
          plot.margin = unit(c(0, 0.1, 0, 0), "cm"),
          legend.key.size = unit(0.2, "cm"),
          legend.text = element_text(color = 'black', size = 10),
          legend.position = 'bottom')
  return(plot)
}

#函数描述：
##cell line 定量分析
#参数：
##correlation_filed_df - 已经过滤的结合的表格
##cell.line_file - cell line的表格
#返回值：
##所有属和种哪些是污染：在样本里的99百分位数高于cell line里的99百分位数
cell_line_quantitive = function(correlation_filed_df,cell.line_file,qtile =0.99){
  cell_line_qtile <- cell.line_file %>%
    group_by(name, rank) %>% 
    summarize(CLrpmm = 10^quantile(log10(rpmm), qtile, na.rm = T), .groups = 'keep')
  correlation_filed_df_qtile <- correlation_filed_df %>%
    group_by(name, rank) %>% 
    summarize(Studyrpmm = 10^quantile(log10(rpmm), qtile, na.rm = T), .groups = 'keep')
    
  cell_line_quantitive_result <- correlation_filed_df_qtile %>% left_join(cell_line_qtile, by = 'name') %>% 
    mutate(cell_line_contaminate = Studyrpmm < CLrpmm )
  return(cell_line_quantitive_result)
}

#函数描述：
##读入原始的kreport
#参数：
##kreport_path - 路径
#返回值：
##原始的kreport的列表

read_raw_kreport = function(kreport_path){
  files <- list.files(path = kreport_path,full.names = T)
  raw_input_list <- list()
  sample_names = file_path_sans_ext(basename(files))
  for (i in 1:length(files)){
    sample_name <- file_path_sans_ext(basename(files[[i]]))
    raw_input_list[[sample_name]] = read.delim(files[[i]], header = T,skip = 2)
    colnames(raw_input_list[[sample_name]])[1] <- '%'
  }
  return(raw_input_list)
}

#函数描述：
##检查删除完之后剩余物种的所有父节点是否存在，如果不存在在后续bracken定量的时候会报错。
#参数：
##cell_line_result - 09文件夹里的结果文件，以样本为单位的。
#返回值：
##新增列，如果是F则代表有问题，如果是T则是无问题
check_parent_node = function(check_file,parent_mapping){
  for (n in nrow(check_file):1){
    #n = nrow(check_file)
    child_id <- check_file$taxID[[n]]
    parent_id <- parent_mapping[as.character(child_id)]
    while (child_id != 1){
     # print(paste0("child_id is ",child_id))
     # print(parent_id)
      parent_row <- which(check_file$taxID == parent_id)
     # print(parent_row)
      if (length(parent_row) == 0){
        check_file$pass_check[n] <- FALSE
        break
      } else {
        check_file$pass_check[n] <- TRUE
        child_id <- parent_id
        parent_id <- parent_mapping[as.character(child_id)]
      }
    }
  }
  check_file$pass_check[1:2] <- TRUE
  return(check_file)
  }


