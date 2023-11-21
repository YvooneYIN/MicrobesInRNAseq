# MicrobesInRNAseq
Find microbial sequences in human RNA sequencing data
## Input:
1. krakenunique report
2. taxDB
3. cell.lines.txt
## Process: 
1. unique kmer filtering
2. total kmer filtering
3. study level filtering
4. correlation spearman filtering
5. cell line filtering
## Output: 
- /Output
  - /Study1
    - /decontamination_statistic.txt
    - /01raw_kreport
      - /sample.txt
        |study|sample|rank|taxid|name|reads|totalkmer|uniq|rpm|rpmm|
        |-----|------|----|-----|----|-----|---------|----|---|----|
        |PRJNA512027|SRR8378431|no rank|0|unclassified|141913|109152660.9|65754615|32290.8963486642|51717565.5976676|
    - /02unique_kmer_filted
      - /sample_unique_kmer_filted.txt
        |%|reads|taxReads|kmers|dup|cov|taxID|rank|taxName|
        |-|-----|--------|-----|---|---|-----|----|-------|
        |3.229|141913|141913|65754615|1.66|2453|0|no rank|unclassified|
    - /03total_kmer_filted
      - /sample_total_kmer_filted.txt
                    |study|sample|rank|taxid|name|reads|totalkmer|uniq|rpm|rpmm|unique_contamination|total_contamination|n()|sample_contamination|
      |-----|------|----|-----|----|-----|---------|----|---|----|--------------------|-------------------|---|--------------------|
      |PRJNA512027|SRR8378431|no rank|0|unclassified|141913|109152660.9|65754615|32290.8963486642|51717565.5976676|FALSE|FALSE|192|FALSE|
    - 04study_level_filted
      - sample_level_filted.txt
      - sample_level_results.txt
      - sample_study_level_filted.txt
    - 05correlation_filted
      - correlation_contamination.txt
      - correlation_filted.txt
      - correlation_summary.txt
      - sample_correlation_spearman_filted.txt 
     - 06correlation_plot
      - correlation plot.pdf
     - 07cell_line_plot
      - cell line plot_genus.pdf
      - cell line plot_species.pdf
     - 08cell_line_quantative
      - cell_line_quantative.txt
     - 09sample_decontaminated_results
      - sample_filted.kreport


     
