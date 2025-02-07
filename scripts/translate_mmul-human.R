library(biomaRt)
# ensembl <- useEnsembl(biomart = "genes")
# ensembl <- useEnsembl(biomart = "genes", dataset = "mmulatta_gene_ensembl", version=111)
# # get gene names
# trans_table <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description', "hsapiens_homolog_ensembl_gene",
# 									"hsapiens_homolog_associated_gene_name"),
# 					 # filters = "with_hsapiens_homolog", values = TRUE,
#                mart = ensembl)
# # View(trans_table)
# # create hs_homolog+ensm gene list
# trans_table$hs_homolog_ensm <- trans_table$hsapiens_homolog_associated_gene_name
# # if empty, fill with ensbl_gene_id
# trans_table$hs_homolog_ensm[trans_table$hs_homolog_ensm == ""] <- trans_table$ensembl_gene_id[trans_table$hs_homolog_ensm == ""]
# # View(trans_table)

# load trans_table as tsv
trans_table <- read.table("data/Macaca_mulatta_genome/trans.table.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# function

# geneList is the named vector of log2Fold
translate_to_human_geneList <- function(geneList, trans_table) {

	geneList_human <- trans_table %>%
		filter(ensembl_gene_id %in% names(geneList)) %>%
		dplyr::select(ensembl_gene_id, hsapiens_homolog_associated_gene_name) %>%
		mutate(log2FoldChange = geneList[ensembl_gene_id]) %>%
		dplyr::select(hsapiens_homolog_associated_gene_name, log2FoldChange) %>%
		mutate(log2FoldChange = as.numeric(log2FoldChange)) %>%
		drop_na()
	# sort
	gene_list_horth <- geneList_human$log2FoldChange
	names(gene_list_horth) <- geneList_human$hsapiens_homolog_associated_gene_name
	return(gene_list_horth)
}


translate_to_human_count_matrix <- function(count_matrix, trans_table ) {
  ensembl_names <- rownames(count_matrix)
  print(length(ensembl_names))
    human_names <- trans_table %>%
      filter(ensembl_gene_id %in% ensembl_names) %>%
      filter(hsapiens_homolog_associated_gene_name != "") %>% # filter out genes without human homolog name
      dplyr::select(ensembl_gene_id, hsapiens_homolog_associated_gene_name) %>%
		# arrrange by order in count matrix
		arrange(match(ensembl_gene_id, ensembl_names))
  print(dim(human_names))

  # print percentage lost # round to 2 decimal places
    print(paste("Percentage of genes lost in translation:", round(100 * (1 - nrow(human_names) / length(ensembl_names)), 2), "%"), )
	  # filter count matrix
    count_matrix_human <- count_matrix[human_names$ensembl_gene_id, ]
	# change to human gene names
	rownames(count_matrix_human) <- human_names$hsapiens_homolog_associated_gene_name
    return(count_matrix_human)
}


get_heatmap_res <- function(res_name, vsd_scaled, dds, p_threshold = 0.05, trans_table){
   # get significant genes
  res <- results(dds, name = res_name, tidy = TRUE)

  sig_genes <- res[res$padj < p_threshold,][["row"]]
  if (length(sig_genes) == 0) {
    text <- paste("No significant genes for", res_name, "at" , p_threshold)
    no_res_plot <- ggplot() +
            annotate("text", x = 0.5, y = 0.5, label = text, size = 10)
    return(no_res_plot)
  }
  # get vsd for significant genes
  vsd_scaled <- vsd_scaled[rownames(vsd_scaled) %in% sig_genes,, drop=FALSE] # avoid transforming to a vector

  # get metadata
  sample_metadata <- colData(dds) %>% as.data.frame()
  vsd_long <- vsd_scaled %>% as.data.frame() %>%
      mutate(Gene = rownames(.)) %>%
    gather(key = "Sample", value = "counts", -Gene) %>%
      left_join(sample_metadata, by = "Sample")
   # change to symbol gene names
	vsd_long <- vsd_long %>%
		left_join(trans_table, by = c("Gene" = "ensembl_gene_id"), multiple = "first")


  # order genes by LFC
  print(res_name)
  resLFC <- lfcShrink(dds, coef=res_name, type="apeglm", parallel = TRUE)
  order <- resLFC %>% as.data.frame() %>%
            mutate(Gene = rownames(.)) %>%
		left_join(trans_table, by = c("Gene" = "ensembl_gene_id"), multiple = "first") %>%
            filter(Gene %in% sig_genes) %>% # keep only significant genes
            arrange(log2FoldChange)
  #reorder factor hs_homolog_ensm
  vsd_long <- vsd_long %>%
      mutate(hs_homolog_ensm = factor(hs_homolog_ensm, levels=unique(order$hs_homolog_ensm)))


    # make heatmap
  heatmap_week <-  vsd_long %>%
    group_by(Animal, hs_homolog_ensm) %>%
    mutate(counts = counts - counts[Week_post_FHT == 0]) %>% # difference from week 0
    filter(Week_post_FHT != 0) %>%
          ggplot(aes(x=Animal, y=hs_homolog_ensm, fill=counts)) +
          facet_nested(~Week_post_FHT + Treatment, scales = "free") +
          geom_tile() +
    				scale_fill_gradient2(midpoint = 0, mid = "gray", low = "blue", high = "red",
							  breaks = c(2,-2), labels = c("Upregulated", "Downregulated"), limits = c(-2,2),
								oob=scales::squish) +
          # geom_text(aes(label = round(counts, 2)), size = 2, color="white") +
          theme(axis.text.x=element_text(angle=65, hjust=1), axis.ticks.y=element_blank()) +
            labs(
                x = "Animals",
                y = "Human ortholog genes", fill=expression(paste(Delta, "expression level")))


  return (heatmap_week + labs(title = res_name))
}


get_vln_plot <- function(dds, res_name, trans_table, p_threshold = 0.05, log2FC_threshold = 1, xlim= c(-20, 20)) {
  data <- results(dds, name = res_name, tidy = TRUE) %>%
          mutate(Gene = row,
                 log10p = -log10(padj),
                 significant = ifelse(padj < p_threshold, "<=0.05", ">0.05" ))
     # translate names
	data <- data %>%
		left_join(trans_table, by = c("Gene" = "ensembl_gene_id"), multiple = "first")

  plot <- data %>%
          ggplot(aes(x = log2FoldChange, y = log10p, color = as.factor(significant))) +
          geom_point() +
          geom_hline(yintercept = -log10(p_threshold), linetype = "dashed") +
          geom_vline(xintercept = log2FC_threshold, linetype = "dashed") +
          geom_vline(xintercept = -log2FC_threshold, linetype = "dashed") +
          geom_text_repel(data = data %>% filter(significant == "<=0.05") %>%
            mutate(abs_log2FC = abs(log2FoldChange)) %>% slice_max(abs_log2FC, n = 15),
                          aes(label = hs_homolog_ensm), size = 5) +
          scale_color_manual(values = c("<=0.05" = "red", ">0.05" = "black")) +
          theme_bw() +
          labs(title = res_name, color="p.adj" ) +
            xlim(xlim)
  return(plot)
}