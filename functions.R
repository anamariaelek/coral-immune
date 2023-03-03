# Plotting

ggplotMA <- function(
  res,
  padj_thr = 0.05,
  sign_col = c("red", "blue"),
  lims_fc = c(NA, NA),
  lims_mean = c(NA, NA),
  trans_mean = "identity",
  title = ""
) {
  res_dt <- as.data.table(res)
  res_dt[log2FoldChange < lims_fc[1], log2FoldChange := lims_fc[1]]
  res_dt[log2FoldChange > lims_fc[2], log2FoldChange := lims_fc[2]]
  res_dt[, dir := ifelse(log2FoldChange > 0, "up", "down")]
  res_dt[is.na(padj), padj := 1]
  res_dt[, sign := padj < padj_thr]
  res_dt[sign == FALSE, dir := NA]
  gp <- ggplot(
    res_dt,
    aes(x = baseMean, y = log2FoldChange)
  )
  if (length(sign_col) == 1) {
    gp <- gp +
    geom_point(aes(colour = sign), size = 0.5) +
    scale_color_manual(
      values = c("FALSE" = "grey", "TRUE" = sign_col),
      name = sprintf("p.adjusted < %s", padj_thr)
    )
  } else if (length(sign_col) > 1) {
    gp <- gp +
    geom_point(aes(colour = dir), size = 0.5) +
    scale_color_manual(
      values = c("up" = sign_col[2], down = sign_col[1]),
      na.value = "grey",
      name = sprintf("p.adjusted < %s", padj_thr)
    )
  }
  gp <- gp +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    scale_y_continuous(
      limits = lims_fc,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_x_continuous(
      limits = lims_mean,
      expand = expansion(mult = c(0.01, 0.01)),
      trans = trans_mean
    ) +
    geom_hline(aes(yintercept = 0), size = 1) +
    labs(
      x = "mean of normalized counts",
      y = "log2 fold change",
      subtitle = sprintf(
        "up=%s; down=%s",
        nrow(res_dt[padj < padj_thr & log2FoldChange > 0]),
        nrow(res_dt[padj < padj_thr & log2FoldChange < 0])
      )
    ) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 20),
      plot.subtitle = element_text(size = 18)
    )
  if (title != "")
    gp <- gp + labs(title = title)
  gp
}
ggplotVolcano <- function(
  res,
  padj_thr = 0.05,
  sign_col = c("red", "blue"),
  lims_fc = c(NA, NA),
  lims_sig = c(NA, NA),
  title = ""
) {
  res_dt <- as.data.table(res)
  res_dt[log2FoldChange < lims_fc[1], log2FoldChange := lims_fc[1]]
  res_dt[log2FoldChange > lims_fc[2], log2FoldChange := lims_fc[2]]
  res_dt[, dir := ifelse(log2FoldChange > 0, "up", "down")]
  res_dt[, minuslog10padj := -1 * log10(padj)]
  res_dt[, sign := padj < padj_thr]
  res_dt[sign == FALSE, dir := NA]
  gp <- ggplot(
    res_dt,
    aes(x = log2FoldChange, y = minuslog10padj)
  )
  if (length(sign_col) == 1) {
    gp <- gp +
    geom_point(aes(colour = sign), size = 0.5) +
    scale_color_manual(
      values = c("FALSE" = "grey", "TRUE" = sign_col),
      name = sprintf("p.adjusted < %s", padj_thr)
    )
  } else if (length(sign_col) > 1) {
    gp <- gp +
    geom_point(aes(colour = dir), size = 0.5) +
    scale_color_manual(
      values = c("up" = sign_col[2], down = sign_col[1]),
      na.value = "grey",
      name = sprintf("p.adjusted < %s", padj_thr)
    )
  }
  gp <- gp +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    scale_x_continuous(
      limits = lims_fc,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
      limits = lims_sig,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      x = "log2 fold change",
      y = "- log10 adjusted p value",
      subtitle = sprintf(
        "up=%s; down=%s",
        nrow(res_dt[padj < padj_thr & log2FoldChange > 0]),
        nrow(res_dt[padj < padj_thr & log2FoldChange < 0])
      )
    ) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 20),
      plot.subtitle = element_text(size = 18)
    )
  if (title != "")
    gp <- gp + labs(title = title)
  gp
}

# GO analysis

require(topGO)

read_eggnog <- function (file){
  # col_names <- c(
  #   "query_name", "seed_eggNOG_ortholog", "seed_ortholog_evalue", 
  #   "seed_ortholog_score", "predicted_gene_name", "GO_terms", "KEGG_KOs", 
  #   "BiGG_reactions", "Annotation_tax_scope", "OGs", "bestOG_evalue_score", 
  #   "COG_cat", "eggNOG_annot"
  # )
  col_names <- c(
    "query_name","seed_eggNOG_ortholog","seed_ortholog_evalue","seed_ortholog_score",
    "best_tax_level","Preferred_name","GOs","EC",
    "KEGG_ko","KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass",
    "BRITE","KEGG_TC","CAZy","BiGG_Reaction"
  )
  fread(cmd=sprintf("grep -v '^#' %s",file), col.names=col_names, sep="\t", select=1:length(col_names))
}

#' @param list_interest character, gene names
#' @param gomap named list of GO annotations for genes, names of list should be gene names
#' @param output_name character, prefix for output file names
#' @param name_geneset character, used to construct file name for output files
#' @param ontology_set character(s) indicating GO ontology to use, `c("BP","CC","MF")` 
#' @param tg_test character, which test to use, one of `c("fisher","t")`, see `statistic` in `?topGO::runTest`
#' @param tg_algorithm character, which algorithm to use, see `algorithm` in `?topGO::runTest`
#' @param printfile logical, whether to save plot and table
#' @param p_adj_method character, multiple correction method to use
topgofun  <- function(
  list_interest, gomap, output_name, name_geneset, ontology_set, tg_test="fisher", tg_algorithm="classic", 
  topnum=20, nodesize=10, printfile=TRUE, p_adj_method="BH", firstSigNodes=10
) {
  
  library(topGO)
  
  # Input 
  list_interest = unique(list_interest)
  genom = names(gomap)
  gesel = factor(as.integer(genom %in% list_interest))
  names(gesel) = genom
  
  # shortened go mappings without empty transcripts
  gomap_nonempty = gomap[lapply(gomap,length)>0]
  
  namepref <- paste0(output_name,".",name_geneset,".topgo",".",tg_test,tg_algorithm)
  # if(printfile){
  #   pdf(file=paste0(namepref,".pdf"),height=4.5,width=4)
  # }
  par(mar=c(5,12,5,2))
  
  topgo_tau_tot = data.frame()
  
  if (length(list_interest[list_interest %in% names(gomap_nonempty)])>1) {
    
    for (ontology_seti in ontology_set) {
      # topGO setup 
      
      GOdata = new(
        "topGOdata", ontology=ontology_seti, allGenes=gesel,
        annot=annFUN.gene2GO, gene2GO=gomap
      )
      
      num_interest_feasible = sum(GOdata@feasible & genom %in% list_interest)
      
      # topGO analysis
      topgo_res = runTest(GOdata, algorithm = tg_algorithm, statistic = tg_test)
      topgo_tau = GenTable(
        GOdata, pval_test = topgo_res, orderBy = "pval_test", 
        topNodes = length(usedGO(object = GOdata))
      )
      topGO::printGraph(
        GOdata, result=topgo_res, firstSigNodes=firstSigNodes, # all the nodes in the graph: length(usedGO(object = GOdata)) -- a mess
        useInfo="all", fn.prefix=paste(namepref,ontology_seti,sep="."), pdfSW=TRUE
      )
      topgo_tau$pval_test = as.numeric(topgo_tau$pval_test)
      topgo_tau$pval_adj  = p.adjust(topgo_tau$pval_test, method=p_adj_method)
      topgo_tau$ontology = ontology_seti
      topgo_tau_tot = rbind(topgo_tau_tot,topgo_tau)
      
      # Output 
      # ploti=barplot(height = rev(head(log(topgo_tau$pval_test,10),topnum)),
      #               names.arg = rev(head(paste(topgo_tau$Term,topgo_tau$GO.ID),topnum)),
      #               xlim=c(0,-5),horiz=T,las=1,col="slategray3",border=NA,
      #               cex.names=0.35,cex.axis=0.6,cex.lab=0.6,cex.sub=0.6,cex.main=0.6,
      #               main=paste(name_geneset,"top GO:",ontology_seti,tg_test,tg_algorithm),
      #               sub =paste("n=",num_interest_feasible,"/",length(list_interest), sep=""),
      #               xlab="log(p)")
      # abline(v=log(0.01,10),lty=2,lwd=0.5,col="pink")
      # abline(v=log(0.05,10),lty=2,lwd=0.5,col="pink")
      # text(x=0,ploti,labels = paste("p =",signif(rev(head(topgo_tau$pval_test,topnum)),3)),
      #      col="red",pos=4,cex=0.35)
    }
    
  }else {
    print("skip, no annotations in interest list!")
  }
  
  if(printfile){
    write.table(
      topgo_tau_tot,
      file=paste(output_name,".",name_geneset,".topgo",".",tg_test,tg_algorithm,".txt",sep=""),
      sep="\t", quote=F, col.names=T, row.names=F, append = F)
    dev.off()
  }
  
  return(topgo_tau_tot)
}

# Revigo

plot_revigo <- function(
  df,
  legend_position = "bottom",
  legend_box = "vertical"
) {

    pdt <- copy(df)
    pdt[, plot_X := as.numeric(as.character(pdt$PC_0))]
    pdt[, plot_Y := as.numeric(as.character(pdt$PC_1))]
    pdt <- pdt[!is.na(plot_X) & !is.na(plot_Y)]
    pdt[, num_annots := as.numeric(as.character(pdt$LogSize))]
    pdt[, log10padj := as.numeric(as.character(pdt$Value)) * -1]
    pdt[, frequency := as.numeric(as.character(pdt$Frequency))]
    pdt[, uniqueness := as.numeric(as.character(pdt$Uniqueness))]
    pdt[, dispensability := as.numeric(as.character(pdt$Dispensability))]
    class(pdt) <- "data.frame"
    
    ex <- pdt[pdt$dispensability < 0.15, ]
    x_range <- max(pdt$plot_X) - min(pdt$plot_X)
    y_range <- max(pdt$plot_Y) - min(pdt$plot_Y)
    p_range <- max(pdt$log10padj)
    
    p1 <- ggplot(pdt) +
        geom_point(
            aes(plot_X, plot_Y, fill = log10padj, size = num_annots),
            shape = 21,
            alpha = 1
        ) +
        scale_fill_gradientn(
          name = "log10\nadjusted pvalue",
          colours = c(RColorBrewer::brewer.pal(4, "Blues")[-1], "#012d66", "#01153d"),
          limits = c(0, p_range)
        ) +
        scale_size(
          name = "number of\nannotations",
          range = c(2, 10)
        ) +
        geom_text_repel(
            data = ex,
            aes(plot_X, plot_Y, label = Name),
            colour = I(alpha("black", 0.85)),
            size = 4
        ) +
        labs(y = "MDS2", x = "MDS1") +
        coord_fixed() +
        xlim(
            min(pdt$plot_X) - x_range / 10,
            max(pdt$plot_X) + x_range / 10
        ) +
        ylim(
            min(pdt$plot_Y) - y_range / 10,
            max(pdt$plot_Y) + y_range / 10
        ) +
        theme(
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            legend.position = legend_position,
            legend.box = legend_box
        )

    p1
}
