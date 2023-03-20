# Plotting

ggplotMA <- function(
  res,
  pval_thr = NULL,
  padj_thr = 0.05,
  lfc_thr = 1,
  sign_col = c("red", "blue"),
  lims_fc = c(NA, NA),
  lims_mean = c(NA, NA),
  trans_mean = "identity",
  title = "",
  xlab = "mean of normalized counts",
  ylab = "log2 fold change",
  raster = TRUE
) {
  res_dt <- as.data.table(res)
  res_dt[log2FoldChange < lims_fc[1], log2FoldChange := lims_fc[1]]
  res_dt[log2FoldChange > lims_fc[2], log2FoldChange := lims_fc[2]]
  res_dt[, dir := ifelse(log2FoldChange > 0, "up", "down")]
  if (!is.null(padj_thr)) {
    res_dt[is.na(padj), padj := 1]
    res_dt[, signif := padj < padj_thr & abs(log2FoldChange) > lfc_thr]
    res_dt[signif == FALSE, dir := NA]
    legend_name <- sprintf("adjusted p value < %s", padj_thr)
    gp <- ggplot(
      res_dt,
      aes(x = baseMean, y = log2FoldChange)) +
      labs(
        x = xlab,
        y = ylab,
        subtitle = sprintf(
          "up=%s; down=%s",
          nrow(res_dt[padj < padj_thr & log2FoldChange > lfc_thr]),
          nrow(res_dt[padj < padj_thr & log2FoldChange < -1 * lfc_thr])
        )
      )
  } else if (!is.null(pval_thr)) {
    res_dt[is.na(pvalue), pvalue := 1]
    res_dt[, signif := pvalue < pval_thr & abs(log2FoldChange) > lfc_thr]
    res_dt[signif == FALSE, dir := NA]
    legend_name <- sprintf("p value < %s", pval_thr)
    gp <- ggplot(
      res_dt,
      aes(x = baseMean, y = log2FoldChange)
    ) +
      labs(
        x = xlab,
        y = ylab,
        subtitle = sprintf(
          "up=%s; down=%s",
          nrow(res_dt[pvalue < pval_thr & log2FoldChange > lfc_thr]),
          nrow(res_dt[pvalue < pval_thr & log2FoldChange < -1 * lfc_thr])
        )
      )
  } else {
    stop("One of padj_thr or pval_thr has to be supplied!")
  }
  if (length(sign_col) == 1) {
    if (raster == TRUE) {
      gp <- gp + ggrastr::geom_point_rast(aes(colour = signif), size = 0.5)
    } else {
      gp <- gp + geom_point(aes(colour = signif), size = 0.5)
    }
    gp <- gp + scale_color_manual(
      values = c("FALSE" = "grey", "TRUE" = sign_col),
      name = legend_name
    )
  } else if (length(sign_col) > 1) {
    if (raster == TRUE) {
      gp <- gp + ggrastr::geom_point_rast(aes(colour = dir), size = 0.5)
    } else {
      gp <- gp + geom_point(aes(colour = dir), size = 0.5)
    }
    gp <- gp + scale_color_manual(
      values = c("up" = sign_col[2], down = sign_col[1]),
      na.value = "grey",
      name = legend_name
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
  pval_thr = NULL,
  padj_thr = 0.05,
  lfc_thr = 1,
  sign_col = c("#b2182b", "#2166ac"),
  lims_fc = c(NA, NA),
  lims_sig = c(NA, NA),
  label_column = NULL,
  shape_column = NULL,
  label_fc_thr = 2,
  title = "",
  xlab = "log2 fold change",
  ylab = "- log10 adjusted p value",
  raster = TRUE
) {

  # parse data
  res_dt <- as.data.table(res)
  res_dt[log2FoldChange < lims_fc[1], log2FoldChange := lims_fc[1]]
  res_dt[log2FoldChange > lims_fc[2], log2FoldChange := lims_fc[2]]
  res_dt[, dir := ifelse(log2FoldChange > 0, "up", "down")]
  if (is.null(shape_column)) {
    res_dt[, shape_column := 1]
  }

  # using adjusted p value as threshold
  if (!is.null(padj_thr)) {
    res_dt[, sign := padj < padj_thr & abs(log2FoldChange) > lfc_thr]
    res_dt[, minuslog10padj := -1 * log10(padj)]
    res_dt[sign == FALSE, dir := NA]
    legend_name <- sprintf("adjusted p value < %s", padj_thr)
    gp <- ggplot(
      res_dt, aes_string(
        x = "log2FoldChange", y = "minuslog10padj",
        shape = shape_column
      )
    ) +
    labs(
      x = xlab,
      y = ylab,
      subtitle = sprintf(
        "up=%s; down=%s",
        nrow(res_dt[padj < padj_thr & log2FoldChange > lfc_thr]),
        nrow(res_dt[padj < padj_thr & log2FoldChange < -1 * lfc_thr])
      )
    )
  # using p value as threshold
  } else if (!is.null(pval_thr)) {
    res_dt[, sign := pvalue < pval_thr & abs(log2FoldChange) > lfc_thr]
    res_dt[, minuslog10pval := -1 * log10(pvalue)]
    res_dt[sign == FALSE, dir := NA]
    legend_name <- sprintf("p value < %s", pval_thr)
    gp <- ggplot(
      res_dt, aes_string(
        x = "log2FoldChange", y = "minuslog10pval",
        shape = shape_column
      )
    ) +
    labs(
      x = xlab,
      y = ylab,
      subtitle = sprintf(
        "up=%s; down=%s",
        nrow(res_dt[pvalue < pval_thr & log2FoldChange > lfc_thr]),
        nrow(res_dt[pvalue < pval_thr & log2FoldChange < -1 * lfc_thr])
      )
    )
  } else {
    stop("One of padj_thr or pval_thr has to be supplied!")
  }

  # using significance for colour
  if (length(sign_col) == 1) {
    if (raster == TRUE) {
      gp <- gp + ggrastr::geom_point_rast(aes(colour = signif), size = 0.5)
    } else {
      gp <- gp + geom_point(aes(colour = signif), size = 0.5)
    }
    gp <- gp + scale_color_manual(
      values = c("FALSE" = "grey", "TRUE" = sign_col),
      name = legend_name
    ) +
    scale_fill_manual(
      values = c("FALSE" = "grey", "TRUE" = sign_col),
      name = legend_name
    )
  # using significance and fold change direction for colour
  } else if (length(sign_col) > 1) {
    if (raster == TRUE) {
      gp <- gp + ggrastr::geom_point_rast(aes(colour = dir), size = 0.5)
    } else {
      gp <- gp + geom_point(aes(colour = dir), size = 0.5)
    }
    gp <- gp + scale_color_manual(
      values = c("up" = sign_col[2], down = sign_col[1]),
      na.value = "grey",
      name = legend_name
    ) + scale_fill_manual(
      values = c("up" = sign_col[2], down = sign_col[1]),
      na.value = "grey",
      name = legend_name
    )
  }

  # highlight genes by different shape
  if (!is.null(shape_column)) {
    shape_levels <- unique(res_dt[[I(shape_column)]])
    shape_levels <- structure(
      c(16, 22, 24, 23, 25)[seq_along(shape_levels)],
      names = shape_levels
    )
    gp <- gp +
      geom_point(
        data = res_dt[!I(shape_column) %in% c("none", "")][
          pvalue < pval_thr & abs(log2FoldChange) > lfc_thr],
        aes_string(shape = shape_column),
        size = 0.9, color = "grey"
      ) +
        scale_shape_manual(values = shape_levels) +
        guides(
          shape = guide_legend(override.aes = list(size = 4)),
          colour = guide_legend(override.aes = list(size = 4))
        )
  } else {
  gp <- gp +
    guides(
      shape = "none",
      colour = guide_legend(override.aes = list(size = 4))
    )
  }

  # plot limits
  gp <- gp +
    scale_x_continuous(
      limits = lims_fc,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
      limits = lims_sig,
      expand = expansion(mult = c(0.01, 0.01))
    )

  # labeling genes
  if (!is.null(label_column)) {
    if (!is.null(padj_thr)) {
      lab_dt_up <- res_dt[
        log2FoldChange > label_fc_thr & padj < padj_thr
      ]
      lab_dt_down <- res_dt[
        log2FoldChange < -1 * label_fc_thr & padj < padj_thr
      ]
    } else if (!is.null(pval_thr)) {
      lab_dt_up <- res_dt[
        log2FoldChange > label_fc_thr & pvalue < pval_thr
      ]
      lab_dt_down <- res_dt[
        log2FoldChange < -1 * label_fc_thr & pvalue < pval_thr
      ]
    }
    gp <- gp +
      geom_text_repel(
        mapping       = aes_string(label = label_column),
        data          = lab_dt_up,
        nudge_x       = 0.25 - lab_dt_up$log2FoldChange,
        segment.size  = 0.2,
        segment.color = "grey50",
        direction     = "y",
        hjust         = 0,
        max.overlaps  = Inf,
        min.segment.length = 0.5
      ) +
      geom_text_repel(
        mapping       = aes_string(label = label_column),
        data          = lab_dt_down,
        nudge_x       = -0.25 - lab_dt_down$log2FoldChange,
        segment.size  = 0.2,
        segment.color = "grey50",
        direction     = "y",
        hjust         = 1,
        max.overlaps  = Inf,
        min.segment.length = 0
      )
  }

  # theme elements
  gp <- gp +
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


# create a transcript to gene dictionary from a GTF annotation file
dictionary_t2g <- function(
  gtf_fn,
  vector_to_fix,
  t2g = TRUE,
  transcript_field = "transcript",
  transcript_id = "transcript_id",
  gene_id = "gene_id",
  return_elements_not_in_gtf = TRUE
) {

  # import gtf
  gene_txgtf <- rtracklayer::import(gtf_fn)

  if (t2g) {
    dic = as.data.frame(
      GenomicRanges::mcols(gene_txgtf[gene_txgtf$type == transcript_field])
    )[, gene_id]
    names(dic) <- as.data.frame(
      GenomicRanges::mcols(gene_txgtf[gene_txgtf$type == transcript_field])
    )[, transcript_id]
  } else {
    dic <- as.data.frame(
      GenomicRanges::mcols(gene_txgtf[gene_txgtf$type == transcript_field])
    )[, transcript_id]
    names(dic) <- as.data.frame(
      GenomicRanges::mcols(gene_txgtf[gene_txgtf$type == transcript_field])
    )[, gene_id]
  }

  # return object

  out <- dic [vector_to_fix]

  # return elements not in GTF dictionary, unaltered
  if (return_elements_not_in_gtf) {
    ixs_to_keep <- is.na(out)
    out[ixs_to_keep] <- vector_to_fix[ixs_to_keep]
    names(out[ixs_to_keep]) <- out[ixs_to_keep]
  }
  
  # return
  out

}