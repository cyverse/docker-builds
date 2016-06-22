bar_compare<-function(gene.name,n1,estimates.file,output.file=paste("Barplot_",gene.name,".pdf",sep=""),group.name){
    estimates<-read.delim(estimates.file,stringsAsFactors=FALSE,comment.char="#")
    estimates.gene<-estimates[estimates[,1]==gene.name,,drop=FALSE]
    
    if (nrow(estimates.gene)==0 | sum(is.na(estimates.gene[1,-(1:2)]))==(ncol(estimates.gene)-2)){
      try("No data for the input gene! \n")
    }
    else{
    estimates.gene1 <- estimates.gene[,-c(1)]
    usage.data<-t(estimates.gene[,-(1:2)])
    
    estimates.gene1.m <- melt(estimates.gene1, measure.vars = names(estimates.gene1)[-1])
    
    te1 <- c(replicate(n1,replicate(nrow(estimates.gene1), unique(rep(group.name[1],n1)))))
    te2 <- c(replicate(n1,replicate(nrow(estimates.gene1), unique(rep(group.name[2],n1)))))
    
    estimates.gene1.m$group <- c(te1,te2)
    
    main=paste("Estimated isoform usages for ",gene.name,"\n in ",sum(complete.cases(usage.data)[1:n1]),"+",sum(complete.cases(usage.data)[-(1:n1)])," samples", sep="")
    
    facetAdjust <- function(x, pos = c("up", "down"), newpage = is.null(vp), vp = NULL)
    {
      # part of print.ggplot
      ggplot2:::set_last_plot(x)
      if(newpage)
        grid.newpage()
      pos <- match.arg(pos)
      p <- ggplot_build(x)
      gtable <- ggplot_gtable(p)
      # finding dimensions
      dims <- apply(p$panel$layout[2:3], 2, max)
      nrow <- dims[1]
      ncol <- dims[2]
      # number of panels in the plot
      panels <- sum(grepl("panel", names(gtable$grobs)))
      space <- ncol * nrow
      # missing panels
      n <- space - panels
      # checking whether modifications are needed
      if(panels != space){
        # indices of panels to fix
        idx <- (space - ncol - n + 1):(space - ncol)
        # copying x-axis of the last existing panel to the chosen panels 
        # in the row above
        gtable$grobs[paste0("axis_b",idx)] <- list(gtable$grobs[[paste0("axis_b",panels)]])
        if(pos == "down"){
          # if pos == down then shifting labels down to the same level as 
          # the x-axis of last panel
          rows <- grep(paste0("axis_b\\-[", idx[1], "-", idx[n], "]"), 
                       gtable$layout$name)
          lastAxis <- grep(paste0("axis_b\\-", panels), gtable$layout$name)
          gtable$layout[rows, c("t","b")] <- gtable$layout[lastAxis, c("t")]
        }
      }
      # again part of print.ggplot, plotting adjusted version
      if(is.null(vp)){
        grid.draw(gtable)
      }
      else{
        if (is.character(vp)) 
          seekViewport(vp)
        else pushViewport(vp)
        grid.draw(gtable)
        upViewport()
      }
      invisible(p)
    }
    
    pdf(file=output.file,height=15,width=25)
    d <- ggplot(estimates.gene1.m, aes(variable, value, fill=group)) + geom_bar(stat = "identity") + facet_wrap(~isoform) + 
      theme_bw() + theme(axis.text = element_text(size = 12), axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5)) +
      labs(x="Isoform",y="Relative abundance") + ggtitle(main) + theme(text = element_text(size=20)) + theme(legend.position="top") + labs(fill="") +
      theme(strip.text.x = element_text(size = 12, colour = "black"))
    facetAdjust(d)
    dev.off()
    }
}

