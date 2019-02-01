gene_decode <- function (dataDir, image, measure,dbs,term,maxterms,prefix,...){

  ## Generic function to:
  ## >> download neurosynth surface maps
  ## >> extract the AIBS gene decoding
  ## >> do some thresholding
  ## >> run enrichment analysis on thresholded genes
  ## >> create some basic plots for the enrichment analysis
  
  ## By:
  ## Richard A.I. Bethlehem
  ## University of Cambridge
  ## ©rb643 2019
  
  ## EXAMPLE USAGE
  # dataDir <- "/Users/Richard/Dropbox/Research/Projects/Other/" #your top output directory
  # measure <- "Mean" #the measure you're using, just a place holder for creating a subdirectory and label
  # image <- "107923" #the image identifier on neurosynth
  # term <- "GO:" #how you wish to shorten the terms using the regexpr function in case the database terms are very long, leave as empty string if irrelevant
  # maxterms <- 5 #how many terms to display in your boxplot
  # prefix <- "GO" #string prefix for output files
  # dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018","GO_Biological_Process_2018")
  ## dbs can be queried using `dbs <- listEnrichrDbs()` #this is a list of databases that Enrichr searches
  
  ## gene_decode(dataDir, image, measure, dbs, term)
  # needs these libraries
  require(enrichR)
  require(ggplot2)
  require(dplyr)
  require(rjson)
  require(viridis)
  require(ggpubr)
  
  # do some folder setup
  folderOut <- paste(dataDir,measure,sep="/")
  folderOutTables <- paste(dataDir,measure,"tables",sep="/")
  
  # check if it exists and if not make it so
  if (!dir.exists(dataDir))(
    dir.create(dataDir)
  )
  
    if (!dir.exists(folderOut))(
    dir.create(folderOut)
  )
  
  if (!dir.exists(folderOutTables))(
    dir.create(folderOutTables)
  )
  
  # set the working directory
  setwd(dataDir)
  
  # load the json from Neurosynth
  data = rjson::fromJSON(file=paste("https://neurovault.org/images/",image,"/gene_expression/json?mask=cortex",sep = ""))

  # I don't like lists so convert to a usable dataframe (there's probably a better way to do this...)
  df <- data.frame(matrix(t(unlist(data$data)), nrow=length(data$data), byrow=T))
  colnames(df) <- c("symbol","page?","name","t","p","p_corr","var_explained","var_sd")
  
  # now make sure they have the correct format again
  df$t <- as.numeric(as.character(df$t))
  df$p <- as.numeric(as.character(df$p))
  df$p_corr <- as.numeric(as.character(df$p_corr))
  df$var_explained <- as.numeric(as.character(df$var_explained))
  df$var_sd <- as.numeric(as.character(df$var_sd))

  # split positive and negative and threshold
  genelist.pos <- df[ which( df$p < 0.05 & df$t >= 0) , ]
  genelist.neg <- df[ which( df$p < 0.05 & df$t <= 0) , ]
  genelist.pos.thres <- df[ which( df$p_corr < 0.05 & df$t >= 0) , ]
  genelist.neg.thres <- df[ which( df$p_corr < 0.05 & df$t <= 0) , ]
  
  # save the tables
  write.table(genelist.pos$symbol,file = paste(folderOut,"/genesample_pos_",measure,".txt",sep = ""),quote=FALSE,row.names = FALSE,col.names = FALSE, sep = "")
  write.table(genelist.neg$symbol,file = paste(folderOut,"/genesample_neg_",measure,".txt",sep = ""),quote=FALSE,row.names = FALSE,col.names = FALSE, sep = "")
  write.table(genelist.pos.thres$symbol,file = paste(folderOut,"/genesample_pos_thres_",measure,".txt",sep = ""),quote=FALSE,row.names = FALSE,col.names = FALSE, sep = "")
  write.table(genelist.neg.thres$symbol,file = paste(folderOut,"/genesample_neg_thres_",measure,".txt",sep = ""),quote=FALSE,row.names = FALSE,col.names = FALSE, sep = "")

  # also do the ontology analysis while we are at it but only if there are enough terms in the list 
  if ( (nrow(genelist.neg.thres)>5) & (nrow(genelist.pos.thres>5))) {
  
  # and run my dear run!!
  enriched.pos <- enrichr((as.character(genelist.pos.thres$symbol)), dbs)
  enriched.neg <- enrichr((as.character(genelist.neg.thres$symbol)), dbs)
    
    # convert them to manageable dataframes and get a GO term shorter ID, 
    # order by p-value write the full file and subset to top 10
    
    # Positive first
    plotlist.pos <- list()
    for (n in 1:length(unique(dbs))){ 
      i <- dbs[n]
      dfG <- enriched.pos[[i]]
      dfG$GO <- substring(as.character(dfG$Term), regexpr(term, as.character(dfG$Term)))
      dfG <- dfG[order(dfG$P.value),]
      write.csv(dfG, file = paste(folderOutTables,"/",prefix,i,"_pos.csv",sep = ""))
      
      # just a small check to make sure we are not introducing additional empty rows
      if (nrow(dfG)>maxterms){
      dfG <- dfG[1:maxterms,]
      } else {
      dfG <- dfG  
      }
      
      pl <- ggplot(dfG, aes(x=reorder(GO,sort(as.numeric(Z.score))), y=Z.score,fill = (Adjusted.P.value))) +
        geom_bar(stat = "identity")  +
        scale_fill_gradientn(colors = viridis_pal()(100), limits=c(0, 1), 
                             na.value = "#FFFFFF") +
        expand_limits(y=c(min(dfG$Z.score*1.5),max(dfG$Z.score*1.5))) + # just to create some space for labels
        geom_text(aes(label=round(Adjusted.P.value,3)),vjust=1) +
        xlab("Term") + ylab("Z-Score") + ggtitle(i) +
        theme_bw() +
        theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
              legend.position = "none") + 
        coord_flip()
      plotlist.pos[[n]] <- pl
    }
    
    posplot <- ggarrange(plotlist = plotlist.pos, ncol = 1, nrow = length(dbs))
    posplot <- annotate_figure(posplot, top = text_grob("Positive Genes",face = "bold", size = 14))
    
    # negative ones
    plotlist.neg <- list()
    for (n in 1:length(unique(dbs))){ 
      i <- dbs[n]
      dfG <- enriched.neg[[i]]
      dfG$GO <- substring(as.character(dfG$Term), regexpr(term, as.character(dfG$Term)))
      dfG <- dfG[order(dfG$P.value),]
      write.csv(dfG, file = paste(folderOutTables,"/",prefix,i,"_neg.csv",sep = ""))
      
      # just a small check to make sure we are not introducing additional empty rows
      if (nrow(dfG)>maxterms){
        dfG <- dfG[1:maxterms,]
      } else {
        dfG <- dfG  
      }
      
      pl <- ggplot(dfG, aes(x=reorder(GO,sort(as.numeric(Z.score))), y=Z.score,fill = (Adjusted.P.value))) +
        geom_bar(stat = "identity")  +
        scale_fill_gradientn(colors = viridis_pal()(100), limits=c(0, 1), 
                             na.value = "#FDE725FF") +
        expand_limits(y=c(min(dfG$Z.score*1.5),max(dfG$Z.score*1.5))) + # just to create some space for labels
        geom_text(aes(label=round(Adjusted.P.value,3)),vjust=1) +
        xlab("Term") + ylab("Z-Score") + ggtitle(i) +
        theme_bw() +
        theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
              legend.position = "none")+ 
        coord_flip()
      plotlist.neg[[n]] <- pl
    }
    
    negplot <- ggarrange(plotlist = plotlist.neg, ncol = 1, nrow = length(dbs))
    negplot <- annotate_figure(negplot, top = text_grob("Negative Genes",face = "bold", size = 14))
    
    figure <- ggarrange(posplot,negplot,ncol = 1, nrow = 2)
    figure <- annotate_figure(figure, top = text_grob("Gene Ontology",face = "bold", size = 16))
    
    pdf(file = paste(folderOut,"/",prefix,"_Overview.pdf",sep = ""), height = 9, width = 6)
      plot(figure)
    dev.off()
  }
  else {
    warning("Not enough genes to conduct gene ontology")
  }
    
}
