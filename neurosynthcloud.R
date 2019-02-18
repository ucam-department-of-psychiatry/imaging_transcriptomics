## Generic function to:
## Plot wordclouds for top and bottom terms of a neurosynth decoding

## By:
## Richard A.I. Bethlehem
## University of Cambridge
## ©rb643 2019

## EXAMPLE USAGE:
## neurosynthcloud(image = 107925,nterms = 20)

neurosynthcloud <- function(image, nterms, ...) {
  require(ggplot2)
  require(dplyr)
  require(rjson)
  require(viridis)
  require(ggpubr)
  require(ggwordcloud)
  
  ## load an image
  #image <- 107925
  data <-
    rjson::fromJSON(file = paste(
      "http://neurosynth.org/api/v2/decode/?neurovault=",
      image,
      sep = ""
    ))
  
  ## reformat the dataframe
  df <- as.data.frame(unlist(data$data$values))
  colnames(df) <- "freq"
  df$word <- rownames(df)
  df$word <- as.character(as.factor(df$word))
  df <- df[order(df$freq,decreasing = TRUE), ]
  
  ## get the top set of terms
  df_top <- df[1:nterms, ]
  df_bottom <- df[(nrow(df) - (nterms-1)):nrow(df), ]
  
  top <-
    ggplot(df_top, aes(
      label = word,
      size = freq,
      color = freq
    )) +
    geom_text_wordcloud_area(area_corr_power = 1,
                             rm_outside = FALSE,
                             shape = "triangle-upright") +
    scale_size_area(max_size = 10) +
    scale_color_gradient2(low = "red", high = "red4",guide = "colourbar",aesthetics = "colour") +
    theme_minimal()
  
  plotdefault1 <- data.frame(freq = seq(min(df_top$freq),max(df_top$freq), length.out=nterms),
                            y = as.factor(1))  
  colorbar1 <- ggplot(plotdefault1,aes(freq,y)) +
    geom_tile(aes(fill=freq)) + 
    scale_fill_gradient2(low = "red", high = "red4") +
    theme_minimal() +
    theme(legend.title = element_blank(),
          legend.position = "none",
          #axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank()) 
  
  bottom <-
    ggplot(df_bottom, aes(
      label = word,
      size = freq,
      color = freq
    )) +
    geom_text_wordcloud_area(area_corr_power = 1,
                             rm_outside = FALSE,
                             shape = "triangle-upright") +
    scale_size_area(max_size = 10) +
    scale_color_gradient2(low = "royalblue", high = "royalblue4",guide = "colourbar",aesthetics = "colour") +
    theme_minimal()
  
  plotdefault2 <- data.frame(freq = seq(min(df_bottom$freq),max(df_bottom$freq), length.out=nterms),
                             y = as.factor(1))  
  colorbar2 <- ggplot(plotdefault2,aes(freq,y)) +
    geom_tile(aes(fill=freq)) + 
    scale_fill_gradient2(low = "royalblue", high = "royalblue4") +
    theme_minimal() +
    theme(legend.title = element_blank(),
          legend.position = "none",
          axis.title.y=element_blank(),
          axis.text.y=element_blank()) 
  
  figure <-
    ggarrange(bottom, top, colorbar2, colorbar1,
              ncol = 2, nrow = 2,
              labels = c("Bottom Terms", "Top Terms","",""),
              widths = c(1,1), heights = c(4,1))
  
  pdf(file = "wordcloud.pdf")
  plot(figure)
  dev.off()
  
}
