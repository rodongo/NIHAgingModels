

piePlot <- function(df, titleText){
  require(ggplot2)
  df$Percent <- rep(NA,nrow(df))
  for(i in 1:nrow(df)){
    df$Percent[i] <- df$Frequency[i]/sum(df$Frequency)*100
  }
  pie <- ggplot(df, aes(x="", y=Frequency, fill=Compartment)) + 
    geom_bar(stat="identity", width=1) + 
    coord_polar("y", start=0) + 
    geom_text(aes(label = paste0(Percent, "%")), position = position_stack(vjust = 0.5)) + 
    # scale_fill_brewer("RdYlBu") + 
    labs(x = NULL, y = NULL, fill = NULL, title = paste0(titleText)) + 
    theme_classic() + theme(axis.line = element_blank(),
                                      axis.text = element_blank(),
                                      axis.ticks = element_blank(),
                                      plot.title = element_text(hjust = 0.5, color = "#666666"))
  return(pie)
}


# Create a basic bar
comps$Percent <- rep(NA,nrow(comps))
for(i in 1:nrow(comps)){
  comps$Percent[i] <- comps$Frequency[i]/sum(comps$Frequency)*100
}
pie = ggplot(comps, aes(x="", y=Frequency, fill=Compartment)) + 
  geom_bar(stat="identity", width=1)

# Convert to pie (polar coordinates) and add labels
pie = pie + coord_polar("y", start=0) + 
  geom_text(aes(label = paste0(Percent, "%")), position = position_stack(vjust = 0.5))

# Add color scale (hex colors)
pie = pie + scale_fill_brewer("Blues") 

# Remove labels and add title
pie = pie + labs(x = NULL, y = NULL, fill = NULL, title = '')

# Tidy up the theme
pie = pie + theme_classic() + theme(axis.line = element_blank(),
                                    axis.text = element_blank(),
                                    axis.ticks = element_blank(),
                                    plot.title = element_text(hjust = 0.5, color = "#666666"))
pie
