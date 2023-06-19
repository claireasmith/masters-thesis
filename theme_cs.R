# Define custom theme for plots
# Claire Smith 17 Nov 2022
# Based off this tutorial: https://rpubs.com/mclaire19/ggplot2-custom-themes

theme_cs <- function(fontsize=12, font = "sans"){ 
  # font <- "serif"   #assign font family up front
  
  theme_classic() %+replace%    #replace elements we want to change
    
    theme(
      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = fontsize,                #set font size
        face = 'italic',            #bold typeface
        hjust = 0,                #left align
        vjust = 2),               #raise slightly
      
      plot.subtitle = element_text(          #subtitle
        family = font,            #font family
        size = fontsize),               #font size
      
      plot.caption = element_text(           #caption
        family = font,            #font family
        size = fontsize,                 #font size
        hjust = 1),               #right align
      
      #grid elements
      panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      axis.ticks = element_blank(),          #strip axis ticks
      
      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = fontsize),               #font size
      
      axis.text = element_text(              #axis text
        family = font,            #axis famuly
        size = fontsize),                #font size
      
      axis.text.x = element_text(            #margin for axis text
        margin=margin(5, b = 10)),
      
      legend.title = element_text(
        family = font,
        size = fontsize),
      
      legend.text = element_text(
        family = font,
        size = fontsize)
      
    )
}
