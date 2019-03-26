Fig6_FUN <- function(data) {
  ###################################################################################
  #Name: Wetlandscape Hydrology - Figure 3 Region Validation
  
  #Date: 14 JAN 2019
  #Purpose: Create figure 3 to plot observed/simulated 
  
  ##################################################################################
  
  df <- data %>%
    filter(day == 0) %>%
    filter(scale %in% c("wetland")) %>%
    dplyr::select(-c('day', 'V1', 'scale'))
  
  df.wide <- spread(df, var, value)
  regional.pca <- prcomp(df.wide[,c(5, 6, 7, 8, 11, 12, 13, 14 )], center = TRUE,scale. = TRUE)
  
  #setwd("~/Wetland_Hydrologic_Capacitance_Model/data")
  #save(regional.pca, df.wide,file =  'PCA2.rdata')
  
  autoplot(regional.pca, data = df.wide, colour = 'region',
           loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3)
  
  

}