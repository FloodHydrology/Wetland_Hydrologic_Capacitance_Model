Fig6_FUN <- function(data) {
###################################################################################
#Name: Wetlandscape Hydrology - Figure 3 Region Validation

#Date: 14 JAN 2019
#Purpose: Create figure 3 to plot observed/simulated 

##################################################################################
df <- data %>%
  filter(day == 0) %>%
  filter(!var %in% c("GW_in")) %>%
  filter(scale %in% c("weltand")) 

%>%
  select(-X)

select(-one_of('day', 'X'))

df.wide <- spread(df, var, value)
regional.pca <- prcomp(df.wide[,c(6, 8, 12, 14)], center = TRUE,scale. = TRUE)

#setwd("~/Wetland_Hydrologic_Capacitance_Model/data")
#save(regional.pca, df.wide,file =  'PCA2.rdata')


#install.packages('ggfortify')
#library(ggfortify)

autoplot(regional.pca, data = df.wide, colour = 'region',
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3)

}