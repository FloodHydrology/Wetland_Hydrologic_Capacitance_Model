Fig5_FUN <- function(data) {
  ###################################################################################
  #Name: Wetlandscape Hydrology - Figure 5 Regional Comparisons
  #Date: 14 JAN 2019
  #Purpose: Create figure 5 to compare regional characteristics
  ##################################################################################
  
  # # 1. Load data ---------------------------------------------------------------------
  # 
  # rawdata <- read.table('/nfs/WHC-data/Figure Generation/delmarva_output.csv', sep=",", header=TRUE)
  # rawdata$value <- as.numeric(as.character(rawdata$value))
  # rawdata$region <- 'DELMARVA'
  # 
  # rawdata2 <- read.table('/nfs/WHC-data/Figure Generation/ppr_output.csv', sep=",", header=TRUE)
  # rawdata2$value <- as.numeric(as.character(rawdata2$value))
  # rawdata2$region <- 'PPR'
  # 
  # rawdata3 <- read.table('/nfs/WHC-data/Figure Generation/ppr_output.csv', sep=",", header=TRUE)
  # rawdata3$value <- as.numeric(as.character(rawdata2$value))
  # rawdata3$region <- 'FLORIDA_PLACEHOLDER'
  # 
  # data <- do.call("rbind", list(rawdata, rawdata2, rawdata3))
  
  
  # Inter-regional comparison of wetland fluxes -------------------------------------------
  levels <- subset(data, scale == 'weltand' & day == 0)   # need to fix typo!
  levels <- subset(data, var == 'qf_in' | var == 'qf_out' | 
                     var == 'SW_in' | var == 'SW_out' |
                     var == 'GW_out')
  
  
  p1 <- ggplot(levels, aes( x = var, y = value, fill = region)) + 
    geom_boxplot(color = 'black', alpha = 0.7, outlier.shape = NA) +
    scale_y_continuous(limits = quantile(levels$value, c(0.1, 0.8)))+
    theme_bw() +
    theme(plot.title = element_text(size = 16, face = 'bold'),
          axis.text = element_text(size = 12, face = 'bold', color = 'black'),
          axis.title = element_text(size = 14, face = 'bold')) +
    labs(x = '', y = 'Annual Flux', title = 'Figure 5a - Wetland Scale') 
  
  setwd("/nfs/WHC-data/Figure Generation/Output")
  ggsave(filename = paste('Fig5_interregion_wetlandscale.jpg',sep = ""), plot = p1, 
         units = 'in', width = 6, height = 4, dpi = 500)
  
  # Inter-regional comparison of catchment fluxes ------------------------------------------
  levels <- subset(data, scale == 'catchment' & day == 0)   
  levels <- subset(data, var == 'sw_out' | var == 'gw_out' | 
                     var == 'sw_out_day' | var == 'gw_out_day')
  
  
 p2 <- ggplot(levels, aes( x = var, y = value, fill = region)) + 
    geom_boxplot(color = 'black', alpha = 0.7, outlier.shape = NA) +
    facet_wrap(~ var, ncol = 2, scales='free',
               strip.position = "top", 
               labeller = as_labeller(c(gw_out = "Groundwater Out", gw_out_day = "GW_OUT_DAY",
                                        sw_out = "Surface Water Out", sw_out_day = "SW_OUT_DAY"))) +
    theme_bw() +
    theme(plot.title = element_text(size = 16, face = 'bold'),
          axis.text = element_text(size = 12, face = 'bold', color = 'black'),
          axis.title.y = element_text(size = 14, face = 'bold'),
          axis.title.x  = element_blank(),
          axis.text.x  = element_blank(),
          legend.position = "bottom",
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text    = element_text(size = 12, face = 'bold')) +
    labs(x = '',y = '', title = 'Figure 5b - Catchment Scale') 
  
  setwd("/nfs/WHC-data/Figure Generation/Output")
  ggsave(filename = paste('Fig5b_interregion_catchmentscale.jpg',sep = ""), plot = p2, 
         units = 'in', width = 4.5, height = 4.5, dpi = 500)
  
  
}



