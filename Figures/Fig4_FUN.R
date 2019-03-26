Fig4_FUN <- function(region = 'delmarva') {
###################################################################################
#Name: Wetlandscape Hydrology - Figure 4 Seasonality in depth and fluxes
#Date: 14 JAN 2019
#Purpose: Create figure 4 to plot seasonal behaviour
##################################################################################

# Step 1: Setup Worskspace ---------------------------------------------------------
# 1. Load data ---------------------------------------------------------------------
rawdata <- data.table::fread(paste('/nfs/WHC-data/Figure Generation/', region,'.csv', sep = ""))
  
rawdata$value <- as.numeric(as.character(rawdata$value))

levels <- subset(rawdata, scale == 'wetland' &        # need to fix typo!
                          (var == 'y_w' | var == 'gw_local' |
                          var == 'spill_out'  | var == 'runoff_in') 
                          )   

stat <- levels %>%
  group_by(day, var) %>%
  summarise(x25th  = quantile(value, probs = 0.25, na.rm = T),
            median = quantile(value, probs = 0.50, na.rm = T),
            x75th  = quantile(value, probs = 0.75, na.rm = T))

stat2 <- melt(stat, by_id = 'day')
stat2$day <- as.numeric(as.character(stat$day))
stat$day <- as.numeric(as.character(stat$day))

# 2. Plot Data ---------------------------------------------------------------------
p1 <- ggplot() + 
  geom_line(data = stat2, aes(x = day, y = value, 
                              group   = variable, 
                              linetype= ifelse(variable=='median', 'dotted', 'solid')),
                              size    = 0.5) +
  geom_ribbon(data = stat, aes(x = day, ymin = x25th, ymax = x75th),
              fill = 'blue', alpha = 0.1) +
  theme_bw() +
  theme(plot.title    = element_text(size = 12, face = 'bold'),
        axis.text     = element_text(size = 8, face = 'bold', color = 'black'),
        axis.title.x  = element_text(size = 10, face = 'bold'),
        axis.title.y  = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text    = element_text(size = 8, face = 'bold')) +
  labs(x = 'Day of Year', title = paste('Figure 4a - ', toupper(region), sep = "" )) +
  facet_wrap(vars(var), 
             nrow = 2, 
             scales = 'free_y',
             strip.position = "left", 
             labeller = as_labeller(c(y_w = "Normalized Water Level", spill_out = "Spill Out (UNIT)",
                                      gw_local = "Groundwater Exchange (UNIT)", runoff_in = "Runoff In (UNIT)")))

setwd("/nfs/WHC-data/Figure Generation/Output")
ggsave(filename = paste('Fig4_',region, '_wetlandscale.jpg',sep = ""), plot = p1, 
       units = 'in', width = 6, height = 4, dpi = 500)



# =================================================================================
# Catchment Scale Values
# =================================================================================
# 1. Load data ---------------------------------------------------------------------
levels <- subset(rawdata, scale == 'catchment' &        # need to fix typo!
                   (var == 'y_w' | var == 'y_wt' |
                    var == 'spill_out'  | var == 'bf_out'))   

                   
stat <- levels %>%
  group_by(day, var) %>%
  summarise(x25th  = quantile(value, probs = 0.25, na.rm = T),
            median = quantile(value, probs = 0.50, na.rm = T),
            x75th  = quantile(value, probs = 0.75, na.rm = T))

stat2 <- melt(stat, by_id = 'day')
stat2$day <- as.numeric(as.character(stat$day))
stat$day <- as.numeric(as.character(stat$day))

# 2. Plot Data ---------------------------------------------------------------------
p2 <- ggplot() +
  geom_line(data = stat2, aes(x = day, y = value, group = variable,
                              linetype=ifelse(variable=='median', 'dotted', 'solid')),
            size = 0.1) +
  geom_ribbon(data = stat, aes(x = day, ymin = x25th, ymax = x75th),
              fill = 'blue', alpha = 0.1) +
  theme_bw() +
  theme(plot.title    = element_text(size = 12, face = 'bold'),
        axis.text     = element_text(size = 8, face = 'bold', color = 'black'),
        axis.title.x  = element_text(size = 10, face = 'bold'),
        axis.title.y  = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text    = element_text(size = 8, face = 'bold')) +
  labs(x = 'Day of Year', title = paste('Figure 4a - ', toupper(region), ' Catchment Scale', sep = "" )) +
  facet_wrap(vars(var), 
             nrow = 3, 
             scales = 'free_y',
             strip.position = "left", 
             labeller = as_labeller(c(y_w = "Normalized Wetland Stage", spill_out = "Spill Out (UNIT)",
                                      bf_out = "Baseflow Out (GW Out) (UNIT)", y_wt = "Normalized Water Table (UNIT)")))

setwd("/nfs/WHC-data/Figure Generation/Output")
ggsave(filename = paste('Fig4_',region , '_catchmentscale.jpg',sep = ""), plot = p2, 
       units = 'in', width = 6, height = 4, dpi = 500)




# # plot climate data
# load('/nfs/WHC-data/Regional_Analysis/PPR/inputs/climate.Rdata')
# a <- rep(1:365, times = 1000)
# b <- data.frame(a, pet.VAR, precip.VAR)

}



