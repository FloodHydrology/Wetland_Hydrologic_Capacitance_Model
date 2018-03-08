# ======================================================================
# Function takes clay and sand percents and interpolates values 
# based on ternary plot from Laio et al, 2009

# REQUIRED INPUTS: function requires clay and sand percents and 
# OUTPUTS: Returns yc in m

# NOTES ================================================================
# x values -> % clay and x0 and x1 are bounds for % clay
# y values -> yc and y0 and y1 are the bounds for yc

get.yc <- function(clay, sand){
  x <- clay
  if (clay >= 65) {
    x0 <- 65; x1 <- 100; y0 <- 5; y1 <- 0.1
  }else if  (clay > 55) {
    x0 <- 55; x1 <- 65; y0 <- 30; y1 <- 5
  }else if  (clay > 43) {
    x0 <- 43; x1 <- 55; y0 <- 60; y1 <- 30
  }else if  (clay > 30) {
    x0 <- 30; x1 <- 43; y0 <- 80; y1 <- 60
  }else if  (sand >= 80) {
    x0 <- 80; x1 <- 100; y0 <- 80; y1 <- 60
    x <- sand
  }else if  (clay > 20) {
    x0 <- 20; x1 <- 30; y0 <- 100; y1 <- 80
  }else if  (sand > 50) {
    x0 <- 50; x1 <- 80; y0 <- 100; y1 <- 80
    x <- sand
  }else{  # force yc to be 100
    x0 <- 0; x1 <- 1; y0 <- 100; y1 <- 0
    x <- 0
  }
  # linear interpolation 
  Yc = y0 + (x - x0) * (y1- y0) / (x1 - x0)
  # and conversion to m from cm
  Yc = Yc / 100
  return(Yc)
}
