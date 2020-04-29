# setup -------------------------------------------------------------------

library(pool); library(tidyverse); library(ineq); library(geosphere);library(TapeR); library(zoo); library(vegan) 

source("pw.R")
source("functions.R")

data(SK.par.lme) # Taper model parameters for spruce in Germany based on BWI3 data obtained using

# data --------------------------------------------------------------------

plot.id <- tbl(KELuser, "plot") %>% pull(id)

data.list <- get_data(plot.id)

# parameters --------------------------------------------------------------

data.param <- calculate_parameters(data = data.list, 
                                   dataType = c("plot", "tree", "disturbance", "deadwood", "deadwood_tree", "core",
                                                "regeneration", "regeneration_subplot", "canopy", "temperature", "mortality"))

# collect everything ------------------------------------------------------

data.all <- collect_data(data.param)

# disconnection -----------------------------------------------------------

poolClose(KELuser)
