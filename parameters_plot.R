pacman::p_load(tidyverse, pool, ineq, geosphere, TapeR, zoo, vegan)

# connection --------------------------------------------------------------


source('pw.R')


# functions ---------------------------------------------------------------


deg2rad <- function(deg) {(deg * pi) / 180}

interquantilerange <- function(x){quantile(x, 0.75, na.rm = T) - quantile(x, 0.25, na.rm = T)}

movingSum <- function(x, windowLength = 30){
  rollapply( x, 
             width = windowLength,
             FUN = sum,
             fill = NA,
             align = "center",
             na.rm = T,
             partial = TRUE)
}

volume <- function(SM, y){(SM * 0.001) / y}

# data --------------------------------------------------------------------


data(SK.par.lme) # Taper model parameters for spruce in Germany based on BWI3 data obtained using

plot <- tbl(KELuser, "plot")
tree <- tbl(KELuser, "tree") %>%  inner_join(., plot, by = c('plot_id' = 'id'))
deadwood <- tbl(KELuser, "deadwood") %>%  inner_join(., plot, by = c('plot_id' = 'id'))
core <- tbl(KELuser, "core") %>% inner_join(., tree, by = c('tree_id' = 'id'))
ring <- tbl(KELuser, "ring") 
regeneration <- tbl(KELuser, "regeneration") %>% inner_join(., plot, by = c("plot_id" = "id"))
regeneration_subplot <- tbl(KELuser, "regeneration_subplot") %>% inner_join(., plot, by = c("plot_id" = "id"))
dist_plot <- tbl(KELuser, "dist_plot") %>% inner_join(., plot %>% select(id, date, plotid), by = c("plot_id" = "id"))
dist_plot_event <- tbl(KELuser, "dist_plot_event") %>% inner_join(., dist_plot, by = c("dist_plot_id" = "id"))
canopy_analysis <- tbl(KELuser, "canopy_analysis") %>% inner_join(., plot, by = c("plot_id" = "id"))


# # tree parameters -------------------------------------------------------


tree %>%
  filter(!onplot %in% 0) %>%
  mutate(ba = pi * dbh_mm ^ 2 / 4 / 1000000) %>%
  collect() %>%
  group_by(date, plotid) %>%
  summarise(plotsize = first(plotsize),
            n_trees_live_500 = length(treeid[dbh_mm >= 500 & status %in% c(1:4)]),
            n_trees_live_700 = length(treeid[dbh_mm >= 700 & status %in% c(1:4)]),
            n_trees_dead_500 = length(treeid[dbh_mm >= 500 & status %in% c(11:23)]),
            n_trees_dead_700 = length(treeid[dbh_mm >= 700 & status %in% c(11:23)]),
            n_trees_dead_60 = length(treeid[dbh_mm >= 60 & status %in% c(11:23)]),
            dbh_quadrmean_dead_60 = sqrt(mean(dbh_mm[dbh_mm >= 60 & status %in% c(11:23)]^2, na.rm = T)),
            dbh_quadrmean_dead_60 = ifelse(dbh_quadrmean_dead_60 %in% "NaN", NA, dbh_quadrmean_dead_60),
            ba_dead_60 = sum(ba[dbh_mm >= 60 & status %in% c(11:23)], na.rm = T),
            n_trees_dead_100 = length(treeid[dbh_mm >= 100 & status %in% c(11:23)]),
            dbh_quadrmean_dead_100 = sqrt(mean(dbh_mm[dbh_mm >= 100 & status %in% c(11:23)]^2, na.rm = T)),
            dbh_quadrmean_dead_100 = ifelse(dbh_quadrmean_dead_100 %in% "NaN", NA, dbh_quadrmean_dead_100),
            ba_dead_100 = sum(ba[dbh_mm >= 100 & status %in% c(11:23)], na.rm = T),
            n_trees_live_60 = length(treeid[dbh_mm >= 60 & status %in% c(1:4)]),
            dbh_quadrmean_live_60 = sqrt(mean(dbh_mm[dbh_mm >= 60 & status %in% c(1:4)]^2, na.rm = T)),
            dbh_quadrmean_live_60 = ifelse(dbh_quadrmean_live_60 %in% "NaN", NA, dbh_quadrmean_live_60),
            dbh_mean_live_60 = mean(dbh_mm[dbh_mm >= 60 & status %in% c(1:4)], na.rm = T),
            dbh_mean_live_60 = ifelse(dbh_mean_live_60 %in% "NaN", NA, dbh_mean_live_60),
            dbh_gini_live_60 = ineq(dbh_mm[dbh_mm >= 60 & status %in% c(1:4)], type = "Gini"),
            dbh_gini_live_60 = ifelse(dbh_gini_live_60 %in% "NaN", NA, dbh_gini_live_60),
            ba_live_60 = sum(ba[dbh_mm >= 60 & status %in% c(1:4)], na.rm = T),
            n_trees_live_100 = length(treeid[dbh_mm >= 100 & status %in% c(1:4)]),
            dbh_quadrmean_live_100 = sqrt(mean(dbh_mm[dbh_mm >= 100 & status %in% c(1:4)]^2, na.rm = T)),
            dbh_quadrmean_live_100 = ifelse(dbh_quadrmean_live_100 %in% "NaN", NA, dbh_quadrmean_live_100),
            dbh_mean_live_100 = mean(dbh_mm[dbh_mm >= 100 & status %in% c(1:4)], na.rm = T),
            dbh_mean_live_100 = ifelse(dbh_mean_live_100 %in% "NaN", NA, dbh_mean_live_100),
            dbh_gini_live_100 = ineq(dbh_mm[dbh_mm >= 100 & status %in% c(1:4)], type = "Gini"),
            dbh_gini_live_100 = ifelse(dbh_gini_live_100 %in% "NaN", NA, dbh_gini_live_100),
            ba_live_100 = sum(ba[dbh_mm >= 100 & status %in% c(1:4)], na.rm = T)) %>%
  mutate(ba_dead_60 = ifelse(ba_dead_60 %in% 0 , NA, ba_dead_60),
         ba_dead_100 = ifelse(ba_dead_100 %in% 0, NA, ba_dead_100),
         ba_live_60 = ifelse(ba_live_60 %in% 0, NA, ba_live_60),
         ba_live_100 = ifelse(ba_live_100 %in% 0, NA, ba_live_100)) %>%
  mutate_at(vars(n_trees_live_500, n_trees_live_700, n_trees_dead_500, n_trees_dead_700, n_trees_dead_60,
                 ba_dead_60, n_trees_dead_100, ba_dead_100, n_trees_live_60, ba_live_60, n_trees_live_100, ba_live_100), 
            funs(.*10000/plotsize)) %>%
  mutate_at(vars(dbh_gini_live_60, dbh_gini_live_100), funs(round(., 2))) %>%
  select(-plotsize) ->
  tree_parameters


# # height_max ------------------------------------------------------------


tree %>%
  group_by(date, plotid) %>%
  summarise(height_max = max(height_m)) %>%
  collect() ->
  height_max


# # aspect_southness ------------------------------------------------------


plot %>%
  collect() %>%
  filter(!aspect %in% NA) %>%
  group_by(date, plotid) %>%
  summarise(aspect_southness = cos(deg2rad(45) - deg2rad(aspect)) + 1) ->
  aspect_southness


# # dominant_species ------------------------------------------------------


tree %>%
  filter(!onplot %in% 0,
         status %in% c(1:4)) %>%
  mutate(ba = pi * dbh_mm ^ 2 / 4 / 1000000) %>%
  group_by(date, plotid, species) %>%
  summarise(ba = sum(ba) * 10000 / min(plotsize)) %>%
  group_by(date, plotid) %>%
  arrange(desc(ba)) %>%
  filter(row_number() == 1) %>%
  mutate(dominant_species = species) %>%
  select(date, plotid, dominant_species) %>%
  collect() ->
  dominant_species


# # volume_dead_lying_decay -----------------------------------------------


deadwood %>%
  group_by(date, plotid, decay) %>%
  summarise(volume_dead_lying_decay = ((pi ^ 2 * sum((dbh_mm * 0.001) ^ 2)) / 800) * 10000) %>%
  filter(!decay %in% 99) %>%
  collect() %>%
  mutate(decay = paste0("volume_dead_lying_decay", decay)) %>%
  spread(decay, volume_dead_lying_decay, fill = 0) ->
  volume_dead_lying_decay


# # volume_dead_lying -----------------------------------------------------


deadwood %>%
  group_by(date, plotid) %>%
  summarise(volume_dead_lying = ((pi ^ 2 * sum((dbh_mm * 0.001) ^ 2)) / 800) * 10000) %>%
  collect() ->
  volume_dead_lying 


# # age (mean, 90quantile, 5oldest, gini, median, iqr) --------------------


core %>% 
  filter(treetype %in% c('0', "m", "x"),    
         growth %in% 1,
         !dbh_mm < 100,
         !corestatus %in% c(2, 3)) %>% 
  select(date, plotid, treeid, subcore, missing_years, id) %>%
  inner_join(.,
             ring,
             by = c('id' = 'core_id')) %>%
  collect() %>%
  filter(missing_years <= 20 | missing_years %in% NA) %>%
  group_by(date, plotid, treeid, subcore) %>%
  summarise(age = sum(n(), min(missing_years), na.rm = T)) %>%
  group_by(treeid) %>%
  arrange(desc(age)) %>%
  filter(row_number() == 1) %>%
  group_by(date, plotid) %>%
  arrange(desc(age)) %>%
  summarise(age_mean = mean(age),
            age_90quantile = quantile(age, 0.90),
            age_5oldest = mean(age[1:5]),
            age_gini = round(ineq(age, type = "Gini"), 2),
            age_median = median(age),
            age_iqr = interquantilerange(age)) ->
  age


# # regeneration_htclass --------------------------------------------------


regeneration %>%
  collect() %>%
  mutate(plotsize = ifelse(plotsize > 1000, 1000, plotsize)) %>%
  group_by(date, plotid, htclass) %>%
  summarise(regeneration_htclass = sum(count) * 10000 / min(plotsize)) %>%
  filter(!htclass %in% c(0, 3, 99)) %>%
  mutate(htclass = as.character(htclass), 
         htclass = ifelse(htclass %in% 1, "50_130", "130_250"),
         htclass = paste("regeneration", htclass, sep = "_"),
         regeneration_htclass = ifelse(regeneration_htclass %in% NA, -1, regeneration_htclass)) %>%
  spread(htclass, regeneration_htclass, fill = 0) %>%
  mutate(regeneration_130_250 = ifelse(regeneration_130_250 %in% -1, NA, regeneration_130_250),
         regeneration_50_130 = ifelse(regeneration_50_130 %in% -1, NA, regeneration_50_130)) ->
  regeneration_htclass


# # regeneration_250_dbh_min ----------------------------------------------


regeneration %>%
  collect() %>%
  mutate(plotsize = ifelse(plotsize > 1000, 1000, plotsize),
         foresttype = ifelse(foresttype %in% "spruce", "regeneration_250_100", "regeneration_250_60")) %>%
  filter(htclass %in% 3) %>%
  group_by(date, plotid, foresttype) %>%
  summarise(regeneration_250_dbh_min = sum(count) * 10000 / min(plotsize)) %>%
  spread(foresttype, regeneration_250_dbh_min) ->
  regeneration_250_dbh_min


# # regeneration_0_50 -----------------------------------------------------


regeneration_subplot %>%
  filter(htclass %in% 0) %>%
  group_by(date, plotid) %>%
  summarise(regeneration_0_50 = sum(count) * 10000 / 20) %>%
  collect() -> 
  regeneration_0_50


# # regeneration_browsing -------------------------------------------------


regeneration_subplot %>%
  filter(htclass %in% 1,
         !browsing %in% c(5, 99)) %>%
  group_by(date, plotid) %>%
  summarise(regeneration_browsing = mean(browsing)) %>%
  collect() -> 
  regeneration_browsing
  

# # nearest_plot ----------------------------------------------------------


plot.df <- plot %>% filter(census %in% 1) %>% collect()

P <- plot.df %>% select(lat, lng) %>%  collect() %>% as.matrix()
p <- distm(P, P)
diag(p) <- NA

plot.key <- plot.df %>% mutate(id = row_number()) %>% select(id, plotid) %>% deframe()

reshape2::melt(as.matrix(p), varnames = c("plotid", "second_plot"), value.name = "nearest_plot") %>%
  filter(!nearest_plot %in% NA) %>%
  mutate(plotid = plot.key[plotid],
         second_plot = plot.key[second_plot]) %>%
  inner_join(., plot.df[ , c("date", "plotid", "foresttype")], by = "plotid") %>%
  group_by(plotid) %>%
  mutate(nearest_plot = ifelse(foresttype %in% "spruce",
                               min(nearest_plot),
                               sort(nearest_plot)[2])) %>%
  filter(row_number() == 1) %>%
  arrange(plotid) %>%
  select(-second_plot, -foresttype) -> 
  nearest_plot

plot.remeasured <- plot %>% filter(!census %in% 1) %>% collect()

plot.remeasured %>% 
  select(date, plotid) %>%
  inner_join(., nearest_plot[, c("plotid", "nearest_plot")], by = "plotid") %>%
  bind_rows(., nearest_plot) -> 
  nearest_plot

# # volume_dead_standing --------------------------------------------------


tree %>%
  filter(!onplot %in% 0,
         status %in% c(11:23),
         !decayht %in% 99) %>%
  collect() %>%
  filter(!dbh_mm %in% NA) %>%
  mutate(
    decayht = case_when(
      decayht == 0 ~ 5,
      decayht == 1 ~ 15,
      decayht == 2 ~ 25,
      decayht == 3 ~ 35,
      decayht == 4 ~ 45,
      decayht == 5 ~ 55)) %>%
  rowwise() %>%
  mutate(volume_snag = E_VOL_AB_HmDm_HT.f(Hm=1.3, Dm=(dbh_mm * 0.1), mHt = (log(dbh_mm * 0.1)-1.08261)^2/0.275541, sHt = 0, par.lme = SK.par.lme, A=0, B=decayht, iDH = "H")$E_VOL) %>%
  group_by(date, plotid) %>%
  summarise(volume_dead_standing_60 = sum(volume_snag[dbh_mm >= 60]) * 10000 / min(plotsize),
            volume_dead_standing_100 = sum(volume_snag[dbh_mm >= 100]) * 10000 / min(plotsize)) ->
  volume_dead_standing


# # biomass + volume ------------------------------------------------------


tree %>%
  filter(!onplot %in% 0,
         status %in% c(1:4)) %>%
  collect() %>%
  left_join(., tree_parameters %>% select(date, plotid, ba_live_60), by = c("date", "plotid")) %>%
  mutate(species_group = case_when(
    species %in% "Fagus sylvatica" ~ "fagus",
    species %in% "Acer pseudoplatanus" ~ "acer",
    species %in% "Acer platanoides" ~ "broadleaves",
    species %in% "Acer" ~ "broadleaves",
    species %in% "Ulmus" ~ "broadleaves",
    species %in% "Sorbus aria" ~ "broadleaves",
    species %in% "Acer obtusifolium" ~ "broadleaves",
    species %in% "Sorbus aucuparia" ~ "broadleaves",
    species %in% "Betula pendula" ~ "broadleaves",
    species %in% "Salix" ~ "broadleaves",
    species %in% "Sambucus racemosa" ~ "broadleaves",
    species %in% "Fraxinus excelsior" ~ "broadleaves",
    species %in% "Tilia cordata" ~ "broadleaves",
    species %in% "Tilia" ~ "broadleaves",
    species %in% "Acer obtusatum" ~ "broadleaves",
    species %in% "Rhamnus" ~ "broadleaves",
    species %in% "Lians" ~ "broadleaves",
    species %in% "Populus tremula" ~ "broadleaves",
    species %in% "Corylus avellana" ~ "broadleaves",
    species %in% "Fraxinus" ~ "broadleaves",
    species %in% "Ulmus glabra" ~ "broadleaves",
    species %in% "Broadleaves" ~ "broadleaves",
    species %in% "Sambucus nigra" ~ "broadleaves",
    species %in% "Salix nigra" ~ "broadleaves",
    species %in% "Laburnum anagyroides" ~ "broadleaves",
    species %in% "Fraxinus ornus" ~ "broadleaves",
    species %in% "Salix caprea" ~ "broadleaves",
    species %in% "Betula" ~ "broadleaves",
    species %in% "Carpinus betulus" ~ "broadleaves",
    species %in% "Larix decidua" ~ "larix",
    species %in% "Abies alba" ~ "abies",
    species %in% "Picea abies" ~ "picea",
    species %in% "Abies" ~ "abies",
    species %in% "Coniferous" ~ "coniferous",
    species %in% "Taxus baccata" ~ "coniferous",
    species %in% "Pinus sylvestris" ~ "pinus",
    species %in% "Pinus cembra" ~ "pinus"
  ),
  TB = case_when(
    species_group %in% "broadleaves" ~ exp(-3.7241 + 2.4069 * log(dbh_mm * 0.1)),
    species_group %in% "acer" ~ exp(-3.7241 + 2.4069 * log(dbh_mm * 0.1)),
    species_group %in% "coniferous" ~ exp(-3.248 + 2.3695 * log(dbh_mm * 0.1) + (-0.0254 * ba_live_60)),
    species_group %in% "fagus" ~ exp(-3.7694 + 2.8003 * log(dbh_mm * 0.1) + (-0.0247 * ba_live_60)),
    species_group %in% "larix" ~ exp(-3.2409 + 2.1412 * log(dbh_mm * 0.1)),
    species_group %in% "picea" ~ exp(-3.3163 + 2.1983 * log(dbh_mm * 0.1)),
    species_group %in% "abies" ~ exp(-3.3163 + 2.1983 * log(dbh_mm * 0.1)),
    species_group %in% "pinus" ~ exp(-3.6641 + 2.1601 * log(dbh_mm * 0.1))
  ),
  FM = case_when(
    species_group %in% "broadleaves" ~ exp(-4.2286 + 1.8625 * log(dbh_mm * 0.1)),
    species_group %in% "acer" ~ exp(-4.0625 + 2.0662 * log(dbh_mm * 0.1)),
    species_group %in% "coniferous" ~ exp(-2.6019 + 2.1097 * log(dbh_mm * 0.1) + (-0.0404 * ba_live_60)),
    species_group %in% "fagus" ~ exp(-4.4813 + 1.9073 * log(dbh_mm * 0.1)),
    species_group %in% "larix" ~ exp(-3.8849 + 1.7502 * log(dbh_mm * 0.1)),
    species_group %in% "picea" ~ exp(-2.1305 + 2.0087 * log(dbh_mm * 0.1) + (-0.0324 * ba_live_60)),
    species_group %in% "abies" ~ exp(-2.1305 + 2.0087 * log(dbh_mm * 0.1) + (-0.0324 * ba_live_60)),
    species_group %in% "pinus" ~ exp(-2.4122 + 1.8683 * log(dbh_mm * 0.1) + (-0.0537 * ba_live_60))
  ),
  RM = case_when(
    species_group %in% "broadleaves" ~ exp(-2.6183 + 2.1353 * log(dbh_mm * 0.1)),
    species_group %in% "acer" ~ exp(-2.6183 + 2.1353 * log(dbh_mm * 0.1)),
    species_group %in% "coniferous" ~ exp(-4.0287 + 2.4957 * log(dbh_mm * 0.1)),
    species_group %in% "fagus" ~ exp(-3.1432 + 2.3794 * log(dbh_mm * 0.1) + (-0.0125 * ba_live_60)),
    species_group %in% "larix" ~ exp(-3.6347 + 2.3038 * log(dbh_mm * 0.1)),
    species_group %in% "picea" ~ exp(-3.7387 + 2.4323 * log(dbh_mm * 0.1)),
    species_group %in% "abies" ~ exp(-3.7387 + 2.4323 * log(dbh_mm * 0.1)),
    species_group %in% "pinus" ~ exp(-3.6347 + 2.3038 * log(dbh_mm * 0.1))
  ),
  SM = case_when(
    species_group %in% "broadleaves" ~ exp(-2.4521 + 2.4115 * log(dbh_mm * 0.1)),
    species_group %in% "acer" ~ exp(-2.4521 + 2.4115 * log(dbh_mm * 0.1)),
    species_group %in% "coniferous" ~ exp(-2.7693 + 2.3761 * log(dbh_mm * 0.1) + (0.0072 * ba_live_60)),
    species_group %in% "fagus" ~ exp(-1.4487 + 2.1661 * log(dbh_mm * 0.1)),
    species_group %in% "larix" ~ exp(-2.4105 + 2.424 * log(dbh_mm * 0.1)),
    species_group %in% "picea" ~ exp(-2.5027 + 2.3404 * log(dbh_mm * 0.1)),
    species_group %in% "abies" ~ exp(-3.2683 + 2.5768 * log(dbh_mm * 0.1)),
    species_group %in% "pinus" ~ exp(-2.3583 + 2.308 * log(dbh_mm * 0.1))
  ),
  volume_live = case_when(
    species %in% "Fagus sylvatica" ~ volume(SM, 0.59), 
    species %in% "Acer pseudoplatanus" ~ volume(SM, 0.51),
    species %in% "Acer platanoides" ~ volume(SM, 0.51),
    species %in% "Acer" ~ volume(SM, 0.51),
    species %in% "Ulmus" ~ volume(SM, 0.55),
    species %in% "Sorbus aria" ~ volume(SM, 0.63),
    species %in% "Acer obtusifolium" ~ volume(SM, 0.51),
    species %in% "Sorbus aucuparia" ~ volume(SM, 0.63),
    species %in% "Betula pendula" ~ volume(SM, 0.53),
    species %in% "Salix" ~ volume(SM, 0.35),
    species %in% "Sambucus racemosa" ~ volume(SM, 0.45),
    species %in% "Fraxinus excelsior" ~ volume(SM, 0.56) ,
    species %in% "Tilia cordata" ~ volume(SM, 0.42),
    species %in% "Tilia" ~ volume(SM, 0.42),
    species %in% "Acer obtusatum" ~ volume(SM, 0.51),
    species %in% "Rhamnus" ~ volume(SM, 0.61),
    species %in% "Lians" ~ volume(SM, 0.54),
    species %in% "Populus tremula" ~ volume(SM, 0.37),
    species %in% "Corylus avellana" ~ volume(SM, 0.52),
    species %in% "Fraxinus" ~ volume(SM, 0.56),
    species %in% "Ulmus glabra" ~ volume(SM, 0.55),
    species %in% "Broadleaves" ~ volume(SM, 0.54),
    species %in% "Sambucus nigra" ~ volume(SM,0.45),
    species %in% "Salix nigra" ~ volume(SM, 0.35),
    species %in% "Laburnum anagyroides" ~ volume(SM, 0.72),
    species %in% "Fraxinus ornus" ~ volume(SM, 0.56),
    species %in% "Salix caprea" ~ volume(SM, 0.35),
    species %in% "Betula" ~ volume(SM, 0.53),
    species %in% "Carpinus betulus" ~ volume(SM, 0.71),
    species %in% "Larix decidua" ~ volume(SM, 0.47),
    species %in% "Abies alba" ~ volume(SM, 0.35),
    species %in% "Picea abies" ~ volume(SM, 0.37),
    species %in% "Abies" ~ volume(SM, 0.35),
    species %in% "Coniferous" ~ volume(SM, 0.41),
    species %in% "Taxus baccata" ~ volume(SM, 0.55),
    species %in% "Pinus sylvestris" ~ volume(SM, 0.42),
    species %in% "Pinus cembra" ~ volume(SM, 0.42)
  ),
  biomass_aboveground = TB + FM + SM,
  biomass_underground = RM) %>%
  group_by(date, plotid) %>%
  summarise(volume_live_60 = sum(volume_live[dbh_mm >= 60], na.rm = T) * 10000 / min(plotsize),
            volume_live_100 = sum(volume_live[dbh_mm >= 100], na.rm = T) * 10000 / min(plotsize),
            biomass_aboveground_60 = sum(biomass_aboveground[dbh_mm >= 60], na.rm = T) * 10000 / min(plotsize),
            biomass_aboveground_100 = sum(biomass_aboveground[dbh_mm >= 100], na.rm = T) * 10000 / min(plotsize),
            biomass_underground_60 = sum(biomass_underground[dbh_mm >= 60], na.rm = T) * 10000 / min(plotsize),
            biomass_underground_100 = sum(biomass_underground[dbh_mm >= 100], na.rm = T) * 10000 / min(plotsize)) %>%
  mutate(volume_live_60 = ifelse(volume_live_60 %in% 0, NA, volume_live_60),
         volume_live_100 = ifelse(volume_live_100 %in% 0, NA, volume_live_100),
         biomass_aboveground_60 = ifelse(biomass_aboveground_60 %in% 0, NA, biomass_aboveground_60),
         biomass_aboveground_100 = ifelse(biomass_aboveground_100 %in% 0, NA, biomass_aboveground_100),
         biomass_underground_60 = ifelse(biomass_underground_60 %in% 0, NA, biomass_underground_60),
         biomass_underground_100 = ifelse(biomass_underground_100 %in% 0, NA, biomass_underground_100)) -> 
  biomass


# # disturbance parameters ------------------------------------------------


dist_plot_event %>%
  collect() %>%
  group_by(date, plotid) %>% 
  complete(year = min(year):max(year), fill = list(kde = 0)) %>% 
  mutate(disturbance_year = year[which.max(kde)],
         disturbance_severity = max(kde),
         disturbance_over_15ca = length(id[kde > 15]),
         disturbance_first_year = min(year[kde > 15]),
         disturbance_first_year = ifelse(disturbance_over_15ca %in% 0, NA, disturbance_first_year),
         disturbance_last_year = max(year[kde > 15]),
         disturbance_last_year = ifelse(disturbance_over_15ca %in% 0, NA, disturbance_last_year),
         disturbance_severity_mean = mean(kde[kde > 15]),
         disturbance_severity_mean = ifelse(disturbance_over_15ca %in% 0, NA, disturbance_severity_mean),
         disturbance_first_severity = kde[year == disturbance_first_year],
         disturbance_last_severity = kde[year == disturbance_last_year],
         ms30 = movingSum(kde),
         disturbance_max_30y_severity = max(ms30),
         disturbance_max_30y_year = median(year[ms30 == disturbance_max_30y_severity]),
         disturbance_1965_severity = mean(kde[year >= 1965 & kde > 15]),
         disturbance_1965_severity = ifelse(disturbance_1965_severity %in% "NaN", NA, disturbance_1965_severity)) %>%
  select(date, plotid, disturbance_year, disturbance_severity, disturbance_over_15ca, disturbance_severity_mean,
         disturbance_first_year, disturbance_first_severity, disturbance_last_year, disturbance_last_severity,
         disturbance_max_30y_year, disturbance_max_30y_severity, disturbance_1965_severity) %>%
  filter(row_number() == 1) -> 
  disturbance_parameters


# # disturbance_recent ----------------------------------------------------


tree %>% 
  select(date, plotid, country, foresttype, location, species, onplot, dbh_mm, growth, decay) %>%
  filter(onplot %in% c(1:3), growth %in% c(1,-1,99)) %>%              
  inner_join(., tbl(KELuser, 'species_fk') %>% select(species = id, sp_group_dist), by = "species") %>%
  inner_join(., tbl(KELuser, 'dist_group'), by = c("country", "foresttype", "location", "sp_group_dist")) %>% 
  inner_join(., tbl(KELuser, 'dist_param') %>% select(dist_param = id, Tdbh_mm = dbh_mm, dbh_ca_f), by = "dist_param") %>% 
  filter(dbh_mm >= Tdbh_mm) %>%
  select(date, plotid, dbh_mm, dbh_ca_f, decay) %>% 
  collect() %>%
  rowwise() %>% 
  mutate(ca = eval(parse(text = dbh_ca_f))) %>%
  group_by(date, plotid) %>%
  summarise(disTot = sum(ca[decay %in% c(-1:3)]),
            rDis1_3N = sum(ca[decay %in% c(1:3)])) %>%
  mutate(disturbance_recent = rDis1_3N/disTot * 100,
         disturbance_recent = ifelse(disturbance_recent %in% "NaN" | disturbance_recent %in% 0, NA, disturbance_recent)) %>%
  select(plotid, date, disturbance_recent) ->
  disturbance_recent


# # disturbance_index -----------------------------------------------------


dist_plot_event %>%
  collect() %>%
  unite(plotid, date, plotid, sep = "-") %>%
  mutate(year = case_when(
    year >= 1800 & year < 1810 ~ 1800,
    year >= 1810 & year < 1820 ~ 1810,
    year >= 1820 & year < 1830 ~ 1820,
    year >= 1830 & year < 1840 ~ 1830,
    year >= 1840 & year < 1850 ~ 1840,
    year >= 1850 & year < 1860 ~ 1850,
    year >= 1860 & year < 1870 ~ 1860,
    year >= 1870 & year < 1880 ~ 1870,
    year >= 1880 & year < 1890 ~ 1880,
    year >= 1890 & year < 1900 ~ 1890,
    year >= 1900 & year < 1910 ~ 1900,
    year >= 1910 & year < 1920 ~ 1910,
    year >= 1920 & year < 1930 ~ 1920,
    year >= 1930 & year < 1940 ~ 1930,
    year >= 1940 & year < 1950 ~ 1940,
    year >= 1950 & year < 1960 ~ 1950,
    year >= 1960 & year < 1970 ~ 1960,
    year >= 1970 & year < 1980 ~ 1970,
    year >= 1980 & year < 1990 ~ 1980
  )) %>%
  group_by(plotid, year) %>%
  summarise(ca = sum(ca_per)) %>%
  filter(!year %in% NA) %>%
  spread(., year, ca, fill = 0) %>%
  ungroup() %>%
  mutate(disturbance_index = diversity(.[ ,c(2:20)])) %>%
  separate(plotid, c("date", "plotid"), sep = "-") %>%
  mutate(date = as.numeric(date),
         disturbance_index = round(disturbance_index, 2)) %>%
  select(date, plotid, disturbance_index) -> 
  disturbance_index


# # openness mean, gini ---------------------------------------------------


canopy_analysis %>%
  filter(parameter %in% "openness") %>%
  collect() %>%
  group_by(date, plotid) %>%
  summarise(openness_mean = mean(value),
            openness_gini = ineq(value, type = "Gini")) %>%
  mutate_at(vars(openness_mean, openness_gini), funs(round(., 2))) ->
  openness


# # temperature data -------------------------------------------------------------


 # setwd("C:/Users/Ondrej_Vostarek/Desktop/MVP/DB/database/data/temperatures")
 # 
 # temp <- read.table("plots_database_temperatures.txt")


# collect everything ----------------------------------------------------


dominant_species %>%
  full_join(., height_max, by = c("date", "plotid")) %>%
  full_join(., volume_dead_lying, by = c("date", "plotid")) %>%
  full_join(., regeneration_0_50, by = c("date", "plotid")) %>%
  full_join(., regeneration_browsing, by = c("date", "plotid")) %>%
  full_join(., nearest_plot, by = c("date", "plotid")) %>%
  full_join(., aspect_southness, by = c("date", "plotid")) %>%
  full_join(., age, by = c("date", "plotid")) %>%  
  full_join(., regeneration_250_dbh_min, by = c("date", "plotid")) %>%
  full_join(., regeneration_htclass, by = c("date", "plotid")) %>%
  full_join(., tree_parameters, by = c("date", "plotid")) %>%
  full_join(., volume_dead_lying_decay, by = c("date", "plotid")) %>%
  full_join(., volume_dead_standing, by = c("date", "plotid")) %>%
  full_join(., biomass, by = c("date", "plotid")) %>%
  full_join(., disturbance_parameters, by = c("date", "plotid")) %>%
  full_join(., disturbance_recent, by = c("date", "plotid")) %>%
  full_join(., disturbance_index, by = c("date", "plotid")) %>%
  full_join(., openness, by = c("date", "plotid")) %>%
  mutate_at(vars(-dominant_species, -disturbance_index, -openness_mean, -openness_gini, -age_gini, -dbh_gini_live_60, -dbh_gini_live_100),
            funs(round(., 0))) -> 
  final_table



# wide to long ------------------------------------------------------------


final_table %>% 
  gather(., parameter, value, dominant_species:openness_gini) %>%
  filter(!value %in% NA) -> final_table


# disconnection -----------------------------------------------------------


poolClose(KELuser)



