deg2rad <- function(deg) {
  #'@description converts degrees into radians
  #'@param deg value in degrees
  
  (deg * pi) / 180
  
}

interquantilerange <- function(x){
  #'@description calculate the range between 0.75 and 0.25 quantile
  #'@param x vector of numerical values
  
  quantile(x, 0.75, na.rm = T) - quantile(x, 0.25, na.rm = T)
  
}

movingSum <- function(x, windowLength = 30){
  #'@description calculate the moving sum of values
  #'@param x vector of numerical values
  #'@param windowLength length of the moving window
  
  rollapply( x, 
             width = windowLength,
             FUN = sum,
             fill = NA,
             align = "center",
             na.rm = T,
             partial = TRUE)
}

volume <- function(SM, y){
  #'@description calculate tree volume in m3 from stem mass and wood density
  #'@param SM stem mass in kg
  #'@param y wood density in g/cm3
  
  (SM * 0.001) / y
  
}

get_data <- function(plot.id){
  #'@description collect necessary data from database and return them as data.list
  #'@param plot.id id of plots as in 'plot' table in database
  
  plot <- tbl(KELuser, "plot") %>% 
    filter(id %in% plot.id) %>% 
    select(plot_id = id, date, plotid, census, country, location, lng, lat, plotsize, foresttype, altitude_m, aspect)
  
  tree <- tbl(KELuser, "tree") %>% 
    select(tree_id = id, plot_id, treeid, onplot, treetype, status, growth, species, dbh_mm, height_m, decayht, decay) %>% 
    inner_join(., plot, by = 'plot_id')
  
  deadwood <- tbl(KELuser, "deadwood") %>%  
    select(plot_id, dbh_mm, decay) %>% 
    inner_join(., plot, by = 'plot_id')
  
  core <- tbl(KELuser, "core") %>% 
    select(core_id = id, tree_id, subcore, missing_years, corestatus) %>% 
    inner_join(., tree, by = 'tree_id')
  
  ring <- tbl(KELuser, "ring") 
  
  regeneration <- tbl(KELuser, "regeneration") %>% 
    select(plot_id, htclass, count) %>%
    inner_join(., plot, by = "plot_id")
  
  regeneration_subplot <- tbl(KELuser, "regeneration_subplot") %>% 
    select(plot_id, htclass, browsing, count) %>%
    inner_join(., plot, by = "plot_id")
  
  dist_plot_event <- tbl(KELuser, "dist_plot_event") %>% 
    inner_join(., tbl(KELuser, "dist_plot"), by = c("dist_plot_id" = "id")) %>%
    inner_join(., plot, by = "plot_id") %>%
    select(id, plot_id, year, ca_per, kde)
  
  canopy_analysis <- tbl(KELuser, "canopy_analysis") %>% inner_join(., plot, by = "plot_id")
 
  data <- list()
  
  data$plot <- collect(plot)
  data$tree <- collect(tree)
  data$deadwood <- collect(deadwood)
  data$core <- collect(core)
  data$ring <- collect(ring)
  data$regeneration <- collect(regeneration)
  data$regeneration_subplot <- collect(regeneration_subplot)
  data$dist_plot_event <- collect(dist_plot_event)
  data$canopy_analysis <- collect(canopy_analysis)
  
  return(data) 
  
}

calculate_parameters <- function(data){
  #'@description calculate the structural parameters
  #'@param data data.list as produced by the get_data function
  
  parameters <- list()
  
  # aspect southness
  
  parameters$aspect_southness <- data$plot %>% filter(!aspect %in% NA) %>% group_by(plot_id) %>%
    summarise(aspect_southness = cos(deg2rad(45) - deg2rad(aspect)) + 1)
  
  # nearest plot
  
  plot.new <- data$plot %>% filter(census %in% 1)
  
  P <- plot.new %>% select(lat, lng) %>%  collect() %>% as.matrix()
  p <- distm(P, P)
  diag(p) <- NA
  
  plot.key <- plot.new %>% mutate(id = row_number()) %>% select(id, plotid) %>% deframe()
  
  plot.nearest <- reshape2::melt(as.matrix(p), varnames = c("plotid", "second_plot"), value.name = "nearest_plot") %>%
    filter(!nearest_plot %in% NA) %>%
    mutate(plotid = plot.key[plotid],
           second_plot = plot.key[second_plot]) %>%
    inner_join(., plot.new[ , c("plot_id", "plotid", "foresttype")], by = "plotid") %>%
    group_by(plotid) %>%
    mutate(nearest_plot = ifelse(foresttype %in% "spruce",
                                 min(nearest_plot),
                                 sort(nearest_plot)[2])) %>%
    filter(row_number() == 1) %>%
    arrange(plotid)

  plot.remeasured <- data$plot %>% filter(!census %in% 1) %>% collect()
  
  parameters$nearest_plot <- plot.remeasured %>% 
    select(plot_id, plotid) %>%
    inner_join(., plot.nearest[, c("plotid", "nearest_plot")], by = "plotid") %>%
    bind_rows(., plot.nearest) %>%
    select(plot_id, nearest_plot)
  
  # tree parameters
  
  parameters$tree_parameters <- data$tree %>%
    filter(!onplot %in% 0) %>%
    mutate(ba = pi * dbh_mm ^ 2 / 4 / 1000000) %>%
    group_by(plot_id) %>%
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
    select(-plotsize)
  
  parameters$height_max <- data$tree %>% group_by(plot_id) %>% summarise(height_max = max(height_m))
  
  parameters$dominant_species <- data$tree %>%
    filter(!onplot %in% 0,
           status %in% c(1:4)) %>%
    mutate(ba = pi * dbh_mm ^ 2 / 4 / 1000000) %>%
    group_by(plot_id, species) %>%
    summarise(ba = sum(ba) * 10000 / min(plotsize)) %>%
    group_by(plot_id) %>%
    arrange(desc(ba)) %>%
    filter(row_number() == 1) %>%
    select(plot_id, dominant_species = species)
  
  # biomass and volume of alive trees
  
  parameters$biomass_volume <- data$tree %>%
    filter(!onplot %in% 0,
           status %in% c(1:4)) %>%
    left_join(., parameters$tree_parameters %>% select(plot_id, ba_live_60), by = "plot_id") %>%
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
    group_by(plot_id) %>%
    summarise(plotsize = first(plotsize),
              volume_live_60 = sum(volume_live[dbh_mm >= 60], na.rm = T),
              volume_live_100 = sum(volume_live[dbh_mm >= 100], na.rm = T),
              biomass_aboveground_60 = sum(biomass_aboveground[dbh_mm >= 60], na.rm = T),
              biomass_aboveground_100 = sum(biomass_aboveground[dbh_mm >= 100], na.rm = T),
              biomass_underground_60 = sum(biomass_underground[dbh_mm >= 60], na.rm = T),
              biomass_underground_100 = sum(biomass_underground[dbh_mm >= 100], na.rm = T)) %>%
    mutate_at(vars(-plot_id), funs(.*10000/plotsize)) %>%
    mutate(volume_live_60 = ifelse(volume_live_60 %in% 0, NA, volume_live_60),
           volume_live_100 = ifelse(volume_live_100 %in% 0, NA, volume_live_100),
           biomass_aboveground_60 = ifelse(biomass_aboveground_60 %in% 0, NA, biomass_aboveground_60),
           biomass_aboveground_100 = ifelse(biomass_aboveground_100 %in% 0, NA, biomass_aboveground_100),
           biomass_underground_60 = ifelse(biomass_underground_60 %in% 0, NA, biomass_underground_60),
           biomass_underground_100 = ifelse(biomass_underground_100 %in% 0, NA, biomass_underground_100)) %>%
    select(-plotsize)

  # volume of dead standing trees
  
  parameters$volume_dead_standing <- data$tree %>%
    filter(!onplot %in% 0,
           status %in% c(11:23),
           !decayht %in% 99,
           !dbh_mm %in% NA) %>%
    mutate(
      decayht = case_when(
        decayht == 0 ~ 5,
        decayht == 1 ~ 15,
        decayht == 2 ~ 25,
        decayht == 3 ~ 35,
        decayht == 4 ~ 45,
        decayht == 5 ~ 55)) %>%
    rowwise() %>%
    mutate(volume_snag = E_VOL_AB_HmDm_HT.f(Hm=1.3, Dm=(dbh_mm * 0.1), 
                                            mHt = (log(dbh_mm * 0.1)-1.08261)^2/0.275541, 
                                            sHt = 0, par.lme = SK.par.lme, A=0, B=decayht, iDH = "H")$E_VOL) %>%
    group_by(plot_id) %>%
    summarise(volume_dead_standing_60 = sum(volume_snag[dbh_mm >= 60]) * 10000 / min(plotsize),
              volume_dead_standing_100 = sum(volume_snag[dbh_mm >= 100]) * 10000 / min(plotsize))
  
  # recent disturbance
  
  parameters$disturbance_recent <- data$tree %>% 
    filter(onplot %in% c(1:3), 
           growth %in% c(1,-1,99)) %>%              
    inner_join(., tbl(KELuser, 'species_fk') %>% select(species = id, sp_group_dist), by = "species", copy = TRUE) %>%
    inner_join(., tbl(KELuser, 'dist_group'), by = c("country", "foresttype", "location", "sp_group_dist"), copy = TRUE) %>% 
    inner_join(., tbl(KELuser, 'dist_param') %>% select(dist_param = id, Tdbh_mm = dbh_mm, dbh_ca_f), by = "dist_param", copy = TRUE) %>% 
    filter(dbh_mm >= Tdbh_mm) %>%
    rowwise() %>% 
    mutate(ca = eval(parse(text = dbh_ca_f))) %>%
    group_by(plot_id) %>%
    summarise(disTot = sum(ca[decay %in% c(-1:3)]),
              rDis1_3N = sum(ca[decay %in% c(1:3)])) %>%
    mutate(disturbance_recent = rDis1_3N/disTot * 100,
           disturbance_recent = ifelse(disturbance_recent %in% "NaN" | disturbance_recent %in% 0, NA, disturbance_recent)) %>%
    select(plot_id, disturbance_recent)
  
  # disturbance parameters
  
  parameters$disturbance_parameters <- data$dist_plot_event %>%
    group_by(plot_id) %>% 
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
    select(plot_id, disturbance_year, disturbance_severity, disturbance_over_15ca, disturbance_severity_mean,
           disturbance_first_year, disturbance_first_severity, disturbance_last_year, disturbance_last_severity,
           disturbance_max_30y_year, disturbance_max_30y_severity, disturbance_1965_severity) %>%
    filter(row_number() == 1)
  
  parameters$disturbance_index <- data$dist_plot_event %>%
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
    filter(!year %in% NA) %>%
    group_by(plot_id, year) %>%
    summarise(ca = sum(ca_per)) %>%
    spread(., year, ca, fill = 0) %>%
    ungroup() %>%
    mutate(disturbance_index = diversity(.[ ,c(2:ncol(.))]),
           disturbance_index = round(disturbance_index, 2)) %>%
    select(plot_id, disturbance_index)  
  
  # volume of laying deadwood
  
  parameters$volume_dead_lying_decay <- data$deadwood %>%
    filter(!decay %in% 99) %>%
    group_by(plot_id, decay) %>%
    summarise(volume_dead_lying_decay = ((pi ^ 2 * sum((dbh_mm * 0.001) ^ 2)) / 800) * 10000) %>%
    mutate(decay = paste0("volume_dead_lying_decay", decay)) %>%
    spread(decay, volume_dead_lying_decay, fill = 0)
  
  parameters$volume_dead_lying <- data$deadwood %>%
    group_by(plot_id) %>%
    summarise(volume_dead_lying = ((pi ^ 2 * sum((dbh_mm * 0.001) ^ 2)) / 800) * 10000)
  
  # age parameters
  
  parameters$age_parameters <- data$core %>% 
    filter(treetype %in% "0" & onplot %in% c(1, 2) | treetype %in% c("m", "x"),    
           growth %in% 1,
           !dbh_mm < 100,
           !corestatus %in% c(2, 3)) %>% 
    inner_join(., data$ring, by = "core_id") %>%
    filter(missing_years <= 20 | missing_years %in% NA) %>%
    group_by(date, plotid, treeid, subcore) %>%                            
    summarise(age = sum(n(), min(missing_years), na.rm = T)) %>%
    group_by(treeid) %>%
    arrange(desc(age)) %>%
    filter(row_number() == 1) %>%
    group_by(plotid) %>%
    mutate(age = ifelse(date %in% max(date), age - (max(date) - min(date)), age),
           date = min(date)) %>%
    inner_join(., data$plot %>% select(plot_id, date, plotid), by = c("date", "plotid")) %>%
    group_by(plot_id) %>%
    arrange(desc(age)) %>%
    summarise(age_mean = mean(age),
              age_90quantile = quantile(age, 0.90),
              age_5oldest = mean(age[1:5]),
              age_gini = round(ineq(age, type = "Gini"), 2),
              age_median = median(age),
              age_iqr = interquantilerange(age))
  
  # regeneration
  
  parameters$regeneration_htclass <- data$regeneration %>%
    filter(!htclass %in% c(0, 3, 99)) %>%
    mutate(plotsize = ifelse(plotsize > 1000, 1000, plotsize)) %>%
    group_by(plot_id, htclass) %>%
    summarise(regeneration_htclass = sum(count) * 10000 / min(plotsize)) %>%
    mutate(htclass = ifelse(htclass %in% 1, "regeneration_50_130", "regeneration_130_250"),
           regeneration_htclass = ifelse(regeneration_htclass %in% NA, -1, regeneration_htclass)) %>%
    spread(htclass, regeneration_htclass, fill = 0) %>%
    mutate(regeneration_130_250 = ifelse(regeneration_130_250 %in% -1, NA, regeneration_130_250),
           regeneration_50_130 = ifelse(regeneration_50_130 %in% -1, NA, regeneration_50_130))
  
  parameters$regeneration_250_dbh_min <- data$regeneration %>%
    filter(htclass %in% 3) %>%
    mutate(plotsize = ifelse(plotsize > 1000, 1000, plotsize),
           foresttype = ifelse(foresttype %in% "spruce", "regeneration_250_100", "regeneration_250_60")) %>%
    group_by(plot_id, foresttype) %>%
    summarise(regeneration_250_dbh_min = sum(count) * 10000 / min(plotsize)) %>%
    spread(foresttype, regeneration_250_dbh_min)
  
  # regeneration_subplot
  
  parameters$regeneration_0_50 <- data$regeneration_subplot %>%
    filter(htclass %in% 0) %>%
    group_by(plot_id) %>%
    summarise(regeneration_0_50 = sum(count) * 10000 / 20)
  
  parameters$regeneration_browsing <- data$regeneration_subplot %>%
    filter(htclass %in% 1,
           !browsing %in% c(5, 99)) %>%
    group_by(plot_id) %>%
    summarise(regeneration_browsing = mean(browsing))
  
  # canopy openness
  
  parameters$openness <- data$canopy_analysis %>%
    filter(parameter %in% "openness") %>%
    group_by(plot_id) %>%
    summarise(openness_mean = mean(value),
              openness_gini = ineq(value, type = "Gini")) %>%
    mutate_at(vars(openness_mean, openness_gini), funs(round(., 2)))
  
  # temperature
  
  load("C:/Users/Ondrej_Vostarek/Desktop/MVP/DB/data/temperatures/model_annual.rda")
  load("C:/Users/Ondrej_Vostarek/Desktop/MVP/DB/data/temperatures/model_vegetation.rda")
  
  plots <- data.list$plot %>% 
    select(plot_id, longitude = lng, latitude = lat, altitude = altitude_m, country, location)
  
  plots$temperature_annual <- stats::predict(model_annual,plots)^2-100
  plots$temperature_vegetation <- stats::predict(model_vegetation,plots)^2-100
  
  parameters$temperature <- plots %>% select(plot_id, temp_mean_year = temperature_annual, temp_mean_vegetseason = temperature_vegetation)
 
  return(parameters)
   
}

collect_data <- function(data){
  #'@description collect all data.list sheets into one data.frame
  #'@param data data.list as produced by the calculate_parameters function
  
  data.all <- data$aspect_southness %>%
    full_join(., data$nearest_plot, by = "plot_id") %>%
    full_join(., data$tree_parameters, by = "plot_id") %>%
    full_join(., data$height_max, by = "plot_id") %>%
    full_join(., data$dominant_species, by = "plot_id") %>%
    full_join(., data$biomass_volume, by = "plot_id") %>%
    full_join(., data$volume_dead_standing, by = "plot_id") %>%
    full_join(., data$disturbance_recent, by = "plot_id") %>%  
    full_join(., data$disturbance_parameters, by = "plot_id") %>%
    full_join(., data$disturbance_index, by = "plot_id") %>%
    full_join(., data$volume_dead_lying_decay, by = "plot_id") %>%
    full_join(., data$volume_dead_lying, by = "plot_id") %>%
    full_join(., data$age_parameters, by = "plot_id") %>%
    full_join(., data$regeneration_htclass, by = "plot_id") %>%
    full_join(., data$regeneration_250_dbh_min, by = "plot_id") %>%
    full_join(., data$regeneration_0_50, by = "plot_id") %>%
    full_join(., data$regeneration_browsing, by = "plot_id") %>%
    full_join(., data$openness, by = "plot_id") %>%
    full_join(., data$temperature, by = "plot_id") %>%
    mutate_at(vars(-dominant_species, -disturbance_index, -openness_mean, 
                   -openness_gini, -age_gini, -dbh_gini_live_60, -dbh_gini_live_100),
              funs(round(., 0))) %>% 
    gather(., parameter, value, aspect_southness:temp_mean_vegetseason) %>%
    filter(!value %in% NA) %>%
    mutate(value = as.character(value))
  
  return(data.all)
  
}