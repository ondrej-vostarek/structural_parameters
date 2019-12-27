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

volume <- function(SM, density_gCm3){
  #'@description calculate tree volume in m3 from stem mass and wood density
  #'@param SM stem mass in kg
  #'@param density_gCm3 wood density in g/cm3
  
  (SM * 0.001) / density_gCm3
  
}

get_data <- function(plot.id){
  #'@description collect necessary data from database and return them as data.list
  #'@param plot.id id of plots as in 'plot' table in database
  
  plot <- tbl(KELuser, "plot") %>% 
    filter(id %in% plot.id) %>% 
    select(plot_id = id, date, plotid, census, country, location, lng, lat, plotsize, foresttype, altitude_m, aspect)
  
  tree <- tbl(KELuser, "tree") %>% 
    select(tree_id = id, plot_id, treeid, onplot, treetype, status, growth, species, dbh_mm, height_m, decayht, decay, decay_wood) %>% 
    inner_join(., plot, by = 'plot_id')
  
  deadwood <- tbl(KELuser, "deadwood") %>% 
    inner_join(., plot, by = 'plot_id') %>%
    select(plot_id, species, dbh_mm, decay)

  deadwood_tree <- tbl(KELuser, "deadwood_tree") %>%
    inner_join(., plot, by = "plot_id") %>%
    filter(volume_m3 > 0,
           !is.na(plotsize)) %>%
    select(plot_id, plotsize, species, decay, volume_m3)

  core <- tbl(KELuser, "core") %>% 
    inner_join(., tree, by = 'tree_id') %>%
    filter(treetype %in% "0" & onplot %in% c(1, 2) | treetype %in% c("m", "x"),    
           growth %in% 1,
           !dbh_mm < 100,
           !corestatus %in% c(2, 3)) %>% 
    inner_join(., tbl(KELuser, "ring"), by = c("id" = "core_id")) %>%
    select(date, plotid, treeid, subcore, missing_years) 

  regeneration <- tbl(KELuser, "regeneration") %>% 
    inner_join(., plot, by = "plot_id") %>%
    filter(!is.na(plotsize)) %>%
    select(plot_id, foresttype, plotsize, htclass, count)
  
  regeneration_subplot <- tbl(KELuser, "regeneration_subplot") %>% 
    inner_join(., plot, by = "plot_id") %>%
    select(plot_id, htclass, browsing, count)

  recent_dist <- tree %>% 
    filter(!onplot %in% c(0, 99), 
           !growth %in% 0,
           decay %in% c(-1:3)) %>%              
    inner_join(., 
               tbl(KELuser, 'species_fk') %>% select(species = id, sp_group_dist), 
               by = "species") %>%
    inner_join(.,
               tbl(KELuser, 'dist_group'), 
               by = c("country", "foresttype", "location", "sp_group_dist")) %>% 
    inner_join(.,
               tbl(KELuser, 'dist_param') %>% select(dist_param = id, Tdbh_mm = dbh_mm, dbh_ca_f), 
               by = "dist_param") %>% 
    filter(dbh_mm >= Tdbh_mm) %>%
    select(plot_id, dbh_mm, decay, dbh_ca_f)
  
  dist_plot_event <- tbl(KELuser, "dist_plot_event") %>% 
    inner_join(., tbl(KELuser, "dist_plot"), by = c("dist_plot_id" = "id")) %>%
    inner_join(., plot, by = "plot_id") %>%
    select(id, plot_id, year, ca_per, kde)
  
  canopy_analysis <- tbl(KELuser, "canopy_analysis") %>% 
    inner_join(., plot, by = "plot_id") %>%
    filter(parameter %in% "openness") %>%
    select(plot_id, value)
 
  wood_density <- tbl(KELuser, "wood_density")
  
  biomass_eq <- tbl(KELuser, "biomass_eq")
    
  data <- list()
  
  data$plot <- collect(plot)
  data$tree <- collect(tree)
  data$deadwood <- collect(deadwood)
  data$deadwood_tree <- collect(deadwood_tree)
  data$core <- collect(core)
  data$regeneration <- collect(regeneration)
  data$regeneration_subplot <- collect(regeneration_subplot)
  data$recent_dist <- collect(recent_dist)
  data$dist_plot_event <- collect(dist_plot_event)
  data$canopy_analysis <- collect(canopy_analysis)
  data$wood_density <- collect(wood_density)
  data$biomass_eq <- collect(biomass_eq)
  
  return(data) 
  
}

calculate_parameters <- function(data, dataType){
  #'@description calculate the structural parameters
  #'@param data data.list as produced by the get_data function
  #'@param dataType reflects the type of data used for calculation of specific parameters
  
  parameters <- list()
  
  for (i in dataType) {
    
    if (i == "plot"){
      
      # aspect southness --------------------------------------------------------
      
      parameters$aspect_southness <- data$plot %>% filter(!aspect %in% NA) %>% group_by(plot_id) %>%
        summarise(aspect_southness = round(cos(deg2rad(45) - deg2rad(aspect)) + 1, 0))
      
      # nearest plot ------------------------------------------------------------
      
      plot.new <- data$plot %>% filter(census %in% 1)
      
      P <- plot.new %>% select(lat, lng) %>%  collect() %>% as.matrix()
      p <- distm(P, P)
      diag(p) <- NA
      
      plot.key <- plot.new %>% mutate(id = row_number()) %>% select(id, plotid) %>% deframe()
      
      plot.nearest <- reshape2::melt(as.matrix(p), varnames = c("plotid", "second_plot"), value.name = "nearest_plot") %>%
        filter(!nearest_plot %in% c("NaN", NA)) %>%
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
        select(plot_id, nearest_plot) %>%
        mutate(nearest_plot = round(nearest_plot, 0))
      
    } else {
      
      if(i == "tree") {
        
        # tree parameters ---------------------------------------------------------
        
        parameters$tree_parameters <- data$tree %>%
          filter(!onplot %in% c(0, 99), 
                 !is.na(dbh_mm)) %>%
          mutate(ba = pi * dbh_mm ^ 2 / 4 / 1000000) %>%
          group_by(plot_id) %>%
          summarise(plotsize = first(plotsize),
                    n_trees_live_500 = length(treeid[dbh_mm >= 500 & status %in% c(1:4)]),
                    n_trees_live_700 = length(treeid[dbh_mm >= 700 & status %in% c(1:4)]),
                    n_trees_dead_500 = length(treeid[dbh_mm >= 500 & status %in% c(11:23)]),
                    n_trees_dead_700 = length(treeid[dbh_mm >= 700 & status %in% c(11:23)]),
                    n_trees_dead_60 = length(treeid[dbh_mm >= 60 & status %in% c(11:23)]),
                    dbh_quadrmean_dead_60 = sqrt(mean(dbh_mm[dbh_mm >= 60 & status %in% c(11:23)]^2)),
                    dbh_quadrmean_dead_60 = ifelse(dbh_quadrmean_dead_60 %in% "NaN", NA, dbh_quadrmean_dead_60),
                    ba_dead_60 = sum(ba[dbh_mm >= 60 & status %in% c(11:23)]),
                    n_trees_dead_100 = length(treeid[dbh_mm >= 100 & status %in% c(11:23)]),
                    dbh_quadrmean_dead_100 = sqrt(mean(dbh_mm[dbh_mm >= 100 & status %in% c(11:23)]^2)),
                    dbh_quadrmean_dead_100 = ifelse(dbh_quadrmean_dead_100 %in% "NaN", NA, dbh_quadrmean_dead_100),
                    ba_dead_100 = sum(ba[dbh_mm >= 100 & status %in% c(11:23)]),
                    n_trees_live_60 = length(treeid[dbh_mm >= 60 & status %in% c(1:4)]),
                    dbh_quadrmean_live_60 = sqrt(mean(dbh_mm[dbh_mm >= 60 & status %in% c(1:4)]^2)),
                    dbh_quadrmean_live_60 = ifelse(dbh_quadrmean_live_60 %in% "NaN", NA, dbh_quadrmean_live_60),
                    dbh_mean_live_60 = mean(dbh_mm[dbh_mm >= 60 & status %in% c(1:4)]),
                    dbh_mean_live_60 = ifelse(dbh_mean_live_60 %in% "NaN", NA, dbh_mean_live_60),
                    dbh_gini_live_60 = ineq(dbh_mm[dbh_mm >= 60 & status %in% c(1:4)], type = "Gini"),
                    dbh_gini_live_60 = ifelse(dbh_gini_live_60 %in% "NaN", NA, dbh_gini_live_60),
                    ba_live_60 = sum(ba[dbh_mm >= 60 & status %in% c(1:4)]),
                    n_trees_live_100 = length(treeid[dbh_mm >= 100 & status %in% c(1:4)]),
                    dbh_quadrmean_live_100 = sqrt(mean(dbh_mm[dbh_mm >= 100 & status %in% c(1:4)]^2)),
                    dbh_quadrmean_live_100 = ifelse(dbh_quadrmean_live_100 %in% "NaN", NA, dbh_quadrmean_live_100),
                    dbh_mean_live_100 = mean(dbh_mm[dbh_mm >= 100 & status %in% c(1:4)]),
                    dbh_mean_live_100 = ifelse(dbh_mean_live_100 %in% "NaN", NA, dbh_mean_live_100),
                    dbh_gini_live_100 = ineq(dbh_mm[dbh_mm >= 100 & status %in% c(1:4)], type = "Gini"),
                    dbh_gini_live_100 = ifelse(dbh_gini_live_100 %in% "NaN", NA, dbh_gini_live_100),
                    ba_live_100 = sum(ba[dbh_mm >= 100 & status %in% c(1:4)])) %>%
          mutate(ba_dead_60 = ifelse(ba_dead_60 %in% 0 , NA, ba_dead_60),
                 ba_dead_100 = ifelse(ba_dead_100 %in% 0, NA, ba_dead_100),
                 ba_live_60 = ifelse(ba_live_60 %in% 0, NA, ba_live_60),
                 ba_live_100 = ifelse(ba_live_100 %in% 0, NA, ba_live_100)) %>%
          mutate_at(vars(n_trees_live_500, n_trees_live_700, n_trees_dead_500, n_trees_dead_700, n_trees_dead_60,
                         ba_dead_60, n_trees_dead_100, ba_dead_100, n_trees_live_60, ba_live_60, n_trees_live_100, ba_live_100), 
                    funs(.*10000/plotsize)) %>%
          mutate_at(vars(dbh_gini_live_60, dbh_gini_live_100), funs(round(., 2))) %>%
          mutate_at(vars(-dbh_gini_live_60, -dbh_gini_live_100), funs(round(., 0))) %>%
          select(-plotsize)
        
        parameters$height_max <- data$tree %>% 
          group_by(plot_id) %>% 
          summarise(height_max = max(height_m, na.rm = T)) %>%
          mutate(height_max = ifelse(height_max %in% -Inf, NA, height_max),
                 height_max = round(height_max, 0))
        
        parameters$dominant_species <- data$tree %>%
          filter(!onplot %in% c(0, 99),
                 !is.na(dbh_mm),
                 status %in% c(1:4)) %>%
          mutate(ba = pi * dbh_mm ^ 2 / 4 / 1000000) %>%
          group_by(plot_id, species) %>%
          summarise(ba = sum(ba) * 10000 / min(plotsize)) %>%
          group_by(plot_id) %>%
          arrange(desc(ba)) %>%
          filter(row_number() == 1) %>%
          select(plot_id, dominant_species = species)
        
        # biomass and volume of alive trees ---------------------------------------
        
        parameters$biomass_volume <- data$tree %>%
          filter(!onplot %in% c(0, 99),
                 status %in% c(1:4),
                 !is.na(dbh_mm),
                 !species %in% "99") %>%
          left_join(., data$wood_density %>% distinct(., species, density_gCm3), by = "species") %>%
          left_join(., data$biomass_eq, by = "species") %>%
          left_join(., parameters$tree_parameters %>% select(plot_id, ba_live_60), by = "plot_id") %>%
          rowwise() %>%
          mutate(TB = eval(parse(text = branches_mass_f)),
                 FM = eval(parse(text = foliage_mass_f)),
                 RM = eval(parse(text = root_mass_f)),
                 SM = eval(parse(text = stem_mass_f)),
                 volume_live = volume(SM, density_gCm3),
                 biomass_aboveground = TB + FM + SM,
                 biomass_underground = RM) %>%
          group_by(plot_id) %>%
          summarise(plotsize = first(plotsize),
                    volume_live_60 = sum(volume_live[dbh_mm >= 60]),
                    volume_live_100 = sum(volume_live[dbh_mm >= 100]),
                    biomass_aboveground_60 = sum(biomass_aboveground[dbh_mm >= 60]),
                    biomass_aboveground_100 = sum(biomass_aboveground[dbh_mm >= 100]),
                    biomass_underground_60 = sum(biomass_underground[dbh_mm >= 60]),
                    biomass_underground_100 = sum(biomass_underground[dbh_mm >= 100])) %>%
          mutate_at(vars(-plot_id, -plotsize), funs(.*10000/plotsize)) %>%
          mutate_at(vars(-plot_id, -plotsize), funs(round(., 0))) %>%
          select(-plotsize)
        
        # volume and biomass of dead standing trees -------------------------------
        
        parameters$volume_dead_standing <- data$tree %>%
          filter(!onplot %in% c(0, 99),
                 status %in% c(11:23),
                 !decayht %in% 99,
                 !dbh_mm %in% NA) %>%
          mutate(decayht = case_when(
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
                    volume_dead_standing_100 = sum(volume_snag[dbh_mm >= 100]) * 10000 / min(plotsize)) %>%
          mutate_at(vars(-plot_id), funs(round(., 0)))
        
        parameters$biomass_dead_standing <- data$tree %>%
          filter(!onplot %in% c(0, 99),
                 status %in% c(11:23),
                 !decayht %in% 99,
                 !dbh_mm %in% NA,
                 !species %in% "99",
                 decay %in% c(1:5) | decay_wood %in% c(1:5)) %>%
          mutate(decayht = case_when(
            decayht == 0 ~ 5,
            decayht == 1 ~ 15,
            decayht == 2 ~ 25,
            decayht == 3 ~ 35,
            decayht == 4 ~ 45,
            decayht == 5 ~ 55),
            decay_class = ifelse(decay_wood %in% 99, round(0.2924953 + 0.7131269 * decay, 0), decay_wood)) %>%
          left_join(., data$wood_density, by = c("species", "decay_class")) %>%
          rowwise() %>%
          mutate(volume_snag = E_VOL_AB_HmDm_HT.f(Hm=1.3, Dm=(dbh_mm * 0.1), 
                                                  mHt = (log(dbh_mm * 0.1)-1.08261)^2/0.275541, 
                                                  sHt = 0, par.lme = SK.par.lme, A=0, B=decayht, iDH = "H")$E_VOL,
                 biomass = volume_snag * (density_gCm3 * relative_density * 1000)) %>%
          group_by(plot_id) %>%
          summarise(biomass_dead_standing_60 = sum(biomass[dbh_mm >= 60]) * 10000 / min(plotsize),
                    biomass_dead_standing_100 = sum(biomass[dbh_mm >= 100]) * 10000 / min(plotsize)) %>%
          mutate_at(vars(-plot_id), funs(round(., 0)))
        
        # recent disturbance ------------------------------------------------------
        
        parameters$disturbance_recent <- data$recent_dist %>% 
          rowwise() %>% 
          mutate(ca = eval(parse(text = dbh_ca_f))) %>%
          group_by(plot_id) %>%
          summarise(disTot = sum(ca[decay %in% c(-1:3)]),
                    rDis1_3N = sum(ca[decay %in% c(1:3)])) %>%
          mutate(disturbance_recent = rDis1_3N/disTot * 100,
                 disturbance_recent = ifelse(disturbance_recent %in% "NaN" | disturbance_recent %in% 0, NA, disturbance_recent),
                 disturbance_recent = round(disturbance_recent, 0)) %>%
          select(plot_id, disturbance_recent)
        
      } else {
        
        if(i == "disturbance"){
          
          # disturbance parameters --------------------------------------------------
          
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
            select(-year, -id, -kde, -ms30, -ca_per) %>%
            filter(row_number() == 1) %>%
            mutate_at(vars(-plot_id), funs(round(., 0)))
          
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
          
        } else {
          
          if(i == "deadwood"){
            
            # volume and biomass of lying deadwood - transects ------------------------
            
            parameters$volume_dead_lying_decay <- data$deadwood %>%
              filter(!decay %in% 99) %>%
              group_by(plot_id, decay) %>%
              summarise(volume_dead_lying_decay = round(((pi ^ 2 * sum((dbh_mm * 0.001) ^ 2)) / 800) * 10000, 0)) %>%
              mutate(decay = paste0("volume_dead_lying_decay", decay)) %>%
              spread(decay, volume_dead_lying_decay, fill = 0) %>%
            
            parameters$volume_dead_lying <- data$deadwood %>%
              group_by(plot_id) %>%
              summarise(volume_dead_lying = round(((pi ^ 2 * sum((dbh_mm * 0.001) ^ 2)) / 800) * 10000, 0))
            
            parameters$biomass_dead_lying <- data$deadwood %>%
              filter(!decay %in% 99 & !species %in% "99") %>%
              left_join(., data$wood_density, by = c("species", "decay" = "decay_class")) %>%
              group_by(plot_id, decay, species) %>%
              summarise(volume = ((pi ^ 2 * sum((dbh_mm * 0.001) ^ 2)) / 800) * 10000,
                        biomass = volume * (first(density_gCm3) * first(relative_density) * 1000)) %>%
              group_by(plot_id) %>%
              summarise(biomass_dead_lying = round(sum(biomass), 0))
            
          } else {
            
            if(i == "deadwood_tree"){
              
              # volume and biomass of lying deadwood - biodiversity ---------------------
              
              parameters$volume_dead_tree_lying_decay <- data$deadwood_tree %>%
                filter(!decay %in% 99) %>%
                group_by(plot_id, decay) %>%
                summarise(volume_dead_tree_lying_decay = round(sum(volume_m3) * 10000 / min(plotsize), 0)) %>%
                mutate(decay = paste0("volume_dead_tree_lying_decay", decay)) %>%
                spread(decay, volume_dead_tree_lying_decay, fill = 0)
              
              parameters$volume_dead_tree_lying <- data$deadwood_tree %>%
                group_by(plot_id) %>%
                summarise(volume_dead_tree_lying = round(sum(volume_m3) * 10000 / min(plotsize), 0))
              
              parameters$biomass_dead_tree_lying <- data$deadwood_tree %>%
                filter(!decay %in% 99 & !species %in% "99") %>%
                left_join(., data$wood_density, by = c("species", "decay" = "decay_class")) %>%
                mutate(biomass = volume_m3 * (density_gCm3 * relative_density * 1000)) %>%
                group_by(plot_id) %>%
                summarise(biomass_dead_tree_lying = round(sum(biomass) * 10000 / min(plotsize), 0))
              
            } else {
              
              if(i == "core"){
                
                # age parameters ----------------------------------------------------------
                
                parameters$age_parameters <- data$core %>% 
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
                            age_iqr = interquantilerange(age)) %>%
                  mutate_at(vars(-plot_id, -age_gini), funs(round(., 0)))
                
              } else {
                
                if(i == "regeneration"){
                  
                  # regeneration ------------------------------------------------------------
                  
                  parameters$regeneration_htclass <- data$regeneration %>%
                    filter(!htclass %in% c(0, 3, 99)) %>%
                    mutate(plotsize = ifelse(plotsize > 1000, 1000, plotsize)) %>%
                    group_by(plot_id, htclass) %>%
                    summarise(regeneration_htclass = round(sum(count) * 10000 / min(plotsize), 0)) %>%
                    mutate(htclass = ifelse(htclass %in% 1, "regeneration_50_130", "regeneration_130_250")) %>%
                    spread(htclass, regeneration_htclass, fill = 0)
                  
                  parameters$regeneration_250_dbh_min <- data$regeneration %>%
                    filter(htclass %in% 3) %>%
                    mutate(plotsize = ifelse(plotsize > 1000, 1000, plotsize),
                           foresttype = ifelse(foresttype %in% "spruce", "regeneration_250_100", "regeneration_250_60")) %>%
                    group_by(plot_id, foresttype) %>%
                    summarise(regeneration_250_dbh_min = round(sum(count) * 10000 / min(plotsize), 0)) %>%
                    spread(foresttype, regeneration_250_dbh_min)
                  
                  
                } else {
                  
                  if(i == "regeneration_subplot"){
                    
                    # regeneration subplots ---------------------------------------------------
                    
                    parameters$regeneration_0_50 <- data$regeneration_subplot %>%
                      filter(htclass %in% 0) %>%
                      group_by(plot_id) %>%
                      summarise(regeneration_0_50 = round(sum(count) * 10000 / 20, 0))
                    
                    parameters$regeneration_browsing <- data$regeneration_subplot %>%
                      filter(htclass %in% 1,
                             !browsing %in% c(5, 99)) %>%
                      group_by(plot_id) %>%
                      summarise(regeneration_browsing = round(mean(browsing), 0))
                    
                    
                  } else {
                    
                    if(i == "canopy"){
                      
                      # canopy openness ---------------------------------------------------------
                      
                      parameters$openness <- data$canopy_analysis %>%
                        group_by(plot_id) %>%
                        summarise(openness_mean = mean(value),
                                  openness_gini = ineq(value, type = "Gini")) %>%
                        mutate_at(vars(openness_mean, openness_gini), funs(round(., 2)))
                      
                      
                    } else {
                      
                      if(i == "temperature"){
                        
                        # temperature -------------------------------------------------------------
                        
                        load("C:/Users/Ondrej_Vostarek/Desktop/MVP/DB/data/temperatures/model_annual.rda")
                        load("C:/Users/Ondrej_Vostarek/Desktop/MVP/DB/data/temperatures/model_vegetation.rda")
                        
                        plots <- data.list$plot %>% 
                          select(plot_id, longitude = lng, latitude = lat, altitude = altitude_m, country, location)
                        
                        plots$temperature_annual <- stats::predict(model_annual,plots)^2-100
                        plots$temperature_vegetation <- stats::predict(model_vegetation,plots)^2-100
                        
                        parameters$temperature <- plots %>% 
                          select(plot_id, temp_mean_year = temperature_annual, temp_mean_vegetseason = temperature_vegetation) %>%
                          mutate_at(vars(-plot_id), funs(round(., 0)))
                        
                        
                      } else {
                        
                        stop("Unknown dataType, possible values are: 'plot', 'tree', 'disturbance', 'deadwood', 'deadwood_tree', 'core', 'regeneration', 'regeneration_subplot', 'canopy', 'temperature'.")
                        
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  return(parameters)
  
}

collect_data <- function(data){
  #'@description collect all data.list sheets into one data.frame
  #'@param data data.list as produced by the calculate_parameters function
  
  data.all <- data.frame(plot_id = NA)
  
  for(i in names(data)){
    
    data.all <- full_join(data.all, data[[i]], by = "plot_id")
    
  }
  
  data.all <- data.all %>%
    gather(., parameter, value, 2:ncol(.)) %>%
    filter(!value %in% NA) %>%
    mutate(value = as.character(value))
  
  return(data.all)
  
}
