library(EcoHydRology) #Lyne and Hollick BFLOW low pass filter a: filter parameter
library(DVstats) #HYSEP and PART
library(FlowScreen) #Eckhardt 2 param a + BFI (a: recession parameter, BFI: maximum long term BFI)
library(gsubfn)
library(zoo)
library(hydroTSM)
library(magrittr)
library(grid)
library(gridExtra)
library(lubridate)
library(foreign)
library(rloadest)
library(readxl)
library(tidyverse)
library(Kendall)

# GENERAL FUNCTIONS -------------------------------------------------------
ReplaceSites <- function(SiteColumn, tribnames, sitenames) {
  for (index in 1:length(tribnames)) {
    SiteColumn <- SiteColumn %>% gsub(sitenames[index], tribnames[index], .)
  }
  return(SiteColumn)
}

rsq <- function (x, y) cor(x, y) ^ 2

# p values
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# Interpolation Function
linear.interp <- function(df, start, end) {
  df <- df[df$Date>=as.Date(start) & df$Date <= as.Date(end),]
  date.seq <- data.frame(Date = seq.Date(as.Date(start), as.Date(end), by = 1))
  df <- merge(df, date.seq, by = 'Date', all = T)
  
  dates <- df$Date
  df$numdate <- as.numeric(df$Date)
  df <- zoo(df)
  df <- data.frame(na.approx(df[,2:ncol(df)], x = df$numdate), stringsAsFactors = FALSE)
  df <- mutate_all(df, as.numeric)
  df$Date <- dates
  df <- select(df, -'numdate')
  df <- select(df, "Date", everything())
  return(df)
}


areas <- read.csv("data/areas.txt")
areas <- areas %>% select("UMASS_WSHD", "AREA") %>% spread(key = UMASS_WSHD, value = AREA)

conversion <- 1000 / 1e6 * 86400 # conversion mg m3 / L s <- kg/day
threshold <- 0.90 #baseflow separation threshold (i.e. 90%)

# CONDUCTIVITY/CHLORIDE ---------------------------------------------------
trib.names <- c("French.Brook", "Malagasco.Brook", "Malden.Brook", "Muddy.Brook", 
                "Waushacum.Brook", "West.Boylston.Brook", "Gates.Brook", "Quinapoxet.River", "Stillwater.River", "Trout.Brook")
site.names <- c(".*MD01.*", '.*MD02.*', '.*MD06.*', '.*MD03.*', '.*MD83.*', '.*MD05.*', '.*MD04.*', '.*MD69.*', '.*MD07.*', ".*M110.*")

chloride.cond <- readRDS("data/ChlorideCond.RDS")
reg.cl <- lm(chloride.cond$Chloride~chloride.cond$`Specific Conductivity`)
cond.to.chloride <- function(conductivity, reg) {
  chloride <- conductivity * reg$coefficients[2] + reg$coefficients[1]
  return(chloride)
}

conductivity <- read_excel("data/Historic Tributary Conductivity.xlsx")
conductivity <- conductivity %>% select(SampleDateTime, Site, Result) %>% mutate(SampleDateTime = as.Date(SampleDateTime))
conductivity$Site <- ReplaceSites(conductivity$Site, trib.names, site.names)
names(conductivity) <- c("Date", "Site", "Conductivity")

# generate chloride data
cl.edit <- function(conductivity, start.date, end.date) {
  chloride <- conductivity
  names(chloride) <- c("Date", "Site", "Chloride")
  chloride$Chloride <- cond.to.chloride(chloride$Chloride, reg.cl)
  chloride$Chloride[chloride$Chloride<0] <- 0
  chloride <- spread(chloride, Site, Chloride)
  date.sequence <- data.frame(Date = seq.Date(as.Date(start.date), as.Date(end.date), by = 1))
  chloride <- merge(chloride, date.sequence, "Date", all.y = T)
  chloride <- gather(chloride, "Site", "Chloride", 2:ncol(chloride))
  
  outliers <- aggregate(Chloride~Site, chloride, function(x){mean(x)*10})
  chloride <- chloride %>% merge(outliers, by = 'Site', all = T) %>% subset(Chloride.x < Chloride.y) %>% 
    select(Date,Site,Chloride.x) %>%
    dplyr::rename('Chloride' = Chloride.x)
  
  chloride <- chloride %>% spread(Site, Chloride) %>% 
    merge(date.sequence, "Date", all.y = T) %>%
    gather("Site", "Chloride", 2:11)
  
  
  return(chloride)
}

# FLOWS -------------------------------------------------------------------

runoff.model <- function(trib.area, stillwater.flow, stillwater.area) {
  trib.flow = stillwater.flow/stillwater.area * trib.area
  return(trib.flow)
}

return.flows <- function(start.date, end.date) {
  flows <- readRDS("data/USGSFlows.RDS")
  flows <- spread(flows, key = Site, value = Flow)
  # compute flows for rest of timeseries 
  trib.names <- c("French Brook", "Malagasco Brook", "Malden Brook", "Muddy Brook", "Waushacum Brook", "West Boylston Brook", "Trout Brook")
  for (trib in trib.names) {
    flows[trib] <- runoff.model(areas[[trib]], flows$`Stillwater River`, areas$`Stillwater River`)
  }
  
  # add non-existing gates flow data
  gates.dates <- flows[is.na(flows$`Gates Brook`),]$Date
  flows[flows$Date %in% gates.dates,]$`Gates Brook` <-
    runoff.model(areas[["Gates Brook"]], subset(flows, Date %in% gates.dates)$`Stillwater River`, areas$`Stillwater River`)
  flows <- linear.interp(flows, start.date, end.date)
  return(flows)
}

# BASEFLOW CONDITIONS
add.baseflow <- function(flows, areas, tribs, sep.filter, start.date, end.date) {
  baseflow <- flows
  for (trib in tribs) {
    if (sep.filter == 'hysep-l') {
      BaseQ <- hysep(Flow = flows[[trib]], Dates = flows$Date, Start = start.date, End = end.date, da = areas[[trib]]* 0.386, select = 'l')$BaseQ
    } else if (sep.filter == 'hysep-f') {
      BaseQ <- hysep(Flow = flows[[trib]], Dates = flows$Date, Start = start.date, End = end.date, da = areas[[trib]]* 0.386, select = 'f')$BaseQ
    } else if (sep.filter == 'hysep-s') {
      BaseQ <- hysep(Flow = flows[[trib]], Dates = flows$Date, Start = start.date, End = end.date, da = areas[[trib]]* 0.386, select = 's')$BaseQ
    } else if (sep.filter == 'part') {
      BaseQ <- part(Flow = flows[[trib]], Dates = flows$Date, Start = start.date, End = end.date, da = areas[[trib]]* 0.386)$BaseQ
    } else if (sep.filter == 'L&H3') {
      BaseQ <- BaseflowSeparation(streamflow = flows[[trib]], filter_parameter = 0.925, passes = 3)$bt
    } else if (sep.filter == 'L&H2') {
      BaseQ <- BaseflowSeparation(streamflow = flows[[trib]], filter_parameter = 0.925, passes = 2)$bt
    } else if (sep.filter == 'L&H1') {
      BaseQ <- BaseflowSeparation(streamflow = flows[[trib]], filter_parameter = 0.925, passes = 1)$bt
    } else if (sep.filter == 'eck') {
      BaseQ <- bf_eckhardt(discharge = flows[[trib]], a = 0.925, BFI = 0.8)
    } else {stop('Improper filter selection')}
    
    baseflow[[trib]] <- BaseQ
    # remove rounding errors when bflow > flow
    baseflow[[trib]][baseflow[[trib]] >= flows[[trib]]] <- flows[[trib]][baseflow[[trib]] >= flows[[trib]]]
  }
  
  runoff <- flows
  runoff[,2:11] <- flows[,2:11] -  baseflow[,2:11]
  tmp <- list("Total" = flows, "Baseflow" = baseflow, "Runoff" = runoff)
  tmp <- bind_rows(tmp, .id = 'Flowtype')
  tmp <-  gather( tmp, "Site", "Flow", 3:ncol(tmp))
  return(tmp)
}

bflow.dom <- function(baseflow, total, chloride, tributary, threshold, start.date, end.date) {
  t.s <- seq.Date(as.Date(start.date), as.Date(end.date), by = 1)
  dates <- total$Date[baseflow[[tributary]] >= total[[tributary]] * threshold]
  
  baseflow.cl <- merge(data.frame(Date=t.s), select(chloride[chloride$Date %in% dates,], Date, tributary), all = T, by  = "Date")
  baseflow.cl[[tributary]][1] <- baseflow.cl[[tributary]][!(is.na(baseflow.cl[[tributary]]))][1]
  baseflow.cl[[tributary]][length(baseflow.cl[[tributary]])] <- baseflow.cl[[tributary]][!(is.na(baseflow.cl[[tributary]]))][length(baseflow.cl[[tributary]][!(is.na(baseflow.cl[[tributary]]))])]
  return(baseflow.cl)
}

separate.cl.interp <- function(chloride, tributary, threshold, start.date, end.date, interpolate = T, remove.dominant.days = F, remove.negatives = T) {
  dates <- total$Date[baseflow[[tributary]] >= total[[tributary]] * threshold]
  
  baseflow.cl <- merge(data.frame(Date=t.s), select(chloride[chloride$Date %in% dates,], Date, tributary), all = T, by  = "Date")
  baseflow.cl[[tributary]][1] <- baseflow.cl[[tributary]][!(is.na(baseflow.cl[[tributary]]))][1]
  baseflow.cl[[tributary]][length(baseflow.cl[[tributary]])] <- baseflow.cl[[tributary]][!(is.na(baseflow.cl[[tributary]]))][length(baseflow.cl[[tributary]][!(is.na(baseflow.cl[[tributary]]))])]
  baseflow.cl <- linear.interp(baseflow.cl, start.date, end.date)
  
  runoff.cl <- baseflow.cl
  runoff.cl[[tributary]] <- ((total[[tributary]] * chloride[[tributary]]) - (baseflow[[tributary]] * baseflow.cl[[tributary]])) / runoff[[tributary]]
  
  runoff.cl[[tributary]][runoff.cl[[tributary]] == -Inf] <- 0
  
  if (remove.dominant.days == T) {
    runoff.cl <- runoff.cl[!(runoff.cl$Date %in% dates),]
  }
  
  if (remove.negatives == T) {
    runoff.cl <- runoff.cl %>% subset(get(tributary) >= 0) %>%
      merge(data.frame(Date = seq.Date(as.Date(start.date), as.Date(end.date), 1)), all = T)
  }
  
  if (interpolate == T) {
    runoff.cl <- runoff.cl %>% linear.interp(start.date, end.date)
  }
  
  
  return(list(runoff.cl = runoff.cl, baseflow.cl = baseflow.cl, Dates = dates))
}

