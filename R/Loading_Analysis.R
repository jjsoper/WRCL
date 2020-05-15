rm(list = ls())

loadest.folder <- "C:/Users/Josh Soper/Documents/Master's Thesis/Research/LOADEST - DCR"
setwd(loadest.folder)
source('R/Data_Init.R')
source('R/LOADEST_Init.R')

# outputs loading dataframes aggregated by month and flow path
load.setup <- function(start.date, end.date, sep.filter = 'hysep-l') {
  loads <- sapply(c("Total", "Baseflow", "Runoff"),function(x) NULL)
  flows <- return.flows(start.date, end.date)
  tribs <- names(flows)[-1]
  chloride <- cl.edit(conductivity, start.date, end.date)
  chloride <- spread(chloride, "Site", "Chloride")
  
  tmp <- merge(gather(flows, 'Site', 'Flow', 2:ncol(flows)), gather(chloride, "Site", 'Chloride', 2:ncol(chloride)), by = c('Date', 'Site'))
  
  
  for(flowsource in c("Total", "Baseflow")) {
    for (tributary in tribs[-10]) {
      source <- return.source(flows = flows, chloride = chloride, source = flowsource, tributary = tributary, sep.filter = sep.filter,
                              start.date = start.date,end.date = end.date)
      tmp <- generate.data(flows = source$Q, wq = source$C)
      load.mod <- return.best.model(data = tmp)
      dat <- load.obs.sim(tributary = tributary,load.mod = load.mod, flow = source$Q, tmp = tmp)
      loads[[flowsource]][[tributary]] <- subset(dat, Data == 'Simulated')
      loads[[flowsource]][[tributary]] <- loads[[flowsource]][[tributary]] %>% transform(Month = month(Date)) %>% 
        transform(Year = year(Date)) %>%
        transform(Flow = Flow * 86400) %>%
        select(-c(Date,Data)) %>% 
        aggregate(.~Month+Year, ., sum)
    }
  }
  
  for (tributary in tribs) {
    loads[["Runoff"]][[tributary]] <- loads$Total[[tributary]]
    loads[["Runoff"]][[tributary]]$Load <- loads$Total[[tributary]]$Load - loads$Baseflow[[tributary]]$Load
    loads[["Runoff"]][[tributary]]$Flow <- loads$Total[[tributary]]$Flow - loads$Baseflow[[tributary]]$Flow
  }
  
  loads <- lapply(loads, bind_rows, .id = "Tributary")
  loads <- loads %>% bind_rows(.id = "Source")%>% transform(Tributary = gsub('[.]', ' ', Tributary))
  return(loads)
}

start.date='2000-01-01'
end.date = '2019-12-31'
sep.filter <- 'hysep-l'
saveRDS(load.setup(start.date,end.date,sep.filter), file = "data/loads.RDS")
loads <- readRDS("data/loads.RDS")
head(loads)