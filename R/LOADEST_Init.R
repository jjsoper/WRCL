rm(list = ls())

loadest.folder <- "C:/Users/Josh Soper/Documents/Master's Thesis/Research/LOADEST - DCR"
setwd(loadest.folder)
source('R/Data_Init.R')

# RETURN LOADEST MODEL ------------------------------------------------------------
return.source <- function(flows, chloride, source, tributary, sep.filter, start.date, end.date) {
  tribs <- names(flows)[-1]
  all.flows <- add.baseflow(flows, areas, tribs, sep.filter = sep.filter, start.date = start.date, end.date = end.date)
  baseflow <- spread(subset(all.flows, Flowtype == 'Baseflow'), 'Site', 'Flow')
  runoff <- spread(subset(all.flows, Flowtype == 'Runoff'), 'Site', 'Flow')
  # set negative runoffs to zero (generally very small numbers)
  runoff[,3:ncol(runoff)][runoff[,3:ncol(runoff)]<0] <- 0
  total <- spread(subset(all.flows, Flowtype == 'Total'), 'Site', 'Flow')
  
  if (source == 'Total') {
    Q <- select(flows, Date, tributary)
    C <- select(chloride, Date, tributary)
  }
  if (source == 'Baseflow') {
    Q <- select(baseflow, Date, tributary)
    C <- bflow.dom(baseflow = baseflow, total = total, chloride = chloride, tributary, threshold, start.date, end.date)
  }
  
  return(list(Q=Q, C=C))
}
generate.data <- function(flows, wq) {
  tmp <- merge(flows,
               wq,
               by = "Date")
  names(tmp)[2:3] <- c("Flow", "Chloride")
  tmp <- tmp[order(tmp$Date),]
  tmp$Load <- tmp$Flow * tmp$Chloride * conversion
  tmp$Data <- "Observed"
  # clean data
  tmp <- na.omit(tmp)
  tmp <- tmp[tmp$Chloride !=0,]
  return(tmp)
}
return.best.model <- function(data) {
  
  best.mod <- selBestModel(constituent = 'Chloride', 
                           data = data, 
                           flow = 'Flow',
                           dates = "Date", 
                           conc.units = 'mg/L', 
                           flow.units = 'cms', 
                           load.units = 'kg')
  return(best.mod)
}

# WORKAROUND FOR "MODEL NUMBER NOT IN RANGE" ERROR
return.load.mod <- function(tmp, no.) {
  if (no. == 1) {load.mod <- loadReg(Chloride~model(1), data = tmp, flow = 'Flow',dates = "Date", conc.units = 'mg/L', flow.units = 'cms', load.units = 'kg')}
  if (no. == 2) {load.mod <- loadReg(Chloride~model(2), data = tmp, flow = 'Flow',dates = "Date", conc.units = 'mg/L', flow.units = 'cms', load.units = 'kg')}
  if (no. == 3) {load.mod <- loadReg(Chloride~model(3), data = tmp, flow = 'Flow',dates = "Date", conc.units = 'mg/L', flow.units = 'cms', load.units = 'kg')}
  if (no. == 4) {load.mod <- loadReg(Chloride~model(4), data = tmp, flow = 'Flow',dates = "Date", conc.units = 'mg/L', flow.units = 'cms', load.units = 'kg')}
  if (no. == 5) {load.mod <- loadReg(Chloride~model(5), data = tmp, flow = 'Flow',dates = "Date", conc.units = 'mg/L', flow.units = 'cms', load.units = 'kg')}
  if (no. == 6) {load.mod <- loadReg(Chloride~model(6), data = tmp, flow = 'Flow',dates = "Date", conc.units = 'mg/L', flow.units = 'cms', load.units = 'kg')}
  if (no. == 7) {load.mod <- loadReg(Chloride~model(7), data = tmp, flow = 'Flow',dates = "Date", conc.units = 'mg/L', flow.units = 'cms', load.units = 'kg')}
  if (no. == 8) {load.mod <- loadReg(Chloride~model(8), data = tmp, flow = 'Flow',dates = "Date", conc.units = 'mg/L', flow.units = 'cms', load.units = 'kg')}
  if (no. == 9) {load.mod <- loadReg(Chloride~model(9), data = tmp, flow = 'Flow',dates = "Date", conc.units = 'mg/L', flow.units = 'cms', load.units = 'kg')}
  if (no. == 10) {load.mod <- loadReg(Chloride~model(10), data = tmp, flow = 'Flow',dates = "Date", conc.units = 'mg/L', flow.units = 'cms', load.units = 'kg')}
  if (no. == 11) {load.mod <- loadReg(Chloride~model(11), data = tmp, flow = 'Flow',dates = "Date", conc.units = 'mg/L', flow.units = 'cms', load.units = 'kg')}
  return(load.mod)
}

load.obs.sim <- function(tributary, load.mod, flow,tmp) {
  pred <- predLoad(load.mod, 
                   data.frame(Date = flow$Date, Flow = flow[[tributary]]), 
                   by = 'day')
  pred <- pred %>% dplyr::rename("Load" = Flux) %>% transform(Data = "Simulated")
  tmp$L95 <- NA
  tmp$U95 <- NA
  dat <- rbind(select(tmp, Date, Flow, Load, L95, U95, Data),
               select(pred, Date, Flow, Load, L95, U95, Data))
  return(dat)
}

return.loadest.models <- function(tributary, flows, chloride, flowsource, start.date, end.date, sep.filter, write.report = F) {
  source <- return.source(flows = flows, chloride = chloride, source = flowsource, tributary = tributary, sep.filter = sep.filter, start.date,end.date)
  tmp <- generate.data(flows = source$Q, wq = source$C)
  load.mod <- return.best.model(data = tmp)
  return(load.mod)
}

# CUSTOM LOADEST PLOTS ------------------------------------------
loadest.plts <- function(dat, load.mod, tributary, sep.filter, start.date, end.date) {
  G2 <- signif(2*(load.mod$lfit$LLR - load.mod$lfit$LLR1), 10)
  pval <- 1 - pchisq(G2, load.mod$lfit$NPAR - 1)
  stats <- data.frame(Rsq = load.mod$lfit$RSQ/100,
                      NSE = loadStats(load.mod, 'load')$outBias[[3]],
                      Bp = loadStats(load.mod, 'load')$outBias[[1]]
  )
  # stats <- stats %>% round(3) %>% t()
  # max.axis <- mean(spread(dat, Data, Load)$Simulated, na.rm = T) + 3*sd(spread(dat, Data, Load)$Simulated, na.rm = T)
  max.axis <- max(dat$Load, na.rm = T)
  
  # timeseries
  plt <- 
    ggplot(dat, aes(x = Date, y = Load, color = Data, linetype = Data, shape = Data)) +
    geom_point(size = 3) +
    geom_line() +
    scale_color_manual(values = c("black", 'brown2')) +
    scale_shape_manual(values = c(1, NA)) +
    scale_linetype_manual(values = c(NA, 1)) +
    scale_y_continuous(limits = c(0, max.axis)) +
    scale_x_date(breaks = 'year', date_labels = '%Y') +
    annotate(label = bquote(R^2 == .(round(stats$Rsq, 3))), geom = 'text', x = as.Date(start.date) + 365, y = 0.95*max.axis) +
    annotate(label = bquote(NSE == .(round(stats$NSE, 3))), geom = 'text', x = as.Date(start.date) + 365, y = 0.85*max.axis) +
    annotate(label = bquote(Bp == .(round(stats$Bp, 2))*'%'), geom = 'text', x = as.Date(start.date) + 365, y = 0.75*max.axis) +
    labs(x ="Date", y = "Chloride Load (kg/day)", 
         subtitle = paste(gsub("[.]", " ",tributary), " (Model No.", load.mod$model.no, ")", ' ', sep.filter, sep = ""),
         color = element_blank(), shape = element_blank(), linetype = element_blank()) +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face = 'bold'),
          panel.grid.major = element_blank(),
          legend.position = 'right',
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(color = 'black', fill = NA),
          strip.background = element_rect(fill="gray95"),
          strip.text = element_text(face = 'bold'),
          panel.spacing.x = unit(6, 'mm'))
  
  # sim vs observed
  pltLComp <- 
    ggplot(spread(select(dat, -c(L95, U95)), Data, Load), aes(x = Observed, y = Simulated)) +
    geom_point(size = 3, shape = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = 'longdash') +
    scale_x_continuous(limits = c(0, max.axis)) +
    scale_y_continuous(limits = c(0, max.axis)) +
    labs(y = "Simulated Chloride Load (kg/day)", x = "Observed Chloride Load (kg/day)") +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face = 'bold'),
          panel.grid.major = element_blank(),
          legend.position = 'bottom',
          panel.grid.minor = element_blank(),
          panel.background = element_rect(color = 'black', fill = NA),
          strip.background = element_rect(fill="gray95"),
          strip.text = element_text(face = 'bold'),
          panel.spacing.x = unit(6, 'mm'))
  
  # NORMALIZED ERRORS PPCC
  residuals.df <- data.frame(
    date = dat$Date[dat$Data == 'Observed'],
    sresid = load.mod$lfit$RESID / sqrt(load.mod$lfit$PARAML[load.mod$lfit$NPAR+1])
  )
  
  Res <- residuals(load.mod$lfit, type="working", suppress.na.action=TRUE)
  if(load.mod$method == "AMLE") {ppcc <- censPPCC.test(as.lcens(Res, censor.codes=load.mod$lfit$CENSFLAG))}
  
  n.samples <- nrow(residuals.df)
  bin.size <- 5 + 3.322 * log10(n.samples)
  scale.x <- ceiling(max(abs(residuals.df$sresid)))
  
  density <-
    ggplot(residuals.df, aes(sresid)) + 
    geom_density(aes(y = ..density..), fill = 'brown2', alpha = 0.7, color = 'brown2') +
    labs(x = "Standardized Residuals", y = "Density") +
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1), size = 1, linetype = 5) +
    scale_x_continuous(limits = c(-scale.x,scale.x), expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face = 'bold'),
          panel.grid.major = element_blank(),
          legend.position = 'bottom',
          panel.grid.minor = element_blank(),
          panel.background = element_rect(color = 'black', fill = NA),
          strip.background = element_rect(fill="gray95"),
          strip.text = element_text(face = 'bold'),
          panel.spacing.x = unit(6, 'mm'))
  
  QQ.plt <-
    ggplot(residuals.df, aes(sample=sresid))+
    stat_qq(shape = 1, size = 3, color = alpha('brown2', alpha = 0.7))+
    geom_abline(slope = 1, intercept = 0) +
    labs(x = "Normal Quantiles", y = "Standardized Residual") +
    scale_x_continuous(limits = c(-3,3), breaks = seq(-3,3,1), expand = c(0,0)) +
    scale_y_continuous(limits = c(-3,3), breaks = seq(-3,3,1), expand = c(0,0)) +
    annotate(label = bquote(italic(PPCC) == .(round(ppcc$statistic, 3))), geom = 'text', x=-1.5, y=2) +
    annotate(label = bquote(SCORR == .(round(load.mod$lfit$SCORR, 3))), geom = 'text', x=-1.5, y=1.5) +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face = 'bold'),
          panel.grid.major = element_blank(),
          legend.position = 'bottom',
          panel.grid.minor = element_blank(),
          panel.background = element_rect(color = 'black', fill = NA),
          strip.background = element_rect(fill="gray95"),
          strip.text = element_text(face = 'bold'),
          panel.spacing.x = unit(6, 'mm'))
  
  loads <- rbind(data.frame(Key = "Simulated",
                            Load = load.mod$lfit$YPRED),
                 data.frame(Key = "Observed",
                            Load = exp(load.mod$lfit$YLCAL)))
  
  boxplt <- ggplot(loads, aes(x = Key, y = Load)) +
    geom_boxplot() +
    scale_y_log10() +
    labs(x =element_blank(), y = "Load (kg/day)", title = element_blank()) +
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face = 'bold'),
          legend.position = 'none',
          panel.background = element_rect(color = 'black', fill = NA),
          strip.background = element_rect(fill="gray95"),
          strip.text = element_text(face = 'bold'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing.x = unit(6, 'mm'))
  
  return(grid.arrange(plt,
                      pltLComp,
                      density,
                      QQ.plt,
                      nrow = 2,
                      widths = rep(1,3),
                      layout_matrix = rbind(c(1,1,1),c(2,3,4))))
  
}

# GENERATE CUSTOM PLOTS AND REPORTS ----------------------------------------------------------
# will generate custom plots and LOADEST reports w/ diagnostic plots

generate.reports <- function(start.date, end.date, filter.list, flowsource.list) {
  
  flows <- return.flows(start.date, end.date)
  tribs <- names(flows)[-1]
  chloride <- cl.edit(conductivity, start.date, end.date)
  chloride <- spread(chloride, "Site", "Chloride")
  
  loadest.models <- list()
  write.report = T
  
  for (sep.filter in filter.list)  {
    for (flowsource in flowsource.list)  {
      for (tributary in tribs) {
        source <- return.source(flows = flows, chloride = chloride, source = flowsource, tributary = tributary, sep.filter = sep.filter,
                                start.date = start.date,end.date = end.date)
        tmp <- generate.data(flows = source$Q, wq = source$C)
        load.mod <- return.best.model(data = tmp)
        # load.mod <- return.load.mod(tmp = tmp, no. = 7)
        
        if (write.report == T) {
          wd <- getwd()
          setwd(file.path(wd,"reports"))
          loadReport(load.mod, paste(tributary, flowsource, sep.filter))
          setwd(wd)
        }
        
        dat <- load.obs.sim(tributary = tributary,load.mod = load.mod, flow = source$Q, tmp = tmp)
        loadest.models[[tributary]] <- loadest.plts(dat = dat, load.mod = load.mod, tributary = tributary, sep.filter = sep.filter, 
                                                    start.date = start.date, end.date = end.date)
        pdf(file.path(loadest.folder, 'figures', paste(flowsource, sep.filter, "Loadest Models (Flow and Time).pdf")), width = 10, height = 7)
        for (plt in loadest.models) {
          grid.draw(plt)
          grid.newpage()
        }
        dev.off()
      }
    }
  }
}
# generate.reports(start.date = '2000-01-01', end.date = '2019-12-31', filter.list = c('hysep-l', 'hysep-f', 'hysep-s'), flowsource.list = c('Total', 'Baseflow'))
# generate.reports(start.date = '2000-01-01', end.date = '2019-12-31', filter.list = c('hysep-l'), flowsource.list = c('Total', 'Baseflow'))
# ERROR METRICS -------------------------------------------------

return.stats <- function(model) {
  stats <- data.frame(Rsq = model$lfit$RSQ/100,
                      NSE = loadStats(model, 'load')$outBias[[3]],
                      Bp = loadStats(model, 'load')$outBias[[1]],
                      SCORR = model$lfit$SCORR,
                      model.no = model$model.no)
  return(stats)
}

model.stats <- function(start.date, end.date, sep.filter, write.clipboard) {
  flows <- return.flows(start.date, end.date)
  tribs <- names(flows)[-1]
  chloride <- cl.edit(conductivity, start.date, end.date)
  chloride <- spread(chloride, "Site", "Chloride")
  
  stats <- list()
  
  for (flowsource in c('Total', 'Baseflow')) {
    loadest.models <- lapply(tribs, return.loadest.models, flows, chloride, flowsource, start.date, end.date, sep.filter)
    stats[[flowsource]] <- lapply(loadest.models, return.stats)
    names(stats[[flowsource]]) <- tribs
    stats[[flowsource]] <- stats[[flowsource]] %>% bind_rows(.id = 'Site') %>% transform(Site = gsub('[.]', ' ', Site))
  }
  
  stats <- bind_rows(stats, .id = 'Model')
  
  if (write.clipboard == T) {
    write.table(stats, 'clipboard', quote = F, sep = '\t', row.names = F)
  } else {write.csv(x = stats, file = file.path('stats', paste(sep.filter,'stats.csv', sep = '_')), quote = F, row.names = F)}
    
}

# model.stats(start.date = '2000-01-01', end.date = '2019-12-31', sep.filter = 'hysep-l', write.clipboard = F)