# function to estimate changepoints
# works for both 1 or multiple variables [vars] 

# -- Note from Rbeast package ---
# If a list of multiple variables is provided for y, the multivariate version 
# of the BEAST algorithm will be invoked to decompose the multiple time series 
# and detect common changepoints altogether. This feature is experimental only 
# and under further development. 

estimate_cp <- function(dataset, sessionid, vars, maxcp = 10, timeinterval = 1,
                        trim = 0){
  
  if(trim == 10){
    # remove first 10 percent of writing process
    dataset2 = dataset %>%
      filter(as.numeric(startClock2) >= (maxend/1000)*.1)
    starttime = first(as.numeric(dataset2$startClock2))
  } else{
    # use full writing process
    dataset2 = dataset
    starttime = 0
  }
  
  if(length(vars) == 1){
    log_segm2 <- dataset2[[vars]]
  } else{
    log_segm2 <- dataset2 %>%
      ungroup() %>%
      dplyr::select(!!vars) 
  }
  
  out = beast(log_segm2, 
              start = starttime, deltat = timeinterval,
              tcp.minmax = c(0, maxcp),
              tseg.min = 10/timeinterval, # minimal 10 secs apart
              season='none',#  'none': trend-only data without seasonality  
              mcmc.seed = 42,
              mcmc.burnin = 1500,
              ci = F,
              print.options = F) 
  
  if (is.na(out$trend$ncp)){
    # no changepoints detected
    output2 <- data.frame(
      cp_nof =  0,
      cp_nof_median = 0)
  } else {
    # at least one changepoint detected
    output2 <- data.frame(
      cp_prob = out$trend$cpPr,
      # Changepoints Credible Interval
      cp_ci_lower = matrix(out$trend$cpCI, ncol = 2)[,1],
      cp_ci_upper = matrix(out$trend$cpCI, ncol = 2)[,2],
      cp_loc = out$trend$cp,
      cp_pos = (out$trend$cp - starttime)/timeinterval,
      cp_change_abrupt = out$trend$cpAbruptChange,
      cp_nof_mean = out$trend$ncp,
      cp_nof_prob = which(out$trend$ncpPr %in% max(out$trend$ncpPr))-1,
      cp_nof_median = out$trend$ncp_median,
      R2 = out$R2,
      RMSE = out$RMSE,
      mlik = out$marg_lik
      ) %>%
      arrange(cp_loc)
    if(length(vars) == 1){
      output2 <- output2  %>%
        # estimated slopes only available if 1 variable included
        arrange(desc(cp_prob)) %>%
        # only keep nof changepoints related to highest probability
        slice(1:max(cp_nof_prob)) %>%
        rowwise() %>%
        mutate(
          slope_before = ifelse(!is.na(cp_pos), 
                                mean(out$trend$slp[0:cp_pos]), NA),
          slope_after  = ifelse(!is.na(cp_pos), 
                                mean(out$trend$slp[cp_pos:length(out$trend$slp)]), 
                                NA),
          slope_after_pr_zero = ifelse(!is.na(cp_pos), 
                                       mean(out$trend$slpSgnZeroPr[cp_pos:length(out$trend$slp)]), 
                                       NA),
          est_Y_before = ifelse(!is.na(cp_pos), 
                                mean(out$trend$Y[0:cp_pos]), NA),
          est_Y_after  = ifelse(!is.na(cp_pos), 
                                mean(out$trend$Y[cp_pos:length(out$trend$slp)]), 
                                NA),
          est_Y = ifelse(!is.na(cp_pos), 
                         mean(out$trend$Y[cp_pos]), 
                         NA)
        ) %>%
        ungroup() 
    }
    output2
  } 
}

write_cp_to_file <- function(vars, maxcp = 10, timeinterval = 1, trim = 0){
  
  trimyes = ifelse(trim > 0, "trim", "")

  models = c("perc_doclength", "dist_process_prod",
            "rel_position", "sum_source_length", 
            "m_pause_length")
  
  input = read_csv(paste0("app/input/log_time_", timeinterval, "sec.csv"),
                   show_col_types = FALSE) 

beast_1 <- input %>%
  group_by(session_id) %>%
  group_modify(~estimate_cp(.x, vars = vars, maxcp = maxcp,
                            timeinterval = timeinterval, trim = trim))

write.csv(beast_1, paste0("app/models5/rbeast_time", timeinterval, trimyes, "_", 
                          paste0(which(models %in% vars), collapse = ""), "_", 
                          maxcp, "cp.csv"), 
          row.names = F)

beast_1

}