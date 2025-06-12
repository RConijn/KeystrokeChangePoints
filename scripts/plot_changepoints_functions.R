# function to plot breakpoints for time input data
plot_time_breakpoints <- function(dataset, model, maxcp = 10, timeinterval = 1,
                             trim = 0, subsetsessions, min_cpprob = 0){

  
  trimyes = ifelse(trim > 0, "trim", "")

  dataplot = read_csv(paste0("models5/rbeast_time", timeinterval, trimyes,
                             "_", model, "_", maxcp, "cp.csv"),
                      show_col_types = FALSE)
  
  dataset %>%
    filter(session_id %in% subsetsessions) %>%
    group_by(participant, session, session_id) %>%
    mutate(
      rev_start = ifelse(final_revision_present == 1, 
                         startClock2[final_revision_start == 1], NA)
    ) %>%
    summarize_by_time(startClock2, .by = paste0(timeinterval, " seconds"),
                      charProduction = last(charProduction),
                      doclengthFull = last(doclengthFull),
                      positionFull = last(positionFull),
                      rev_start = last(rev_start),
    ) %>%
    ggplot() +
    # positionline
    geom_line(aes(x = startClock2, y = positionFull, 
                  color = "Position"), linetype = "dashed") +
    # predicted changepoints
    geom_vline(data = dataplot %>% filter(session_id %in% subsetsessions) %>% 
                 group_by(session_id) %>%
                 filter(cp_prob > min_cpprob) %>%
                 select(cp_loc, session_id), 
               mapping = aes(xintercept = cp_loc, 
                             color = factor("EstimatedBreakpoint"))) +
    # confidence interval of changepoints (shaded area)
    geom_rect(data = dataplot %>% filter(session_id %in% subsetsessions) %>% 
                group_by(session_id) %>%
                filter(cp_prob > min_cpprob) %>%
                mutate(
                  cp_ci_lower = as.POSIXct(cp_ci_lower, 
                                           origin = "1970-01-01"),
                  cp_ci_upper = as.POSIXct(cp_ci_upper, 
                                           origin = "1970-01-01")) %>%
                select(cp_ci_upper, cp_ci_lower, session_id),
              aes(xmin=cp_ci_lower, xmax=cp_ci_upper, ymin=-Inf, ymax=Inf), 
              alpha = .2, fill = "grey", inherit.aes = FALSE) +
    # processline
    geom_line(aes(x = startClock2, y = charProduction, color = "Process")) +
    # productline
    geom_line(aes(x = startClock2, y = doclengthFull, color = "Product")) +
    # actual revision end
    geom_vline(aes(xintercept = rev_start, color = factor("RevisionStart")),
               linetype = "dashed", linewidth = 0.8) +
    theme(text = element_text(size=16)) +
    theme_bw() +
    scale_colour_manual(name = "", 
                        breaks = c("Process", "Product", "Position", 
                                   "RevisionStart","EstimatedBreakpoint"),
                        values=c(Process="blue", Product="darkgreen", 
                                 Position="green", RevisionStart = "red",
                                 EstimatedBreakpoint = "darkgrey")) +

    ylab("Number of characters") +
    xlab("Time") +
    facet_wrap(~session_id, scales = "free") 
}

plot_accuracy_breakpoints <- function(preproc = "time", var = "correct_cp", 
                                      models = c(1,2,3,4,5,123,12345)){
  
  dataplot = read_csv(paste0("output/", preproc, "_diff_mod5all.csv"),
                      show_col_types = FALSE) %>%
    mutate( AnnotationRule = ordered(as.character(model),
                            levels = c("1","2","3","4","5","123", "12345"),
                            labels = c("1. Document length", 
                                       "2. Distance process/product",
                                       "3. Relative position",
                                       "4. Source use",
                                       "5. Pause duration",
                                       "Rule 1-3",
                                       "Rule 1-5")),
            maxcp = ordered(as.character(maxcp), 
                            levels = c("1", "3", "5","10", "20")),
            timeinterval = ordered(as.character(timeinterval),
                                   levels = c("1", "5", "10"),
                                   labels = c("Per 1 sec", "Per 5 sec", "Per 10 sec")),  
            trim = ordered(as.character(trim), 
                           labels = c("No trimming", "10% trimming"))) %>%
    filter(model %in% models)
  
  y_labels_scale = waiver()
  
  #change y axis based on selected var
  if(var == "correct_cp"){
    varname = "Percentage of correct changepoints"
    var = paste0(var, "/rev_present")
    y_labels_scale = scales::percent
  }else if(var == "close_cp"){
    varname = "Percentage of correct changepoints within 10 sec"
    var = paste0(var, "/rev_present")
    y_labels_scale = scales::percent
  }else if(var == "cp_in_ci"){
    varname = "Percentage of correct changepoints within 95% credible interval"
    var = paste0(var, "/rev_present")   
    y_labels_scale = scales::percent
  }else if(var == "mean_diff"){
    varname = "Mean difference of change points (sec)"
  }else if(var == "median_diff"){
    varname = "Median difference of change points (sec)"
  }else if(var == "sd_diff"){
    varname = "Standard deviation difference of change points (sec)"
  }else{
    varname = "unknown"
  }
  
  # plot the differences
  dataplot  %>%
    ggplot() +
    geom_bar(aes_string(x = "maxcp", y = var, color = "AnnotationRule", 
                        fill = "AnnotationRule"),
             stat = "identity", position = position_dodge()) +
    facet_grid(trim ~ timeinterval) +
    theme_bw() +
    theme(text = element_text(size=16)) +
    scale_y_continuous(labels = y_labels_scale) +
    ylab(varname) +
    xlab("Maximum number of change points") 

}

plot_final_breakpoint <- function(dataset, timeinterval = 5,
                                  subsetsessions, min_cpprob = 0){
  
  beast_123_pick <- read_csv("output/beast_123_pick.csv", 
                             show_col_types = FALSE)
  beast_pick <- read_csv("output/beast_pick.csv", show_col_types = FALSE)
  
  dataset %>%
  filter(session_id %in% subsetsessions) %>%
  group_by(participant, session, session_id) %>%
  mutate(
    rev_start = ifelse(final_revision_present == 1, 
                       startClock2[final_revision_start == 1], NA)
  ) %>%
  summarize_by_time(startClock2, .by = paste0(timeinterval, " seconds"),
                    charProduction = last(charProduction),
                    doclengthFull = last(doclengthFull),
                    positionFull = last(positionFull),
                    rev_start = last(rev_start),
  ) %>%
  ggplot() +
  # positionline
  geom_line(aes(x = startClock2, y = positionFull, 
                color = "Position"), linetype = "dashed") +
  # predicted segments
  geom_vline(data = beast_123_pick %>% filter(session_id %in% subsetsessions) %>% 
               group_by(session_id) %>%
               select(cp_mean), 
             mapping = aes(xintercept = cp_mean, 
                           color = "EstimatedBreakpoints")) +
  # picked segment
  geom_vline(data = beast_pick %>% filter(session_id %in% subsetsessions) %>% 
               group_by(session_id) %>%
               select(selected_cp3), 
             mapping = aes(xintercept = selected_cp3, 
                           color = "PickedBreakpoint"), linewidth = 0.8, 
             linetype = "dotdash") +
  # processline
  geom_line(aes(x = startClock2, y = charProduction, color = "Process")) +
  # productline
  geom_line(aes(x = startClock2, y = doclengthFull, color = "Product")) +
  # actual revision end
  geom_vline(aes(xintercept = rev_start, color = "RevisionStart"),
             linewidth = 0.8, 
             linetype = "dashed") +
  theme_bw() +
  theme(text = element_text(size=16)) +
  scale_colour_manual(name = "", 
                      breaks = c("Process", "Product", "Position", 
                                 "RevisionStart","EstimatedBreakpoints", 
                                 "PickedBreakpoint"),
                      values=c(Process="blue", Product="darkgreen", 
                               Position="green", RevisionStart = "red",
                               EstimatedBreakpoints = "darkgrey",
                               PickedBreakpoint = "purple")) +
  ylab("Number of characters") +
  xlab("Time") +
  facet_wrap(~session_id, scales = "free") 
}

