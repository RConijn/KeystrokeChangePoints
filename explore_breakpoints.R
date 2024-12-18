# explore all breakpoints as estimated in rbeast3
library(Rbeast)
library(tidyverse)
options(scipen = 9999)
options(digits.secs=3)

library(lubridate)
library(timetk)


# all vars:
# 1) doclengthFull percentage flattens
# 2) difference process - product becomes larger
# 3) relative position mouse to start
# 4) source usage (focus) decreases --> check if source use > 0
# 5) average pause length increases

# load scripts with plot functions
source("app/scripts/plot_changepoints_functions.R")

################################################################################
# plot breakpoints per ppt (see also dashboard)

# plot for  e.g. model 1
log_filt2 <- read_csv("app/input/log_filt2.csv", show_col_types = FALSE)
plot_time_breakpoints(log_filt2, model = 1, maxcp = 5, timeinterval = 1,
                      trim = 0, subsetsessions = 1:4, min_cpprob = 0.0)

################################################################################
# calculate the distance from the manual annotation
distance_from_cp <- function(model, maxcp = 10, timeinterval = 1,
                             trim = 0){
  
  
  trimyes = ifelse(trim > 0, "trim", "")
  
  input = read_csv(paste0("app/input/log_time_", timeinterval, "sec.csv"),
                   show_col_types = FALSE)
  
  estimate_cp = read_csv(paste0("app/models4/rbeast_time", timeinterval, trimyes,
                             "_", model, "_", maxcp, "cp.csv"),
                      show_col_types = FALSE)
  
  log_rev <- input %>%
    group_by(session_id) %>%
    summarize(
      rev_start = ifelse(sum(rev_start, na.rm = T) > 0,
                         first(startClock2[rev_start == 1]), NA),
      max_end = as.numeric(last(startClock2)),
      n = n()) %>%
    mutate(
      rev_present = as.numeric(!is.na(rev_start))
    )

  # check how much overlap
  beast_check <- estimate_cp %>% 
    arrange(session_id, cp_loc) %>%
    right_join(log_rev, by = join_by(session_id)) %>%
    mutate(
      correct_cp = rev_start == cp_loc,
      diff = abs(rev_start - cp_loc),
      close_cp = abs(rev_start - cp_loc) <= 10,
      cp_in_ci = rev_start >= cp_ci_lower & rev_start <= cp_ci_upper
    ) %>%
    group_by(session_id) %>%
    slice(which.min(abs(diff))) %>%
    group_by(pre_proc = "time", timeinterval = timeinterval, trim = trim,
             model = model, maxcp = maxcp) 
}

# run for 1 example
lr2 <- distance_from_cp(model = 123, maxcp = 10, timeinterval = 10, trim = 10)

models = c(1,2,3,4,5)
maxcps = c(1,3,5,10,20)
timeintervals = c(1,5,10)
trims = c(0,10)

FullGrid <- expand.grid(models, maxcps, timeintervals, trims)
colnames(FullGrid) <- c("models","maxcps","timeintervals", "trims")

# run for the full sample
output <- as_tibble(do.call(rbind,
                  mapply(distance_from_cp, FullGrid$models, FullGrid$maxcps, 
                         FullGrid$timeintervals, FullGrid$trims, 
                         SIMPLIFY=FALSE)))


outputsum <- output %>%
  group_by(model, trim, timeinterval, session_id) %>%
  arrange( .by_group = TRUE) %>%
  mutate(bf = (mlik) / (lag(mlik)),
         bflog = exp(mlik) / exp(first(mlik)),
         k = maxcp *2,
         adj_r2 = 1 - (((1 - R2) * (n - 1)) / (n - k - 1))) %>%
  group_by(timeinterval, trim, model, maxcp) %>%
  summarize(
    rev_present = sum(rev_present == 1, na.rm = T),
    correct_cp = sum(correct_cp, na.rm = T),
    close_cp = sum(close_cp, na.rm = T),
    cp_in_ci = sum(cp_in_ci, na.rm = T),
    median_diff = median(diff, na.rm = T),
    mean_diff = mean(diff, na.rm = T),
    sd_diff = sd(diff, na.rm = T),
    mean_r2 = mean(R2, na.rm = T),
    sd_r2 = sd(R2,na.rm = T),
    mean_adjr2 = mean(adj_r2, na.rm = T),
    sd_adjr2 = sd(adj_r2,na.rm = T),
    mean_rmse = mean(RMSE, na.rm = T),
    sd_rmse = sd(RMSE,na.rm = T),
    mean_bf = mean(bf, na.rm = T),
    sd_bf = sd(bf, na.rm = T),
    mean_mlik = mean(mlik, na.rm = T),
    sd_mlik = sd(mlik, na.rm = T)
  ) 
  
write_csv(outputsum, "app/output/time_diff_mod4.csv")

###############################################################################
## plot outcomes for timeinterval 5, trim = 0
outputsum  %>%
  filter(timeinterval == 5, trim == 0, ) %>%
  pivot_longer(cols = c("correct_cp", "close_cp")) %>%
  mutate( percent = value / 138,
          name = ordered(as.character(name),
                         levels = c("correct_cp", "close_cp"),
                         labels = c("Correct", "Within 10sec")),
          AnnotationRule = ordered(as.character(model),
                                   levels = c("1","2","3","4","5"),
                                   labels = c("1. Document length", 
                                              "2. Distance process/product",
                                              "3. Relative position",
                                              "4. Source use",
                                              "5. Pause duration")),
          maxcp = ordered(as.character(maxcp), 
                          levels = c("1", "3", "5","10", "20"))) %>%
  ggplot() +
  geom_bar(aes_string(x = "maxcp", y = "percent", color = "AnnotationRule", 
                      fill = "AnnotationRule"),
           stat = "identity", position ="identity", alpha = 0.5) +
  facet_grid(.~AnnotationRule, 
             labeller = label_wrap_gen()) +
  theme_bw() +
  theme(text = element_text(size=14)) +
  scale_y_continuous(labels = scales::percent) +
  ylab("Change points (correct / within 10s)") +
  xlab("Maximum number of change points") 

################################################################################

# run for Multivariate models (cannot output RMSE/R2/mlik)
models = c(123,12345)
maxcps = c(1,3,5,10,20)
timeintervals = c(1,5,10)
trims = c(0,10)

FullGrid <- expand.grid(models, maxcps, timeintervals, trims)
colnames(FullGrid) <- c("models","maxcps","timeintervals", "trims")


output <- as_tibble(do.call(rbind,
                            mapply(distance_from_cp, FullGrid$models, FullGrid$maxcps, 
                                   FullGrid$timeintervals, FullGrid$trims, 
                                   SIMPLIFY=FALSE)))

outputsummv <- output %>%
  group_by(model, trim, timeinterval, session_id) %>%
  group_by(timeinterval, trim, model, maxcp) %>%
  summarize(
    rev_present = sum(rev_present == 1, na.rm = T),
    correct_cp = sum(correct_cp, na.rm = T),
    close_cp = sum(close_cp, na.rm = T),
    cp_in_ci = sum(cp_in_ci, na.rm = T),
    median_diff = median(diff, na.rm = T),
    mean_diff = mean(diff, na.rm = T),
    sd_diff = sd(diff, na.rm = T),
  ) 

write_csv(outputsummv, "app/output/time_diff_mod4mv.csv")

# combine all outcomes in one
outputsumall <- outputsum %>% rbind(outputsummv)

write_csv(outputsumall, "app/output/time_diff_mod4all.csv")


################################################################################
# plot the differences (also see dashboard)
plot_accuracy_breakpoints(var= "close_cp")
plot_accuracy_breakpoints(var= "correct_cp")
plot_accuracy_breakpoints(var= "cp_in_ci")
plot_accuracy_breakpoints(var= "mean_diff")
plot_accuracy_breakpoints(var= "median_diff")

# combine univariate models
beast_12345 <- mutate(beast_1, model = "perc_doclength") %>% 
  rbind(mutate(beast_2, model = "dist_process_prod")) %>%
  rbind(mutate(beast_3, model = "rel_position")) %>%
  rbind(mutate(beast_4, model = "source_use")) %>%
  rbind(mutate(beast_5, model = "m_pause_length")) 




