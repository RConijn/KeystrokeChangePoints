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
  
  estimate_cp = read_csv(paste0("app/models5/rbeast_time", timeinterval, trimyes,
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
lr2 <- distance_from_cp(model = 123, maxcp = 10, timeinterval = 5, trim = 0)

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
    # CrI around mean
    ci_lower_diff = quantile(diff, probs = c(0.025, 0.975), na.rm = TRUE)[1],
    ci_upper_diff = quantile(diff, probs = c(0.025, 0.975), na.rm = TRUE)[2],
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
  
write_csv(outputsum, "app/output/time_diff_mod5.csv")

outputsumforpaper <- outputsum %>%
  filter(trim == 0, timeinterval == 5) %>%
  mutate(percentage_correct = round(correct_cp / rev_present *100,1),
         percentage_CI = round(cp_in_ci / rev_present *100,1),
         CI_diff = paste0("[", round(ci_lower_diff), ";",
                          round(ci_upper_diff), "]"),
         median_diff = round(median_diff),
         mean_mlik = round(mean_mlik),
         mean_adjr2 = round(mean_adjr2, 3)) %>%
  select(model, maxcp, percentage_correct, percentage_CI, median_diff,
         CI_diff, mean_mlik, mean_adjr2)

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
                                              "4. Length source use",
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
models = c(123)
maxcps = c(1,3,5,10,20)
timeintervals = c(1,5,10)
trims = c(0,10)

FullGrid <- expand.grid(models, maxcps, timeintervals, trims)
colnames(FullGrid) <- c("models","maxcps","timeintervals", "trims")


outputmv <- as_tibble(do.call(rbind,
                            mapply(distance_from_cp, FullGrid$models, FullGrid$maxcps, 
                                   FullGrid$timeintervals, FullGrid$trims, 
                                   SIMPLIFY=FALSE)))

outputsummv <-  outputmv %>%
  group_by(timeinterval, trim, model, maxcp) %>%
  summarize(
    rev_present = sum(rev_present == 1, na.rm = T),
    correct_cp = sum(correct_cp, na.rm = T),
    close_cp = sum(close_cp, na.rm = T),
    cp_in_ci = sum(cp_in_ci, na.rm = T),
    median_diff = median(diff, na.rm = T),
    mean_diff = mean(diff, na.rm = T),
    sd_diff = sd(diff, na.rm = T),
    mean_mlik = mean(mlik, na.rm = T),
    sd_mlik = sd(mlik, na.rm = T),
    # CrI around mean
    ci_lower_diff = quantile(diff, probs = c(0.025, 0.975), na.rm = TRUE)[1],
    ci_upper_diff = quantile(diff, probs = c(0.025, 0.975), na.rm = TRUE)[2]
  ) 

outputsummvforpaper <- outputsummv %>%
  filter(trim == 0, timeinterval == 5) %>%
  mutate(percentage_correct = round(correct_cp / rev_present *100,1),
         percentage_CI = round(cp_in_ci / rev_present *100,1),
         CI_diff = paste0("[", round(ci_lower_diff), ";",
                          round(ci_upper_diff), "]"),
         median_diff = round(median_diff)) %>%
  select(model, maxcp, percentage_correct, percentage_CI, median_diff,
         CI_diff, mean_mlik)

write_csv(outputsummv, "app/output/time_diff_mod5mv.csv")

# combine all outcomes in one
outputsumall <- outputsum %>% rbind(outputsummv)

write_csv(outputsumall, "app/output/time_diff_mod5all.csv")

################################################################################

# Check overlap between selected changepoints across models
cps1 = read_csv("app/models5/rbeast_time5_1_10cp.csv",
                       show_col_types = FALSE)
cps2 = read_csv("app/models5/rbeast_time5_2_10cp.csv",
                show_col_types = FALSE)
cps3 = read_csv("app/models5/rbeast_time5_3_10cp.csv",
                show_col_types = FALSE)
cps123 = read_csv("app/models5/rbeast_time5_123_10cp.csv",
                show_col_types = FALSE)

# identify number of overlapping changepoints per model
allcps <- mutate(cps1, model = "model1") %>% 
  rbind(mutate(cps2, model = "model2")) %>%
  rbind(mutate(cps3, model = "model3"))  %>%
  select(session_id:cp_loc, model) %>%
  rbind(cps123 %>% 
          select(session_id:cp_loc) %>%
          mutate(model = "model123")) %>% 
  # only select changepoints with higher probability
  filter(cp_prob > 0.7) 

allcpsoverlap <- allcps %>%
group_by(session_id) %>%
  arrange(session_id,cp_loc) %>%
  mutate(
    # overlapping cp's will be combined
    ci_non_overlap = 1- ifelse(row_number() == 1, FALSE,
                               # overlap in confidence interval
                               (lag(cp_ci_upper) > cp_ci_lower) |
                                 # overlap within 10 secs
                                 ((cp_loc - 10) < lag(cp_loc))),
    cumsum = cumsum(ci_non_overlap)
  ) %>% 
  group_by(session_id, cumsum) %>%
  summarize(
    overlap = n(),
    n_models = n_distinct(model),
    cp_mean = mean(cp_loc, na.rm = T),
    cp_median = median(cp_loc, na.rm = T),
    model1 = sum(model == "model1"),
    model2 = sum(model == "model2"),
    model3 = sum(model == "model3"),
    model123 = sum(model == "model123"),
    )

# get mean/sd changepoint detected per model
allcp_m_sd <- allcps %>%
  group_by(session_id, model) %>%
  summarize(
    n = n()
  ) %>%
  group_by(model) %>%
  summarize(mean(n), sd(n))

# get overlap across model pairs
allcp_sum <- allcpsoverlap %>%
  ungroup() %>%
  summarise(
    n1 = sum(model1 > 0),
    n2 = sum(model2 > 0),
    n3 = sum(model3 > 0),
    n123 = sum(model123 > 0),
    perc_model1_model2 = sum(model1 > 0 & model2 > 0) / (n1+n2),
    perc_model1_model3 = sum(model1 > 0 & model3 > 0) / (n1+n3),
    perc_model1_model123 = sum(model1 > 0 & model123 > 0) / (n1+n123),
    perc_model2_model3 = sum(model1 > 0 & model3 > 0) / (n2+n3),
    perc_model2_model123 = sum(model1 > 0 & model123 > 0) / (n2+n123),
    perc_model3_model123 = sum(model1 > 0 & model123 > 0) / (n3+n123),
  ) %>%
  select(starts_with("perc_")) %>%
  pivot_longer(cols = starts_with("perc_"), 
               names_to = "pair", 
               values_to = "value") %>%
  mutate(value = round(value,2)) %>%
  separate(pair, into = c("prefix", "model", "to"), sep = "_") %>%
  pivot_wider(names_from = to, values_from = value) %>%
  full_join(allcp_m_sd)


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




