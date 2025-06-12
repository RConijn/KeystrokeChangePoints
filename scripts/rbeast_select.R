library(tidyverse)
library(ggplot2)
library(timetk)
library(ggdist)
theme_set(theme_ggdist())

# read models
beast_1 <- read_csv("app/models5/rbeast_time5_1_10cp.csv")
beast_2 <- read_csv("app/models5/rbeast_time5_2_10cp.csv")
beast_3 <- read_csv("app/models5/rbeast_time5_3_10cp.csv")
log_time_sec <- read_csv("app/input/log_time_5sec.csv")

# get manual annotated revision end
log_rev1 <- log_time_sec %>%
  group_by(session_id) %>%
  mutate(
    max_end = as.numeric(last(startClock2))      
  )


log_rev <- log_rev1 %>%
  summarize(
    rev_start = ifelse(sum(rev_start, na.rm = T) > 0,
                       first(startClock2[rev_start == 1]), NA),
    maxend = first(max_end)) %>%
  mutate(
    rev_present = as.numeric(!is.na(rev_start))
  )

################################################################################
##### Pick changepoint per model ###############################################
# pick based on model 1
# 1) product line flattens
# slope decreases or stays flat (decreasing slope /flattest slope)
# LENIENT RULE: slope after is smaller than slope before
beast_1_pick <- beast_1 %>%
  left_join(select(log_rev1, session_id, perc_doclength, 
                   startClock2, session, max_end) %>%
              mutate(startClock2 = as.numeric(startClock2)),
            by = join_by(cp_loc == startClock2, session_id)) %>%
  group_by(session_id) %>%
  filter(
    # rule 1. after first 1/3rd of the process
    cp_loc > 0.33*max(max_end),
    # rule 2. majority of final product written
    perc_doclength > .7,
    # rule 3. probability of cp relatively high
    cp_prob > .7,
    # rule 4. relatively narrow CI
    cp_ci_upper - cp_ci_lower <= 60) %>%
  mutate(
    productline_flat1 = rank(abs(slope_after)),
    productline_flat2 = rank(-abs(slope_after_pr_zero)),
    productline_dec = slope_before > slope_after) 

# pick based on model 2
# 2) difference process - product becomes larger
# sudden increase in distance prod/process (deletion of a larger chunk)
# steep increase in slope
beast_2_pick <- beast_2 %>%
  group_by(session_id) %>%
  left_join(select(log_rev1, session_id, perc_doclength, 
                   startClock2, session, max_end) %>%
              mutate(startClock2 = as.numeric(startClock2)),
            by = join_by(cp_loc == startClock2, session_id)) %>%
  group_by(session_id) %>%
  filter(
    # rule 1. after first 1/3rd of the process
    cp_loc > 0.33*max(max_end),
    # rule 2. majority of final product written
    perc_doclength > .7,
    # rule 3. probability of cp relatively high
    cp_prob > .7,
    # rule 4. relatively narrow CI
    cp_ci_upper - cp_ci_lower <= 60) %>%
  mutate(
    largest_procprod_diff1 = rank(-(slope_after-lag(slope_after))),
    largest_procprod_diff2 = rank(-cp_change_abrupt)) 

# pick based on model 3
# 3) relative position mouse to start
# steep drop in Y (lead est_Y of prev. segment is a lot smaller)
beast_3_pick <- beast_3 %>%
  group_by(session_id) %>%
  left_join(select(log_rev1, session_id, perc_doclength, 
                   startClock2, session, max_end) %>%
              mutate(startClock2 = as.numeric(startClock2)),
            by = join_by(cp_loc == startClock2, session_id)) %>%
  group_by(session_id) %>%
  filter(
    # rule 1. after first 1/3rd of the process
    cp_loc > 0.33*max(max_end),
    # rule 2. majority of final product written (max process length > 70%)
    perc_doclength > .7,
    # rule 3. probability of cp relatively high
    cp_prob > .7,
    # rule 4. relatively narrow CI
    cp_ci_upper - cp_ci_lower <= 60) %>%
  mutate(
    largest_dec_cursor = rank(cp_change_abrupt)) 


################################################################################
##### PICK changepoint on all models combined ##################################

# calculate position of abrupt change
# 1. big increase in distance process vs. product (100 characters or more)
# 2. big move towards start of document (rel position drop 30% or more)
big_jumps <- log_time_sec %>% 
  group_by(session_id) %>%
  mutate(big_increase_dist = dist_process_prod - lag(dist_process_prod),
         big_drop_cursor = rel_position - lag(rel_position)) %>%
  filter(big_increase_dist > 100 | big_drop_cursor < -0.30) %>%
  mutate(startClock2 = as.numeric(startClock2)) %>%
  select(session_id, startClock2, big_increase_dist, big_drop_cursor) %>%
  nest() 

# combine univariate models
beast_123 <- mutate(beast_1_pick, model = "perc_doclength") %>% 
  rbind(mutate(beast_2_pick, model = "dist_process_prod")) %>%
  rbind(mutate(beast_3_pick, model = "rel_position"))  %>%
  group_by(session_id) %>%
  # add info about manual annotated final revision
  right_join(log_rev, by = join_by(session_id)) %>%
  group_by(session_id, model) %>%
  mutate(
    diff_cp = cp_ci_lower - lag(cp_ci_upper)
  )



################################################################################
# apply general rules to create subset of change point candidates
beast_123_pick <- beast_123 %>%
  left_join(big_jumps) %>%
  group_by(session_id) %>%
  rowwise()%>%
  mutate(
    # check if big jump is within cp_ci
    which_cp = list(which(data$startClock2 > cp_ci_lower & 
                            data$startClock2 < cp_ci_upper)),
    in_cp = !length(which_cp) == 0,
    drop_cursor = ifelse(in_cp, 
                         min(data$big_drop_cursor[which_cp],
                             na.rm = T), NA),
    increase_dist = ifelse(in_cp, 
                          max(data$big_increase_dist[which_cp],
                              na.rm = T), NA)) %>%
  filter(
    # filter the breakpoints that adhere to the rules (per criteria)
    (productline_flat1 <= 3 &  productline_flat1 != 0)| 
    (productline_flat2 <= 3 &  productline_flat2 != 0)|
      productline_dec |
    (largest_procprod_diff2 <= 3 & largest_procprod_diff2 != 0 & 
       increase_dist > 100 & !is.na(increase_dist)) | 
    (largest_dec_cursor <= 3 & largest_dec_cursor != 0 & 
       drop_cursor < -0.3 & !is.na(drop_cursor))
    ) %>% 
  arrange(session_id, cp_loc) %>%
  select(session_id:cp_loc, perc_doclength, productline_flat1:rev_start, in_cp,
         session, increase_dist, drop_cursor, perc_doclength) %>%
  group_by(session_id) %>%
  mutate(
    # rule 5. overlapping cp's will be combined
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
    cp_ci_lower = min(cp_ci_lower, na.rm = T),
    cp_ci_upper = max(cp_ci_upper, na.rm = T),
    cp_prob = mean(cp_prob, na.rm = T),
    cp_mean = mean(cp_loc, na.rm = T),
    cp_median = median(cp_loc, na.rm = T),
    cps = paste0(cp_loc, collapse = ","),
    productline_flat1 = sum(productline_flat1, na.rm = T),
    productline_flat2 = sum(productline_flat2, na.rm = T),
    largest_procprod_diff2 = sum(largest_procprod_diff2, na.rm = T),
    largest_dec_cursor = sum(largest_dec_cursor, na.rm = T),
    drop_cursor = ifelse(!all(is.na((drop_cursor))), 
                                min(drop_cursor, na.rm = T), NA),
    increase_dist = ifelse(!all(is.na((increase_dist))), 
                           max(increase_dist, na.rm = T), NA),
    perc_doclength = ifelse(!all(is.na((perc_doclength))), 
                            mean(perc_doclength, na.rm = T), NA),
    rev_start = first(rev_start),
    in_cp = first(in_cp),
    session = first(session)
            )

write.csv(beast_123_pick, "app/output/beast_123_pick.csv",row.names = F)
  
################################################################################
# check accuracy after initial filtering step ##################################

# check number of points left
beast_123_pick_sum <- beast_123_pick %>%
  filter(!is.na(rev_start)) %>%
  group_by(session_id) %>%
  summarize(n = n(), 
            # check how many more than one rule
            n_3 = sum(n_models == 3),
            n_2 = sum(n_models == 2),
            n_1 = sum(n_models == 1),
            rev_start = ifelse((sum(rev_start, na.rm = T)) > 0, 
                               "Revision phase", "No revision phase"))

# plot points left per session
beast_123_pick_sum %>%
  mutate( session_id2 = factor(session_id, 
                  levels = as.character(session_id)[order(n,session_id)])) %>%
  pivot_longer(n_3:n_1, names_to = "Models", values_to = "changepoints") %>%
  filter(changepoints != 0) %>%
  mutate(Models = substr(Models, 3,3)) %>%
  ggplot() +
  geom_bar(aes(y = changepoints, x = session_id2, fill = Models),
           stat = "Identity") +
  facet_grid(rev_start~Models, scales = "free_y") +
  theme_bw() +
  theme(text = element_text(size=14),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_discrete("Session") +
  ylab("Number of changepoints") +
  coord_flip()

# plot points left
beast_123_pick_sum %>%
  pivot_longer(n_3:n_1, names_to = "Models", values_to = "changepoints") %>%
  filter(changepoints != 0) %>%
  mutate(Models = substr(Models, 3,3)) %>%
  ggplot() +
  geom_bar(aes(x = n, fill = Models)) +
  theme_bw() +
  theme(text = element_text(size=14)) +
  scale_x_continuous("Number of changepoints", breaks = 1:24) +
  ylab("Number of sessions")

table(beast_123_pick_sum$n)
sum(beast_123_pick_sum$n)
mean(beast_123_pick_sum$n)
sd(beast_123_pick_sum$n)

# check how much overlap with remaining cps
beast_check <- beast_123_pick %>% 
  arrange(session_id, cp_mean) %>%
  mutate(
    correct_cp = rev_start == cp_mean,
    diff = abs(rev_start - cp_mean),
    close_cp = abs(rev_start - cp_mean) <= 10,
    cp_in_ci = rev_start >= cp_ci_lower & rev_start <= cp_ci_upper
  ) %>%
  group_by(session_id) %>%
  slice(which.min(abs(diff))) 

# check overall accuracy
beast_check_sum <- beast_check %>%
  group_by(session) %>%
  summarize(
    n = n(),
    correct_cp = sum(correct_cp, na.rm = T),
    close_cp = sum(close_cp, na.rm = T),
    cp_in_ci = sum(cp_in_ci, na.rm = T),
    median_diff = median(diff, na.rm = T),
    mean_diff = mean(diff, na.rm = T),
    sd_diff = sd(diff, na.rm = T)
  ) 


################################################################################
# Rule-based to select ONE specific cp_location ################################

beast_pickt <- beast_123_pick %>%
  group_by(session_id, session) %>%
  mutate(drop_cursor = ifelse(is.na(drop_cursor), 0, drop_cursor),
         increase_dist = ifelse(is.na(increase_dist), 0, increase_dist)) %>%
  mutate(
    # if only 1 model left, pick that one
    selected_cp2 = ifelse(n() == 1 | sum(drop_cursor) == 0, 
                          first(cp_mean), NA),
    perc_perc_doc = perc_doclength / max(perc_doclength,na.rm =T)) %>%
  filter(perc_perc_doc > 0.80  | (n_models == 3 & perc_perc_doc > 0.70)
         | !is.na(selected_cp2),
        # distance decreases a lot - indicates weir copy-paste behavior
         increase_dist > -1000 | !is.na(selected_cp2))


beast_pick <- beast_pickt %>%
  summarize(
    selected_cp3 = case_when(
      any(!is.na(selected_cp2)) ~ first(selected_cp2),
      n() == 1  ~ first(cp_mean),
      # if 2/3 models & drop: pick the FIRST cp that shows this 
      any(n_models >1 & drop_cursor < -0.30 & largest_dec_cursor != 0) ~ {
        valid_indices <- which(n_models >1 & drop_cursor < -0.30 & 
                                 largest_dec_cursor != 0)
        if (length(valid_indices) > 0) cp_mean[min(valid_indices)] else NA_real_
      },
      any(n_models >= 1 & increase_dist > 50 & largest_dec_cursor != 0) ~ {
        valid_indices <- which(n_models >= 1 & increase_dist > 50 
                               & largest_dec_cursor != 0)
        if (length(valid_indices) > 0) cp_mean[min(valid_indices)] else NA_real_
      },
      any(n_models > 1  & largest_dec_cursor < 0) ~ {
        valid_indices <- which(n_models >1 & largest_dec_cursor != 0)
        if (length(valid_indices) > 0) cp_mean[min(valid_indices)] else NA_real_
      },
      any(largest_dec_cursor != 0) ~ {
        valid_indices <- which(largest_dec_cursor != 0)
        if (length(valid_indices) > 0) cp_mean[min(valid_indices)] else NA_real_
      },
      any(n_models > 1) ~ {
        valid_indices <- which(n_models > 1)
        if (length(valid_indices) > 0) cp_mean[min(valid_indices)] else NA_real_
      },
      
      TRUE ~ first(cp_mean)
    )
  ) %>%
  right_join(log_rev, by = "session_id") %>%
  mutate(dif = rev_start - selected_cp3)

write.csv(beast_pick, "app/output/beast_pick.csv",row.names = F)

