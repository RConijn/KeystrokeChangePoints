library(tidyverse)
library(tsibble)
library(lubridate)
library(timetk)
options(scipen = 9999)
options(digits.secs=3)

log_plantra <- read_csv("PlanTra_Data/data_out/plantra_annotated.csv") %>%
  mutate(session = ifelse(session == 1, "plantra-pre", "plantra-post"))
log_lift <- read_csv("LIFT_Data/data_out/lift_annotated.csv")

# create extra variables for model
log_filt2 <- log_plantra %>%
  bind_rows(log_lift) %>%
  group_by(participant, session) %>%
  mutate(
    rel_dist_leadingedge = (doclengthFull - positionFull) /max(doclengthFull),
    first_time = min(startTime),
    char_per_doclength = charProduction / doclengthFull,
    # set first keystroke to time = 0
    startClock2 = as_datetime((startTime-first_time)/1000),
    # session_id per session
    session_id = cur_group_id(),
    session_id2 = cur_group_id(),
    final_revision_start = ifelse(is.na(final_revision_start), 0,
                                  final_revision_start)) %>%
  # remove synchronous keystrokes
  filter(startTime != lag(startTime) | row_number() == 1) %>%
  mutate(index = row_number()) %>%
  ungroup()

write.csv(log_filt2, "app/input/log_filt2.csv", row.names = F)

################################################################################
# create different input files (transformed log files)
source("app/scripts/transform_input_functions.R")

# aggregate per x seconds
log_time_1sec  <- create_timelog(log_filt2, 1)
log_time_5sec  <- create_timelog(log_filt2, 5)
log_time_10sec <- create_timelog(log_filt2, 10)

# aggregate per x % of total session 
log_time_1perc <- create_perclog(log_filt2, 1)
log_time_0.1perc <- create_perclog(log_filt2, 0.1)


