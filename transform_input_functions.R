################################################################################
# create time log (per x sec)

create_timelog <- function(data, timeinterval = 1) {

log_time_sec <- data %>%
  group_by(participant, session, session_id) %>%
  mutate(
    pausetime_99 = ifelse(pauseTime < quantile(pauseTime, .99),
                          pauseTime, NA)) %>%
  group_by(participant, session, session_id,
           isna = is.na(pausetime_99)) %>%
  mutate(perc_doclength = doclengthFull / last(doclengthFull),
         dist_process_prod = charProduction - doclengthFull,
         rel_position = positionFull / doclengthFull,
         max_end = ifelse(isna,  NA, max(endTime, na.rm = T)),
         source_use = ifelse(isna,  NA, cumsum(type == "focus")),
         # (only include pauses in 99% percentile)
         m_pause_length = ifelse(isna,  NA,
                                 cummean(pausetime_99))) %>%
  ungroup() %>%
  dplyr::select(-isna) %>%
  group_by(participant, session, session_id) %>%
  fill(max_end:m_pause_length) %>%
  # summarize per x sec of the time elapsed
  summarize_by_time(startClock2, .by = paste0(timeinterval, " seconds"),
                    perc_doclength = last(perc_doclength),
                    dist_process_prod = last(dist_process_prod),
                    rel_position = last(rel_position),
                    source_use = last(source_use),
                    m_pause_length = last(m_pause_length),
                    maxend = last(max_end),
                    rev_start = sum(final_revision_start)
  ) %>%
  # flesh out data such that every row is x sec (equally-spaced time intervals)
  pad_by_time(.date_var = startClock2) %>%
  mutate(rev_start = ifelse(is.na(rev_start), 0, rev_start)) %>%
  fill(perc_doclength:maxend) 


write.csv(log_time_sec, paste0("app/input/log_time_", timeinterval, 
                                "sec.csv"), row.names = F)
log_time_sec

}

################################################################################
# create time log (per % of length of total session)
library(padr)

create_perclog <- function(data, percinterval = 1) {
  
  percentages = data %>%
    group_by(participant, session, session_id) %>%
    reframe(perc_time_elapsed = seq(0,1, percinterval/100))
  
log_perc_time <- data %>%
  group_by(participant, session, session_id) %>%
  mutate(
    pausetime_99 = ifelse(pauseTime < quantile(pauseTime, .99),
                          pauseTime, NA)) %>%
  group_by(participant, session, session_id,
           isna = is.na(pausetime_99)) %>%
  mutate(perc_doclength = doclengthFull / last(doclengthFull),
         dist_process_prod = charProduction - doclengthFull,
         rel_position = positionFull / doclengthFull,
         max_end = ifelse(isna,  NA, max(endTime, na.rm = T)),
         source_use = ifelse(isna,  NA, cumsum(type == "focus")),
         # (only include pauses in 99% percentile)
         m_pause_length = ifelse(isna,  NA,
                                 cummean(pausetime_99)),
         # calculate for every x percent of time elapsed
         perc_time_elapsed = round(as.numeric(startClock2) /
           as.numeric(max(endClock)), digits = -log10(percinterval/100))) %>%
  ungroup() %>%
  dplyr::select(-isna) %>%
  group_by(participant, session, session_id) %>%
  fill(max_end:m_pause_length) %>%
  group_by(participant, session, session_id, perc_time_elapsed) %>%
  # summarize per 1% of the time elapsed
  summarize(
    perc_doclength = last(perc_doclength),
    dist_process_prod = last(dist_process_prod),
    rel_position = last(rel_position),
    source_use = last(source_use),
    m_pause_length = last(m_pause_length),
    maxend = first(max_end)
  ) %>%
  # flesh out data such that every row is x percent (equally-spaced intervals)
  right_join(percentages) %>%
  arrange(participant, session, session_id, perc_time_elapsed) %>%
  fill(perc_doclength:maxend) %>%
   mutate(
     time = round(perc_time_elapsed*maxend))   

write.csv(log_perc_time, paste0("app/input/log_perc_", percinterval, "perc.csv"),
                                row.names = F)
log_perc_time

}