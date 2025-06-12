library(Hmisc)

# explore indicators
log_filt2 <- read.csv("app/input/log_filt2.csv")

log_filt2_indicators <- log_filt2 %>%
  group_by(participant, session, session_id) %>%
  mutate(
    pausetime_99 = ifelse(pauseTime < quantile(pauseTime, .99),
                          log1p(pauseTime), NA),
    sourcetime = ifelse(!is.na(lag(endTime)) & type == "focus",
                        (startTime - lag(endTime)), 0)) %>%
  group_by(participant, session, session_id,
           isna = is.na(pausetime_99)) %>%
  mutate(perc_doclength = doclengthFull / last(doclengthFull),
         dist_process_prod = charProduction - doclengthFull,
         rel_position = positionFull / doclengthFull,
         max_end = ifelse(isna,  NA, max(endTime, na.rm = T)),
         sum_source_length = ifelse(isna,  NA, log1p(cumsum(sourcetime))),
         # (only include pauses in 99% percentile)
         m_pause_length = ifelse(isna,  NA,
                                 (cummean(pausetime_99)))) %>%
  ungroup() %>%
  dplyr::select(-isna) %>%
  group_by(participant, session, session_id) %>%
  fill(max_end:m_pause_length) %>%
  select(perc_doclength, dist_process_prod, rel_position,sum_source_length,
         m_pause_length) %>%
  ungroup()

cor <- rcorr(as.matrix(log_filt2_indicators  %>%
                         select(perc_doclength:m_pause_length)),
             type = "pearson")
options(digits = 16)
print(cor)
