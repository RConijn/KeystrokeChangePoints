# RBEAST
library(Rbeast)
library(tidyverse)
options(scipen = 9999)
options(digits.secs=3)

library(lubridate)
library(timetk)


# load sample data
log_time_10sec <- read_csv("app/input/log_time_10sec.csv")

###############################################################################
# test for one ppt (use for debugging)

log_segm2 <- log_time_10sec %>%
  filter(session_id == 63) %>%
  # with trimming of 10 %
  # (don't include the first 10% of the writing phase in cp detection)
  filter(as.numeric(startClock2) >= (maxend/1000)*.1) %>%
  ungroup() 

out = beast(log_segm2$perc_doclength, 
            start = first(as.numeric(log_segm2$startClock2)), 
            deltat = 10,
            tcp.minmax = c(0, 1),
            season='none',#  'none': trend-only data without seasonality  
            mcmc.seed = 42,
            ci = F) # symmetric credible interval suffices

vars = c("perc_doclength", "dist_process_prod", "rel_position",
         "source_use","m_pause_length")

outmv = log_segm2 %>% select(all_of(vars)) %>%
  beast( 
            start = first(as.numeric(log_segm2$startClock2)), 10,
            tcp.minmax = c(0, 5),
            tseg.min = 10,
            season='none',#  'none': trend-only data without seasonality  
            mcmc.seed = 42,
            ci = F,
            print.options = F) 

            
print(out)                   
plot(out, vars  = c('y','t','tcp','error'))

###############################################################################
# load function
source("app/scripts/estimate_changepoints_function.R")

vars = c("perc_doclength", "dist_process_prod", "rel_position",
         "source_use","m_pause_length")
maxcps = c(3,5,10,20)
timeintervals = c(1,5,10)
trims = c(0,10)

FullGrid <- expand.grid(vars, maxcps, timeintervals, trims)
colnames(FullGrid) <- c("vars","maxcps","timeintervals", "trims")

# write all univariate model outputs to file
mapply(write_cp_to_file, FullGrid$vars, FullGrid$maxcps, 
                                   FullGrid$timeintervals, FullGrid$trims)


# OR ONE BY ONE
###########################
# 1) product line flattens
# run changepoint detection on model 1 (perc_doclength) with trim, 10sec
beast_1 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 3, timeinterval = 10, trim = 10)
beast_2 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 5, timeinterval = 10, trim = 10)
beast_3 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 10, timeinterval = 10, trim = 10)
beast_4 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 20, timeinterval = 10, trim = 10)
beast_5 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 1, timeinterval = 10, trim = 10)

# run changepoint detection on model 1 (perc_doclength) without trim, 10sec
beast_1 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 3, timeinterval = 10, trim = 0)
beast_2 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 5, timeinterval = 10, trim = 0)
beast_3 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 10, timeinterval = 10, trim = 0)
beast_4 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 20, timeinterval = 10, trim = 0)
beast_5 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 1, timeinterval = 10, trim = 0)

# run changepoint detection on model 1 (perc_doclength) with trim, 5sec
beast_1 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 3, timeinterval = 5, trim = 10)
beast_2 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 5, timeinterval = 5, trim = 10)
beast_3 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 10, timeinterval = 5, trim = 10)
beast_4 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 20, timeinterval = 5, trim = 10)
beast_5 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 1, timeinterval = 5, trim = 10)

# run changepoint detection on model 1 (perc_doclength) without trim, 5sec
beast_1 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 3, timeinterval = 5, trim = 0)
beast_2 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 5, timeinterval = 5, trim = 0)
beast_3 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 10, timeinterval = 5, trim = 0)
beast_4 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 20, timeinterval = 5, trim = 0)

# run changepoint detection on model 1 (perc_doclength) with trim, 1sec
beast_1 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 3, timeinterval = 1, trim = 10)
beast_2 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 5, timeinterval = 1, trim = 10)
beast_3 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 10, timeinterval = 1, trim = 10)
beast_4 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 20, timeinterval = 1, trim = 10)

# run changepoint detection on model 1 (perc_doclength) without trim, 1sec
beast_1 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 3, timeinterval = 1, trim = 0)
beast_2 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 5, timeinterval = 1, trim = 0)
beast_3 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 10, timeinterval = 1, trim = 0)
beast_4 <- write_cp_to_file(vars = c("perc_doclength"), 
                            maxcp = 20, timeinterval = 1, trim = 0)


###########################
# 2) difference process - product becomes larger
# run changepoint detection on model 2 (dist_process_prod) with trim, 10sec
beast_1 <- write_cp_to_file(vars = c("dist_process_prod"), 
                            maxcp = 3, timeinterval = 10, trim = 10)
beast_2 <- write_cp_to_file(vars = c("dist_process_prod"), 
                            maxcp = 5, timeinterval = 10, trim = 10)
beast_3 <- write_cp_to_file(vars = c("dist_process_prod"), 
                            maxcp = 10, timeinterval = 10, trim = 10)
beast_4 <- write_cp_to_file(vars = c("dist_process_prod"), 
                            maxcp = 20, timeinterval = 10, trim = 10)

# run changepoint detection on model 2 (dist_process_prod) without trim, 10sec
beast_1 <- write_cp_to_file(vars = c("dist_process_prod"), 
                            maxcp = 3, timeinterval = 10, trim = 0)
beast_2 <- write_cp_to_file(vars = c("dist_process_prod"), 
                            maxcp = 5, timeinterval = 10, trim = 0)
beast_3 <- write_cp_to_file(vars = c("dist_process_prod"), 
                            maxcp = 10, timeinterval = 10, trim = 0)
beast_4 <- write_cp_to_file(vars = c("dist_process_prod"), 
                            maxcp = 20, timeinterval = 10, trim = 0)

# run changepoint detection on model 2 (dist_process_prod) with trim, 5sec
beast_1 <- write_cp_to_file(vars = c("dist_process_prod"), 
                            maxcp = 3, timeinterval = 5, trim = 10)
beast_2 <- write_cp_to_file(vars = c("dist_process_prod"), 
                            maxcp = 5, timeinterval = 5, trim = 10)
beast_3 <- write_cp_to_file(vars = c("dist_process_prod"), 
                            maxcp = 10, timeinterval = 5, trim = 10)
beast_4 <- write_cp_to_file(vars = c("dist_process_prod"), 
                            maxcp = 20, timeinterval = 5, trim = 10)

# run changepoint detection on model 2 (dist_process_prod) without trim, 5sec
beast_1 <- write_cp_to_file(vars = c("dist_process_prod"), 
                            maxcp = 3, timeinterval = 5, trim = 0)
beast_2 <- write_cp_to_file(vars = c("dist_process_prod"), 
                            maxcp = 5, timeinterval = 5, trim = 0)
beast_3 <- write_cp_to_file(vars = c("dist_process_prod"), 
                            maxcp = 10, timeinterval = 5, trim = 0)
beast_4 <- write_cp_to_file(vars = c("dist_process_prod"), 
                            maxcp = 20, timeinterval = 5, trim = 0)

# run changepoint detection on model 2 (dist_process_prod) with trim, 1sec
beast_1 <- write_cp_to_file(vars = c("dist_process_prod"), 
                            maxcp = 3, timeinterval = 1, trim = 10)
beast_2 <- write_cp_to_file(vars = c("dist_process_prod"), 
                            maxcp = 5, timeinterval = 1, trim = 10)
beast_3 <- write_cp_to_file(vars = c("dist_process_prod"), 
                            maxcp = 10, timeinterval = 1, trim = 10)
beast_4 <- write_cp_to_file(vars = c("dist_process_prod"), 
                            maxcp = 20, timeinterval = 1, trim = 10)

# run changepoint detection on model 2 (dist_process_prod) without trim, 1sec
beast_1 <- write_cp_to_file(vars = c("dist_process_prod"), 
                            maxcp = 3, timeinterval = 1, trim = 0)
beast_2 <- write_cp_to_file(vars = c("dist_process_prod"), 
                            maxcp = 5, timeinterval = 1, trim = 0)
beast_3 <- write_cp_to_file(vars = c("dist_process_prod"), 
                            maxcp = 10, timeinterval = 1, trim = 0)
beast_4 <- write_cp_to_file(vars = c("dist_process_prod"), 
                            maxcp = 20, timeinterval = 1, trim = 0)

###########################
# 3) relative position mouse to start
# run changepoint detection on model 3 (rel_position) with trim, 10sec
beast_1 <- write_cp_to_file(vars = c("rel_position"), 
                            maxcp = 3, timeinterval = 10, trim = 10)
beast_2 <- write_cp_to_file(vars = c("rel_position"), 
                            maxcp = 5, timeinterval = 10, trim = 10)
beast_3 <- write_cp_to_file(vars = c("rel_position"), 
                            maxcp = 10, timeinterval = 10, trim = 10)
beast_4 <- write_cp_to_file(vars = c("rel_position"), 
                            maxcp = 20, timeinterval = 10, trim = 10)

# run changepoint detection on model 3 (rel_position) without trim, 10sec
beast_1 <- write_cp_to_file(vars = c("rel_position"), 
                            maxcp = 3, timeinterval = 10, trim = 0)
beast_2 <- write_cp_to_file(vars = c("rel_position"), 
                            maxcp = 5, timeinterval = 10, trim = 0)
beast_3 <- write_cp_to_file(vars = c("rel_position"), 
                            maxcp = 10, timeinterval = 10, trim = 0)
beast_4 <- write_cp_to_file(vars = c("rel_position"), 
                            maxcp = 20, timeinterval = 10, trim = 0)

# run changepoint detection on model 3 (rel_position) with trim, 5sec
beast_1 <- write_cp_to_file(vars = c("rel_position"), 
                            maxcp = 3, timeinterval = 5, trim = 10)
beast_2 <- write_cp_to_file(vars = c("rel_position"), 
                            maxcp = 5, timeinterval = 5, trim = 10)
beast_3 <- write_cp_to_file(vars = c("rel_position"), 
                            maxcp = 10, timeinterval = 5, trim = 10)
beast_4 <- write_cp_to_file(vars = c("rel_position"), 
                            maxcp = 20, timeinterval = 5, trim = 10)

# run changepoint detection on model 3 (rel_position) without trim, 5sec
beast_1 <- write_cp_to_file(vars = c("rel_position"), 
                            maxcp = 3, timeinterval = 5, trim = 0)
beast_2 <- write_cp_to_file(vars = c("rel_position"), 
                            maxcp = 5, timeinterval = 5, trim = 0)
beast_3 <- write_cp_to_file(vars = c("rel_position"), 
                            maxcp = 10, timeinterval = 5, trim = 0)
beast_4 <- write_cp_to_file(vars = c("rel_position"), 
                            maxcp = 20, timeinterval = 5, trim = 0)

# run changepoint detection on model 3 (rel_position) with trim, 1sec
beast_1 <- write_cp_to_file(vars = c("rel_position"), 
                            maxcp = 3, timeinterval = 1, trim = 10)
beast_2 <- write_cp_to_file(vars = c("rel_position"), 
                            maxcp = 5, timeinterval = 1, trim = 10)
beast_3 <- write_cp_to_file(vars = c("rel_position"), 
                            maxcp = 10, timeinterval = 1, trim = 10)
beast_4 <- write_cp_to_file(vars = c("rel_position"), 
                            maxcp = 20, timeinterval = 1, trim = 10)

# run changepoint detection on model 3 (rel_position) without trim, 1sec
beast_1 <- write_cp_to_file(vars = c("rel_position"), 
                            maxcp = 3, timeinterval = 1, trim = 0)
beast_2 <- write_cp_to_file(vars = c("rel_position"), 
                            maxcp = 5, timeinterval = 1, trim = 0)
beast_3 <- write_cp_to_file(vars = c("rel_position"), 
                            maxcp = 10, timeinterval = 1, trim = 0)
beast_4 <- write_cp_to_file(vars = c("rel_position"), 
                            maxcp = 20, timeinterval = 1, trim = 0)

###########################
# 4) source usage (focus) decreases
# run changepoint detection on model 4 (source_use) with trim, 10sec
beast_1 <- write_cp_to_file(vars = c("source_use"), 
                            maxcp = 3, timeinterval = 10, trim = 10)
beast_2 <- write_cp_to_file(vars = c("source_use"), 
                            maxcp = 5, timeinterval = 10, trim = 10)
beast_3 <- write_cp_to_file(vars = c("source_use"), 
                            maxcp = 10, timeinterval = 10, trim = 10)
beast_4 <- write_cp_to_file(vars = c("source_use"), 
                            maxcp = 20, timeinterval = 10, trim = 10)

# run changepoint detection on model 4 (source_use) without trim, 10sec
beast_1 <- write_cp_to_file(vars = c("source_use"), 
                            maxcp = 3, timeinterval = 10, trim = 0)
beast_2 <- write_cp_to_file(vars = c("source_use"), 
                            maxcp = 5, timeinterval = 10, trim = 0)
beast_3 <- write_cp_to_file(vars = c("source_use"), 
                            maxcp = 10, timeinterval = 10, trim = 0)
beast_4 <- write_cp_to_file(vars = c("source_use"), 
                            maxcp = 20, timeinterval = 10, trim = 0)

# run changepoint detection on model 4 (source_use) with trim, 5sec
beast_1 <- write_cp_to_file(vars = c("source_use"), 
                            maxcp = 3, timeinterval = 5, trim = 10)
beast_2 <- write_cp_to_file(vars = c("source_use"), 
                            maxcp = 5, timeinterval = 5, trim = 10)
beast_3 <- write_cp_to_file(vars = c("source_use"), 
                            maxcp = 10, timeinterval = 5, trim = 10)
beast_4 <- write_cp_to_file(vars = c("source_use"), 
                            maxcp = 20, timeinterval = 5, trim = 10)

# run changepoint detection on model 4 (source_use) without trim, 5sec
beast_1 <- write_cp_to_file(vars = c("source_use"), 
                            maxcp = 3, timeinterval = 5, trim = 0)
beast_2 <- write_cp_to_file(vars = c("source_use"), 
                            maxcp = 5, timeinterval = 5, trim = 0)
beast_3 <- write_cp_to_file(vars = c("source_use"), 
                            maxcp = 10, timeinterval = 5, trim = 0)
beast_4 <- write_cp_to_file(vars = c("source_use"), 
                            maxcp = 20, timeinterval = 5, trim = 0)

# run changepoint detection on model 4 (source_use) with trim, 1sec
beast_1 <- write_cp_to_file(vars = c("source_use"), 
                            maxcp = 3, timeinterval = 1, trim = 10)
beast_2 <- write_cp_to_file(vars = c("source_use"), 
                            maxcp = 5, timeinterval = 1, trim = 10)
beast_3 <- write_cp_to_file(vars = c("source_use"), 
                            maxcp = 10, timeinterval = 1, trim = 10)
beast_4 <- write_cp_to_file(vars = c("source_use"), 
                            maxcp = 20, timeinterval = 1, trim = 10)

# run changepoint detection on model 4 (source_use) without trim, 1sec
beast_1 <- write_cp_to_file(vars = c("source_use"), 
                            maxcp = 3, timeinterval = 1, trim = 0)
beast_2 <- write_cp_to_file(vars = c("source_use"), 
                            maxcp = 5, timeinterval = 1, trim = 0)
beast_3 <- write_cp_to_file(vars = c("source_use"), 
                            maxcp = 10, timeinterval = 1, trim = 0)
beast_4 <- write_cp_to_file(vars = c("source_use"), 
                            maxcp = 20, timeinterval = 1, trim = 0)

###########################
# 5) average pause length increases
# run changepoint detection on model 5 (m_pause_length) with trim, 10sec
beast_1 <- write_cp_to_file(vars = c("m_pause_length"), 
                            maxcp = 3, timeinterval = 10, trim = 10)
beast_2 <- write_cp_to_file(vars = c("m_pause_length"), 
                            maxcp = 5, timeinterval = 10, trim = 10)
beast_3 <- write_cp_to_file(vars = c("m_pause_length"), 
                            maxcp = 10, timeinterval = 10, trim = 10)
beast_4 <- write_cp_to_file(vars = c("m_pause_length"), 
                            maxcp = 20, timeinterval = 10, trim = 10)

# run changepoint detection on model 5 (m_pause_length) without trim, 10sec
beast_1 <- write_cp_to_file(vars = c("m_pause_length"), 
                            maxcp = 3, timeinterval = 10, trim = 0)
beast_2 <- write_cp_to_file(vars = c("m_pause_length"), 
                            maxcp = 5, timeinterval = 10, trim = 0)
beast_3 <- write_cp_to_file(vars = c("m_pause_length"), 
                            maxcp = 10, timeinterval = 10, trim = 0)
beast_4 <- write_cp_to_file(vars = c("m_pause_length"), 
                            maxcp = 20, timeinterval = 10, trim = 0)

# run changepoint detection on model 5 (m_pause_length) with trim, 5sec
beast_1 <- write_cp_to_file(vars = c("m_pause_length"), 
                            maxcp = 3, timeinterval = 5, trim = 10)
beast_2 <- write_cp_to_file(vars = c("m_pause_length"), 
                            maxcp = 5, timeinterval = 5, trim = 10)
beast_3 <- write_cp_to_file(vars = c("m_pause_length"), 
                            maxcp = 10, timeinterval = 5, trim = 10)
beast_4 <- write_cp_to_file(vars = c("m_pause_length"), 
                            maxcp = 20, timeinterval = 5, trim = 10)

# run changepoint detection on model 5 (m_pause_length) without trim, 5sec
beast_1 <- write_cp_to_file(vars = c("m_pause_length"), 
                            maxcp = 3, timeinterval = 5, trim = 0)
beast_2 <- write_cp_to_file(vars = c("m_pause_length"), 
                            maxcp = 5, timeinterval = 5, trim = 0)
beast_3 <- write_cp_to_file(vars = c("m_pause_length"), 
                            maxcp = 10, timeinterval = 5, trim = 0)
beast_4 <- write_cp_to_file(vars = c("m_pause_length"), 
                            maxcp = 20, timeinterval = 5, trim = 0)

# run changepoint detection on model 5 (m_pause_length) with trim, 1sec
beast_1 <- write_cp_to_file(vars = c("m_pause_length"), 
                            maxcp = 3, timeinterval = 1, trim = 10)
beast_2 <- write_cp_to_file(vars = c("m_pause_length"), 
                            maxcp = 5, timeinterval = 1, trim = 10)
beast_3 <- write_cp_to_file(vars = c("m_pause_length"), 
                            maxcp = 10, timeinterval = 1, trim = 10)
beast_4 <- write_cp_to_file(vars = c("m_pause_length"), 
                            maxcp = 20, timeinterval = 1, trim = 10)

# run changepoint detection on model 5 (m_pause_length) without trim, 1sec
beast_1 <- write_cp_to_file(vars = c("m_pause_length"), 
                            maxcp = 3, timeinterval = 1, trim = 0)
beast_2 <- write_cp_to_file(vars = c("m_pause_length"), 
                            maxcp = 5, timeinterval = 1, trim = 0)
beast_3 <- write_cp_to_file(vars = c("m_pause_length"), 
                            maxcp = 10, timeinterval = 1, trim = 0)
beast_4 <- write_cp_to_file(vars = c("m_pause_length"), 
                            maxcp = 20, timeinterval = 1, trim = 0)

##############################
# First 3 combined (exploratory)
beast_1 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 3, timeinterval = 10, trim = 10)
beast_2 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 5, timeinterval = 10, trim = 10)
beast_3 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 10, timeinterval = 10, trim = 10)
beast_4 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 20, timeinterval = 10, trim = 10)
beast_5 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 1, timeinterval = 10, trim = 10)


# run changepoint detection without trim, 10sec
beast_1 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 3, timeinterval = 10, trim = 0)
beast_2 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 5, timeinterval = 10, trim = 0)
beast_3 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 10, timeinterval = 10, trim = 0)
beast_4 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 20, timeinterval = 10, trim = 0)
beast_5 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 1, timeinterval = 10, trim = 0)

# run changepoint detection with trim, 5sec
beast_1 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 3, timeinterval = 5, trim = 10)
beast_2 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 5, timeinterval = 5, trim = 10)
beast_3 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 10, timeinterval = 5, trim = 10)
beast_4 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 20, timeinterval = 5, trim = 10)
beast_5 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 1, timeinterval = 5, trim = 10)



# run changepoint detection without trim, 5sec
beast_1 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 3, timeinterval = 5, trim = 0)
beast_2 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 5, timeinterval = 5, trim = 0)
beast_3 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 10, timeinterval = 5, trim = 0)
beast_4 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 20, timeinterval = 5, trim = 0)
beast_5 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 1, timeinterval = 5, trim = 0)

# run changepoint detection with trim, 1sec
beast_1 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 3, timeinterval = 1, trim = 10)
beast_2 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 5, timeinterval = 1, trim = 10)
beast_3 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 10, timeinterval = 1, trim = 10)
beast_4 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 20, timeinterval = 1, trim = 10)
beast_5 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 1, timeinterval = 1, trim = 10)

# run changepoint detection without trim, 1sec
beast_1 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 3, timeinterval = 1, trim = 0)
beast_2 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 5, timeinterval = 1, trim = 0)
beast_3 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 10, timeinterval = 1, trim = 0)
beast_4 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 20, timeinterval = 1, trim = 0)
beast_5 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position" ), 
                            maxcp = 1, timeinterval = 1, trim = 0)

##############################
# all 5 combined (exploratory)
beast_1 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 3, timeinterval = 10, trim = 10)
beast_2 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 5, timeinterval = 10, trim = 10)
beast_3 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 10, timeinterval = 10, trim = 10)
beast_4 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 20, timeinterval = 10, trim = 10)
beast_5 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 1, timeinterval = 10, trim = 10)


# run changepoint detection without trim, 10sec
beast_1 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 3, timeinterval = 10, trim = 0)
beast_2 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 5, timeinterval = 10, trim = 0)
beast_3 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 10, timeinterval = 10, trim = 0)
beast_4 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 20, timeinterval = 10, trim = 0)
beast_5 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 1, timeinterval = 10, trim = 0)

# run changepoint detection with trim, 5sec
beast_1 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 3, timeinterval = 5, trim = 10)
beast_2 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 5, timeinterval = 5, trim = 10)
beast_3 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 10, timeinterval = 5, trim = 10)
beast_4 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 20, timeinterval = 5, trim = 10)
beast_5 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 1, timeinterval = 5, trim = 10)

###FROM HERE
# run changepoint detection without trim, 5sec
beast_1 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 3, timeinterval = 5, trim = 0)
beast_2 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 5, timeinterval = 5, trim = 0)
beast_3 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 10, timeinterval = 5, trim = 0)
beast_4 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 20, timeinterval = 5, trim = 0)
beast_5 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 1, timeinterval = 5, trim = 0)

# run changepoint detection with trim, 1sec
beast_1 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 3, timeinterval = 1, trim = 10)
beast_2 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 5, timeinterval = 1, trim = 10)
beast_3 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 10, timeinterval = 1, trim = 10)
beast_4 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 20, timeinterval = 1, trim = 10)
beast_5 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 1, timeinterval = 1, trim = 10)

# run changepoint detection without trim, 1sec
beast_1 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 3, timeinterval = 1, trim = 0)
beast_2 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 5, timeinterval = 1, trim = 0)
beast_3 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 10, timeinterval = 1, trim = 0)
beast_4 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 20, timeinterval = 1, trim = 0)
beast_5 <- write_cp_to_file(vars = c("perc_doclength",  "dist_process_prod", "rel_position", "source_use", "m_pause_length" ), 
                            maxcp = 1, timeinterval = 1, trim = 0)





