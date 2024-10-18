Code to estimate and select change points in keystroke log data
For more details see Rianne Conijn, Alessandra Rossetti, Nina Vandermeulen, Luuk Van Waes (under review). Phase to Phase: Towards an Automated Procedure to Identify Phases in Writing Processes Using Keystroke Data.

[Data / models not included - please contact the first author for more information]

#create_input_tables.R
- Runs create_timelog() and create_perclog() from transform_input_functions.R
- Creates aggregated keystroke logs, per x seconds or per x percent of total time

#rbeast_detect.R
- Runs estimate_changepoints_function.R
- Apply of function to multiple vars, CPs, timeintervals, trimming
- Expects a pre-processed dataset, where rows - equal spaced in time (can be done via create_input_tables.R)

#explore_breakpoints.R
- Runs on BEAST models (output of rbeast_detect.R)
- plots all changepoints using plot_changepoints_functions.R
- Calculates the distance between model changepoint and manual annotated changepoint (for closest changepoint)

#rbeast_select.R
- Runs on BEAST models (output of rbeast_detect.R)
- picks best fitting changepoint for revision phase
- outputs performance (distance from manual annotated changepoint)
- plots best & worst predicted sessions

#app.R
- runs writing phases dashboard app 
