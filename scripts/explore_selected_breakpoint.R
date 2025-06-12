# Explore selected change point
library(tidyverse)
library(ggplot2)
library(timetk)
library(ggdist)
library(caret)
library(rpart.plot)
library(MLmetrics)
theme_set(theme_ggdist())

# read in final break points
beast_123_pick <- read_csv("output/beast_123_pick.csv", 
                           show_col_types = FALSE)
beast_pick <- read_csv("output/beast_pick.csv", show_col_types = FALSE)
log_filt2 <- read_csv("app/input/log_filt2.csv")
log_time_sec <- read_csv("app/input/log_time_5sec.csv")


################################################################################
# check accuracy selected change points
# per dataset / condition (not reported in paper)
beast_pick_sum <- beast_pick %>%
  filter(!is.na(selected_cp3)) %>%
  group_by(session) %>%
  summarize(
    m_seconds_off = mean(abs(dif), na.rm = T),
    meadian_seconds_off = median(abs(dif), na.rm = T),
    sd_seconds_off = sd(abs(dif), na.rm = T),
    perc_correct = sum(dif == 0, na.rm = T) / sum(!is.na(dif), na.rm = T),
    perc_almost_correct = sum(abs(dif) < 10, na.rm = T) / sum(!is.na(dif), na.rm = T),
    perc_almost60_correct = sum(abs(dif) < 60, na.rm = T) / sum(!is.na(dif), na.rm = T),
    time_min = mean(maxend)/ 60) 

# overall accuracy
total_beast_pick <- beast_pick %>%
  filter(!is.na(selected_cp3)) %>%
  group_by("Total") %>%
  summarize(
    m_seconds_off = mean(abs(dif), na.rm = T),
    meadian_seconds_off = median(abs(dif), na.rm = T),
    sd_seconds_off = sd(abs(dif), na.rm = T),
    perc_correct = sum(dif == 0, na.rm = T) / sum(!is.na(dif), na.rm = T),
    perc_almost_correct = sum(abs(dif) < 10, na.rm = T) / sum(!is.na(dif), na.rm = T),
    perc_almost60_correct = sum(abs(dif) < 60, na.rm = T) / sum(!is.na(dif), na.rm = T),
    time_min = mean(maxend)/ 60) 

# check with session characteristics
log_filt2_sum <- log_filt2 %>%
  group_by(session_id) %>%
  summarize(
    nchars = max(charProduction),
    nprod = max(doclengthFull),
    totaltime_min = (max(endTime)- min(startTime))/1000/60
  )

beast_pick_session <- beast_pick %>%
  filter(!is.na(rev_start)) %>%
  left_join(log_filt2_sum) %>%
  mutate(dif = abs(dif),
         correct = dif == 0)

# test session characteristics
t.test(nchars ~ correct, data = beast_pick_session)
t.test(totaltime_min ~ correct, data = beast_pick_session)
t.test(nprod ~ correct, data = beast_pick_session)

# test how often the annotation is really far away
sum(beast_pick_session$dif >60)/nrow(beast_pick_session)
sum(beast_pick_session$dif >300)/nrow(beast_pick_session)


################################################################################
# plot worst/best breakpoints per ppt
source("app/scripts/plot_changepoints_functions.R")


worst_sessions <- beast_pick %>%
  ungroup() %>%
  arrange(desc((dif))) %>%
  slice(1:12) 

best_sessions <- beast_pick %>%
  ungroup() %>%
  arrange((abs(dif))) %>%
  slice(1:12) 

# plot worst/best sessions
plot_final_breakpoint(log_filt2, subsetsessions = best_sessions$session_id)
plot_final_breakpoint(log_filt2, subsetsessions = worst_sessions$session_id)

# specific selection of worst sessions for paper
plot_final_breakpoint(log_filt2, subsetsessions = c(55,82,77,135))



################################################################################
# look into sessions without manual annotated change point

selected_cps <- beast_123_pick %>%
 full_join(beast_pick, by = join_by(session, session_id, rev_start)) 

no_cp <- selected_cps %>%
  filter(is.na(rev_start))

# no sessions
no_cp2 <- log_time_sec %>%
  group_by(session_id) %>%
  summarize(rev_start = sum(rev_start)) %>%
  filter(rev_start == 0)

# session 11 no change point found(!)

# summarize no change points
sum_no_cp <- no_cp %>%
  group_by(session, session_id) %>%
  summarize(n_cp = n(),
            model1 = sum(productline_flat1 != 0),
            model2 = sum(largest_procprod_diff2 != 0),
            model3 = sum(largest_dec_cursor != 0),
            overlap2 = sum(overlap ==2),
            overlap3 = sum(overlap ==3),
            cp_prob = first(cp_prob[cp_mean == selected_cp3]))

mean(sum_no_cp$n_cp)
sd(sum_no_cp$n_cp)

## Rule based on beast 123 pick to identify rules for why no cp
beast_123_sum <- selected_cps %>%
  group_by(session_id) %>%
  summarize(
    n_cp = ifelse(is.na(first(cp_mean)), 0, n()),
    m_cploc = mean(cp_mean/maxend),
    sd_cploc = ifelse(n_cp == 1, 0, sd(cp_mean/maxend)),
    model1 = sum(productline_flat1 != 0),
    model2 = sum(largest_procprod_diff2 != 0),
    model3 = sum(largest_dec_cursor != 0),
    overlap2 = sum(overlap ==2),
    overlap3 = sum(overlap ==3),
    cp_prob = mean(cp_prob),
    perc_doclength = mean(perc_doclength),
    rev_present = first(rev_present)
  ) %>%
  mutate(
    rev_present = as.factor(ifelse(rev_present == 0, "NoSecondPhase", 
                                  "SecondPhase"))
  ) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .))) %>%
  select(-session_id)


# test Machine Learning Models

set.seed(42)
train <- caret::createDataPartition(beast_123_sum$rev_present, p = 0.7, 
                                    list = FALSE)
data.trn <- beast_123_sum[train,]
data.tst <- beast_123_sum[-train,]

# fix unbalanced dataset using SMOTE
library(performanceEstimation)
data.trn.smote <- smote(rev_present ~. ,data.trn, perc.under = 20)

# normal decision tree
beast.tree = train(rev_present ~ ., 
                   data=data.trn.smote, 
                   method="rpart", # decision tree
                   na.action  = na.pass,
                   metric = "F",
                   trControl = trainControl(method = "cv", number = 10,
                                            summaryFunction = prSummary,
                                            classProbs = TRUE),
                   tuneLength = 30)

beast.tree
plot(beast.tree$finalModel, uniform=TRUE,
     main="Classification Tree")
text(beast.tree$finalModel,  all=TRUE, cex=.7)
print(varImp(beast.tree))
plot(varImp(beast.tree))


rpart.plot(beast.tree$finalModel, fallen.leaves = FALSE)

preds <- predict(beast.tree, data.tst)
confusionMatrix(preds, data.tst$rev_present, mode = "everything")

# random forest model
set.seed(42)
beast.rf = train(rev_present ~ ., 
                   data=data.trn.smote, 
                   method="rf", # decision tree
                   na.action  = na.pass,
                   metric = "F",
                   trControl = trainControl(method = "cv", number = 10,
                                            summaryFunction = prSummary,
                                            classProbs = TRUE))

plot(beast.rf)
beast.rf$results
print(varImp(beast.rf))
plot(varImp(beast.rf))
predsrf <- predict(beast.rf, data.tst)
confusionMatrix(predsrf, data.tst$rev_present, mode = "everything")

# still not a very good model
