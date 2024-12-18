# create shiny dashboard
library(shiny)
library(tidyverse)
library(timetk)
library(shinyWidgets)
library(shinythemes)
library(shinyBS)
library(ggdist)
theme_set(theme_ggdist())

# get data
models = c("1", "2", "3", "4", "5", "123", "12345")
log_filt2 <- read_csv("input/log_filt2.csv", show_col_types = FALSE)

# load author names
source("initialize_anom.R")

# Define UI for the application
ui <- 
  shinyUI(
    navbarPage(title = "Change Point Detection",
               header = tags$head(
                 tags$meta(name="author", content=author),
                 tags$meta(name="creation_date", content=Sys.Date())
               ),
               tabPanel("Detect change points",
                    sidebarPanel(
                        radioButtons("trim", "Trimming:",
                                     c("No trimming" = 0,
                                       "Trimming first 10 percent" = 10), 
                                     selected = 0),
                        radioButtons("preprocessing", "Log file summarized by:",
                                     c("Every 1 second" = "time1",
                                       "Every 5 seconds" = "time5",
                                       "Every 10 seconds" = "time10")),
                        selectInput("model", "Annotation rule:",
                                    choices = models,  selected = "1"),
                        textInput("participant", 
                                  HTML(paste0("Participant number(s) [1,",
                                              max(log_filt2$session_id), 
                                              "] <br/>  (comma delimited):" )), 
                                  value = "1,2,3,4", 
                                  placeholder = "Type participant number" ),
                        shinyWidgets::sliderTextInput("nof_cp",
                                              "Max. number of change points:",
                                              choices=c(1, 3, 5, 10, 20),
                                              selected=5, grid = T),
                        sliderInput("mincp_prob", 
                                    "Minimum probability of change points:",
                                    min = 0, max = 1,  value = 0),
                        bsTooltip("mincp_prob", 
                                  "Provide the minimum probability for each of the change points.",
                                  placement = "bottom", trigger = "hover",
                                  options = NULL)),
                    mainPanel(plotOutput("plot", height = 800)),
                    ),
               tabPanel("Accuracy change point detection",
                    sidebarPanel(
                      radioButtons("varplot", "Plot accuracy variable:",
                                   c("Correct change point" = "correct_cp",
                                     "Correct change point within 10sec" = 
                                       "close_cp",
                                     "Correct change point within 95% CI" = 
                                       "cp_in_ci",
                                     "Mean difference change point (sec)" =
                                       "mean_diff",
                                     "Median difference change point (sec)" =
                                       "median_diff",
                                     "Mean adjusted R2 (Note: not available for multivariate models)" = "mean_adjr2")),
                      selectInput("model2", "Annotation rule:",
                                  choices = models,  selected = models,
                                  multiple = T),
                              ),
                    mainPanel( plotOutput("plot_sum",  height = 500))
                    ),
               tabPanel("Select change point",
                        sidebarPanel(
                          radioButtons("trim3", "Trimming:",
                                       c("No trimming" = 0), 
                                       selected = 0),
                          radioButtons("preprocessing3", "Log file summarized by:",
                                       c("Every 5 seconds" = "time5")),
                          textInput("participant3", 
                                    HTML(paste0("Participant number(s) [1,",
                                                max(log_filt2$session_id), 
                                                "] <br/>  (comma delimited):" )), 
                                    value = "1,2,3,4", 
                                    placeholder = "Type participant number" ),
                          radioButtons("nof_cp3",  "Max. number of change points:",
                                      c(10))
                          ),
                        
                        mainPanel(plotOutput("plot_select", height = 800)),
               ),
                tabPanel("About", 
                         fluidPage(
                  p(HTML(paste("By", "<strong>", author, "</strong>", " ~ ", 
                               "<strong>", Sys.Date(), "</strong>"))),
                  p("This app can be used to detect change points in keystroke log data."),
                  p(paste("For more details see", authors, "(under review). Phase to Phase: 
                    Towards an Automated Procedure to Identify Phases in 
                    Writing Processes Using Keystroke Data.")))
                  )
               )
    )
  

    

# Define server logic
server <- function(input, output, session) {
    
  output$plot <- renderPlot( { 
    subsetsessions <-  as.numeric(unlist(strsplit(input$participant,",")))
    timeinterval = as.numeric(gsub("time", "", input$preprocessing))
    
    source("scripts/plot_changepoints_functions.R")
    
    plot_time_breakpoints(log_filt2, model = input$model, 
                          maxcp = input$nof_cp, timeinterval = timeinterval,
                          trim = input$trim, subsetsessions = subsetsessions,
                          min_cpprob = input$mincp_prob)
  })
  output$plot_sum <- renderPlot( { 
    
    source("scripts/plot_changepoints_functions.R")
    plot_accuracy_breakpoints(preproc = "time", var= input$varplot,
                              model = input$model2)
  })
  output$plot_select<- renderPlot( { 
    subsetsessions <-  as.numeric(unlist(strsplit(input$participant3,",")))
    
    source("scripts/plot_changepoints_functions.R")
    beast_123_pick <- read_csv("output/beast_123_pick.csv", 
                               show_col_types = FALSE)
    beast_pick <- read_csv("output/beast_pick.csv", show_col_types = FALSE)
    
    plot_final_breakpoint(log_filt2, subsetsessions = subsetsessions,
                          beast_123_pick = beast_123_pick,
                          beast_pick = beast_pick)
  })
}

# Run the application
shinyApp(ui = ui, server = server)

