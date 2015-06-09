shinyUI(fluidPage(
  titlePanel("EMMA online - End-member modelling analysis online"),
  
  sidebarLayout(position = "left",
                sidebarPanel(tabsetPanel(
                               tabPanel(title = "Data",
                                        fileInput(inputId = "file.user", 
                                                  label = strong("Data set"),
                                                  accept="text/plain"),
                                        tags$hr(),
                                        checkboxInput('header', 'Header', TRUE),
                                        checkboxInput('ID', 'ID (col 1)', TRUE),
                                        checkboxInput('Classes', 'Classes (row 1)', TRUE),
                                        radioButtons('sep', 'Separator',
                                                     c(Space = " ",
                                                       Comma = ",",
                                                       Semicolon = ";",
                                                       Tab = "\t"),
                                                     " "),
                                        tags$hr(),
                                        radioButtons('show', 'Show', selected = "curves",
                                                     c(Curves = "curves",
                                                       Map = "map",
                                                       EMMA = "EMMA"))
                               ),
                               tabPanel(title = "EMMA",
                                        sliderInput("q", "Number of EM (q):", 
                                                    min = 2, max = 10, value = 2),
                                        sliderInput("lw", "Weight transf. limit (lw):", 
                                                    min = 0, max = 0.5, value = 0),
                                        selectInput("rotation", "Rotation type:", 
                                                    choices = c("Varimax", "other", "yet_other"))
                               )
                )
                ),
                
                mainPanel(tabsetPanel(tabPanel("Tabular", tableOutput("tabular")),
                                      tabPanel("Graphical", 
                                               plotOutput(outputId = "graphical"))
                )
                )
  )
))