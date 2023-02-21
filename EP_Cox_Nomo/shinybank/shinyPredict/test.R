library(shinydashboard)

ui <- dashboardPage(
  dashboardHeader(
    title = "Basic dashboard", dropdownMenuOutput("msgOutput"),
    dropdownMenu(
      type = "message",
      messageItem(from = "Finance update", message = "We are on threshold")
    )
  ),
  dashboardSidebar(
    sliderInput("bins", "Numbers of Breaks", 1, 100, 50),
    sidebarMenu(
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
        menuSubItem("Dashboard Finance", tabName = "finance", icon = icon("list-alt")),
        menuSubItem("Dashboard Sales", tabName = "sales", icon = icon("shop")),
      menuItem("Detailed Analysis"),
      menuItem("Raw Data")
    )
  ),
  dashboardBody(
    # Boxes need to be put in a row (or column)
    fluidRow(
      box(plotOutput("plot", height = 250)),
      box(
        title = "Controls",
        sliderInput("slider", "Number of observations:", 1, 100, 50)
      )
    )
  )
)

server <- function(input, output) {
  set.seed(123)
  histdata <- rnorm(1500)

  output$plot <- renderPlot({
    data <- histdata[seq_len(input$slider)]
    hist(data)
  })
}

shinyApp(ui, server)
