#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(leaflet)

ui <- fluidPage(leafletOutput("mymap", height = 800))

server <- function(input, output, session) {
  
  rtp <- readRDS("./Data/rtp.RDS")
  
  m <- mapview(rtp,
               method = "ngb", 
               na.color = rgb(0, 0, 255, max = 255, alpha = 0),
               query.type = "click",
               trim = TRUE,
               legend = FALSE,
               map.types =  "Esri.WorldImagery",
               alpha.regions = 0,
               lwd=2,
               color="red")
  
  output$mymap <- renderLeaflet({m@map})
}

shinyApp(ui, server)
