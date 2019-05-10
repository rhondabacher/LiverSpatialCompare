    

source("global.R")
server <- function(input, output, session) {
    

  output$choose_category <- renderUI({
     selectInput("inputCat", "KEGG", as.list(allKEGG), selected = "KEGG: Biosynthesis of amino acids")
   })
   
   
    plotInput <- reactive({
        makeCorPlots(input$inputCat, .1, .90, input$sigThresh)
    })
  
  
 output$genePlot <- renderPlot({
        plotInput()
    }, height=600, width=1200)
    

    
  output$downloadPlot <- downloadHandler(
      filename = function() {
      paste0("categoryPlot_", input$inputCat, ".pdf")
      },
      content = function(file) {
          pdf(file, height=5, width=10)
           makeCorPlots(input$inputCat, .1, .90, input$sigThresh)
      dev.off()
      }
  )

   
}