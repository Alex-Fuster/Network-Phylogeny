# Script to generate web page to graphically see the model equation with different parameters value
# Tutorial : https://shiny.rstudio.com/tutorial/written-tutorial/lesson1/

library(shiny)
#runExample("01_hello")


# Define UI for app that draws a histogram ----
ui <- fluidPage(

  # App title ----
  titlePanel("Macroevolution model parameters"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: Slider for the number of bins ----

			# Positive interaction parameters

      HTML("<font size='4'><b>Positive interaction parameters<br></b></font>"),
			HTML("<font size='4'><b>Establishment probability<br></b></font>"),

			sliderInput(inputId = "a_u",
                  label = "a_u - Shape of the exponential decay of the colonization - interaction relationship",
                  min = 0,
                  max = 10,
                  step = 0.01,
                  value = 0.45),

			sliderInput(inputId = "u_0",
                  label = "u_0 - Asymptotic establishment probability with infinite facilitation interactions",
                  min = 0,
                  max = 1,
                  step = 0.01,
                  value = 1),

			sliderInput(inputId = "u_1",
                  label = "u_1 - Establishment probability with absence of facilitation interactions",
                  min = -10,
                  max = 0,
                  step = 0.01,
                  value = -1),

			HTML("<font size='4'><b>Extinction probability<br></b></font>"),

		  sliderInput(inputId = "a_epos",
					        label = "a_epos - Shape of the exponential decay of the positive extinction - interaction relationship",
					        min = 0,
					        max = 5,
					        step = 0.01,
					        value = 1.2),

			sliderInput(inputId = "e_1pos",
						      label = "Positive interactions - 1pos",
						      min = 0,
						      max = 10,
						      step = 0.01,
						      value = 5.19),

      sliderInput(inputId = "e_0pos",
                  label = "e_0pos - Asymptotic extinction probability with infinite positive interactions",
                  min = 0,
                  max = 1,
									step = 0.005,
                  value = 0.075),

			# Negative interaction parameters

			HTML('<br><hr style="width:100%; size=3; border-top:1px solid #60605a;"></hr>'),
      HTML("<font size='4'><b>Negative interaction parameters<br></b></font>"),
			HTML("<font size='4'><b>Establishment probability<br></b></font>"),

			sliderInput(inputId = "a_uneg",
			            label = "a_uneg - Shape of the exponential decay of the colonization - interaction relationship",
			            min = 0,
			            max = 5,
			            step = 0.005,
			            value = 0.075),

			sliderInput(inputId = "u_0neg",
			            label = "u_0neg - Asymptotic establishment probability with infinite competition interactions",
			            min = 0,
			            max = 1,
			            step = 0.005,
			            value = 0.075),

			sliderInput(inputId = "u_1neg",
			            label = "u_1neg - Establishment probability with absence of competition interaction",
			            min = 0,
			            max = 10,
			            step = 0.01,
			            value = 2),

			HTML("<font size='4'><b>Extinction probability<br></b></font>"),

			sliderInput(inputId = "a_eneg",
                  label = "a_eneg - Shape of the exponential decay of the negative extinction - interaction relationship",
                  min = 0,
                  max = 5,
                  step = 0.005,
                  value = 0.025),

			sliderInput(inputId = "e_0neg",
                  label = "e_0neg - Asymptotic extinction probability with infinite negative interactions",
                  min = 0,
                  max = 1,
                  step = 0.01,
                  value = 0.5),

			sliderInput(inputId = "e_1neg",
                  label = "e_1neg - Extinction probability with absence of interactions",
                  min = 0,
                  max = 10,
                  step = 0.01,
                  value = 1),

			HTML('<br><hr style="width:100%; size=3; border-top:1px solid #60605a;"></hr>'),
      HTML("<font size='4'><b>Curve of niche trait inheritance<br></b></font>"),

			sliderInput(inputId = "curve_width",
                  label = "Facteur jouant sur le SD de la courbe normale",
                  min = 0,
                  max = 2,
                  step = 0.01,
                  value = 2)
    ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Histogram ----
      plotOutput(outputId = "distPlot_positive", height = "800px"),
      plotOutput(outputId = "distPlot_negative", height = "800px"),
			plotOutput(outputId = "sd_Plot")
    )
  )
)


# Define server logic required to draw a histogram ----
server <- function(input, output) {

  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot

  output$distPlot_positive <- renderPlot({

      # Parameters
      pars = list()

      # General parameters

      # Extinction probability, positive interactions communities
      pars$a_epos = input$a_epos  # Shape of the exponential decay of the positive extinction - interaction relationship
      pars$e_0pos = input$e_0pos # Asymptotic extinction probability with infinite positive interactions
      pars$e_1pos = input$e_1pos # 1 - pars$e_0pos

      # Establishment probability, positive interactions communities
      pars$a_u = input$a_u # Shape of the exponential decay of the colonization - interaction relationship
      pars$u_0 = input$u_0 # Asymptotic establishment probability with infinite facilitation interactions
      pars$u_1 = input$u_1 # Establishment probability with absence of facilitation interactions
      pars$d = 0.5 # Decrease speed of the establishment probability

      # Plot settings
      spe = matrix(0,1,60)
      for(i in 1:60){
        spe[1,i] <- (pars$u_0 + pars$u_1 * exp(-pars$a_u * i))#/(1 + exp(pars$d * (i - pars$I_max)))
      }

      ext = matrix(0,1,60)
      for(i in 1:60){
        ext[1,i] <- pars$e_0pos + pars$e_1pos*exp(-pars$a_epos* i)
      }


    	# Probiblities for facilitation interactions
			plot(spe[1,], type = "l", xlab = "Number of interaction", ylim = c(0,1), ylab = "Probability", frame.plot = FALSE, main = "Facilitation", xlim = c(0,50))
			lines(ext[1,], col = rgb(0, 158, 115, max = 255))
			legend("topright",
						legend = c("establishment", "extinction"),
						col = c("black", rgb(0, 158, 115, max = 255)),
						bty = "n",
						lty = 1
			)


    })

    output$distPlot_negative <- renderPlot({

      # Parameters
      pars = list()

			# General parameters

      # Extinction probability, negative interactions communities
      pars$a_eneg = input$a_eneg # Shape of the exponential decay of the negative extinction - interaction relationship
      pars$e_0neg = input$e_0neg # Asymptotic extinction probability with infinite negative interactions
      pars$e_1neg = input$e_1neg  # Extinction probability with absence of interactions

      # Establishment probability, negative interactions communities
      pars$a_uneg = input$a_uneg # Shape of the exponential decay of the colonization - interaction relationship
      pars$u_0neg = input$u_0neg # Asymptotic establishment probability with infinite competition interactions
      pars$u_1neg = input$u_1neg # Establishment probability with absence of competition interaction

      # Plot settings
      spe_neg = matrix(0,1,60)
      for(i in 1:60){
        spe_neg[1,i] <- pars$u_0neg + pars$u_1neg*exp(-pars$a_uneg * i)
      }

      extneg = matrix(0,1,60)
      for(i in 1:60){
        extneg[1,i] <- pars$e_0neg * (1 - exp(-pars$a_eneg * i))
      }


			# Probiblities for competition interactions
			plot(spe_neg[1,], type = "l", xlab = "Number of interaction", ylim = c(0,1),
						ylab = "Probability", frame.plot = FALSE, main = "Competition",
						xlim = c(0,50))
			lines(extneg[1,], col = rgb(0, 158, 115, max = 255))
			legend("topright",
				legend = c("establishment", "extinction"),
				col = c("black", rgb(0, 158, 115, max = 255)),
				bty = "n",
				lty = 1
			)

		})

		output$sd_Plot <- renderPlot({

			# Parameters
    	pars = list()

			pars$av_r = 0.1 # Half range of the niche of the first species
			pars$sd = input$curve_width*pars$av_r + 0.0001 # Standard deviation of the normal distribution used to calculate the niche optimum trait of a new species

			n <- 0.5

			# Plot
			curve(dnorm(x, n, pars$sd, log = FALSE), xlab = "trait n of the new species", ylab = "PDF")
			legend("topleft",
					legend = c(paste0("ancestor n = ", round(n, digit = 4)), paste0("sd = ", pars$sd)),
					bty = "n"
					)
			abline(v=n, col=rgb(0, 158, 115, max = 255))
		})

}

#source("code/simulations/parameters.R") ; shinyApp(ui = ui, server = server)
shinyApp(ui = ui, server = server)
