shinyItemplot <- function(){

    ret <- list(

        ui = shiny::pageWithSidebar(

            # Application title
            shiny::headerPanel("Item plots in mirt"),

            shiny::sidebarPanel(

                shiny::h5('Select an internal mirt item class, the type of plot to display, whether the plot
                   inputs should be the classical or modern IRT parameterizations, and whether the plot
                   should be unidimensional or multidimensional. See ?mirt for more details.'),

                shiny::selectInput(inputId = "itemclass",
                            label = "Class of mirt item:",
                            choices = c('dich' = 'dich',
                                        'graded' = 'graded',
                                        'nominal' = 'nominal',
                                        'gpcm' = 'gpcm',
                                        'partcomp' = 'partcomp',
                                        'ideal' = 'ideal'),
                            selected = 'dich'),

                shiny::selectInput(inputId = "plottype",
                            label = "Type of plot to display:",
                            choices = c('trace', 'info', 'score', 'infocontour', 'SE', 'infoSE', 'tracecontour'),
                            selected = 'trace'),

                shiny::numericInput("theta_lim_low", "Theta lower range:", -6,
                             min = -35, max = 35),

                shiny::numericInput("theta_lim_high", "Theta upper range:", 6,
                             min = -35, max = 35),

                shiny::checkboxInput(inputId = "classical",
                              label = "Use traditional IRT parameterization inputs? \n
                              Only applicable to class dich, graded, gpcm, and nominal.",
                              value = FALSE),

                shiny::checkboxInput(inputId = "nfact",
                              label = "Multidimensional?",
                              value = FALSE),

                #-------------------------------------------------------------------------

                shiny::conditionalPanel(condition = "input.classical == false",

                                        shiny::sliderInput(inputId = "a1par",
                                             label = "a1 value:",
                                             min = -3, max = 3, value = 1, step = 0.2),

                                        shiny::conditionalPanel(condition = "input.nfact == true",
                                                                shiny::sliderInput(inputId = "a2par",
                                                              label = "a2 value:",
                                                              min = -3, max = 3, value = 1, step = 0.2)),

                                        shiny::conditionalPanel(condition = "input.itemclass == 'dich'",
                                                                shiny::sliderInput(inputId = "dpar",
                                                              label = "d value:",
                                                              min = -5, max = 5, value = 0, step = 0.25),
                                                              shiny::sliderInput(inputId = "gpar",
                                                              label = "g value:",
                                                              min = 0, max = 1, value = 0, step = 0.05),
                                                              shiny::sliderInput(inputId = "upar",
                                                              label = "u value:",
                                                              min = 0, max = 1, value = 1, step = 0.05)
                                 ),

                                 shiny::conditionalPanel(condition = "input.itemclass == 'ideal'",
                                                         shiny::sliderInput(inputId = "idpar",
                                                              label = "d value:",
                                                              min = -5, max = 5, value = 0, step = 0.25)
                                 ),

                                 shiny::conditionalPanel(condition = "input.itemclass != 'dich'",
                                                         shiny::conditionalPanel(condition = "input.itemclass == 'gpcm' ||
                                                                   input.itemclass == 'nominal'",
                                                                                 shiny::sliderInput(inputId = "d0par",
                                                                               label = "d0 value (default fixed at 0):",
                                                                               min = -5, max = 5, value = 0, step = 0.25)
                                                  ),

                                                  shiny::conditionalPanel(condition = "input.itemclass != 'ideal'",
                                                                          shiny::sliderInput(inputId = "d1par",
                                                                  label = "d1 value:",
                                                                  min = -5, max = 5, value = 1, step = 0.25),
                                                                  shiny::sliderInput(inputId = "d2par",
                                                                  label = "d2 value:",
                                                                  min = -5, max = 5, value = 0, step = 0.25),
                                                                  shiny::conditionalPanel(condition = "input.itemclass != 'partcomp'",
                                                                                          shiny::sliderInput(inputId = "d3par",
                                                                                   label = "d3 value:",
                                                                                   min = -5, max = 5, value = -1, step = 0.25)
                                                      )
                                                  )
                ),

                shiny::conditionalPanel(condition = "input.itemclass == 'nominal'",
                                        shiny::sliderInput(inputId = "ak0par",
                                                              label = "ak0 value (default fixed at 0):",
                                                              min = -3, max = 3, value = 0, step = 0.2),
                                        shiny::sliderInput(inputId = "ak1par",
                                                              label = "ak0 value:",
                                                              min = -3, max = 3, value = 1, step = 0.2),
                                        shiny::sliderInput(inputId = "ak2par",
                                                              label = "ak2 value:",
                                                              min = -3, max = 3, value = 2, step = 0.2),
                                        shiny::sliderInput(inputId = "ak3par",
                                                              label = "ak3 value (default fixed at (ncat-1)):",
                                                              min = -3, max = 3, value = 3, step = 0.2)
                                 )
                ),

                shiny::conditionalPanel(condition = "input.classical == true",

                                        shiny::conditionalPanel(condition = "input.itemclass == 'dich'",
                                                                shiny::sliderInput(inputId = "a1parc",
                                                              label = "a value:",
                                                              min = -3, max = 3, value = 1, step = 0.2),
                                                              shiny::sliderInput(inputId = "bpar",
                                                              label = "b value:",
                                                              min = -5, max = 5, value = 0, step = 0.25),
                                                              shiny::sliderInput(inputId = "gparc",
                                                              label = "g value:",
                                                              min = 0, max = 1, value = 0, step = 0.05),
                                                              shiny::sliderInput(inputId = "uparc",
                                                              label = "u value:",
                                                              min = 0, max = 1, value = 1, step = 0.05)
                                 ),

                                 shiny::conditionalPanel(condition = "input.itemclass == 'graded'",
                                                         shiny::sliderInput(inputId = "a1pard",
                                                              label = "a value:",
                                                              min = -3, max = 3, value = 1, step = 0.2),
                                                         shiny::sliderInput(inputId = "bpar1d",
                                                              label = "b1 value:",
                                                              min = -5, max = 5, value = -1, step = 0.25),
                                                         shiny::sliderInput(inputId = "bpar2d",
                                                              label = "b2 value:",
                                                              min = -5, max = 5, value = 0, step = 0.25),
                                                         shiny::sliderInput(inputId = "bpar3d",
                                                              label = "b3 value:",
                                                              min = -5, max = 5, value = 1, step = 0.25)
                                 ),

                                 shiny::conditionalPanel(condition = "input.itemclass == 'gpcm'",
                                                         shiny::sliderInput(inputId = "a1pare",
                                                              label = "a value:",
                                                              min = -3, max = 3, value = 1, step = 0.2),
                                                         shiny::sliderInput(inputId = "bpar1e",
                                                              label = "b1 value:",
                                                              min = -5, max = 5, value = -1, step = 0.25),
                                                         shiny::sliderInput(inputId = "bpar2e",
                                                              label = "b2 value:",
                                                              min = -5, max = 5, value = 0, step = 0.25),
                                                         shiny::sliderInput(inputId = "bpar3e",
                                                              label = "b3 value:",
                                                              min = -5, max = 5, value = 1, step = 0.25)
                                 ),

                                 shiny::conditionalPanel(condition = "input.itemclass == 'nominal'",
                                                         shiny::sliderInput(inputId = "a1parf",
                                                              label = "a1 value:",
                                                              min = -3, max = 3, value = -1.4, step = 0.2),
                                                         shiny::sliderInput(inputId = "a2parf",
                                                              label = "a2 value:",
                                                              min = -3, max = 3, value = -0.4, step = 0.2),
                                                         shiny::sliderInput(inputId = "a3parf",
                                                              label = "a3 value:",
                                                              min = -3, max = 3, value = 0.4, step = 0.2),
                                                         shiny::sliderInput(inputId = "a4parf",
                                                              label = "a4 value:",
                                                              min = -3, max = 3, value = 1.4, step = 0.2),
                                                         shiny::sliderInput(inputId = "bpar1f",
                                                              label = "b1 value:",
                                                              min = -5, max = 5, value = -1, step = 0.25),
                                                         shiny::sliderInput(inputId = "bpar2f",
                                                              label = "b2 value:",
                                                              min = -5, max = 5, value = 0, step = 0.25),
                                                         shiny::sliderInput(inputId = "bpar3f",
                                                              label = "b3 value:",
                                                              min = -5, max = 5, value = 1, step = 0.25),
                                                         shiny::sliderInput(inputId = "bpar4f",
                                                              label = "b4 value:",
                                                              min = -5, max = 5, value = 1.5, step = 0.25)
                                 )
                )
                ),

            shiny::mainPanel(
                shiny::verbatimTextOutput("coefs"),
                shiny::plotOutput(outputId = "main_plot", height = "700px", width = "700px")
            )

                ),

        server = function(input, output) {

            genmod <- function(input){
                set.seed(1234)
                itemclass <- c(input$itemclass, input$itemclass)
                itemtype <- switch(input$itemclass,
                                   dich='2PL',
                                   graded='graded',
                                   gpcm='gpcm',
                                   nominal='nominal',
                                   nestlogit='2PLNRM',
                                   partcomp='PC2PL',
                                   nestlogit='2PLNRM',
                                   ideal='ideal')
                nominal <- NULL
                model <- 1
                if(input$nfact) model <- 2
                if(model == 2 && input$plottype == 'infoSE')
                    stop('infoSE only available for single dimensional models')
                a <- matrix(1,2)
                d <- matrix(0,2)
                if(input$itemclass == 'graded'){
                    d <- matrix(c(1,0,-1), 2, 3, byrow=TRUE)
                } else if(input$itemclass == 'gpcm'){
                    d <- matrix(c(0,1,0,-1), 2, 4, byrow=TRUE)
                } else if(input$itemclass == 'nominal'){
                    nominal <- matrix(c(0,1,2,3), 2, 4, byrow=TRUE)
                    d <- matrix(c(0,1,0,-1), 2, 4, byrow=TRUE)
                } else if(input$itemclass == 'nestlogit'){
                    nominal <- matrix(c(0,1,2), 2, 3, byrow=TRUE)
                    d <- matrix(c(0,0,1,-1), 2, 4, byrow=TRUE)
                } else if(input$itemclass == 'partcomp'){
                    if(model != 2) stop('partcomp models require more than 1 dimension')
                    if(input$plottype == 'info' || input$plottype == 'infocontour')
                        stop('information based plots not currently supported for partcomp items')
                    a <- matrix(c(1,1), 2, 2, byrow=TRUE)
                    d <- matrix(c(1,1,1,NA), 2, 2, byrow=TRUE)
                    itemtype[2] <- '2PL'
                    itemclass[2] <- 'dich'
                    model <- mirt.model('F1 = 1,2
                                        F2 = 1', quiet=TRUE)
                }
                dat <- simdata(a=a, d=d, N=100,
                               itemtype=itemclass, nominal=nominal)
                sv <- mirt(dat, model, itemtype=itemtype, pars = 'values', key=c(1, NA),
                           technical=list(message=FALSE))
                sv$est <- FALSE
                mod <- mirt(dat, model, itemtype=itemtype, pars=sv, key=c(1, NA),
                            technical=list(message=FALSE))
                par <- mod@pars[[1]]@par
                if(input$classical){
                    if(itemclass[1L] == 'dich'){
                        par <- c(input$a1parc, input$bpar, logit(input$gparc), logit(input$uparc))
                    } else if(itemclass[1L] == 'graded'){
                        par <- c(input$a1pard, input$bpar1d, input$bpar2d, input$bpar3d)
                    } else if(itemclass[1L] == 'gpcm'){
                        par <- c(input$a1pare, input$bpar1e, input$bpar2e, input$bpar3e)
                    } else if(itemclass[1L] == 'nominal'){
                        par <- c(input$a1parf, input$a2parf, input$a3parf, input$a4parf,
                                 input$bpar1f, input$bpar2f, input$bpar3f, input$bpar4f)
                    } else {
                        stop('Classical parameterization not available for selected item class')
                    }
                    par <- traditional2mirt(x=par, cls=itemclass[1L], ncat=mod@pars[[1]]@ncat)
                } else {
                    par[names(par) == 'a1'] <- input$a1par
                    par[names(par) == 'a2'] <- input$a2par
                    par[names(par) == 'd'] <- input$dpar
                    par[names(par) == 'g'] <- logit(input$gpar)
                    par[names(par) == 'u'] <- logit(input$upar)
                    par[names(par) == 'd0'] <- input$d0par
                    par[names(par) == 'd1'] <- input$d1par
                    par[names(par) == 'd2'] <- input$d2par
                    par[names(par) == 'd3'] <- input$d3par
                    par[names(par) == 'ak0'] <- input$ak0par
                    par[names(par) == 'ak1'] <- input$ak1par
                    par[names(par) == 'ak2'] <- input$ak2par
                    par[names(par) == 'ak3'] <- input$ak3par
                }
                mod@pars[[1]]@par <- par
                mod
                }

            output$main_plot <- shiny::renderPlot({
                mod <- genmod(input)
                print(itemplot(mod, 1, type=input$plottype, rotate = 'none',
                               theta_lim=c(input$theta_lim_low, input$theta_lim_high)))
            })

            output$coefs <- shiny::renderPrint({
                mod <- genmod(input)
                cat('Item parameters: \n\n')
                print(coef(mod, rotate = 'none')[[1L]])
                if(mod@nfact == 1L && !is(mod@pars[[1L]], 'nestlogit')){
                    cat('\n\nItem parameters (traditional IRT metric): \n\n')
                    print(coef(mod, IRTpars = TRUE)[[1L]])
                }
            })

        }
    )

    return(ret)
}
