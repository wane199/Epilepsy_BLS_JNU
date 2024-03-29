ui = bootstrapPage(fluidPage(theme = shinytheme("superhero"), # cyborg/darkly/slate/superhero
                             titlePanel('TLE Relapse Dynamic Nomogram from JNU'),
                             sidebarLayout(sidebarPanel(uiOutput('manySliders'),
                                                        checkboxInput('trans', 'Alpha blending (transparency)', value = TRUE),
                                                        actionButton('add', 'Predict'),
                                                        br(), br(), 
                                                        helpText('Press Quit to exit the application'),
                                                        actionButton('quit', 'Quit'),
                                                        br(),
                                                        img(src = "https://ts1.cn.mm.bing.net/th/id/R-C.c80600d38debc68a12b4b566886c8216?rik=bTkNEfTXK0fisg&riu=http%3a%2f%2fpicture.swwy.com%2fY2UzZDljYTQxNjhmNDI.jpg&ehk=WYS7zLiw1qw9kNUCW14LEMFnE2n0sOPMwjkmxBh71%2fs%3d&risl=&pid=ImgRaw&r=0&sres=1&sresct=1",
                                                            height = 120, width = 240),
                                                        img(src = "https://github.com/wane199/Epilepsy_BLS_JNU/blob/main/EP_Cox_Nomo/shiny/DynNomapp/999logo.png?raw=true",
                                                            height = 100, width = 240)
                             ),
                             mainPanel(tabsetPanel(id = 'tabs',
                                                   tabPanel('Relapse-free plot', plotOutput('plot')),
                                                   tabPanel('Predicted Relapse-free', plotlyOutput('plot2')),
                                                   tabPanel('Numerical Summary', verbatimTextOutput('data.pred')),
                                                   tabPanel('Model Summary', verbatimTextOutput('summary'))
                             )
                             )
                             )))
