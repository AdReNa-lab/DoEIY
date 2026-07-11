# Resources and References Shiny Module

resources_ui <- function(id) {
  ns <- NS(id)
  tagList(
    box(
      id = NULL, title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE, headerBorder = FALSE, width = 12,
      tabBox(
        id = ns("Resources_tabs_id"), title = NULL, side = "left", width = "100%",
        tabPanel("DoEIY Manual",
                 tags$iframe(style="height:600px; width:100%", src="User_Guide_v1.pdf")
        ),
        tabPanel("References",
                 p(HTML("
                   <p style='margin: 0 0 16px 0;'>This app's design generators were build on foundational work by the researchers and package authors listed below.</p>
                   <ul>
                   <li><strong>Box-Behnken</strong><br />Box, G. E. P., &amp; Behnken, D. W. (1960). Some New Three Level Designs for the Study of Quantitative Variables. <em>Technometrics</em>, 2(4), 455-475. DOI: <a href='https://www.tandfonline.com/
                                                              doi/abs/10.1080/00401706.1960.10489912'> 10.1080/00401706.1960.10489912</a>.</li>
                   <li><strong>Central Composite</strong><br />Lenth, R. V. (2008). <code>rsm</code>: Response-Surface Analysis (R package). DOI: <a href='https://doi.org/10.32614/
                                                              CRAN.package.rsm'>10.32614/CRAN.package. rsm</a>.</li>
                   <li><strong>D-Optimal</strong><br />AlgDesign: Algorithmic Experimental Design (R package). DOI: <a href='https://doi.org/10.32614/
                                                              CRAN.package.AlgDesign'>10.32614/ CRAN.package.AlgDesign</a>.</li>
                   <li><strong>Fractional Factorial</strong> <br />Gr&ouml;mping, U. (2025). <code>FrF2</code>: Fractional Factorial Designs with 2-Level Factors (R package). DOI: <a href='https://doi.org/10.32614/
                                                              CRAN.package.FrF2'>10.32614/ CRAN.package.FrF2</a>.</li>
                   <li><strong>Latin Hypercube</strong> <br />Carnell, R. (2024). <code>lhs</code> : Latin Hypercube Samples (R package). DOI: <a href='https://doi.org/10.32614/
                                                              CRAN.package.lhs'>10.32614/ CRAN.package.lhs</a>.</li>
                   <li><strong>Plackett&ndash;Burman</strong> <br />Plackett, R. L., &amp; Burman, J. P. (1946). The Design of Optimum Multifactorial Experiments. <em>Biometrika</em>, 33(4), 305-325. DOI: <a href='https://
                                                              academic.oup.com/biomet/article-abstract/33/
                                                              4/305/225377?redirectedFrom=fulltext'>10.1093/ biomet/33.4.305</a>.</li>
                   </ul>
                 "))
        )
      )
    )
  )
}

resources_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    # Stateless module
  })
}
