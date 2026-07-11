# Analyze the Design Shiny Module

analyze_design_ui <- function(id) {
  ns <- NS(id)
  tagList(
    box(
      id = NULL, title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE, headerBorder = FALSE, width = 12,
      tabBox(
        id = ns("Analysis_results_tabs_id"), title = NULL, side = "left", width = "100%",
        tabPanel(
          "Create Model",
          p("Note: Load a design for analysis on the Enter/Edit Results tab."),
          uiOutput(ns("analysis_table_ui_id"))
        ),
        tabPanel(
          "Model Fitting Results",
          uiOutput(ns("ANOVA_output_ui_id"))
        ),
        tabPanel(
          "Design Diagnostics",
          uiOutput(ns("Design_diagnostics_output_ui_id"))
        )
      )
    )
  )
}

analyze_design_server <- function(id, design_reactive, factor_types_reactive, 
                                  factor_names_reactive, factor_blocks_reactive, responses_reactive,
                                  design_matrix_reactive, model_formula_reactive, model_aov_reactive) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Internal reactive value for manual model builder
    selected_factors_for_model_reactive <- reactiveVal(value = NULL)
    
    # Triggered when design is created, loaded, or changed
    observeEvent(design_reactive(), {
      if (!(is.null(design_reactive()))) {
        design <- design_reactive()
        factor_types <- factor_types_reactive()
        
        types_vec <- as.character(factor_types[1, ])
        names(types_vec) <- colnames(factor_types)
        
        factor_names <- names(types_vec)[sapply(types_vec, function(x) !x %in% c("Response", "Block", "Metadata"))]
        factor_blocks <- names(types_vec)[sapply(types_vec, function(x) x == "Block")]
        responses <- names(types_vec)[sapply(types_vec, function(x) x == "Response")]
        
        factor_names_reactive(factor_names)
        factor_blocks_reactive(factor_blocks)
        responses_reactive(responses)
        
        # Reset selections when design changes
        selected_factors_for_model_reactive(c(factor_names, factor_blocks))
        
        output$analysis_table_ui_id <- renderUI({
          div(
            p("Choose a default model from the dropdown. Factors to include in the model can be manually adjusted as needed. "),
            fluidRow(
              column(
                4,
                selectInput(inputId = ns("select_model_id"), label = "Select Default Model (optional):", 
                            choices = c("Main Effects", "Main Effects + Interactions", "Response Surface"), 
                            selected = "Main Effects", multiple = FALSE, width = "100%")
              ),
              column(
                2,
                actionButton(inputId = ns("apply_selected_model_id"), label = "Apply", width = "100%", class = "Blue-Button", style = "margin-bottom: 0px; margin-top: 20px;")
              )
            ),
            fluidRow(
              column(
                3,
                selectInput(inputId = ns("available_factors_for_model_id"), label = "Select Factor(s):", 
                            choices = c(factor_names, factor_blocks), size = 8, selected = NULL, multiple = TRUE, selectize = FALSE, width = "100%")
              ),
              column(
                1,
                div(style = "height: 25px;"),
                fluidRow(
                  actionButton(inputId = ns("add_factor_to_model_id"), label = "Add", width = "100%", class = "Blue-Button", style = "margin-bottom: 0px; margin-top: 5px;")
                ),
                fluidRow(
                  actionButton(inputId = ns("remove_factor_to_model_id"), label = "Remove", width = "100%", class = "Blue-Button", style = "margin-bottom: 0px; margin-top: 5px;")
                ),
                fluidRow(
                  actionButton(inputId = ns("cross_factors_for_model_id"), label = "Cross", width = "100%", class = "Blue-Button", style = "margin-bottom: 0px; margin-top: 5px;")
                )
              ),
              column(
                3,
                selectInput(inputId = ns("selected_factors_for_model_id"), label = "Factor(s) for Model:", 
                            choices = selected_factors_for_model_reactive(), size = 8, selected = NULL, multiple = TRUE, selectize = FALSE, width = "100%")
              ),
              column(
                3,
                selectInput(inputId = ns("selected_responses_for_model_id"), label = "Select Response:", 
                            choices = responses, selected = if (length(responses) > 0) responses[1] else NULL, multiple = FALSE, width = "100%")
              )
            ),
            fluidRow(
              column(4,
                     offset = 8,
                     div(
                       align = "right",
                       actionButton(inputId = ns("fit_model_button_id"), label = "Analyze", width = "200px", class = "Blue-Button", style = "margin-bottom: 0px; margin-top: 5px;")
                     )
              )
            )
          )
        })
      } else {
        output$analysis_table_ui_id <- renderUI({
          p("Note: Load a design for analysis on the Enter/Edit Results tab.")
        })
      }
    })
    
    # Apply model preset
    observeEvent(input$apply_selected_model_id, {
      req(input$select_model_id)
      selected_model <- input$select_model_id
      
      factor_names <- factor_names_reactive()
      factor_blocks <- factor_blocks_reactive()
      
      if (selected_model == "Main Effects") {
        selected_factors_for_model_reactive(c(factor_names, factor_blocks))
      } else if (selected_model == "Main Effects + Interactions") {
        main_effects <- c(factor_names, factor_blocks)
        interactions <- if (length(factor_names) > 1) combn(factor_names, 2, FUN = paste0, collapse = " * ") else c()
        selected_factors_for_model_reactive(c(main_effects, interactions))
      } else if (selected_model == "Response Surface") {
        main_effects <- c(factor_names, factor_blocks)
        interactions <- if (length(factor_names) > 1) combn(factor_names, 2, FUN = paste0, collapse = " * ") else c()
        quadratics <- unname(sapply(factor_names, function(x) paste(x, x, sep = " * ")))
        selected_factors_for_model_reactive(c(main_effects, interactions, quadratics))
      }
      
      updateSelectInput(session, "selected_factors_for_model_id", choices = selected_factors_for_model_reactive())
    })
    
    # Add terms manually
    observeEvent(input$add_factor_to_model_id, {
      req(input$available_factors_for_model_id)
      factor_to_add <- input$available_factors_for_model_id
      selected_factors <- selected_factors_for_model_reactive()
      
      updated_selected_factors <- unique(c(selected_factors, factor_to_add))
      selected_factors_for_model_reactive(updated_selected_factors)
      updateSelectInput(session, "selected_factors_for_model_id", choices = selected_factors_for_model_reactive())
    })
    
    # Remove terms manually
    observeEvent(input$remove_factor_to_model_id, {
      req(input$selected_factors_for_model_id)
      factor_to_remove <- input$selected_factors_for_model_id
      selected_factors <- selected_factors_for_model_reactive()
      
      remaining_factors <- setdiff(selected_factors, factor_to_remove)
      selected_factors_for_model_reactive(remaining_factors)
      updateSelectInput(session, "selected_factors_for_model_id", choices = selected_factors_for_model_reactive())
    })
    
    # Cross terms manually
    observeEvent(input$cross_factors_for_model_id, {
      req(input$available_factors_for_model_id, input$selected_factors_for_model_id)
      selected_factors_left <- input$available_factors_for_model_id
      factor_blocks <- factor_blocks_reactive()
      
      if (any(selected_factors_left %in% factor_blocks)) {
        showNotification(strong("Interactions involving blocks are not currently supported."), duration = 15, closeButton = TRUE, type = "error")
      }
      req(!any(selected_factors_left %in% factor_blocks), cancelOutput = TRUE)
      
      selected_factors_right <- input$selected_factors_for_model_id
      factor_combinations <- expand.grid(selected_factors_left, selected_factors_right)
      crossed_factors <- apply(factor_combinations, 1, function(row) paste(row, collapse = " * "))
      
      existing_selected_factors <- selected_factors_for_model_reactive()
      all_selected_factors <- unique(c(existing_selected_factors, crossed_factors))
      selected_factors_for_model_reactive(all_selected_factors)
      updateSelectInput(session, "selected_factors_for_model_id", choices = selected_factors_for_model_reactive())
    })
    
    # Analyze and fit the model
    observeEvent(input$fit_model_button_id, {
      req(selected_factors_for_model_reactive(), input$selected_responses_for_model_id)
      selected_factors_all <- selected_factors_for_model_reactive()
      selected_response <- input$selected_responses_for_model_id
      
      design <- design_reactive()
      factor_types <- factor_types_reactive()
      
      types_vec <- as.character(factor_types[1, ])
      names(types_vec) <- colnames(factor_types)
      
      # Decoupled validation of the response column values
      response_data <- design[[selected_response]]
      val_res <- validate_response_values(response_data, selected_response)
      if (!val_res$valid) {
        showNotification(strong(val_res$message), duration = 15, type = "error")
      }
      # Stop execution if validation fails
      validate(need(val_res$valid, val_res$message))
      
      design_matrix <- design
      for (col in colnames(design_matrix)) {
        role <- types_vec[[col]]
        if (role == "Categorical" || role == "Block") {
          design_matrix[[col]] <- as.factor(design_matrix[[col]])
        } else if (role == "Continuous" || role == "Discrete") {
          design_matrix[[col]] <- as.numeric(design_matrix[[col]])
        }
      }
      
      design_matrix_reactive(design_matrix)
      
      transformed_factors <- sapply(selected_factors_all, transform_term)
      model_formula <- as.formula(paste0(selected_response, " ~ ", paste0(transformed_factors, collapse = " + ")))
      model_formula_reactive(model_formula)
      
      model_aov <- aov(model_formula, design_matrix)
      model_aov_reactive(model_aov)
      
      # ANOVA results formatting
      anova_results <- broom::tidy(model_aov)
      colnames(anova_results) <- c("Factor", "Degrees of Freedom", "Sum of Squares", "Mean Square", "F-Statistic", "P-Value")
      
      design_correlations <- model.matrix(model_formula, data = design_matrix)
      design_correlations <- design_correlations[, colnames(design_correlations) != "(Intercept)", drop = FALSE]
      colnames(design_correlations) <- update_factor_names(colnames(design_correlations), design_matrix)
      
      if (ncol(design_correlations) <= 1) {
        cor_data <- NULL
      } else {
        model_correlations <- cor(design_correlations, method = c("pearson"))
        cor_data <- as.data.frame(as.table(model_correlations))
      }
      
      estimates <- summary.lm(model_aov)
      estimates_df <- broom::tidy(estimates)
      colnames(estimates_df) <- c("Factor", "Estimate", "Standard Error", "t-Statistic", "P-Value")
      estimates_df$Factor <- gsub("\\(Intercept\\)", "Intercept", estimates_df$Factor)
      estimates_df$Factor <- update_factor_names(estimates_df$Factor, design_matrix)
      
      alias_complete <- alias(model_aov)$Complete
      if (is.null(alias_complete)) {
        alias_structure <- "Model terms are not aliased."
      } else {
        alias_structure <- data.frame(alias_complete)
        colnames(alias_structure) <- update_factor_names(colnames(alias_structure), design_matrix)
        rownames(alias_structure) <- update_factor_names(rownames(alias_structure), design_matrix)
        names(alias_structure) <- gsub("X.Intercept.", "Intercept", names(alias_structure))
      }
      
      rsq <- signif(estimates$r.squared, 3)
      adjusted_rsq <- signif(estimates$adj.r.squared, 3)
      
      f_statistic <- "N/A"
      model_p_value <- "N/A"
      if (!is.null(estimates$fstatistic)) {
        f_statistic <- signif(unname(estimates$fstatistic[1]), 3)
        model_p_value <- signif(pf(estimates$fstatistic[1], estimates$fstatistic[2], estimates$fstatistic[3], lower.tail = FALSE), 3)
      }
      
      all_data <- data.frame("actual_values" = design_matrix[[selected_response]])
      all_data$predicted_values <- predict(model_aov)
      
      # Renders ANOVA outputs UI
      output$ANOVA_output_ui_id <- renderUI({
        div(
          fluidRow(
            column(6,
                   h4(strong("Predicted vs. Actual")),
                   box(id = NULL, title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE, headerBorder = FALSE, width = 12,
                       plotOutput(ns("pred_vs_actual_output_id"))
                   )
            ),
            column(6,
                   h4(strong("Model Fit Metrics")),
                   box(id = NULL, title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE, headerBorder = FALSE, width = 12,
                       p(HTML(paste0(strong("R-Squared: "), rsq))),
                       p(HTML(paste0(strong("Adjusted R-Squared: "), adjusted_rsq))),
                       if (!is.null(estimates$fstatistic)) {
                         tagList(
                           p(HTML(paste0(strong("F-statistic: "), f_statistic, ", degrees of freedom: ", unname(estimates$fstatistic[2]), ", ", unname(estimates$fstatistic[3])))),
                           p(HTML(paste0(strong("p-value: "), model_p_value)))
                         )
                       } else {
                         p(strong("F-statistic / P-value not available (insufficient degrees of freedom)"))
                       }
                   )
            )
          ),
          fluidRow(
            column(12,
                   h4(strong("Model Summary")),
                   box(id = NULL, title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE, headerBorder = FALSE, width = 12,
                       tableOutput(ns("ANOVA_results_output_id"))
                   )
            ),
            column(12,
                   h4(strong("Estimates")),
                   box(id = NULL, title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE, headerBorder = FALSE, width = 12,
                       tableOutput(ns("Estimates_output_id"))
                   )
            )
          )
        )
      })
      
      min_value <- as.numeric(min(c(all_data$actual_values, all_data$predicted_values), na.rm = TRUE))
      max_value <- as.numeric(max(c(all_data$actual_values, all_data$predicted_values), na.rm = TRUE))
      
      output$pred_vs_actual_output_id <- renderPlot({
        ggplot(all_data, aes(x = as.numeric(actual_values), y = as.numeric(predicted_values))) +
          geom_point(color = "#0097a9") +
          geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black") +
          labs(x = "Actual", y = "Predicted", title = "Parity Plot") +
          xlim(min_value, max_value) +
          ylim(min_value, max_value) +
          coord_fixed(ratio = 1)
      })
      
      output$ANOVA_results_output_id <- renderTable({
        col_names <- colnames(anova_results)
        data_rounded <- as.data.frame(lapply(anova_results, function(x) {
          if (is.numeric(x)) {
            as.character(signif(x, digits = 4))
          } else {
            x
          }
        }))
        
        data_rounded$significance <- sapply(data_rounded$`P-Value`, function(p) {
          p_num <- suppressWarnings(as.numeric(p))
          if (!is.na(p_num) && p_num < 0.001) "***"
          else if (!is.na(p_num) && p_num < 0.01) "**"
          else if (!is.na(p_num) && p_num < 0.05) "*"
          else ""
        })
        
        data_rounded$`P-Value` <- sapply(data_rounded$`P-Value`, function(p) {
          p_val <- suppressWarnings(as.numeric(p))
          if (!is.na(p_val)) {
            if (p_val < 0.001) "<0.001"
            else as.character(signif(p_val, digits = 4))
          } else {
            p
          }
        })
        
        colnames(data_rounded) <- c(col_names, "")
        data_rounded
      })
      
      output$Estimates_output_id <- renderTable({
        col_names <- colnames(estimates_df)
        data_rounded <- as.data.frame(lapply(estimates_df, function(x) {
          if (is.numeric(x)) {
            as.character(signif(x, digits = 4))
          } else {
            x
          }
        }))
        
        data_rounded$significance <- sapply(data_rounded$`P-Value`, function(p) {
          p_num <- suppressWarnings(as.numeric(p))
          if (!is.na(p_num) && p_num < 0.001) "***"
          else if (!is.na(p_num) && p_num < 0.01) "**"
          else if (!is.na(p_num) && p_num < 0.05) "*"
          else ""
        })
        
        data_rounded$`P-Value` <- sapply(data_rounded$`P-Value`, function(p) {
          p_val <- suppressWarnings(as.numeric(p))
          if (!is.na(p_val)) {
            if (p_val < 0.001) "<0.001"
            else as.character(signif(p_val, digits = 4))
          } else {
            p
          }
        })
        
        colnames(data_rounded) <- c(col_names, "")
        data_rounded
      })
      
      # Renders diagnostics outputs
      output$Design_diagnostics_output_ui_id <- renderUI({
        div(
          fluidRow(
            column(12,
                   h3(strong("Alias Structure:")),
                   uiOutput(ns("alias_table_output_id"))
            )
          ),
          fluidRow(
            column(12,
                   h3(strong("Correlations:")),
                   uiOutput(ns("correlations_output_ui_id"))
            )
          )
        )
      })
      
      output$alias_table_output_id <- renderUI({
        if (is.data.frame(alias_structure)) {
          tableOutput(ns("alias_table"))
        } else {
          verbatimTextOutput(ns("alias_text"))
        }
      })
      
      output$alias_table <- renderTable({ alias_structure }, rownames = TRUE)
      output$alias_text <- renderText({ as.character(alias_structure) })
      
      output$correlations_output_ui_id <- renderUI({
        if (is.null(cor_data)) {
          verbatimTextOutput(ns("correlations_text_output_id"))
        } else {
          plotlyOutput(ns("correlations_plot_output_id"))
        }
      })
      
      output$correlations_text_output_id <- renderText({
        "Correlation matrix not available; model contains only one factor."
      })
      
      output$correlations_plot_output_id <- renderPlotly({
        req(cor_data)
        plot_ly(
          data = cor_data, x = ~Var1, y = ~Var2, z = ~Freq,
          type = "heatmap", colors = colorRamp(c("blue", "white", "red")),
          text = ~ paste("Correlation: ", round(Freq, 2)),
          hoverinfo = "text+x+y", zmin = -1, zmax = 1,
          width = 600, height = 600
        ) %>% layout(
          xaxis = list(title = "", tickangle = -45, autorange = TRUE, side = "top"),
          yaxis = list(title = "", autorange = "reversed", scaleanchor = "x", scaleratio = 1),
          autosize = FALSE
        )
      })
      
      updateTabsetPanel(session, "Analysis_results_tabs_id", selected = "Model Fitting Results")
    })
  })
}
