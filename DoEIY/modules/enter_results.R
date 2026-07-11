# Enter and Edit Results Shiny Module

enter_results_ui <- function(id) {
  ns <- NS(id)
  tagList(
    box(
      id = NULL, title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE, headerBorder = FALSE, width = 12,
      fluidRow(
        column(4, 
               uiOutput(ns("augment_design_btn_ui_id"))     
        ),
        column(4,
               div(
                 align = "center",
                 downloadButton(outputId = ns("save_design_id"), "Save", class = "Blue-Button", style = "margin-bottom: 0px; margin-top: 5px; width: 200px;")
               )
        ),
        column(
          4,
          div(
            align = "center",
            actionButton(inputId = ns("load_design_id"), label = "Load", width = "200px", class = "Blue-Button", style = "margin-bottom: 0px; margin-top: 5px;", icon = icon("upload"))
          )
        )
      ),
      hr(),
      uiOutput(outputId = ns("design_table_UI_id"))
    )
  )
}

enter_results_server <- function(id, design_reactive, factor_types_reactive, factor_data_reactive) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Internal reactives for enter_results module
    augment_factor_data_reactive <- reactiveVal(NULL)
    augment_design_model_components <- reactiveVal(value = NULL)
    augment_original_factor_data_reactive <- reactiveVal(NULL)
    
    Augment_min_runs_reactive <- reactiveVal(NULL)
    Augment_max_runs_reactive <- reactiveVal(NULL)
    
    # Observe design creation or change to render the table and Augment buttons
    observe({
      req(design_reactive())
      
      output$augment_design_btn_ui_id <- renderUI({
        div(
          align = "center",
          actionButton(inputId = ns("augment_design_btn_id"), label = "Augment Design", width = "200px", class = "Blue-Button", style = "margin-bottom: 0px; margin-top: 5px;")
        )
      })
      
      output$design_table_UI_id <- renderUI({
        fluidPage(
          fluidRow(
            column(
              4,
              div(
                align = "center",
                actionButton(inputId = ns("remove_column_id"), label = "Remove Column", width = "200px", class = "Blue-Button", style = "margin-bottom: 10px; margin-top: 0px;")
              )
            ),
            column(
              4,
              div(
                align = "center",
                actionButton(inputId = ns("add_column_id"), label = "Add Column", width = "200px", class = "Blue-Button", style = "margin-bottom: 10px; margin-top: 0px;")
              )
            ),
            column(
              4,
              div(
                align = "center",
                actionButton(inputId = ns("update_column_id"), label = "Update Column", width = "200px", class = "Blue-Button", style = "margin-bottom: 10px; margin-top: 0px;")
              )
            )
          ),
          rHandsontableOutput(outputId = ns("design_table_id"))
        )
      })
    })
    
    output$design_table_id <- renderRHandsontable({
      req(design_reactive())
      rhandsontable(design_reactive(), selectCallback = TRUE, width = "1050px", height = "350px", stretchH = "all") %>%
        hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE)
    })
    
    # Save design trigger
    output$save_design_id <- downloadHandler(
      filename = function() {
        "DoEIY_Design.csv"
      },
      content = function(file) {
        req(input$design_table_id, factor_types_reactive())
        # The first row of the CSV is the factor types, followed by the design matrix itself
        write.csv(rbind(factor_types_reactive(), hot_to_r(input$design_table_id)), file, row.names = FALSE)
      }
    )
    
    # Load design triggers upload dialog modal
    observeEvent(input$load_design_id, {
      showModal(modalDialog(
        title = h3(strong("Select Design File to Load"), style = "margin: 0px; padding: 0px; color: #0097a9; margin-bottom: 5px;"),
        fluidRow(
          column(
            10,
            fileInput(inputId = ns("select_load_file_id"), label = NULL, width = "100%", accept = c(".csv", "text/csv"))
          ),
          column(
            2,
            actionButton(inputId = ns("load_file_id"), label = "Upload", class = "Blue-Button", width = "100%", style = "height: 34px; margin-bottom: 0px;")
          )
        ),
        easyClose = FALSE,
        footer = tagList(
          column(2, offset = 10, actionButton(inputId = ns("close_modal_id"), strong("Cancel"), width = "100%", class = "Grey-Button", style = "height: 34px; margin: 5px;"))
        )
      ))
    })
    
    observeEvent(input$load_file_id, {
      req(input$select_load_file_id)
      file_path <- input$select_load_file_id[[1, "datapath"]]
      design <- read.csv(file_path, header = TRUE, sep = ",")
      
      colnames(design) <- gsub("\\.", " ", colnames(design))
      design[is.na(design)] <- ""
      
      factor_types <- design[1, ]
      design <- design[-1, ]
      rownames(design) <- NULL
      
      # Enforce types (continuous/discrete/responses are numeric, blocks/categorical are character)
      for (col_name in colnames(design)) {
        col_type <- factor_types[[col_name]]
        if (col_type %in% c("Continuous", "Discrete", "Response")) {
          design[[col_name]] <- as.numeric(design[[col_name]])
        } else if (col_type %in% c("Block", "Categorical", "Metadata")) {
          design[[col_name]] <- as.character(design[[col_name]])
        }
      }
      
      # Extract metadata/factor details if loading from a previously saved file
      types_vec <- as.character(factor_types[1, ])
      names(types_vec) <- colnames(factor_types)
      factor_names <- names(types_vec)[!types_vec %in% c("Response", "Block", "Metadata")]
      
      # Reconstruct factor_data mock frame
      if (length(factor_names) > 0) {
        reconstructed_factor_df <- data.frame(
          `Factor Names` = factor_names,
          `Factor Type` = factor(types_vec[factor_names], levels = c("Continuous", "Discrete", "Categorical")),
          `Level Values` = sapply(factor_names, function(fn) paste(sort(unique(design[[fn]])), collapse = ",")),
          stringsAsFactors = FALSE, check.names = FALSE
        )
        factor_data_reactive(reconstructed_factor_df)
      }
      
      factor_types_reactive(factor_types)
      design_reactive(design)
      
      removeModal()
      showNotification(strong("Design file loaded successfully!"), duration = 10, type = "message")
    })
    
    # Close modal trigger
    observeEvent(input$close_modal_id, { removeModal() })
    
    # Remove columns triggers dialog
    observeEvent(input$remove_column_id, {
      req(input$design_table_id)
      design <- hot_to_r(input$design_table_id)
      columns <- colnames(design)
      
      showModal(modalDialog(
        title = h3(strong("Select columns to delete"), style = "margin: 0px; padding: 0px; color: #0097a9; margin-bottom: 5px;"),
        fluidRow(
          column(
            6,
            pickerInput(inputId = ns("remove_selected_columns_id"), label = NULL, choices = columns, multiple = TRUE, width = "100%")
          ),
          column(
            3,
            actionButton(inputId = ns("delete_columns_id"), label = "Delete", class = "Blue-Button", width = "100%", style = "margin-bottom: 0px;")
          )
        ),
        easyClose = FALSE,
        footer = tagList(
          column(2, offset = 10, actionButton(inputId = ns("close_modal_id"), strong("Cancel"), width = "100%", class = "Grey-Button", style = "height: 34px; margin: 5px;"))
        )
      ))
    })
    
    observeEvent(input$delete_columns_id, {
      req(input$remove_selected_columns_id)
      selected_columns <- input$remove_selected_columns_id
      design <- hot_to_r(input$design_table_id)
      factor_types <- factor_types_reactive()
      
      design <- design[, !(names(design) %in% selected_columns), drop = FALSE]
      factor_types <- factor_types[, !(names(factor_types) %in% selected_columns), drop = FALSE]
      
      factor_types_reactive(factor_types)
      design_reactive(design)
      removeModal()
    })
    
    # Add columns triggers dialog
    observeEvent(input$add_column_id, {
      showModal(modalDialog(
        title = h3(strong("Create a new column"), style = "margin: 0px; padding: 0px; color: #0097a9; margin-bottom: 5px;"),
        fluidRow(
          column(4, textInput(inputId = ns("add_column_name_id"), label = "Column Name:")),
          column(4, selectInput(inputId = ns("add_column_type_id"), label = "Column Type:", choices = c("Response", "Block", "Metadata"))),
          column(4, selectInput(inputId = ns("add_column_position_id"), label = "Insert Position:", choices = c("Append to end", colnames(design_reactive()))))
        ),
        easyClose = FALSE,
        footer = tagList(
          actionButton(ns("create_column_id"), "Create", class = "Blue-Button", style = "height: 34px; margin: 5px;"),
          actionButton(ns("close_modal_id"), strong("Cancel"), class = "Grey-Button", style = "height: 34px; margin: 5px;")
        )
      ))
    })
    
    observeEvent(input$create_column_id, {
      req(input$add_column_name_id, input$add_column_type_id, input$add_column_position_id)
      column_name <- input$add_column_name_id
      column_type <- input$add_column_type_id
      insert_position_name <- input$add_column_position_id
      
      design <- hot_to_r(input$design_table_id)
      factor_types <- factor_types_reactive()
      
      if (column_name %in% colnames(design)) {
        showNotification("Column name must be unique. Duplicate columns are not allowed.", type = "error")
        return(NULL)
      }
      
      new_col <- rep("", nrow(design))
      design$Temp_column <- new_col
      factor_types$Temp_column <- column_type
      
      col_names_new <- colnames(design)
      
      if (insert_position_name == "Append to end") {
        design_updated <- design
        factor_types_updated <- factor_types
        colnames_updated <- c(col_names_new[-length(col_names_new)], column_name)
      } else {
        insert_position <- which(col_names_new == insert_position_name)
        design_updated <- data.frame(
          design[, 1:(insert_position - 1), drop = FALSE],
          Temp_column = new_col,
          design[, insert_position:(ncol(design) - 1), drop = FALSE]
        )
        factor_types_updated <- data.frame(
          factor_types[, 1:(insert_position - 1), drop = FALSE],
          New_Column = column_type,
          factor_types[, insert_position:(ncol(factor_types) - 1), drop = FALSE]
        )
        colnames_updated <- c(col_names_new[1:(insert_position - 1)], column_name, col_names_new[insert_position:(length(col_names_new) - 1)])
      }
      
      colnames(design_updated) <- colnames_updated
      colnames(factor_types_updated) <- colnames_updated
      
      design_updated$Temp_column <- NULL
      factor_types_updated$Temp_column <- NULL
      
      design_reactive(design_updated)
      factor_types_reactive(factor_types_updated)
      removeModal()
    })
    
    # Update column trigger dialog
    observeEvent(input$update_column_id, {
      req(input$design_table_id)
      columns <- colnames(hot_to_r(input$design_table_id))
      
      showModal(modalDialog(
        title = h3(strong("Update Column Name and Type"), style = "margin: 0px; padding: 0px; color: #0097a9; margin-bottom: 5px;"),
        fluidRow(
          column(4, selectInput(inputId = ns("column_to_update_id"), label = "Column to Update:", choices = columns)),
          column(4, textInput(inputId = ns("update_column_name_id"), label = "New Column Name:")),
          column(4, selectInput(inputId = ns("update_column_type_id"), label = "New Column Type:", choices = c("Response", "Block", "Metadata")))
        ),
        easyClose = FALSE,
        footer = tagList(
          actionButton(ns("save_column_updates_id"), "Save", class = "Blue-Button", style = "height: 34px; margin: 5px;"),
          actionButton(ns("close_modal_id"), strong("Cancel"), class = "Grey-Button", style = "height: 34px; margin: 5px;")
        )
      ))
    })
    
    observeEvent(input$column_to_update_id, {
      req(input$column_to_update_id)
      updateTextInput(session, "update_column_name_id", value = input$column_to_update_id)
      current_type <- as.character(factor_types_reactive()[[input$column_to_update_id]])
      if (current_type %in% c("Response", "Block", "Metadata")) {
        updateSelectInput(session, "update_column_type_id", selected = current_type)
      }
    })
    
    observeEvent(input$save_column_updates_id, {
      req(input$column_to_update_id, input$update_column_name_id, input$update_column_type_id)
      column_name <- input$update_column_name_id
      column_type <- input$update_column_type_id
      
      design <- hot_to_r(input$design_table_id)
      factor_types <- factor_types_reactive()
      
      if (column_name != input$column_to_update_id && column_name %in% colnames(design)) {
        showNotification("Column name must be unique. Duplicate columns are not allowed.", type = "error")
        return(NULL)
      }
      
      factor_types[[input$column_to_update_id]] <- column_type
      names(design)[names(design) == input$column_to_update_id] <- column_name
      names(factor_types)[names(factor_types) == input$column_to_update_id] <- column_name
      
      design_reactive(design)
      factor_types_reactive(factor_types)
      removeModal()
    })
    
    # Watches handsontable updates and writes them back to design_reactive
    observe({
      req(input$design_table_id)
      design_reactive(hot_to_r(input$design_table_id))
    })
    
    # Opens design augmentation modal
    observeEvent(input$augment_design_btn_id, {
      req(factor_types_reactive(), design_reactive())
      current_factor_types <- factor_types_reactive()
      current_design <- design_reactive()
      
      types_vec <- as.character(current_factor_types[1, ])
      names(types_vec) <- colnames(current_factor_types)
      
      keep_types <- c("Continuous", "Discrete", "Categorical")
      kept_names <- names(types_vec)[types_vec %in% keep_types]
      
      design_subset <- current_design[, kept_names, drop = FALSE]
      factor_names <- kept_names
      factor_types <- types_vec[kept_names]
      
      level_values <- sapply(factor_names, function(fn) {
        col <- design_subset[[fn]]
        if (factor_types[[fn]] %in% c("Continuous", "Discrete")) {
          vals <- sort(unique(suppressWarnings(as.numeric(col))))
        } else {
          vals <- sort(unique(as.character(col)))
        }
        vals <- vals[vals != "" & !is.na(vals)]
        paste(vals, collapse = ",")
      })
      
      factor_df <- data.frame(
        Factor = factor_names,
        Type = factor(factor_types, levels = c("Continuous", "Discrete", "Categorical")),
        Level_Values = level_values,
        stringsAsFactors = FALSE
      )
      
      colnames(factor_df) <- c("Factor Names", "Factor Type", "Level Values")
      
      augment_original_factor_data_reactive(factor_df)
      augment_factor_data_reactive(factor_df)
      
      # Initialize model components list
      augment_design_model_components(factor_names)
      
      showModal(modalDialog(
        title = h3(strong("Review Factor Details for Augmentation"), style = "margin: 0px; padding: 0px; color: #0097a9; margin-bottom: 5px;"),
        fluidRow(
          column(12, p("Verify the factor details for the existing design. All factors in the table below must be configured correctly. Factor names cannot contain spaces or special characters.")),
          column(12, div(style = "height: 220px;", rHandsontableOutput(ns("augment_factor_table_id"), height = "100%"))),
          uiOutput(ns("augment_warning_message_id"))
        ),
        easyClose = FALSE,
        footer = tagList(
          actionButton(ns("augment_open_model_modal_id"), "Next", class = "Blue-Button", style = "height: 34px; margin: 5px;"),
          actionButton(ns("close_modal_id"), strong("Cancel"), class = "Grey-Button", style = "height: 34px; margin: 5px;")
        )
      ))
      
      output$augment_factor_table_id <- renderRHandsontable({
        rhandsontable(factor_df, selectCallback = TRUE, width = 1150, height = 250) %>%
          hot_cols(colWidths = c(200, 200, 700)) %>%
          hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE)
      })
    })
    
    # Opens model building dialog for Augmentation
    observeEvent(input$augment_open_model_modal_id, {
      req(input$augment_factor_table_id)
      current_factor_data <- hot_to_r(input$augment_factor_table_id)
      augment_factor_data_reactive(current_factor_data)
      
      warning_msg <- c()
      for (i in seq_len(dim(current_factor_data)[1])) {
        warning_msg_temp <- c()
        factor_name <- current_factor_data$`Factor Names`[i]
        
        if (is.na(factor_name) || grepl("^\\s*$", factor_name) || factor_name == "") {
          warning_msg_temp <- c(warning_msg_temp, paste0("factor name is missing"))
        }
        if (is.na(factor_name) || grepl("^\\s*$", factor_name) || factor_name == "" || factor_name != make.names(factor_name)) {
          warning_msg_temp <- c(warning_msg_temp, "factor name cannot contain spaces or special characters")
        }
        
        factor_type <- current_factor_data$`Factor Type`[i]
        if (is.na(as.character(factor_type))) {
          warning_msg_temp <- c(warning_msg_temp, paste0("factor type is missing"))
        }
        
        level_values <- current_factor_data$`Level Values`[i]
        if (is.na(level_values) || grepl("^\\s*$", level_values) || level_values == "") {
          warning_msg_temp <- c(warning_msg_temp, paste0("level values are missing"))
        } else {
          level_values_vector <- gsub("^\\s+|\\s+$", "", unlist(strsplit(level_values, ",")))
          if (factor_type == "Continuous") {
            if (!all(!is.na(as.numeric(level_values_vector)))) {
              warning_msg_temp <- c(warning_msg_temp, paste0("the level values for a Continuous factor must be numeric"))
            }
          }
        }
        
        if (length(warning_msg_temp) != 0) {
          warning_msg_temp <- paste0("Row ", i, ": ", paste0(warning_msg_temp, collapse = "; "), ".")
          warning_msg <- c(warning_msg, warning_msg_temp)
        }
      }
      
      if (length(warning_msg) != 0) {
        output$augment_warning_message_id <- renderUI({
          lapply(seq_len(length(warning_msg)), function(j) {
            p(strong(warning_msg[j]), class = "text-red")
          })
        })
      }
      req(length(warning_msg) == 0, cancelOutput = TRUE)
      
      removeModal()
      showAugmentModelModal()
    })
    
    showAugmentModelModal <- function() {
      showModal(modalDialog(
        title = h3(strong("Build Augmentation Model"), style = "margin: 0px; padding: 0px; color: #0097a9; margin-bottom: 5px;"),
        fluidRow(
          column(12, p("Provide the components of the model that you wish to fit to the augmented design. Select additional interaction or quadratic terms from the dropdown menus below.")),
          column(6,
                 h4(strong("Allowed Model Terms:")),
                 selectInput(inputId = ns("selected_factors_augment_DOD_id"), label = NULL, choices = augment_design_model_components(), selected = NULL, multiple = TRUE, selectize = FALSE, width = "100%", size = 10)
          ),
          column(6,
                 br(),
                 fluidRow(column(10, offset = 1, actionButton(inputId = ns("augment_add_interactions_id"), "Add Interactions", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;"))),
                 br(),
                 fluidRow(column(10, offset = 1, actionButton(inputId = ns("augment_add_quadratics_id"), "Add Quadratics", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;"))),
                 br(),
                 fluidRow(column(10, offset = 1, actionButton(inputId = ns("augment_remove_interactions_id"), "Remove Selected Terms", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;")))
          )
        ),
        easyClose = FALSE,
        footer = tagList(
          fluidRow(
            column(2, offset = 6, actionButton(inputId = ns("augment_reopen_factors_modal_id"), "Back", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;")),
            column(2, actionButton(inputId = ns("augment_open_dod_modal_id"), "Next", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;")),
            column(2, actionButton(inputId = ns("close_modal_id"), strong("Cancel"), width = "100%", class = "Grey-Button", style = "height: 34px; margin: 5px;"))
          )
        )
      ))
    }
    
    # Augment modal back button
    observeEvent(input$augment_reopen_factors_modal_id, {
      removeModal()
      showModal(modalDialog(
        title = h3(strong("Review Factor Details for Augmentation"), style = "margin: 0px; padding: 0px; color: #0097a9; margin-bottom: 5px;"),
        fluidRow(
          column(12, p("Verify the factor details for the existing design. All factors in the table below must be configured correctly. Factor names cannot contain spaces or special characters.")),
          column(12, div(style = "height: 220px;", rHandsontableOutput(ns("augment_factor_table_id"), height = "100%"))),
          uiOutput(ns("augment_warning_message_id"))
        ),
        easyClose = FALSE,
        footer = tagList(
          actionButton(ns("augment_open_model_modal_id"), "Next", class = "Blue-Button", style = "height: 34px; margin: 5px;"),
          actionButton(ns("close_modal_id"), strong("Cancel"), class = "Grey-Button", style = "height: 34px; margin: 5px;")
        )
      ))
      output$augment_factor_table_id <- renderRHandsontable({
        rhandsontable(augment_factor_data_reactive(), selectCallback = TRUE, width = 1150, height = 250) %>%
          hot_cols(colWidths = c(200, 200, 700)) %>%
          hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE)
      })
    })
    
    # Model components selection observers for Augment modal
    observeEvent(input$augment_add_interactions_id, {
      factor_data <- augment_factor_data_reactive()
      factor_names <- as.character(factor_data$`Factor Names`)
      req(length(factor_names) > 1, cancelOutput = TRUE)
      
      showModal(modalDialog(
        title = h3(strong("Add Interaction Term")),
        selectInput(inputId = ns("augment_interaction_factor1_id"), label = "Factor 1", choices = factor_names),
        selectInput(inputId = ns("augment_interaction_factor2_id"), label = "Factor 2", choices = factor_names),
        easyClose = FALSE,
        footer = tagList(
          actionButton(inputId = ns("augment_interactions_add_btn_id"), "Add", class = "Blue-Button", style = "height:34px; margin:5px;"),
          actionButton(inputId = ns("augment_model_select_back_btn_id"), "Back", class = "Grey-Button", style = "height:34px; margin:5px;")
        )
      ))
    })
    
    observeEvent(input$augment_interactions_add_btn_id, {
      req(input$augment_interaction_factor1_id, input$augment_interaction_factor2_id)
      f1 <- input$augment_interaction_factor1_id
      f2 <- input$augment_interaction_factor2_id
      
      if (f1 == f2) {
        showNotification("Interaction terms must be between different factors.", type = "warning")
      } else {
        new_interaction <- paste(sort(c(f1, f2)), collapse = ":")
        currently_selected_factors <- augment_design_model_components()
        updated_selected_factors <- unique(c(currently_selected_factors, new_interaction))
        term_types <- sapply(updated_selected_factors, classify_terms)
        terms_df <- data.frame(term = updated_selected_factors, type = term_types, stringsAsFactors = FALSE)
        sorted_terms_df <- terms_df[order(terms_df$type, terms_df$term), ]
        sorted_terms <- sorted_terms_df$term
        augment_design_model_components(sorted_terms)
      }
      
      removeModal()
      showAugmentModelModal()
    })
    
    observeEvent(input$augment_add_quadratics_id, {
      factor_data <- augment_factor_data_reactive()
      continuous_discrete_factors <- as.character(factor_data$`Factor Names`[factor_data$`Factor Type` %in% c("Continuous", "Discrete")])
      
      if (length(continuous_discrete_factors) == 0) {
        showNotification("Quadratic terms are only supported for Continuous or Discrete factors.", type = "warning")
      } else {
        showModal(modalDialog(
          title = h3(strong("Add Quadratic Term")),
          selectInput(inputId = ns("augment_quadratic_factor_id"), label = "Factor", choices = continuous_discrete_factors),
          easyClose = FALSE,
          footer = tagList(
            actionButton(inputId = ns("augment_quadratics_add_btn_id"), "Add", class = "Blue-Button", style = "height:34px; margin:5px;"),
            actionButton(inputId = ns("augment_model_select_back_btn_id"), "Back", class = "Grey-Button", style = "height:34px; margin:5px;")
          )
        ))
      }
    })
    
    observeEvent(input$augment_quadratics_add_btn_id, {
      req(input$augment_quadratic_factor_id)
      f <- input$augment_quadratic_factor_id
      new_quadratic <- paste0(f, ":", f)
      
      currently_selected_factors <- augment_design_model_components()
      updated_selected_factors <- unique(c(currently_selected_factors, new_quadratic))
      term_types <- sapply(updated_selected_factors, classify_terms)
      terms_df <- data.frame(term = updated_selected_factors, type = term_types, stringsAsFactors = FALSE)
      sorted_terms_df <- terms_df[order(terms_df$type, terms_df$term), ]
      sorted_terms <- sorted_terms_df$term
      augment_design_model_components(sorted_terms)
      
      removeModal()
      showAugmentModelModal()
    })
    
    observeEvent(input$augment_remove_interactions_id, {
      req(input$selected_factors_augment_DOD_id)
      terms_to_remove <- input$selected_factors_augment_DOD_id
      currently_selected_factors <- augment_design_model_components()
      updated_selected_factors <- setdiff(currently_selected_factors, terms_to_remove)
      
      augment_design_model_components(updated_selected_factors)
      updateSelectInput(session, inputId = "selected_factors_augment_DOD_id", choices = updated_selected_factors)
    })
    
    observeEvent(input$augment_model_select_back_btn_id, {
      removeModal()
      showAugmentModelModal()
    })
    
    # Opens Augmentation Runs spec Modal
    observeEvent(input$augment_open_dod_modal_id, {
      factor_data <- augment_factor_data_reactive()
      components <- augment_design_model_components()
      
      current_design <- design_reactive()
      factor_names <- as.character(factor_data$`Factor Names`)
      existing_design <- current_design[, colnames(current_design) %in% factor_names, drop = FALSE]
      
      factor_types <- setNames(as.character(factor_data$`Factor Type`), factor_names)
      level_values <- setNames(
        lapply(strsplit(as.character(factor_data$`Level Values`), ","), function(x) trimws(x)),
        factor_names
      )
      
      # Rank calculation logic
      var_names <- trimws(unique(all.vars(as.formula(paste("~", paste(sapply(components, fix_component), collapse = " + "))))))
      existing_design_subset <- existing_design[, intersect(var_names, colnames(existing_design)), drop = FALSE]
      
      # Generate candidate set & models parameters count p
      candidate_list <- lapply(var_names, function(f) {
        ftype <- factor_types[[f]]
        levs  <- level_values[[f]]
        deg   <- factor_degree(f, components)
        if (ftype == "Continuous") {
          seq(-1, 1, length.out = deg + 1)
        } else if (ftype == "Discrete") {
          sort(as.numeric(levs))
        } else if (ftype == "Categorical") {
          factor(levs, levels = levs)
        }
      })
      candidate_set <- do.call(expand.grid, candidate_list)
      colnames(candidate_set) <- var_names
      
      formula <- as.formula(paste("~", paste(sapply(components, fix_component), collapse = " + ")))
      X_full <- model.matrix(formula, data = candidate_set)
      p <- ncol(X_full)
      
      # Pad existing design for model.matrix
      typed_existing <- data.frame(row_id = seq_len(nrow(existing_design_subset)), stringsAsFactors = FALSE)
      for (f in var_names) {
        ftype <- factor_types[[f]]
        levs  <- level_values[[f]]
        if (f %in% colnames(existing_design_subset)) {
          col <- existing_design_subset[[f]]
          if (ftype %in% c("Continuous", "Discrete")) {
            typed_existing[[f]] <- as.numeric(col)
          } else if (ftype == "Categorical") {
            typed_existing[[f]] <- factor(as.character(col), levels = levs)
          }
        } else {
          if (ftype %in% c("Continuous", "Discrete")) {
            typed_existing[[f]] <- rep(0, nrow(existing_design_subset))
          } else if (ftype == "Categorical") {
            typed_existing[[f]] <- factor(rep(levs[1], nrow(existing_design_subset)), levels = levs)
          }
        }
      }
      
      # Scale continuous/discrete values to [-1,1]
      for (f in intersect(names(typed_existing), names(factor_types))) {
        ftype <- factor_types[[f]]
        levs  <- level_values[[f]]
        col   <- typed_existing[[f]]
        if (ftype %in% c("Continuous", "Discrete")) {
          x <- suppressWarnings(as.numeric(col))
          xmin <- min(x, na.rm = TRUE)
          xmax <- max(x, na.rm = TRUE)
          if (xmin == xmax) {
            typed_existing[[f]] <- rep(0, length(x))
          } else {
            typed_existing[[f]] <- 2 * (x - xmin) / (xmax - xmin) - 1
          }
        } else if (ftype == "Categorical") {
          typed_existing[[f]] <- factor(match(as.character(col), levs), levels = seq_len(length(levs)))
        }
      }
      typed_existing$row_id <- NULL
      
      X_exist <- model.matrix(formula, data = typed_existing)
      rank_existing <- qr(X_exist)$rank
      
      required_new_runs <- max(0, p - rank_existing)
      max_runs <- p * 3 # Reasonable ceiling for custom entries
      
      Augment_min_runs_reactive(required_new_runs)
      Augment_max_runs_reactive(max_runs)
      
      removeModal()
      
      showModal(modalDialog(
        title = h3(strong("Augmented Design Details"), style = "margin: 0px; padding: 0px; color: #0097a9; margin-bottom: 5px;"),
        fluidRow(
          column(6,
                 p(HTML(paste0(strong("Existing Runs: "), nrow(existing_design)))),
                 DT::dataTableOutput(ns("summary_factors_table_id"))
          ),
          column(6,
                 fluidRow(
                   uiOutput(ns("Augment_min_runs_ui_id")),
                   numericInput(inputId = ns("Augment_num_runs_id"), label = "Number of Additional Runs to Add:", value = required_new_runs, min = required_new_runs, max = max_runs, width = "100%"),
                   materialSwitch(inputId = ns("Augment_randomize_checkbox_id"), label = "Randomize Added Runs:", value = TRUE, width = "100%", status = "primary")
                 )
          )
        ),
        easyClose = FALSE,
        footer = tagList(
          fluidRow(
            column(2, offset = 6, actionButton(inputId = ns("augment_open_dod_modal_back_btn_id"), "Back", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;")),
            column(2, actionButton(inputId = ns("augment_design_id"), "Augment", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;")),
            column(2, actionButton(inputId = ns("close_modal_id"), strong("Cancel"), width = "100%", class = "Grey-Button", style = "height: 34px; margin: 5px;"))
          )
        )
      ))
      
      output$Augment_min_runs_ui_id <- renderUI(p(strong(paste0("Minimum Additional Runs Required: ", required_new_runs))))
      
      output$summary_factors_table_id <- DT::renderDataTable(
        factor_data,
        escape = FALSE, selection = "none", server = FALSE,
        extensions = c("FixedColumns", "FixedHeader", "Scroller"),
        options = list(
          dom = "t",
          pageLength = 10, paging = TRUE,
          headerCallback = DT::JS("function(thead) {", "  $(thead).css('font-size', '1em');", "}")
        ),
        rownames = FALSE, filter = "none"
      )
    })
    
    observeEvent(input$augment_open_dod_modal_back_btn_id, {
      removeModal()
      showAugmentModelModal()
    })
    
    # Executes Augmentation operation
    observeEvent(input$augment_design_id, {
      factor_data <- augment_factor_data_reactive()
      original_factor_types <- factor_types_reactive()
      n_factors <- dim(factor_data)[1]
      factor_names <- factor_data$`Factor Names`
      
      randomize <- input$Augment_randomize_checkbox_id
      model_components <- augment_design_model_components()
      
      existing_design_full <- design_reactive()
      existing_design <- existing_design_full[, colnames(existing_design_full) %in% factor_names, drop = FALSE]
      other_columns <- existing_design_full[, !(colnames(existing_design_full) %in% factor_names), drop = FALSE]
      
      val_res <- validate_num_runs(input$Augment_num_runs_id, Augment_min_runs_reactive(), Augment_max_runs_reactive())
      if (!val_res$valid) {
        showNotification(val_res$message, duration = 20, closeButton = TRUE)
      }
      req(val_res$valid, cancelOutput = TRUE)
      num_runs <- as.integer(round(input$Augment_num_runs_id))
      
      factor_types <- setNames(as.character(factor_data$`Factor Type`), factor_names)
      level_values <- setNames(
        lapply(strsplit(as.character(factor_data$`Level Values`), ","), function(x) trimws(x)),
        factor_names
      )
      
      design <- D_Optimal_Augment(existing_design, level_values, model_components, factor_types, num_runs, randomize)
      req(!is.null(design))
      
      # Re-pad existing non-factor metadata columns
      n_old <- nrow(other_columns)
      n_new <- nrow(design) - n_old
      
      filler <- as.data.frame(
        matrix("", nrow = n_new, ncol = ncol(other_columns)),
        stringsAsFactors = FALSE
      )
      colnames(filler) <- colnames(other_columns)
      
      other_extended <- rbind(other_columns, filler)
      design <- cbind(design, other_extended)
      
      other_factor_types <- original_factor_types[colnames(other_columns)]
      factor_types <- c(factor_types, other_factor_types)
      factor_types_df <- as.data.frame(factor_types)
      colnames(factor_types_df) <- colnames(design)
      
      # Convert columns to correct types
      for (col_name in colnames(design)) {
        col_type <- factor_types_df[[col_name]]
        if (col_type %in% c("Continuous", "Discrete", "Response")) {
          design[[col_name]] <- as.numeric(design[[col_name]])
        } else if (col_type %in% c("Block", "Categorical", "Metadata")) {
          design[[col_name]] <- as.character(design[[col_name]])
        }
      }
      
      design_reactive(design)
      factor_types_reactive(factor_types_df)
      removeModal()
      
      showNotification(strong("Design augmented successfully!"), duration = 10, type = "message")
    })
  })
}
