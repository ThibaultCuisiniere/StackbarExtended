#' @title Clustered stackbar coupled with DESeq2
#'
#' @description Create a plot of taxa relative abundance, clustered by phylogeny with the possibility to highlight the differentially abundant features.
#' 
#' @param ps_object A [phyloseq-class] object. Must have a [sample_data], [otu_table], and [tax_table] components.
#' @param exp_group A column in the [sample_data] containing the experimental group information.
#' @param subset_group Default NULL. Groups among the 'exp_group' column from the [sample_data] object to subset.
#' @param sample_name Name of the column in [sample_data] containing the unique sample identifier.
#' @param main_level  Default 'Phylum'. Level present in the [tax_table] to which taxa will be clustered.
#' @param sub_level Default 'Family'. Level present in the [tax_table] to which taxa will be plotted and analyzed.
#' @param threshold Default 1, % threshold to regroup taxa with lower relative abundance into the 'other' groups
#' @param n_phy Default 4, number of main_level to plot. Same number of colors must be given in the 'hues' parameter.
#' @param mean_group Default FALSE. Whether or not agglomerate samples belonging to the same exp_group.
#' @param hues Color used to represent main_level, should be the same number than n_phy parameter. See [colorRampPalette].
#' @param color_bias Define the gradiant of shades among the colors. See [colorRampPalette].
#' @param n_row Default 1. Define the number of row of the graph. See [facet_wrap]
#' @param n_col Default NULL. Define the number of column of the graph. See [facet_wrap]
#' @param text_size Control the text size of the graph. See [ggplot2].
#' @param legend_size Control the legend text size of the graph. see ggplot2. See [ggplot2].
#' @param x_axis_size Control the x axis text size of the graph. see ggplot2. See [ggplot2].
#' @param differential_analysis Default TRUE. Whether or not use DESeq2 to do a differential abundance analysis. See [DESeq2-package].
#' @param mult_comp Default False. Whether or not perform  differential abundance analysis pairwise comparisons. 
#' @param selected_comparisons Default NULL. A list of vectors of length 2 to subset differential abundance analysis pairwise comparisons. Example: "list(c("0", "7"), c("0", "21"), c("0", "90"))"
#' @param test Default 'Wald'. See [DESeq2-package].
#' @param fdr_threshold Default '0.05'. Threshold to which taxa bellow are considered significant.
#' @param sig_lab Default TRUE. Whether to add stars after taxa name reflecting statistical significance. Only available for 2 groups comparison.
#' @param fitType Default "parametric". See [DESeq2-package].
#' @param sfType Default "ratio". See [DESeq2-package].
#' @param betaPrior Default FALSE.See [DESeq2-package].
#' @param reduced Default FALSE. See [DESeq2-package]. For example: "~1".
#' @param quiet Default TRUE. See [DESeq2-package].
#' @param minReplicatesForReplace Default '7'. See [DESeq2-package].
#' @param modelMatrixType Default "standard". See [DESeq2-package].
#' @param useT Default FALSE. See [DESeq2-package].
#' @param minmu See [DESeq2-package].
#' @param parallel Default FALSE. See [DESeq2-package].
#' 
#'
#'
#' @return A [list] containing  [ggplot2] '$plot' object and, if differential_analysis = TRUE, [data.frame] from \code{\link{results}} DESeq2 function output of significant features '$significant_table_main' and '$significant_table_sub' .
#'
#' @examples
#' 
#' \dontrun{
#' 
#'data(ps)
#'
#'my_plot <- plot_microbiota(
#'  ps_object = ps,
#'  exp_group = 'timepoint',
#'  sample_name = 'SampleID',
#'  hues = c("Purples", "Blues", "Greens", "Oranges"),
#'  differential_analysis = T,
#'  sig_lab = T,
#'  fdr_threshold = 0.05
#' )
#'
#'print(my_plot$plot)
#'print(my_plot$significant_table_main)
#'print(my_plot$significant_table_sub)
#' 
#'
#'  }
#'
#' @export
#' 
#' @seealso [phyloseq]
#' @seealso [DESeq2-package]
#' @seealso [ggplot2]
#' 

#' @author Thibault Cuisiniere,  Manuela M. Santos
#'
#'
#' @import phyloseq
#' @import ggplot2
#' @import ggtext
#' @import dplyr
#' @import reshape2
#' @import RColorBrewer
#' @import DESeq2



plot_microbiota <- function(ps_object = ps,
                                exp_group = 'group',
                                subset_group = NULL,
                                sample_name = 'SampleID',
                                main_level = 'Phylum',
                                sub_level = 'Family',
                                threshold = 1,
                                n_phy = 4,
                                mean_group = FALSE,
                                hues = c("Oranges", "Greens", "Blues", "Purples"),
                                color_bias = 2,
                                n_row = 1,
                                n_col = NULL,
                                text_size = 9,
                                legend_size = 7,
                                x_axis_size = 8,
                                differential_analysis = FALSE,
                                mult_comp = FALSE,
                                selected_comparisons = NULL,
                                test = c("Wald", "LRT")[1],
                                fdr_threshold = 0.05,
                                sig_lab = FALSE,
                                fitType = c("parametric", "local", "mean", "glmGamPoi")[1],
                                sfType = c("ratio", "poscounts", "iterate")[1],
                                betaPrior = FALSE,
                                reduced = FALSE,
                                quiet = TRUE,
                                minReplicatesForReplace = 7,
                                modelMatrixType = c("standard", "expanded")[1],
                                useT = FALSE,
                                minmu = if (fitType == "glmGamPoi") 1e-06 else 0.5,
                                parallel = FALSE
                                ) {
  
  
  
  # Validate inputs
  if (!("phyloseq" %in% class(ps_object)))
    stop("ps_object must be a phyloseq-class object.")
  if (!main_level %in% colnames(tax_table(ps_object)))
    stop("main_level column not found in the tax_table.")
  if (!sub_level %in% colnames(tax_table(ps_object)))
    stop("sub_level column not found in the tax_table.")
  if (!exp_group %in% names(sample_data(ps_object)))
    stop("exp_group column not found in the sample_data.")
  
  
  #Assure taxa are rows
  if (!taxa_are_rows(ps_object)) {
    ps_object <- t(ps_object)
  }
  
  #Subset the samples belonging to the groups selected
  if (is.null(subset_group) == F) {
    keep_samples = as.character(get_variable(ps_object, exp_group)) %in% subset_group
    ps_object = prune_samples(keep_samples, ps_object)
    
  }
  
  #transform to relative abundance
  ps_prop <-
    transform_sample_counts(ps_object, function(OTU)
      ((OTU / sum(OTU)) * 100))
  
  
  #store TAX, OTU and metadata tables
  otu <-  as.data.frame(otu_table(ps_prop))
  tax <-  as.data.frame(tax_table(ps_prop))
  meta <- data.frame(sample_data(ps_prop))
  
  #clean up OTU absent and taxa not present
  
  otu <- otu %>% filter_all(any_vars(. != 0))
  
  tax <- subset(tax, rownames(tax) %in% rownames(otu))
  
  #Agglomerate at the sub level defined
  ps_f <- tax_glom(ps_prop, taxrank = sub_level , NArm = FALSE)
  #Store OTU and TAX
  otu_f <-  as.data.frame(otu_table(ps_f))
  tax_f <-  as.data.frame(tax_table(ps_f))
  
  #Bind TAX and OTU
  
  otu_tax_f <- cbind(tax_f, otu_f)
  
  
  #create a vector, if True, the ASV mean among the samples is higher than the threshold
  
  #keep the position of the column storing the abundance of the taxa in the otu_tax_f object
  position <- ncol(tax) + 1:(ncol(otu_tax_f) - ncol(tax))
  
  message(paste0('\n', length(position), ' samples are analyzed \n'))
  
  #Create 2 vectors in the dataframe :
  #1 storing the abundance mean of the taxa
  #2 T or F this taxa is above the define threshold of filtering
  
  otu_tax_f <- otu_tax_f %>%
    rowwise() %>%
    mutate(Mean = mean(c_across(all_of(position)))) %>%
    mutate(high_abundance = case_when(Mean > threshold  ~ TRUE,
                                      TRUE  ~ FALSE))
  
  #Store top x of the main_level phylogeny
  
  topx <- otu_tax_f %>%
    dplyr::group_by(!!as.name(main_level)) %>%
    dplyr::summarise(sum_top = sum(c_across(all_of(position)))) %>%
    dplyr::arrange(desc(sum_top)) %>%
    dplyr::slice_head(n = n_phy)
  
  
  #Create a vector, TRUE when the phylum belongs to the top X, the ones we want to plot
  otu_tax_f <- mutate(otu_tax_f,
                      selected_top = !!as.name(main_level) %in% pull(topx[, main_level]))
  
  #replace NA by "Unknown" 
  otu_tax_f[is.na(otu_tax_f)] <- "Unknown"
  
  #Add a column 'plot_taxa' to create the legend of the graphs :
  
  otu_tax_f <- otu_tax_f %>%
    mutate(
      plot_taxa = case_when(
        high_abundance == TRUE  &
          selected_top == TRUE &
          !!as.name(sub_level) == "Unknown"  ~
          paste0("Unknown ",!!as.name(main_level)),
        high_abundance == TRUE  &
          selected_top == TRUE  &
          !!as.name(sub_level) != "Unknown"  ~ !!as.name(sub_level),
        high_abundance == FALSE &
          selected_top == TRUE  ~ paste0("Others ", !!as.name(main_level)),
        selected_top   == FALSE ~ paste0("Others ")
      )
    )
  
  
  if (nrow(topx) != length(hues)) {
    message(
      'the number of colors choosen (',
      length(hues),
      ') is different from the defined number of features to plot (',
      nrow(topx),
      ')'
    )
  }
  
  
  #initialize 
  df <-
    as_tibble(matrix(nrow = 0, ncol = length(colnames(otu_tax_f))),
              .name_repair = ~ colnames(otu_tax_f))
  main_level_col <- c()
  
  #loop through selected main_level to add color
  i <- 1
  
  for (i in 1:nrow(topx)) {
    top_loop <- pull(topx[, main_level])[i]
    # print(top_loop)
    
    # Subset the df with unique features and in function of the Means abundance
    
    temp <- otu_tax_f %>%
      dplyr::filter(!!as.name(main_level) == top_loop) %>%
      dplyr::arrange(desc(Mean))
    
    message(
      'Among the ',
      main_level,
      ' ' ,
      top_loop,
      ' ' ,
      length(unique(temp$plot_taxa)) ,
      ' features at ',
      sub_level,
      ' level will be plotted'
    )
    
    #define the color palette
    getPalette = colorRampPalette(brewer.pal(length(unique(temp$plot_taxa)), hues[i]), bias = color_bias)
    col <-
      getPalette(length(unique(temp$plot_taxa)) + 1)
    
    #store the darkestcolor of the palette to plot main_level legend graph
    if (length(col) >= 3) {
      main_level_col[i] <-   col[length(col)-1]
    } else {
      main_level_col[i] <-   col[length(col)]
    }  
    
    #loop among the sub_level feature and add color to it
    u <- 1
    for (u in 1:length(unique(temp$plot_taxa))) {
      print(unique(temp$plot_taxa)[u])
      
      #Add the color
      temp[which(temp$plot_taxa ==  unique(temp$plot_taxa)[u]), 'MyColors'] <-
        col[u + 1]
      
    }
    
    #create the dataframe storing the colors information
    df <- rbind(df, temp)
    
  }
  
  
  #Select features that don't belong to high abundance levels and assigned the black color to them
  unselescted <- otu_tax_f %>%
    filter(!(!!as.name(main_level) %in% pull(topx[, main_level])))
  
  unselescted$MyColors <- '#000000'
  
  #Merge the 2 dataframes
  df <- rbind(df, unselescted)
  
  #Order df to define the order that will be plotted, first the selected main_level, then the main_level in alphabetical order, then the mean abundance
  
  df <- df %>% ungroup %>%
    dplyr::arrange(desc(selected_top), !!as.name(main_level) , desc(Mean))
  
  
  #check if no features have been lost
  nrow(df) == nrow(otu_tax_f)
  
  #factor to keep order in the plot
  df$plot_taxa <-
    factor(df$plot_taxa, levels = unique(df$plot_taxa))
  
  
  #Transform df into long format by sample
  df_long <- melt(
    df,
    id = c("plot_taxa", "MyColors", main_level),
    measure.vars = meta[, sample_name],
    variable.name = sample_name
  )
  
  
  #Add exp_group to df_long
  df_long <- left_join(df_long, meta, by = sample_name)
  
  
  #Replace low abundance features at level X by "Others "
  
  df_long[,main_level] <- ifelse(df_long[,main_level] %in% pull(topx[, main_level]), 
                                 df_long[,main_level], 
                                 "Others ")
  
  
  if (differential_analysis &&
      length(unique(meta[[exp_group]])) != 2) {
    print(
      "The number of experimental group to test is not equal to 2, no statistical significance will appear in the legend"
    )
    mult_comp <- T
    
  }
  
  
  
  # Differential abundance analysis on sub_level using DESeq2 : 2 groups only
  if (differential_analysis &&
      length(unique(meta[[exp_group]])) == 2) {
    fam_glom <- tax_glom(ps_object, taxrank = sub_level)
    
    #remove OTU with 0 count across the dataset
    
    if (sum(rowSums(otu_table(fam_glom)) == 0) > 0) {
      fam_glom <-
        prune_taxa(rowSums(otu_table(fam_glom)) > 0, fam_glom)
      
    }
    
    diagdds = phyloseq_to_deseq2(fam_glom,  formula(paste("~", exp_group)))
    
    # Run DESeq2 analysis
    if (test == "Wald"){
      
      
      
      diag = DESeq(diagdds, test = test, fitType = fitType, sfType = sfType, 
                   betaPrior = betaPrior, quiet= quiet, 
                   minReplicatesForReplace = minReplicatesForReplace,
                   modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                   parallel = parallel)
      
    } else {
      
      diag = DESeq(diagdds, test = test, fitType = fitType, sfType = sfType, 
                   betaPrior = betaPrior, reduced = formula(reduced), quiet= quiet, 
                   minReplicatesForReplace = minReplicatesForReplace,
                   modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                   parallel = parallel)
      
    }
    
    # Get the differentially abundant features
    results <- results(diag)
    
    # Extract features with adjusted p-value below the threshold
    significant_features <- subset(results, padj < fdr_threshold)
    
    #Check if significant features were detected
    if (nrow(significant_features) == 0) {
      message("no significant feature detected at ", sub_level, " level.")
    } else {
      #We have significantly different taxa, we link them to their sub_level
      
      
      significant_features <-
        as.data.frame(cbind(significant_features,
                            tax_table(fam_glom)[rownames(tax_table(fam_glom)) %in% rownames(significant_features)]))
      
      #store the sig table to return it 
      significant_features_sub <- significant_features
      
      # Add information to df_long to mark differentially abundant features
      df_long$differential_abundance <- FALSE
      df_long$differential_abundance[df_long$plot_taxa %in% significant_features[, sub_level]] <-
        TRUE
      
      #Add stars to legends if sig_lab is defined as TRUE
      significant_features$stars <- ""
      
      if (sig_lab == T) {
        significant_features$stars <- symnum(
          significant_features$`padj`,
          symbols   = c("***", "**", "*", ""),
          cutpoints = c(0,  .001, .01, .05, 1),
          corr      = FALSE
        )
        
        #Add stars to the name of the taxa
        star_vec <-
          significant_features$stars[match(df_long$plot_taxa , significant_features[, sub_level])]
        star_vec[is.na(star_vec)]  <- ""
        df_long$plot_taxa <- paste0(df_long$plot_taxa, " ", star_vec)
        
        
      }
      
      #put significant features in bold
      df_long$legend_label <-
        ifelse(
          df_long$differential_abundance,
          paste0("<b>", df_long$plot_taxa, "</b>"),
          as.character(df_long$plot_taxa)
        )
      df_long$plot_taxa <- df_long$legend_label
      df_long$plot_taxa <-
        factor(df_long$plot_taxa, levels = unique(df_long$plot_taxa))
      
    }
    
  }
  
  
  # Differential abundance analysis on main_level using DESeq2 : 2 groups only
  if (differential_analysis &&
      length(unique(meta[[exp_group]])) == 2) {
    main_glom <- tax_glom(ps_object, taxrank = main_level)
    
    #remove OTU with 0 count across the dataset
    
    if (sum(rowSums(otu_table(main_glom)) == 0) > 0) {
      main_glom <-
        prune_taxa(rowSums(otu_table(main_glom)) > 0, main_glom)
      
    }
    
    
    diagdds_main = phyloseq_to_deseq2(main_glom,  formula(paste("~", exp_group)))
    
    # Run DESeq2 analysis
    if (test == "Wald"){
      
      diag_main = DESeq(diagdds_main, test = test, fitType = fitType, sfType = sfType, 
                        betaPrior = betaPrior, quiet= quiet, 
                        minReplicatesForReplace = minReplicatesForReplace,
                        modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                        parallel = parallel)
      
    }else {
      
      diag_main = DESeq(diagdds_main, test = test, fitType = fitType, sfType = sfType, 
                        betaPrior = betaPrior, reduced = formula(reduced), quiet= quiet, 
                        minReplicatesForReplace = minReplicatesForReplace,
                        modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                        parallel = parallel)
      
    }
    
    # Get the differentially abundant features
    results_main <- results(diag_main)
    
    # Extract features with adjusted p-value below the threshold
    significant_features_main <- subset(results_main, padj < fdr_threshold)
    
    #Check if significant features were detected
    if (nrow(significant_features_main) == 0) {
      message("no significant feature detected at ", main_level, " level.")
    } else {
      #We have significantly different taxa, we link them to their sub_level
      
      significant_features_main <-
        as.data.frame(cbind(significant_features_main,
                            tax_table(main_glom)[rownames(tax_table(main_glom)) %in% rownames(significant_features_main)]))
      

      # Add information to df_long to mark differentially abundant features
      df_long$differential_abundance_main <- FALSE
      df_long$differential_abundance_main[df_long[,main_level] %in% significant_features_main[, main_level]] <-
        TRUE
      
      #store the initial main_level names
      df_long$legend_label_main <- df_long[,main_level]
      
      #Add stars to legends if sig_lab is defined as TRUE
      significant_features_main$stars <- ""
      
      if (sig_lab == T) {
        significant_features_main$stars <- symnum(
          significant_features_main$`padj`,
          symbols   = c("***", "**", "*", ""),
          cutpoints = c(0,  .001, .01, .05, 1),
          corr      = FALSE
        )
        
        #Add stars to the name of the taxa
        star_vec_main <-
          significant_features_main$stars[match(df_long[,main_level], significant_features_main[, main_level])]
        star_vec_main[is.na(star_vec_main)]  <- ""
        df_long[,main_level] <- paste0(df_long[,main_level], " ", star_vec_main)
        
        #remove the stars column from the datafram 
        significant_features_main <-  significant_features_main[ ,-c(ncol(significant_features_main)) ]
        
      }
      
      #put significant features in bold
      
      
      df_long[,main_level] <-
        ifelse(
          df_long$differential_abundance_main,
          paste0("<b>", df_long[,main_level], "</b>"),
          as.character(df_long[,main_level])
        )
      #df_long[,main_level] <- df_long$legend_label_main
      df_long[,main_level] <-
        factor(df_long[,main_level], levels = unique(df_long[,main_level]))
      
    }
    
  }
  
  
  # Differential abundance analysis on sub_level using DESeq2 : multiple comparisons
  if (differential_analysis && mult_comp == T) {
    
    fam_glom <- tax_glom(ps_object, taxrank = sub_level)
    
    #remove OTU with 0 count across the dataset
    
    if (sum(rowSums(otu_table(fam_glom)) == 0) > 0) {
      fam_glom <-
        prune_taxa(rowSums(otu_table(fam_glom)) > 0, fam_glom)
      
    }
    
    diagdds = phyloseq_to_deseq2(fam_glom,  formula(paste("~", exp_group)))
    
    # Run DESeq2 analysis
    if (test == "Wald"){
      
      diag = DESeq(diagdds, test = test, fitType = fitType, sfType = sfType, 
                   betaPrior = betaPrior, quiet= quiet, 
                   minReplicatesForReplace = minReplicatesForReplace,
                   modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                   parallel = parallel)
      
    } else {
      
      diag = DESeq(diagdds, test = test, fitType = fitType, sfType = sfType, 
                   betaPrior = betaPrior, reduced = formula(reduced), quiet= quiet, 
                   minReplicatesForReplace = minReplicatesForReplace,
                   modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                   parallel = parallel)
      
    }
    
    comparisons <- combn(levels(meta[[exp_group]]), 2, simplify = FALSE)
    
    if (is.null(selected_comparisons) == F)  {
      
      # Filter only selected comparisons
      comparisons <- Filter(function(cmp) {
        any(sapply(selected_comparisons, function(sel) all(sel == cmp)))
      }, comparisons)
    }
    
    
    # Generate names for the results based on comparisons
    comparison_names <- sapply(comparisons, function(cmp) paste(cmp, collapse = "_vs_"))
    
    # Using correctly formatted contrasts based on results names
    results_list <- list()
    results_list <- lapply(comparisons, function(cmp) {
      contrast_vec <- c(exp_group, cmp[1], cmp[2])  # cmp[1] is the numerator, cmp[2] is the denominator
      tryCatch({
        res <- results(diag, contrast = contrast_vec)
        return(res)
      }, error = function(e) {
        message("Failed to compute results for comparison: ", paste(cmp, collapse = " vs "), "\nError: ", e$message)
        return(NULL)  
      })
    })
    
    # Set names based on comparisons
    results_list <- setNames(results_list, comparison_names)
    
    
    # Initialize an empty list to store significant features for each comparison
    significant_features_sub <- list()
    
    # Iterate over the results_list to process each comparison
    for (i in seq_along(results_list)) {
      res <- results_list[[i]]
      comparison_name <- names(results_list)[i]
      
      # Apply the significance threshold
      significant_features <- subset(res, padj < fdr_threshold)
      
      # Check if significant features were detected
      if (nrow(significant_features) == 0) {
        message("No significant feature detected for comparison ", comparison_name ," at the ", sub_level, " level.")
      } else {
        # Link significant features with their taxonomy information
        significant_features <- cbind(
          significant_features,
          tax_table(fam_glom)[rownames(tax_table(fam_glom)) %in% rownames(significant_features), , drop = FALSE]
        )
        
        # Convert to data frame and add to the list
        significant_features <- as.data.frame(significant_features)
        significant_features_sub[[comparison_name]] <- significant_features
      }
    }
    
    
    
    
    
  }
  
  # Differential abundance analysis on main_level using DESeq2 : multiple comparisons
  if (differential_analysis && mult_comp == T) {
    
    main_glom <- tax_glom(ps_object, taxrank = main_level)
    
    #remove OTU with 0 count across the dataset
    
    if (sum(rowSums(otu_table(main_glom)) == 0) > 0) {
      fam_glom <-
        prune_taxa(rowSums(otu_table(main_glom)) > 0, main_glom)
      
    }
    
    diagdds = phyloseq_to_deseq2(main_glom,  formula(paste("~", exp_group)))
    
    # Run DESeq2 analysis
    if (test == "Wald"){
      
      diag = DESeq(diagdds, test = test, fitType = fitType, sfType = sfType, 
                   betaPrior = betaPrior, quiet= quiet, 
                   minReplicatesForReplace = minReplicatesForReplace,
                   modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                   parallel = parallel)
      
    } else {
      
      diag = DESeq(diagdds, test = test, fitType = fitType, sfType = sfType, 
                   betaPrior = betaPrior, reduced = formula(reduced), quiet= quiet, 
                   minReplicatesForReplace = minReplicatesForReplace,
                   modelMatrixType = modelMatrixType, useT = useT, minmu = minmu,
                   parallel = parallel)
      
    }
    
    
    comparisons <- combn(levels(meta[[exp_group]]), 2, simplify = FALSE)
    
    if (is.null(selected_comparisons) == F)  {
      
      # Filter only selected comparisons
      comparisons <- Filter(function(cmp) {
        any(sapply(selected_comparisons, function(sel) all(sel == cmp)))
      }, comparisons)
    }
    
    # Generate names for the results based on comparisons
    comparison_names <- sapply(comparisons, function(cmp) paste(cmp, collapse = "_vs_"))
    
    # Using correctly formatted contrasts based on results names
    results_list <- list()
    results_list <- lapply(comparisons, function(cmp) {
      contrast_vec <- c(exp_group, cmp[1], cmp[2])  # cmp[1] is the numerator, cmp[2] is the denominator
      tryCatch({
        res <- results(diag, contrast = contrast_vec)
        return(res)
      }, error = function(e) {
        message("Failed to compute results for comparison: ", paste(cmp, collapse = " vs "), "\nError: ", e$message)
        return(NULL)  
      })
    })
    
    # Set names based on comparisons
    results_list <- setNames(results_list, comparison_names)
    
    
    # Initialize an empty list to store significant features for each comparison
    significant_features_main <- list()
    
    # Iterate over the results_list to process each comparison
    for (i in seq_along(results_list)) {
      res <- results_list[[i]]
      comparison_name <- names(results_list)[i]
      
      # Apply the significance threshold
      significant_features <- subset(res, padj < fdr_threshold)
      
      # Check if significant features were detected
      if (nrow(significant_features) == 0) {
        message("No significant feature detected for comparison ", comparison_name," at the ", main_level, " level.")
      } else {
        # Link significant features with their taxonomy information
        significant_features <- cbind(
          significant_features,
          tax_table(main_glom)[rownames(tax_table(main_glom)) %in% rownames(significant_features), , drop = FALSE]
        )
        
        # Convert to data frame and add to the list
        significant_features <- as.data.frame(significant_features)
        significant_features_main[[comparison_name]] <- significant_features
      }
    }
    
    
    
    
    
  }
  
  
  #Prepare the color vector
  MyColors <- df_long$MyColors
  names(MyColors) <- df_long$plot_taxa
  
  MyColors2 <- unique(df_long$MyColors)
  names(MyColors2) <- unique(df_long$plot_taxa)
  
  #add the black color to "Other_main_level"
  main_level_col[length(main_level_col)+1] <- '#000000'
  
  #Order the colors
  df_long[,main_level] <- factor(df_long[,main_level], levels = unique(df_long[,main_level]))
  
  vec1 <- unique(df_long[,main_level])
  vec2 <- c(pull(topx[,main_level]), paste0("Others"))
  core_text_vec1 <- gsub("(<[^>]*>|\\*|\\s+$)", "", vec1)
  core_text_vec1 <- trimws(core_text_vec1)
  order_index <- match( core_text_vec1, vec2)
  main_level_col <- main_level_col[order_index]
  names(main_level_col) <- as.character(vec1)
  
  #plot
  
if(mean_group == F) {  
  
  
  p <-
    ggplot(df_long, aes(
      x = !!as.name(sample_name),
      y = value,
      fill = plot_taxa
    )) +
    geom_bar(stat = "identity", width = 0.85) +
    ylab("Relative abundance (%)\n") +
    guides(fill = guide_legend(reverse = FALSE, title = sub_level, order = 2)) +
    theme(
      line = element_line(colour = "black", linewidth = .5),
      text = element_text(size = 9),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linewidth = .5)
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      legend.position = "right",
      legend.box = "vertical",
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = x_axis_size
      ),
      text = element_text(size = text_size),
      legend.text = element_markdown(size = legend_size)
    ) +
    geom_bar(aes(alpha = df_long[, main_level]), stat = "identity", show.legend = TRUE) +
    scale_alpha_manual(
      values = rep(1, length(unique(df_long[, main_level]))),
      guide = guide_legend(order = 1,override.aes = list(fill = main_level_col)),
      name = main_level
    ) +
    scale_fill_manual('plot_taxa', values = MyColors2)
  p <- p + facet_wrap( ~ df_long[[exp_group]] , scales  = "free_x", nrow = n_row, ncol = n_col)
  
  p
  
} else {
  
  
  df_long <- df_long %>% 
    group_by(!!as.name(exp_group), plot_taxa, !!as.name(main_level)) %>%
    reframe(
      n = n(),
      sum = sum(as.double(value)))
  df_long <- data.frame(df_long)
  
i<-1
for (i in 1: length(unique(df_long[,exp_group]))) {
  
  df_long[,"sum"][df_long[,exp_group] == unique(df_long[,exp_group] )[i]] <- df_long[,"sum"][df_long[,exp_group] == unique(df_long[,exp_group] )[i]]/count(meta[,exp_group] == unique(meta[,exp_group])[i])
  
}

colnames(df_long)[colnames(df_long) == "sum"] <- "value"  
  
  
  
  #plot
  p <-
    ggplot(df_long, aes(
      x = !!as.name(exp_group),
      y = value,
      fill = plot_taxa
    )) +
    geom_bar(stat = "identity", width = 0.85) +
    ylab("Relative abundance (%)\n") +
    guides(fill = guide_legend(reverse = FALSE, title = sub_level, order = 2)) +
    theme(
      line = element_line(colour = "black", linewidth = .5),
      text = element_text(size = 9),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linewidth = .5)
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      legend.position = "right",
      legend.box = "vertical",
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = x_axis_size
      ),
      text = element_text(size = text_size),
      legend.text = element_markdown(size = legend_size)
    ) +
    geom_bar(aes(alpha = df_long[, main_level]), stat = "identity", show.legend = TRUE) +
    scale_alpha_manual(
      values = rep(1, length(unique(df_long[, main_level]))),
      guide = guide_legend(order = 1,override.aes = list(fill = main_level_col)),
      name = main_level
    ) +
    scale_fill_manual('plot_taxa', values = MyColors2)
  p <- p + facet_wrap( ~ df_long[[exp_group]] , scales  = "free_x", nrow = n_row, ncol = n_col)
  
  p
  
  
  
}
  
  
  
  if (differential_analysis == T ) {
    
    return(
      list(
        significant_table_main = significant_features_main,
        significant_table_sub = significant_features_sub,
        plot = p
      )
    )
    
  } else {
    return(list(plot = p))
  }
  
}  


