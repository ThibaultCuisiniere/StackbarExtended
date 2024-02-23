#' Clustered stackbar coupled with DESeq2
#'
#' Create a plot of taxa relative abundance, clustered by phylogeny with the possibility to highlight the differentially abundant features.
#' 
#' @param ps_object A [phyloseq-class] object. Must have a [sample_data], [otu_table], and [tax_table] components.
#' @param exp_group A column in the [sample_data] containing the experimental group information.
#' @param subset_group Default NULL. Groups among the 'exp_group' column from the [sample_data] object to subset.
#' @param sample_name Name of the column in [sample_data] containing the unique sample identifier.
#' @param main_level  Default 'Phylum'. Level present in the [tax_table] to which taxa will be clustered.
#' @param sub_level Default 'Family'. Level present in the [tax_table] to which taxa will be plotted and analyzed.
#' @param threshold Default 1%, threshold to regroup taxa with lower relative abundance into the 'other' groups
#' @param n_phy Default 4, number of main_level to plot. Same number of colors must be given in the 'hues' parameter.
#' @param hues Color used to represent main_level, should be the same number than n_phy parameter. See [colorRampPalette].
#' @param color_bias Define the gradiant of shades among the colors. See [colorRampPalette].
#' @param text_size Control the text size of the graph. See [ggplot2].
#' @param legend_size Control the legend text size of the graph. see ggplot2. See [ggplot2].
#' @param x_axis_size Control the x axis text size of the graph. see ggplot2. See [ggplot2].
#' @param differential_analysis Default TRUE. Whether or not use DESeq2 to do a differential abundance analysis. See [DESeq2-package].
#' @param test Default 'Wald'. See [DESeq2-package].
#' @param sig_lab Default TRUE. Wheter add stars after taxa name reflecting statistical significance. 
#' @param fitType Default "parametric". See [DESeq2-package].
#' @param sfType Default "ratio". See [DESeq2-package].
#' @param betaPrior Default FALSE.See [DESeq2-package].
#' @param reduced Default FALSE. See [DESeq2-package].
#' @param quiet Default FALSE. See [DESeq2-package].
#' @param minReplicatesForReplace Default '7'. See [DESeq2-package].
#' @param modelMatrixType Default "standard". See [DESeq2-package].
#' @param useT Default FALSE. See [DESeq2-package].
#' @param minmu See [DESeq2-package].
#' @param parallel Default FALSE. See [DESeq2-package].
#' @param BPPARAM See [DESeq2-package].
#' @param ... additionnal parameters passed into the \code{\link{phyloseq_to_deseq2}} functions, see [Phyloseq].
#' 
#'
#'
#' @return A [list] containing  [ggplot2] '$plot' object and, if differential_analysis = TRUE, a [dataframe] from \code{\link{results}} DESeq2 function output of significant features '$significant_table'.
#'
#' @examples
#' 
#' \dontrun{
#' data(GlobalPatterns)
#' ps_unfiltered <- GlobalPatterns
#' 
#' my_plot <- plot_gut_microbiota(ps_object = ps_unfiltered,
#' exp_group = 'SampleType',
#' subset_group = c("Feces", "Skin"),
#' sample_name = 'X.SampleID',
#' differential_analysis = T,
#' sig_lab = T,
#' fdr_threshold = 0.8)
#' 
#' print(my_plot$plot)
#' print(my_plot$significant_table)
#'
#'}
#'
#' @export
#' 
#' @seealso [phyloseq]
#' @seealso [DESeq2-package]
#' @seealso [ggplot2]
#' 

#' @author Thibault Cuisiniere, Marco Constante, Manuela M. Santos
#'
#'
#' @import phyloseq
#' @import ggplot2
#' @import ggtext
#' @import tibble
#' @import tidyverse
#' @import reshape2
#' @import RColorBrewer
#' @import forcats
#' @import DESeq2
#'


plot_gut_microbiota <- function(ps_object = ps_unfiltered,
                                exp_group = 'group',
                                subset_group = NULL,
                                sample_name = 'SampleID',
                                main_level = 'Phylum',
                                sub_level = 'Family',
                                threshold = 1,
                                n_phy = 4,
                                hues = c("Oranges", "Greens", "Blues", "Purples"),
                                color_bias = 2,
                                text_size = 9,
                                legend_size = 7,
                                x_axis_size = 8,
                                differential_analysis = T,
                                test = c("Wald", "LRT")[1],
                                fitType = c("parametric", "local", "mean", "glmGamPoi")[1],
                                sfType = c("ratio", "poscounts", "iterate")[1],
                                betaPrior = FALSE,
                                reduced = FALSE,
                                quiet = FALSE,
                                minReplicatesForReplace = 7,
                                modelMatrixType = c("standard", "expanded")[1],
                                useT = FALSE,
                                minmu = if (fitType == "glmGamPoi") 1e-06 else 0.5,
                                parallel = FALSE,
                                BPPARAM = bpparam(),
                                fdr_threshold = 0.05,
                                sig_lab = T,
                                ...) {
  
  
  
  # Validate inputs
  if (!("phyloseq" %in% class(ps_object))) stop("ps_object must be a phyloseq-class object.")
  if (!main_level %in% colnames(tax_table(ps_object))) stop("main_level column not found in the tax_table")
  if (!sub_level %in% colnames(tax_table(ps_object))) stop("sub_level column not found in the tax_table")
  if (!exp_group %in% names(sample_data(ps_object))) stop("exp_group column not found in the sample_data.")
  
  
  #Assure taxa are row
  if (!taxa_are_rows(ps_object)) {
    ps_object <- t(ps_object)
  }

  #Subset the samples belonging to the groups selected
  if (is.null(subset_group) == F) {
    keep_samples = as.character(get_variable(ps_object, exp_group)) %in% subset_group
    ps_object = prune_samples(keep_samples, ps_object)

  }

  #transform in proportional
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

  ##keep the position of the column storing the abundance of the taxa in the otu_tax_f object
  position <- ncol(tax) + 1:(ncol(otu_tax_f) - ncol(tax))

  message(paste0('\n', length(position), ' samples are analyzed \n'))

  #Create 2 vectors in the dataframe :
  #1 storing the abundance mean of the taxa
  #2 T or F this taxa is above the define thrshold of filtering

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


  #Create a vector, TRUE when the phylum belong to the top X, the ones we want to plot
  otu_tax_f <- mutate(otu_tax_f,
                      selected_top = !!as.name(main_level) %in% pull(topx[, main_level]))

  #Create a column to create the legend of the graphs :
  #1 higher than the threshold and in the main_level selected -> name of the sub_level
  #1.1 if sub_level is NA -> Other_main_level_name

  #2 belong to the main_level selected, but sub_level bellow the threshold -> other_phylum name

  #3 all the others -> "Other"

  #replace NA


  otu_tax_f <- otu_tax_f %>%
    mutate(
      plot_taxa = case_when(
        high_abundance == TRUE  &
          selected_top == TRUE & !!as.name(sub_level) == "Unkwown"  ~
          paste0("Others_", !!as.name(main_level)),
        high_abundance == TRUE  &
          selected_top == TRUE  &
          !!as.name(sub_level) != "Unkwown"  ~ !!as.name(sub_level),
        high_abundance == FALSE &
          selected_top == TRUE  ~ paste0("Others_",!!as.name(main_level)),
        selected_top   == FALSE ~ paste0("Others"),
        high_abundance == TRUE  &
          selected_top == TRUE &
          is.na(!!as.name(sub_level)) == TRUE  ~ paste0("Unknown_",!!as.name(main_level))
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


  #initialize the dataframe
  df <-
    as_tibble(matrix(nrow = 0, ncol = length(colnames(otu_tax_f))),
              .name_repair = ~ colnames(otu_tax_f))

  #loop through selected phylum to add color
  i <- 1

  for (i in 1:nrow(topx)) {
    top_loop <- pull(topx[, main_level])[i]
    # print(top_loop)

    # Subset the df with unique features and in function of the Means abundance

    temp <- otu_tax_f %>%
      dplyr::filter(!!as.name(main_level) == top_loop) %>%
      dplyr::arrange(desc(Mean))

    message(
      'Among the ', main_level, ' ' ,top_loop,' ' ,
      length(unique(temp$plot_taxa)) ,
      ' features at ',sub_level,' level will be plotted'
    )

    #define the color palette
    getPalette = colorRampPalette(brewer.pal(length(unique(temp$plot_taxa)), hues[i]), bias = color_bias)
    col <-
      getPalette(length(unique(temp$plot_taxa)) + 1)


    #loop among the low level feature and add color to it
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


  #Selected features that don't belong to high abundances level and assigned the black color to them
  unselescted <- otu_tax_f %>%
    filter(!(!!as.name(main_level) %in% pull(topx[, main_level])))

  unselescted$MyColors <- '#000000'

  #Merge the 2 dataframe
  df <- rbind(df, unselescted)

  #order df to define the order that will be plotted, first the selected main_level, then the main_level in alphabetical order, then the mean

  df <- df %>% ungroup %>%
    dplyr::arrange(desc(selected_top),!!as.name(main_level) , desc(Mean))


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

  
  if (differential_analysis &&
      length(unique(meta[[exp_group]])) != 2) {
    
    print("The number of experimental group to test is not equal to 2, no differential abundance testing will be perform")
    differential_analysi <- F
    
  }
  

  # Differential abundance analysis using DESeq2 : 2 groups only
  if (differential_analysis &&
      length(unique(meta[[exp_group]])) == 2) {
    fam_glom <- tax_glom(ps_object, taxrank = sub_level)

    #ps_object is already filtered (or not but contains only 2 groups) with user defined group.

    #remove OTU with 0 count accros the dataset

    if (sum(rowSums(otu_table(ps_object)) == 0) > 0) {
      ps_object <-
        prune_taxa(rowSums(otu_table(ps_object)) > 0, ps_object)

    }

    # Coerce count data to matrix of integers
    countData = round(as(otu_table(fam_glom), "matrix"), digits = 0)
    colData = data.frame(sample_data(fam_glom))
    # Create the DESeq data set, dds.

    diagdds = phyloseq_to_deseq2(fam_glom,  formula(paste("~", exp_group)), ...)

    # Run DESeq2 analysis
    diag = DESeq(diagdds, test = test, fitType = fitType)

    # Get differentially abundant features
    results <- results(diag)

    # Extract features with adjusted p-value below a threshold
    significant_features <- subset(results, padj < fdr_threshold)

    #We have significant ASV, we want to link them to the level we are interested (family)

    significant_features <- as.data.frame(cbind(significant_features,
                                                tax_table(fam_glom)[rownames(tax_table(fam_glom)) %in% rownames(significant_features)]))

    # Add information to df_long to mark differentially abundant features
    df_long$differential_abundance <- FALSE
    df_long$differential_abundance[df_long$plot_taxa %in% significant_features[, sub_level]] <-
      TRUE

    #Add stars to legends if sig_lab is defined as TRUE
    significant_features$stars <-""

    if(sig_lab == T) {

      significant_features$stars <- symnum(significant_features$`padj`,
                                           symbols   = c("***","**","*",""),
                                           cutpoints = c(0,  .001,.01,.05,1),
                                           corr      = FALSE
      )

      #Add stars to the name of the taxa
      star_vec <- significant_features$stars[match(df_long$plot_taxa , significant_features[, sub_level])]
      star_vec[is.na(star_vec)]  <-""
      df_long$plot_taxa <- paste0(df_long$plot_taxa, " ", star_vec )


    }

    #put in bold significant features
    df_long$legend_label <-
      ifelse(
        df_long$differential_abundance,
        paste0("<b>", df_long$plot_taxa, "</b>"),
        as.character(df_long$plot_taxa)
      )
    df_long$plot_taxa <-   df_long$legend_label
    df_long$plot_taxa <-
      factor(df_long$plot_taxa, levels = unique(df_long$plot_taxa))

  }


  #Prepare the color vector
  MyColors <- df_long$MyColors
  names(MyColors) <- df_long$plot_taxa
  unique(df_long$plot_taxa)
  unique(MyColors)

  MyColors2 <- unique(df_long$MyColors)
  names(MyColors2) <- unique(df_long$plot_taxa)


  p <-
    ggplot(df_long, aes(
      x = !!as.name(sample_name),
      y = value,
      fill = plot_taxa
    )) +
    geom_bar(stat = "identity", width = 0.85) +
    #    xlab("\nMouse") +
    ylab("Relative abundance (%)\n") +
    #    guides(fill=FALSE) +
    #    theme(legend.position="bottom") +
    #    guides(fill=guide_legend(nrow=2))
    guides(fill = guide_legend(reverse = FALSE, title = "")) +
    theme(
      line = element_line(colour = "black", linewidth = .5),
      text = element_text(size = 9),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      #     plot.title = element_text(size = 6, family="Arial"),
      panel.border = element_blank(),
      #axis.text.x = element_text(angle=90, vjust=1),
      axis.line = element_line(colour = "black", linewidth = .5)
    )
  p <- p + theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = x_axis_size
      ),
      text = element_text(size = text_size),
      # legend.text=element_text(size=legend_size),
      legend.text = element_markdown(size = legend_size)
    )
  #p<-p+ scale_fill_manual(values=MyColors)
  p <- p + scale_fill_manual('plot_taxa', values = MyColors2)

  p <- p + facet_wrap( ~ df_long[[exp_group]] , scales  = "free_x")

  p




  if (differential_analysis == T) {
    return(list(significant_table = significant_features, plot = p))

  } else {
    return(list(plot = p))
  }

}
