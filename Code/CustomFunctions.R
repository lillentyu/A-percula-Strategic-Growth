

##################################################################
###################### plotPCA modifications ######################
##################################################################

### Updated 'plotPCA' (from DESeq2) to save all PCs for plasticity measures

plotPCA_allPCs <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object)) # Variance estimates for each row (column) in a matrix.
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))] # orders genes from largest -> smallest and pulls ntop (default 500) genes
  pca <- prcomp(t(assay(object)[select, ]), center = TRUE, scale. = FALSE) # performs PCA on ntop (default 500) from dataset
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(group = group, intgroup.df, name = colnames(object), pca$x)
  if (returnData) {
    attr(d, "percentVar") <- percentVar # add [1:2] if only want 1st two PC variances
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
    geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 
                                                        100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 
                                                                                                            100), "% variance")) + coord_fixed()
}


################################################################################
################### Sample plasticity based on PC distances ####################
############## Written by Colleen B. Bove (colleenbove@gmail.com) ##############
########################## Laste update: 10 Jan 2023 ###########################
################################################################################

## Required packages to run:
# dplyr


## PC distances can be calculated from prcomp() objects or from plotPCA()


## To run the function, enter the following objects:
# PCAplast(pca = XXX, # the PCA dataframe containing the PCA eigenvalues
#          data = XXX, # the condition/treatment data corresponding to samples
#          sample_ID = "XXX", # the name of column that provide unique ID per sample (if blank, will pull rownames for this)
#          num_pca =  "XXX", # the number of PCAs to include in analysis (default is 'all', but you can specify another number with a minimum of 2 PCAs)
#          control_col = "XXX", # what the 'treatment' column is called
#          control_lvl = "XXX", # control level of the treatment. If blank, a control mean per control level is assumed
#          group = "XXX") # the grouping column (i.e., colony). If blank, will assume control level grouping only!



########################################################################################

PCAplast <- function(pca, data, sample_ID = NA, num_pca = "all", control_col, control_lvl = "none", group = NA) {
  
  # rename the user input info
  pca_df <- pca
  data_df <- data
  control_name <- control_col
  control_lvl <- control_lvl
  group_col <- group
  
  # ifelse statement to pull PCs from dataframe or prcomp objects
  if(class(pca_df) == "prcomp"){
    pca_dist <- pca_df$x # grab PC distances from prcomp() object
  } else {
    pca_dist <- pca_df # grab PC distances from data.frame object
  }
  
  
  # check for correct number of PCAs provided
  if(class(num_pca) == "numeric") {
    if(num_pca < 2) { # will throw error if too few PCAs requested
      stop("please select more than 2 PCs to calculate distance")
    } 
  }
  
  if(class(num_pca) == "numeric") {
    if(num_pca > (data.frame(pca_dist) %>% dplyr::select(starts_with("PC")) %>% ncol())) { # will throw error if too many PCAs requested
      stop(paste(num_pca, "PCs requested for calculation, but only", (pca_dist %>% dplyr::select(starts_with("PC")) %>% ncol()), "PCs available. Select appropriate number of PCAs for calculation."))
    } 
  }
  
  
  # oder the dataframe to ensure correct pairing after calculating distance 
  if(sample_ID %in% colnames(data_df)){
    data_df <- data_df[match(row.names(pca_dist), data_df[[sample_ID]]),]
  } else {
    data_df <- data_df[order(as.numeric(row.names(data_df))),]
  }
  
  
  # combine the datasets
  dist_df <- cbind(data_df, pca_dist) 
  
  
  # if there is no control level, modify so function pulls all levels
  if(control_lvl == "none") {
    control_lvl = unique(dist_df[[control_name]])
  } else {
    control_lvl = control_lvl
  }
  
  
  # make dataframe of control grouping only
  if(!is.na(group_col)) { # calculate mean per grouping ID (if provided)
    mean_control <- dist_df %>%
      filter(dist_df[[control_name]] == list(control_lvl)[[1]]) %>% 
      rename_with(tolower) %>% # renames all pc's with lowercase 'PC' (just to differentiate from all sample PCs)
      dplyr::select(colnames((dist_df %>% rename_with(tolower))[tolower(group_col)]), starts_with("pc")) %>%
      group_by_at(vars(tolower(group_col))) %>% 
      summarise_if(is.numeric, mean)
    
    
    # add the control PCA values to treatment samples per grouping
    dist_df2 <- left_join(dist_df, mean_control, by.x = group_col, by.y = tolower(group_col))
    
  } else { # calculate mean per control treatment
    mean_control <- dist_df %>%
      #filter(dist_df[[control_name]] == list(control_lvl)[[1]]) %>% 
      rename_with(tolower) %>% # renames all pc's with lowercase 'PC' (just to differentiate from all sample PCs)
      dplyr::select(colnames((dist_df %>% rename_with(tolower))[tolower(control_name)]), starts_with("pc")) %>% # select just the PCs 
      dplyr::select(-pco2) %>% 
      #group_by() 
      group_by_at(vars(tolower(control_name))) %>% 
      summarise_if(is.numeric, mean)
    
    # add the control PCA values to all samples 
    dist_df2 <- merge(dist_df, mean_control, by.x = control_col, by.y = tolower(control_col), all.x = TRUE)
    
  }
  
  
  
  if(sample_ID %in% colnames(data_df)){
    data_df <- data_df[match(row.names(pca_dist), data_df[[sample_ID]]),]
  } else {
    data_df <- data_df[order(as.numeric(row.names(data_df))),]
  }
  
  
  # again, reorder data
  if(sample_ID %in% colnames(data_df)){
    dist_df2 <- dist_df2[order(dist_df2[[sample_ID]]),]
  } else {
    rownames(dist_df2) <- rownames(data_df)
    dist_df2 <- dist_df2[order(as.numeric(row.names(dist_df2))),]
  }
  
  
  ### Calculate sample (PCA) distances from control (pca) using all PCAs
  # make dataframe to populate with pca distances
  full_calc_dist <- data.frame(control_name = dist_df2[control_name])
  
  if(num_pca == "all") {
    ## forloop that will calculate distances between control and sample for all PCs (n will be total number)
    for(n in 1:(dist_df %>% dplyr::select(starts_with("PC")) %>% ncol())){
      # makes the PCA column name for control (lowercase) and sample (uppercase)
      PC_col <- paste0("PC", n)
      pc_col <- paste0("pc", n)
      
      # pulls the PC column for control (lowercase) and sample (uppercase)
      PCx <- dist_df2[PC_col]
      pcx <- dist_df2[pc_col]
      
      pca_calc_dist <- data.frame((PCx - pcx)^2) # calculates the distance between 2 PCs
      full_calc_dist <- cbind(full_calc_dist, pca_calc_dist) # add that distance to running dataframe
    }
  } else {
    ## forloop that will calculate distances between control and sample for SPECIFIED # of PCs (n will be total number)
    for(n in 1:as.numeric(num_pca)){
      # makes the PC column name for control (lowercase) and sample (uppercase)
      PC_col <- paste0("PC", n)
      pc_col <- paste0("pc", n)
      
      # pulls the PC column for control (lowercase) and sample (uppercase)
      PCx <- dist_df2[PC_col]
      pcx <- dist_df2[pc_col]
      
      pca_calc_dist <- data.frame((PCx - pcx)^2) # calculates the distance between 2 PCs
      full_calc_dist <- cbind(full_calc_dist, pca_calc_dist) # add that distance to running dataframe
    }
  }
  
  
  ## final distance calculation (adds all PCA distances and takes squareroot)
  distance <- full_calc_dist %>% 
    mutate(dis_sum = rowSums(across(where(is.numeric)))) %>% 
    mutate(dist = sqrt(dis_sum)) %>% 
    dplyr::select(matches("dist"))
  
  
  ## combine the calculated distance with the metadata and remove controls for final dataframe
  dist_df <- data_df %>% 
    bind_cols(distance) %>% 
    filter(!is.na(dist)) 
  
  ## removes the control levels
  if(length(control_lvl) > 1){
    dist_df <- dist_df
  } else {
    dist_df <- dist_df %>%
      filter(dist_df[[control_name]] != control_lvl)  %>% 
      droplevels()
  }
  
}


#################################################################
###### Session information from last update

# R version 3.6.3 (2020-02-29)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS  10.16
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggplot2_3.3.6 dplyr_1.0.8  
# 
# loaded via a namespace (and not attached):
#   [1] rstudioapi_0.13  magrittr_2.0.3   munsell_0.5.0    tidyselect_1.1.1 colorspace_2.0-3 R6_2.5.1        
# [7] rlang_1.0.4      fansi_1.0.3      tools_3.6.3      grid_3.6.3       gtable_0.3.0     xfun_0.29       
# [13] tinytex_0.40     utf8_1.2.2       cli_3.3.0        DBI_1.1.3        withr_2.5.0      ellipsis_0.3.2  
# [19] digest_0.6.29    assertthat_0.2.1 tibble_3.1.7     lifecycle_1.0.1  crayon_1.5.1     farver_2.1.1    
# [25] purrr_0.3.4      vctrs_0.4.1      glue_1.6.2       labeling_0.4.2   compiler_3.6.3   pillar_1.8.0    
# [31] scales_1.2.0     generics_0.1.3   pkgconfig_2.0.3





################################################################################

### Custom functions for WGCNA pipeline

## Function to identify and select significant modules from WGCNA 
sigMods <- function(matrix, alpha = 0.0500001) {
  sig_mat <- (matrix < alpha)
  ind <- which(sig_mat == 'TRUE', arr.ind = TRUE)
  moduleCols <- gsub("ME", "", rownames(sig_mat)[ind[,"row"]])
  whichTrait <- colnames(sig_mat)[ind[,"col"]]
  out <- data.frame(moduleCols = as.data.frame(moduleCols),
                    whichTrait = as.data.frame(whichTrait))
  return(out)
}


## Function to do all the module/trait plots and GO MWU file creations
ModWGCNA <- function(wgcna_data, traitMods, traits, geneModuleMembership, MEs, vsdWG, corrplotpath, heatplotpath, goMWUpath) {
  
  ModWGCNA_out <- list()
  
  ########## Gene Significance
  
  # Create blank dataframes to populate with below forloop
  GSPvalue_df <- data.frame(matrix(nrow = ncol(datt), ncol = 0))
  GTSvalue_df <- data.frame(matrix(nrow = ncol(datt), ncol = 0))
  
  # Forloop to calculate gene trait significance and gene significance p values for each significant module/trait pairing
  for (t in 1:nrow(traitMods)) { #
    #t = 1
    ## trait selection/specification
    Trait <- traitMods$whichTrait[t] # pulls trait name
    selTrait <- as.data.frame(traits[,Trait]) # selects the trait data that matches the trait name
    names(selTrait) <- Trait # renames the selected trait data with trait name
    
    ## calculate trait correlation p values
    geneTraitSignificance <- as.data.frame(cor(datt, selTrait, use = "p")) # calculates correlation between GE data and trait data
    GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples)) # Calculates Student asymptotic p-value
    names(geneTraitSignificance) <- names(selTrait) # rename geneTraitSignificance column with trait name
    names(GSPvalue) = paste("p.GS.", names(selTrait), sep = "") # rename GSPvalue colmn with trait name
    
    ## add the new trait columns to the existing dataframes
    GTSvalue_df <- cbind(GTSvalue_df, geneTraitSignificance)
    GSPvalue_df <- cbind(GSPvalue_df, GSPvalue)
    
  }
  
  ModWGCNA_out[["GTSig"]] <- GTSvalue_df
  ModWGCNA_out[["GSPvalue"]] <- GSPvalue_df
  
  
  ########## Correlation Plots
  
  # Create blank list to populate with correlation plots below
  plots <- list() # create an empty list to populate with plots
  
  # Forloop to create correlation plots per trait/module pairing (also saves as png)
  for (r in 1:nrow(traitMods)) {
    
    module <- traitMods$moduleCols[r]
    selectTrait <- traitMods$whichTrait[r]
    
    column <- match(module, modNames)
    moduleGenes <- (moduleColors == module)
    combo <- paste(module, selectTrait, sep = "; ")
    
    x <- abs(geneModuleMembership[moduleGenes, column])
    y <- abs(GTSvalue_df[moduleGenes, selectTrait])
    
    # plot the correlations
    plot <- ggplot() +
      theme_classic() +
      geom_point(aes_string(x = x, y = y), shape = 1, colour = module, size = 2) +
      xlab(paste(module, "module membership")) +
      ylab(paste("GS for", traitMods$whichTrait[r])) +
      stat_cor(aes_string(x = abs(geneModuleMembership[moduleGenes, column]), y = abs(GTSvalue_df[moduleGenes, selectTrait])), method = "pearson")
    ggsave(paste(corrplotpath, combo, ".png", sep = ""), width = 5, height = 4.3)
    
    plots[[combo]] <- plot
    
  }
  
  ModWGCNA_out[["CorrPlots"]] <- plots  
  
  
  ########## Module Heatmaps
  
  # Create blank list to populate with correlation plots below
  plots2 <- list() # create an empty list to populate with plots
  
  # Forloop to cycle through trait/module pairings to produce expression heatmaps and eigengene bar plots (also saves png)
  for (r in 1:nrow(traitMods)) {
    
    module <- traitMods$moduleCols[r] # pulls module colour
    selectTrait <- traitMods$whichTrait[r] # identified trait name
    combo2 <- paste(module, selectTrait, sep = "; ") # create name of module/trait pair
    module2 <- paste("ME", module, sep="") # name of module with 'ME' added for pulling correct data
    
    MEs$sample <- rownames(MEs) # name column for sample ID (from rownames) 
    ME <- MEs[, c("sample", module2)] # subset MEs for selected module only
    names <- rownames(arrange(ME, ME[,2])) #rearrange ME in ascending order and extract names
    datExpr <- cbind(wgcna_data, ME[,2]) # add mean expression to expression dataframe
    ME <- arrange(ME, ME[,2]) # order expression dataframe in ascending order
    
    ## Heatmap of isogroup expression per sample
    heatmap <- t(scale(datExpr[, moduleColors == module])) %>% 
      as.data.frame() %>%
      rownames_to_column("isogroup") %>%
      pivot_longer(-c(isogroup), names_to = "samples", values_to = "expr") %>%
      mutate(samples = fct_relevel(samples, names)) %>%
      ggplot(aes(x = samples, y = isogroup, fill = expr)) + 
      ggtitle(paste(selectTrait, "; ", module, " module", sep = "")) +
      theme_cowplot() +
      theme(axis.text.y = element_blank(), axis.text.x = element_blank(), legend.position = "left", legend.title = element_blank(), axis.title = element_blank(), legend.text=element_text(size = 6)) +
      scale_fill_gradient2(low = "#5e3c99", high = "#e66101", mid = "gray95") +
      geom_tile()
    
    ## Barplot of eigengene expression per sample
    barplot <- ggplot(data = ME) +
      theme_cowplot() +
      geom_hline(aes_string(yintercept = 0)) +
      theme(text = element_text(size = 10), axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1)) +
      geom_bar(aes_string(x = factor(ME$sample, levels = c(names)), y = ME[,2]), stat = "identity", fill = module, colour = "black", size = 0.2) +
      labs(x = "", y = "eigengene expression")
    
    ## Combine and save heatmap and bar plot
    heatbar <- cowplot::plot_grid(heatmap, barplot, nrow = 2, rel_heights = c(1, 1))
    ggsave(paste(heatplotpath, combo2, ".png", sep = ""), width = 8, height = 6)
    
    plots2[[combo2]] <- heatbar
    
  }
  
  ModWGCNA_out[["Heatmaps"]] <- plots2 
  
  
  ########## GO and KOG Files
  
  
  
  MEs <- MEs[,-length(MEs)] # removed the last column of the dataframe that should correspond with the sample IDs
  
  # calculating module memberships for all genes for all modules
  allkME <- as.data.frame(signedKME(datt, MEs))
  names(allkME) <- gsub("kME", "", names(allkME))
  
  # Forloop to cycle through trait/module pairings to produce GO ouput CSVs
  for (r in 1:nrow(traitMods)) {
    r = 1
    module <- traitMods$moduleCols[r] # pulls module colour
    selectTrait <- traitMods$whichTrait[r] # identified trait name
    
    # Saving data for Fisher-MWU combo test (GO_MWU)
    inModuleBinary <- as.numeric(moduleColors == module)
    combo3 <- data.frame("gene" = row.names(vsdWG), "Fish_kME" = allkME[,module] * inModuleBinary)
    write.csv(combo3, file = paste(goMWUpath, selectTrait, "_",module, ".csv", sep = ""), row.names = FALSE, quote = FALSE)
    
    # saving it also as a binary version
    combo4 <- combo3 %>% 
      mutate(binary = ifelse(Fish_kME > 0, 1, 0)) %>% 
      dplyr::select(-Fish_kME)
    write.csv(combo4, file = paste(goMWUpath, module, "_binary", ".csv", sep = ""), row.names = FALSE, quote = FALSE)
  }
  
  return(ModWGCNA_out)
  
}

