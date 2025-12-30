# Large part of the code regarding calculation and illustration in Figs.2, 4, and 5 is adapted "Tutorial for R microeco package (v1.16.0)"
# which can be accessed from: https://chiliubio.github.io/microeco_tutorial/index.html 

#Input data ####
library(microeco)
setwd('D:\\XX')
  
# Environmental variables
ENV.raw <- read.csv('XX.csv',row.names = 'Sites')

#Proka
  
  Proka.raw <- read.csv('xx.csv');rownames(Proka.raw) <- paste0('OTU_',c(1:nrow(Proka.raw)))
  # OTU table
  Proka.otu <- Proka.raw[,-c(1:8)];head(Proka.otu)
  # Taxa table
  Proka_taxonomy_table <- Proka.raw[,c(1:8)];head(Proka_taxonomy_table)
  
  Proka.tree <- read.tree('xx.tre')
  Proka.fasta <- seqinr::read.fasta("xx.fasta")

  library(magrittr)
  Proka_taxonomy_table %<>% tidy_taxonomy #rm NA, and tidy the format
  head(Proka_taxonomy_table)
  
#Euka
  
  Euka.raw <- read.csv('xx.csv');rownames(Euka.raw) <- paste0('OTU_',c(1:nrow(Euka.raw)))
  # OTU table
  Euka.otu <- Euka.raw[,-c(1:8)];head(Euka.otu)
  # Taxa table
  Euka_taxonomy_table <- Euka.raw[,c(1:8)];head(Euka_taxonomy_table)
  
  Euka.tree <- read.tree('xx.tre')
  Euka.fasta <- seqinr::read.fasta("xx.fasta")
  
  library(magrittr)
  Euka_taxonomy_table %<>% tidy_taxonomy #rm NA, and tidy the format
  head(Euka_taxonomy_table)
  
# Define microtable
  Proka.mt <- microtable$new(sample_table = ENV.raw
                       , otu_table = Proka.otu
                       , tax_table = Proka_taxonomy_table
                       , rep_fasta = Proka.fasta
                       , phylo_tree = Proka.tree
  )
  Proka.raw <- clone(Proka.mt) # make a safety insurance
  
  Euka.mt <- microtable$new(sample_table = ENV.raw
                             , otu_table = Euka.otu
                             , tax_table = Euka_taxonomy_table[Euka_taxonomy_table$Kingdom != 'k__Animalia',] # The rest are microeukaryotes
                             , rep_fasta = Euka.fasta
                             , phylo_tree = Euka.tree
  )
  Euka.raw <- clone(Euka.mt) # make a safety insurance
  
  #trim the data, very important!
  Proka.mt$tidy_dataset()
  Euka.mt$tidy_dataset()

  # define data for low-water period (LWP; January, 2024) and high-water period (HWP; August, 2024)
    #Jan (LWP)
    Proka.mt_Jan <- clone(Proka.mt)
    Proka.mt_Jan$sample_table <- subset(Proka.mt_Jan$sample_table,Month =='Jan')
    Proka.mt_Jan$tidy_dataset() 
    
    Euka.mt_Jan <- clone(Euka.mt)
    Euka.mt_Jan$sample_table <- subset(Euka.mt_Jan$sample_table,Month =='Jan')
    Euka.mt_Jan$tidy_dataset() 
  
    #Aug (HWP)
    Proka.mt_Aug <- clone(Proka.mt)
    Proka.mt_Aug$sample_table <- subset(Proka.mt_Aug$sample_table,Month =='Aug')
    Proka.mt_Aug$tidy_dataset() 
    
    Euka.mt_Aug <- clone(Euka.mt)
    Euka.mt_Aug$sample_table <- subset(Euka.mt_Aug$sample_table,Month =='Aug')
    Euka.mt_Aug$tidy_dataset() 
  
# Fig. S2: PCA of Environmental Variables ####
  plot_pca <- function(mt) {
    env_data <- mt$sample_table[, sapply(mt$sample_table, is.numeric)] # Select numeric columns
    group_info <- mt$sample_table$Month
    
    # PCA (using rda on scaled data)
    library(vegan)
    pca_result <- rda(scale(env_data))
    
    scores_site <- scores(pca_result, display = "sites")
    scores_env <- scores(pca_result, display = "species")
    
    var_expl <- round(100 * summary(pca_result)$cont$importance[2, 1:2], 1)
    
    # Plot setup
    par(mar = c(5, 5, 4, 6) + 0.1, xpd = TRUE)
    plot(scores_site, type = "n",
         xlab = paste0("PC1 (", var_expl[1], "%)"), 
         ylab = paste0("PC2 (", var_expl[2], "%)"),
         xlim = c(-2, 2), ylim = c(-3, 3), main = "Environmental PCA")
    
    # Colors
    col_vec <- ifelse(group_info == 'Jan', '#E69F00', '#AA0D91')
    
    # Points and Ellipses
    points(scores_site, pch = 17, col = col_vec, cex = 1.8, lwd = 1.5)
    ordiellipse(pca_result, groups = group_info, conf = 0.95, col = 'gray', lwd = 2, lty = 2, label = TRUE)
    
    # Env Arrows
    arrows(0, 0, scores_env[,1], scores_env[,2], length = 0.1, lwd = 1.8, col = "black")
    text(scores_env[,1], scores_env[,2], labels = rownames(scores_env), col = "black", cex = 1.1)
    
    legend("topright", legend = c('Jan', 'Aug'), pch = 17, col = c('#E69F00', '#AA0D91'), pt.cex = 1.5, bty = "n")
  }
  plot_pca(Proka.mt) # Using Proka sample_table as representative for Env
  
  
# Fig. 2: Difference in Composition ####
  Diff.compo <- function(mt, title_group){
    # A: Bar Plot
    tax_level <- ifelse(title_group == 'Proka', 'Phylum', 'Kingdom')
    t_abund <- trans_abund$new(dataset = mt, taxrank = tax_level)
    print(t_abund$plot_bar(facet = 'Month', xtext_keep = FALSE))
    
    # B: NMDS
    mt$cal_betadiv(unifrac = FALSE)
    beta_div <- mt$beta_diversity$bray
    nmds_res <- metaMDS(beta_div, trace = 0)
    
    # Plot NMDS
    par(mfrow=c(1,1))
    cols <- ifelse(mt$sample_table$Month == 'Jan', "#E69F00", "#AA0D91")
    plot(nmds_res, type = "n", main = paste("NMDS -", title_group))
    points(nmds_res, col = cols, pch = 1)
    ordiellipse(nmds_res, groups = mt$sample_table$Month, conf = 0.95, col = 'gray', label = TRUE)
    text(min(nmds_res$points[,1]), min(nmds_res$points[,2]), paste('Stress =', round(nmds_res$stress, 2)))
    
    # ANOSIM
    print(anosim(beta_div, mt$sample_table$Month))
    
    # C: Sankey Plot (Genus level)
    t2 <- trans_abund$new(dataset = mt, taxrank = "Genus", ntaxa = 10, groupmean = "Month")
    p_sankey <- t2$plot_bar(use_alluvium = TRUE, clustering = FALSE, legend_text_italic = FALSE,
                            color_values = RColorBrewer::brewer.pal(12, "Paired")) +
      theme_classic() + ylab("Relative Abundance") + ggtitle(paste("Dynamics:", title_group))
    print(p_sankey)
  }
  
  Diff.compo(Proka.mt, 'Proka')
  Diff.compo(Euka.mt, 'Euka')
  
# Fig. 3: Alpha Diversity Patterns ####
  Alpha.plot <- function(mt){
    mt$cal_alphadiv()
    df <- mt$alpha_diversity
    df$Month <- factor(mt$sample_table$Month, levels = c("Jan", "Aug"))
    
    metrics <- c("Observed", "Chao1", "Shannon", "Pielou")
    cols <- c("#E69F00", "#AA0D91") # Jan, Aug
    
    par(mfrow = c(2, 4))
    # Boxplots
    for(m in metrics){
      boxplot(df[[m]] ~ df$Month, col = adjustcolor(cols, 0.3), border = cols, 
              main = m, names = c("LWP", "HWP"), outline = FALSE)
      pval <- wilcox.test(df[[m]] ~ df$Month)$p.value
      text(1.5, max(df[[m]]), paste("p =", format.pval(pval, digits = 3)), cex = 0.8)
    }
    # Stripcharts (for PPT)
    for(m in metrics){
      plot(jitter(as.numeric(df$Month), 0.25), df[[m]], col = cols[as.numeric(df$Month)],
           pch = 1, xaxt = 'n', xlab = "Season", ylab = m, main = m)
    }
    par(mfrow = c(1, 1))
  }
  Alpha.plot(Proka.mt)
  Alpha.plot(Euka.mt)
  
# Fig. S3: Eukaryotic Richness by Kingdom ####
  plot_kingdom_richness <- function(mt) {
    # Convert to binary for richness
    mt_bin <- clone(mt)
    mt_bin$otu_table[mt_bin$otu_table > 0] <- 1
    mt_bin$cal_abund(rel = FALSE)
    
    richness_tab <- mt_bin$taxa_abund$Kingdom
    months <- mt_bin$sample_table$Month
    
    par(mfrow=c(3,3))
    for (i in 1:nrow(richness_tab)){
      vals <- as.numeric(richness_tab[i,])
      res <- wilcox.test(vals ~ months)
      boxplot(vals ~ months, border=c("#AA0D91","#E69F00"),
              col='white',
              ylab = rownames(richness_tab)[i], 
              main = paste("p =", round(res$p.value, 3)))
    }
    par(mfrow=c(1,1))
  }
  plot_kingdom_richness(Euka.mt)
  
  
#Fig.4: Community assembly analysis####
  #Fig.4a-d: Neutral model ####
  # Code from Chen et al., (2019) Microbiome. 10.1186/s40168-019-0749-8
  # https://github.com/Weidong-Chen-Microbial-Ecology/Stochastic-assembly-of-river-microeukaryotes/blob/master/Neutral%20community%20model.r
  Neutral.plot <- function(mt){
    # Fits the neutral model from Sloan et al. 2006 to an OTU table and returns several fitting statistics.
    # Alternatively, will return predicted occurrence frequencies for each OTU based on their abundance in the metacommunity
    
    # Install the following packages if they haven't been available in your computer yet 
    library(Hmisc)
    library(minpack.lm)
    library(stats4)
    
    # Using Non-linear least squares (NLS) to calculate R2:
    # spp: A community table with taxa as rows and samples as columns
    spp <- mt$otu_table
    spp <- t(spp)
    
    N <- mean(apply(spp, 1, sum))
    p.m <- apply(spp, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
    
    spp.bi <- 1*(spp > 0)
    freq <- apply(spp.bi, 2, mean)
    freq <- freq[freq != 0]
    
    C <- merge(p, freq, by=0)
    C <- C[order(C[,2]),]
    C <- as.data.frame(C)
    C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
    
    p <- C.0[,2]
    freq <- C.0[,3]
    names(p) <- C.0[,1]
    names(freq) <- C.0[,1]
    
    d = 1/N
    
    # Fit model parameter m (or Nm) using Non-linear least squares (NLS)
    m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 - p), lower.tail=FALSE), start=list(m=0.1))
    m.fit  # get the m value
    m.ci <- confint(m.fit, 'm', level=0.95)
    
    freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 - p), lower.tail=FALSE)
    pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
    
    Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
    Rsqr  # get the R2 value
    
    # Optional: write 3 files: p.csv, freq.csv and freq.pred.csv
    # write.csv(p, file = "j:/iue/aehg11/neutr_model/p.csv")
    # write.csv(freq, file = "j:/iue/aehg11/neutr_model/freq.csv")
    # write.csv(freq.pred, file = "j:/iue/aehg11/neutr_model/freq.pred.csv")
    
    # Drawing the figure using grid package:
    # p is the mean relative abundance
    # freq is occurrence frequency
    # freq.pred is predicted occurrence frequency
    bacnlsALL <- data.frame(p, freq, freq.pred, pred.ci[,2:3])
    
    inter.col <- rep('black', nrow(bacnlsALL))
    inter.col[bacnlsALL$freq <= bacnlsALL$Lower] <- '#F2CACA'  # define the color of below points
    inter.col[bacnlsALL$freq >= bacnlsALL$Upper] <- '#BFE5E5'  # define the color of up points
    
    library(grid)
    grid.newpage()
    pushViewport(viewport(h=0.6, w=0.6))
    pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0, 1.02), extension=c(0.02, 0)))
    
    grid.rect()
    grid.points(log10(bacnlsALL$p), bacnlsALL$freq, pch=17, gp=gpar(col=inter.col, cex=0.7))
    grid.yaxis()
    grid.xaxis()
    
    grid.lines(log10(bacnlsALL$p), bacnlsALL$freq.pred, gp=gpar(col='blue', lwd=2), default='native')
    grid.lines(log10(bacnlsALL$p), bacnlsALL$Lower, gp=gpar(col='blue', lwd=2, lty=2), default='native') 
    grid.lines(log10(bacnlsALL$p), bacnlsALL$Upper, gp=gpar(col='blue', lwd=2, lty=2), default='native')  
    
    grid.text(y=unit(0, 'npc') - unit(2.5, 'lines'), label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
    grid.text(x=unit(0, 'npc') - unit(3, 'lines'), label='Frequency of Occurance', gp=gpar(fontface=2), rot=90) 
    
    # Function to draw text
    draw.text <- function(just, i, j) {
      grid.text(paste("Rsqr=", round(Rsqr, 3), "\n", "Nm=", round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
    }
    
    x <- unit(1:4/5, "npc")
    y <- unit(1:4/5, "npc")
    draw.text(c("centre", "bottom"), 4, 1)
    
  }
  
  Neutral.plot(Proka.mt_Jan)
  Neutral.plot(Proka.mt_Aug)
  
  Neutral.plot(Euka.mt_Jan)
  Neutral.plot(Euka.mt_Aug)
  #Fig.4e-f：Community assembly####
  # code from Liu et al., (2021)
  # https://chiliubio.github.io/microeco_tutorial/index.html
  
  # test polygenietic signal
  # generate trans_nullmodel object
  # as an example, we only use high abundance OTU with mean relative abundance > 0.0005
  t1 <- trans_nullmodel$new(mt, filter_thres = 0.001,add_data = ENV.raw)
  
  # use pH as the test variable
  t1$cal_mantel_corr(use_env = "pH") # Check if phylogenetic groups respond differently to pH
  
  # plot the mantel correlogram
  t1$plot_mantel_corr()
  
    #betaNRI -  ‘basal’ phylogenetic turnover####
    t1$cal_ses_betampd(runs = 999, abundance.weighted = TRUE) # null model run 99 times for the example
    # return t1$res_ses_betampd
    
    # add betaNRI matrix to beta_diversity list
    mt$beta_diversity[["betaNRI"]] <- t1$res_ses_betampd
    # create trans_beta class, use measure "betaNRI"
    t2 <- trans_beta$new(dataset = mt, group = "Month", measure = "betaNRI")
    
    # transform the distance for each group
    t2$cal_group_distance()
    # see the help document for more methods, e.g. "anova" and "KW_dunn"
    t2$cal_group_distance_diff(method = "wilcox")
    
    # plot the results
    g1 <- t2$plot_group_distance(add = "mean")
    g1 + geom_hline(yintercept = -2, linetype = 2) + geom_hline(yintercept = 2, linetype = 2)
    t2$res_group_distance_diff
    
    #betaNTI - phylogenetic terminal turnover ####
    # null model run 500 times
    t1$cal_ses_betamntd(runs = 999, abundance.weighted = TRUE, null.model = "taxa.labels")
    # return t1$res_ses_betamntd
    
    # Bray-Curtis-based Raup-Crick
    # -  assess whether the compositional turnover was governed primarily by drift
    t1$cal_rcbray(runs = 999)
    
    # use betaNTI and rcbray to evaluate processes
    t1$cal_process(use_betamntd = TRUE, group = "Month") # for each group
    
  t1$res_process
    
  write.table(t1$res_process,'clipboard',sep='\t')
    
    
  #Fig.4g-j: Meta-community, distance decay relationship ####
  #Standardized estimates for Eukaryotes and Bacteria overall ####
  fit_model <- function(Euka.mt, Proka.mt, month_name, Y) {
    kingdom <- clone(Euka.mt)
    
    # Calculate distance matrices and convert to vector 
    Euka.dis <- as.vector(vegdist(t(log(kingdom$otu_table+1)), method = 'bray'))
    Proka.dis <- as.vector(vegdist(t(log(Proka.mt$otu_table+1)), method = 'bray'))
    
    # Environmental distance (assume columns 6+ are environmental factors)
    env_data <- Proka.mt$sample_table[,-c(1:5)]
    Env.dis <- as.vector(vegdist(scale(env_data, scale = TRUE, center = TRUE), method = 'euclidean'))
    
    # Geographic distance (assume columns 2,3 are coordinates)
    coordinates.num <- apply(Proka.mt$sample_table[,c(2,3)], 2, as.numeric)
    Geo.dis <- as.vector(vegdist(scale(coordinates.num, scale = TRUE, center = TRUE), method = 'euclidean'))
    
    lm_data <- data.frame(
      Proka.dis = Proka.dis,
      Euka.dis = Euka.dis,
      Env.dis = Env.dis,
      Geo.dis = Geo.dis
    )
    
    # standardize the data
    your_data_scaled <- lm_data
    your_data_scaled$Env.dis <- scale(lm_data$Env.dis)
    your_data_scaled$Geo.dis <- scale(-1 * lm_data$Geo.dis)
    your_data_scaled$Proka.dis <- scale(lm_data$Proka.dis)
    your_data_scaled$Euka.dis <- scale(lm_data$Euka.dis)
    
    # modelling
    # Note: BEINF distribution requires dependent variable in [0,1], no scaling
    if(Y == 'Proka'){
      # Dependent: Proka (original 0-1), Independent: Euka (scaled)
      your_data_scaled$Proka.dis <- lm_data$Proka.dis 
      mod.beta <- gamlss(Proka.dis ~ Euka.dis + Env.dis + Geo.dis, family = BEINF, data = your_data_scaled, trace = FALSE)
    }
    
    if(Y == 'Euka'){
      # Dependent: Euka (original 0-1), Independent: Proka (scaled)
      your_data_scaled$Euka.dis <- lm_data$Euka.dis
      mod.beta <- gamlss(Euka.dis ~ Proka.dis + Env.dis + Geo.dis, family = BEINF, data = your_data_scaled, trace = FALSE)
    }
    
    # results
    summ <- summary(mod.beta, what = "mu")
    vars_to_keep <- c("Proka.dis", "Euka.dis", "Env.dis", "Geo.dis")
    stand.coefs <- summ[rownames(summ) %in% vars_to_keep, , drop = FALSE]
    
    return(list(stand.coefs = stand.coefs, month = month_name, target = Y))
  }
  
  results_Jan_Euka <- fit_model(Euka.mt_Jan, Proka.mt_Jan, "Jan", Y='Euka')
  results_Jan_Proka <- fit_model(Euka.mt_Jan, Proka.mt_Jan, "Jan", Y='Proka')
  
  results_Aug_Euka <- fit_model(Euka.mt_Aug, Proka.mt_Aug, "Aug", Y='Euka')
  results_Aug_Proka <- fit_model(Euka.mt_Aug, Proka.mt_Aug, "Aug", Y='Proka')
  
    #Convert the data for plot ####
    prepare_plot_data <- function(results) {
      coefs <- results$stand.coefs
      month <- results$month
      target <- results$target # Dependent variable
      
      estimates <- coefs[, "Estimate"]
      se <- coefs[, "Std. Error"]
      p_vals <- coefs[, "Pr(>|t|)"]
      
      lower <- estimates - 1.96 * se
      upper <- estimates + 1.96 * se
      
      sig <- ifelse(p_vals < 0.001, "***", ifelse(p_vals < 0.01, "**", ifelse(p_vals < 0.05, "*", "")))
      
      data.frame(
        Month = month,
        Target = target,
        Variable = rownames(coefs),
        Estimate = estimates,
        Lower = lower,
        Upper = upper,
        P_value = p_vals,
        Significance = sig
      )
    }
    
    plot_data <- rbind(
      prepare_plot_data(results_Jan_Euka),
      prepare_plot_data(results_Jan_Proka),
      prepare_plot_data(results_Aug_Euka),
      prepare_plot_data(results_Aug_Proka)
    )
    
    #Plot ####
    par(mfrow = c(1, 2), mar = c(4, 1, 3, 1))
    
    # Define colors
    colors <- c(
      "Euka.dis" = "#E69F00", 
      "Proka.dis" = "#E69F00",
      "Env.dis" = "#56B4E9", 
      "Geo.dis" = "#009E73"
    )
    
    # Loop
    for (month in c("Jan", "Aug")) {
      
      # Create empty canvas (adjust range to accommodate both models)
      plot(NA, xlim = c(-0.8, 0.8), ylim = c(0.5, 7.5), 
           xlab = "Standardized Estimates", ylab = "", 
           main = month, axes = FALSE)
      axis(1)
      abline(v = 0, lty = 2, col = "gray50")
      
      # Plot two sections:
      # Upper: Target = Euka (explaining Euka variation)
      # Lower: Target = Proka (explaining Proka variation)
      targets <- c("Euka", "Proka")
      base_y_positions <- c(4.5, 0.5) # The initial heights for two models
      
      for (t_idx in 1:2) {
        target_name <- targets[t_idx]
        base_y <- base_y_positions[t_idx]
        
        # Sub title
        text(-0.8, base_y + 3, paste0("Response: ", target_name), pos = 4, font = 2, cex = 0.9)
        
        # sub data to plot
        sub_data <- plot_data[plot_data$Month == month & plot_data$Target == target_name, ]
        
        # Define variable order
        if(target_name == "Euka") {
          vars_to_plot <- c("Geo.dis", "Env.dis", "Proka.dis")
        } else {
          vars_to_plot <- c("Geo.dis", "Env.dis", "Euka.dis")
        }
        
        # Loop to plot
        for (i in seq_along(vars_to_plot)) {
          var_name <- vars_to_plot[i]
          var_data <- sub_data[sub_data$Variable == var_name, ]
          
          if(nrow(var_data) > 0) {
            y_pos <- base_y + (i - 1)
            
            # error bars
            arrows(var_data$Lower, y_pos, var_data$Upper, y_pos, 
                   length = 0.05, angle = 90, code = 3, 
                   col = colors[var_name], lwd = 2)
            
            # points
            points(var_data$Estimate, y_pos, pch = 15, 
                   col = colors[var_name], cex = 1.5)
            
            # siginificant marks
            if (var_data$Significance != "") {
              # Auto-adjust marker position (outside error bars)
              text_x <- ifelse(var_data$Estimate > 0, var_data$Upper + 0.1, var_data$Lower - 0.1)
              text(text_x, y_pos, var_data$Significance, col = "black", cex = 1)
            }
            
          }
        }
      }
    }
    
    # Add legend
    legend("bottomright", legend = c("Biotic", "Env", "Geo"), 
           col = c("#E69F00", "#56B4E9", "#009E73"), pch = 15, bty = "n")
    
    
    
  #Standardized estimates for Eukaryotic kingdoms ####
    Contribution_to_Euka <- function(Euka.mt = Euka.mt
                                     ,Proka.mt= Proka.mt
                                     ,kingdom_name = kingdom_name){
      #Calculating distance matrices ####
      library(vegan)
    
      kingdom <- clone(Euka.mt)
      kingdom$tax_table <- kingdom$tax_table[kingdom$tax_table$Kingdom==kingdom_name,]
      kingdom$tidy_dataset()
    
      Euka.dis <- vegdist(t(log(kingdom$otu_table+1)),method = 'bray')
      Proka.dis <- vegdist(t(log(Proka.mt$otu_table+1)),method = 'bray')

      # the fist 5 col is not environmental variables
      Env.dis <- vegdist(scale(Proka.mt$sample_table[,-c(1:5)],scale = TRUE,center = TRUE)
                         ,method = 'euclidean')
      
      
      # the 2nd and 3rd col are coordinates
      coordinates <- Proka.mt$sample_table[,c(2,3)] 
      # transform coordinates to numeric
      coordinates.num <- apply(coordinates, 2, as.numeric)
      Geo.dis <- vegdist(scale(coordinates.num,scale = TRUE,center = TRUE)
                         ,method = 'euclidean')
      
      #scale the independent variables
      lm_data <- data.frame(Euka.dis = Euka.dis
                            ,Proka.dis = Proka.dis
                            ,Env.dis = Env.dis
                            ,Geo.dis = Geo.dis
      )
      
      your_data_scaled <- lm_data
      your_data_scaled$Euka.dis <- as.vector(lm_data$Euka.dis)
      your_data_scaled$Proka.dis <- scale(as.vector(lm_data$Proka.dis),scale = TRUE,center = TRUE)
      your_data_scaled$Env.dis <- scale(as.vector(lm_data$Env.dis),scale = TRUE,center = TRUE)
      your_data_scaled$Geo.dis <- scale(as.vector(lm_data$Geo.dis),scale = TRUE,center = TRUE)
      
      #Gamlss modeling ####
      
      #check VIF
      library(car)
      mod <- lm(Euka.dis ~ Proka.dis + Env.dis + Geo.dis
                ,data=your_data_scaled)
      (vif(mod))
      
      # when data is [0,1], BEINF distribution is used
      library(gamlss)
      mod.beta <- gamlss(Euka.dis ~  Proka.dis + Env.dis + Geo.dis
                         ,family = BEINF
                         ,data=your_data_scaled)
      
      # validate the model
      wp(mod.beta)
      
      stand.coefs <- summary(mod.beta, what = "mu")[c(2,3,4),]
      return(stand.coefs = stand.coefs)
    }
    
    #Standardized estimates of each distance for eukaryotic kingdoms ####
    kingdom_results_Jan <- list()
    kingdom_results_Aug <- list()
    kingdom.list <- c("k__Amoebozoa","k__Chloroplastida","k__Cryptophyceae"
                      ,"k__Discoba","k__Fungi")
    for (taxa in kingdom.list) {
      kingdom_results_Jan[[taxa]] <- Contribution_to_Euka(
        Euka.mt = Euka.mt_Jan,
        Proka.mt = Proka.mt_Jan,
        kingdom_name = taxa
      )
      
      kingdom_results_Aug[[taxa]] <- Contribution_to_Euka(
        Euka.mt = Euka.mt_Aug,
        Proka.mt = Proka.mt_Aug,
        kingdom_name = taxa
      )
      
      list_metacom_results <- list()
      list_metacom_results[['Jan']]<-kingdom_results_Jan
      list_metacom_results[['Aug']]<-kingdom_results_Aug
    }
    list_metacom_results
    
    #Turn the result list into matrix ####
    library(reshape2)
    library(dplyr)
    
    extract_data <- function(results_list) {
      coefs_list <- lapply(results_list, function(x) x)
      
      # Convert coefficients to long format
      df_coefs <- melt(coefs_list) %>%
        rename(
          Kingdom = L1,
          Variable = Var1,
          Statistic = Var2,
          Value = value
        )
      
      return(df_coefs)
    }
    
    # Extract data for each month
    jan_data <- extract_data(kingdom_results_Jan)
    jan_data$Month <- "Jan"
    
    aug_data <- extract_data(kingdom_results_Aug)
    aug_data$Month <- "Aug"
    
    # Combine all months' data
    df_long <- rbind(jan_data, aug_data)
    
    # Convert to wide format matrix (month+Kingdom as rows, variable+statistic as columns)
    matrix_wide <- dcast(df_long, 
                         Month + Kingdom ~ Variable + Statistic, 
                         value.var = "Value")
    head(matrix_wide)
    
    #Plot the estimate +- Std and p value ####
    library(tidyr)
    tidy_data <- pivot_longer(
      data = matrix_wide,
      cols = -c(Month, Kingdom),
      names_to = c("Variable", "Statistic"),
      names_sep = "_",
      values_to = "Value"
    )
    tidy_data
    
    
    # 95% Confident Intervals
    Estimate = (tidy_data$Value)[tidy_data$Statistic=='Estimate']
    Std.Error = (tidy_data$Value)[tidy_data$Statistic=='Std. Error']
    p_value = (tidy_data$Value)[tidy_data$Statistic=='Pr(>|t|)']
    
    lower <- (Estimate - 1.96 * Std.Error)
    upper <- (Estimate + 1.96 * Std.Error)
    
    # mark significance
    significance <- ifelse(p_value < 0.001, "***",
                           ifelse(p_value < 0.01, "**",
                                  ifelse(p_value < 0.05, "*", "")))

    #kingdom
    name.f <-NULL;kingdom <- c(  for(kingdom in unique(tidy_data$Kingdom)){
      name <- rep(kingdom,3)
      name.f <- c(name.f,name)
    });name.f
    kingdom <- rep( name.f, 2)
    
    Month <- c(rep('Aug',3*length(unique(tidy_data$Kingdom)))
               ,rep('Jan',3*length(unique(tidy_data$Kingdom)))
    )
    
    dis.name <- rep(c('Proka.dis','Env.dis','Geo.dis'),10)
    
    matrix <- data.frame(
      Month = Month
      ,kingdom = kingdom
      ,Variable = dis.name
      ,Estimate = Estimate
      ,lower = lower
      ,upper = upper
      ,sig = significance
      ,p_value = p_value
    )
    matrix
    
    
    par(mfrow = c(1, 2))
    
    data=matrix
    # define colors
    colors <- c("#E69F00", "#56B4E9", "#009E73")
    names(colors) <- c("Proka.dis", "Env.dis", "Geo.dis")
    
    # Loop to plot
    for (month in c("Aug", "Jan")
         ) {
      
      month_data <- data[data$Month == month, ]
      
      # Empty canvas
      plot(NA, xlim = c(-0.8, 0.8), ylim = c(0, 16), 
           xlab = "Estimates", ylab = "", 
           main = month, axes = FALSE)
      axis(1)
      
      # Add y-axis labels and dividers
      kingdoms <- unique(month_data$kingdom)
      for (i in seq_along(kingdoms)) {
        kingdom <- kingdoms[i]
        y_pos <- 16 - i * 2.5
        text(-0.05, y_pos + 1, kingdom, pos = 2, xpd = TRUE, font = 1)
        
        # Plot three variables for current Kingdom
        vars <- rev(c("Proka.dis", "Env.dis", "Geo.dis"))
        for (j in seq_along(vars)) {
          var <- vars[j]
          y_pos_var <- y_pos + (j - 1)
          
          # points
          point_data <- month_data[month_data$kingdom == kingdom & month_data$Variable == var, ]
          
          # error bars and points
          if (nrow(point_data) > 0) {
            arrows(point_data$lower, y_pos_var, point_data$upper, y_pos_var, 
                   length = 0.06, angle = 90, code = 3, col = colors[var],lwd = 2)
            points(point_data$Estimate, y_pos_var, pch = 15, col = colors[var],cex=1.5)
            
            # mark significance
            if (!is.na(point_data$sig) && point_data$sig != "") {
              text_x <- ifelse(point_data$Estimate > 0, point_data$upper + 0.1, point_data$lower - 0.1)
              text(text_x, y_pos_var, point_data$sig, col = "black")
            }
          }
        }
      }
      
      # Add reference line
      abline(v = 0, lty = 2, col = "gray50")
      legend('bottomright'
             ,c("Bac-dis", "Env-dis", "Geo-dis")
             ,col=c("#E69F00", "#56B4E9", "#009E73")
             ,pch=15)
    }
    
     
#Fig.5: Biotic interaction####
  #Fig.5a-b: Co-occurence Network ####
    #3.1 data preparation ####
    rownames(Proka.mt$sample_table); nrow(Proka.mt$sample_table)
    rownames(Euka.mt$sample_table) 
    
    #Rename the OTU, then combine the Euka and Proka
    OTU.Proka.new_names <- paste0('Proka_',rownames(Proka.mt$otu_table));length(rownames(Proka.mt$otu_table))
      rownames(Proka.mt$otu_table) <- OTU.Proka.new_names
      rownames(Proka.mt$tax_table) <- OTU.Proka.new_names
      rownames(Proka.mt$tax_table)
      
    OTU.Euka.new_names <- paste0('Euka_',rownames(Euka.mt$otu_table));length(rownames(Euka.mt$otu_table))
      rownames(Euka.mt$otu_table) <- OTU.Euka.new_names
      rownames(Euka.mt$tax_table) <- OTU.Euka.new_names
      rownames(Euka.mt$tax_table)
  
    # Combine the Proka and Euka
    Pro.Eu_otu_table <- rbind(Proka.mt$otu_table,Euka.mt$otu_table)
    
      #log transform the abundance to reduce the effect of extreme values
      Pro.Eu_otu_table <- log(Pro.Eu_otu_table+1)
    Pro.Eu_tax_table <- rbind(Proka.mt$tax_table,Euka.mt$tax_table)
    
    colnames(Pro.Eu_tax_table)
    
    Combine.Pro.Eu.mt <- microtable$new(sample_table = ENV.raw
                              , otu_table = Pro.Eu_otu_table
                              , tax_table = Pro.Eu_tax_table
    )
    
    head(Combine.Pro.Eu.mt$tax_table)
  
    nrow(Combine.Pro.Eu.mt$tax_table[Combine.Pro.Eu.mt$tax_table$Domain == 'd__Eukaryota',])
    nrow(Combine.Pro.Eu.mt$tax_table[Combine.Pro.Eu.mt$tax_table$Domain == 'd__Bacteria',])
    
    Combine.Pro.Eu.mt$tidy_dataset()
    nrow(Combine.Pro.Eu.mt$tax_table)
    Pro.Euka.raw <- clone(Combine.Pro.Eu.mt)
    Combine.Pro.Eu.mt$save_table('Combined Proka and Euka') # save the dataset
    
    # Season
    Jan_mt <- clone(Combine.Pro.Eu.mt)
    Jan_mt$sample_table <- subset(Jan_mt$sample_table,Month =='Jan')
    Jan_mt$tidy_dataset() # remove other sites
    
    Aug_mt <- clone(Combine.Pro.Eu.mt)
    Aug_mt$sample_table <- subset(Aug_mt$sample_table,Month =='Aug')
    Aug_mt$tidy_dataset() # remove other sites
  
    #3.2 CHOOSE ONE ####
    # Here we need to choose one period to run the rest codes
    mt.network <- clone(Jan_mt);title.month <- 'Jan'
    mt.network <- clone(Aug_mt);title.month <- 'Aug'
    
    # First filter out taxa appearing in less than 20% of sites
    mt.network$filter_taxa(freq = 0.2   # occur more than 20% sites
                           ,rel_abund = 0.00001  # relative abundance more than 0.0001%
    )
    mt.network$tidy_dataset()
    print(c('There are',nrow(mt.network$otu_table),'Taxa remain in the data set'))
    
    #Construct network; require igraph package
      Co_occur <- trans_network$new(dataset = mt.network
                                    , cor_method = "spearman"    #or pearson
                                    ,filter_thres = 0.0001 # more than 0.01% Abundance
      )
      
      # remove some connections
      Co_occur$cal_network(
        COR_p_thres = 0.01       # p value = 0.01
        ,COR_cut = 0.8            # correlation cut = 0.8 
      )          
      Co_occur$res_network
      
      #positive and negative links
      Co_occur$cal_sum_links()
      Co_occur$res_sum_links_pos;sum(Co_occur$res_sum_links_pos)
      Co_occur$res_sum_links_neg;sum(Co_occur$res_sum_links_neg)
      
      # invoke igraph cluster_fast_greedy function for this undirected network
      Co_occur$cal_module(method = "cluster_fast_greedy")
    
    #Network attributes
      Co_occur$cal_network_attr()
      Co_occur$res_network_attr
      write.table(Co_occur$res_network_attr,'clipboard',sep='\t')
      
      Co_occur$get_adjacency_matrix()
      Co_occur$res_adjacency_matrix
      
      Co_occur$get_edge_table()
      Co_occur$res_edge_table[1:5,]
      Co_occur$get_node_table()
      Co_occur$res_node_table[1:5,] # Attributes of each node, including zi and pi
      
      Co_occur$plot_taxa_roles(use_type = 1,add_label = TRUE)
      
      taxa_roles <- Co_occur$res_node_table[
        Co_occur$res_node_table$taxa_roles != 'Peripheral nodes',
        c('z','p','taxa_roles','module','Abundance','Kingdom','Phylum','Class','Genus')
      ]
      taxa_roles <- taxa_roles[is.na(taxa_roles$z)==FALSE,]
      taxa_roles
      sort(rownames(taxa_roles))
  
  #Fig.5c-d: Bipartite network which connect to which ####
  Who_connect_whom <- function(network = Co_occur){
    library(dplyr)
    library(reshape2)
    
    # Extract edge data
    edges_from_to <- network$res_edge_table[, 1:2] %>% 
      setNames(c("from", "to"))
    
    # Create group mapping vector
    group_map <- with(Combine.Pro.Eu.mt$tax_table, {
      ifelse(Domain == 'd__Bacteria', as.character(Phylum),
             ifelse(Domain == 'd__Eukaryota', as.character(Kingdom), NA))
    })
    names(group_map) <- rownames(Combine.Pro.Eu.mt$tax_table)
    
    # Replace nodes with group labels
    edges_grouped <- edges_from_to %>%
      mutate(
        from_group = group_map[as.character(from)],
        to_group = group_map[as.character(to)]
      ) %>%
      na.omit()  # Remove invalid groups
    
    # Generate group connection matrix
    matrix_raw <- table(edges_grouped$from_group, edges_grouped$to_group)
    
    # Normalize
    matrix_normalized <- matrix_raw / sum(matrix_raw)
    rownames(matrix_normalized);rownames(matrix_normalized)
    
    # Extract different matrices
    rnames <- rownames(matrix_normalized)
    cnames <- colnames(matrix_normalized)
    kingdom_rows <- grepl("^k__", rnames)  # Kingdom 
    phylum_rows <- grepl("^p__", rnames)   # Phylum 
    kingdom_cols <- grepl("^k__", cnames)  # Kingdom 
    phylum_cols <- grepl("^p__", cnames)   # Phylum 
    
    # 1. Extract Kingdom-Kingdom submatrix
    kk_mat <- matrix_normalized[kingdom_rows, kingdom_cols, drop = FALSE]
    
    # 2.  Extract Kingdom-Phylum submatrix
    # Includes two parts: Kingdom rows vs Phylum columns + Phylum rows vs Kingdom columns
    kp_mat <- matrix_normalized[phylum_rows, kingdom_cols, drop = FALSE]  # Phylum行 vs Kingdom列

    # 3. Extract Phylum-Phylum submatrix
    pp_mat <- matrix_normalized[phylum_rows, phylum_cols, drop = FALSE]
    
    #3.6 bipartite_net between microeukaryote kingdoms and bacterial phyla ####
    library(bipartite)
    plotweb(kp_mat
            ,labsize = 2
            ,col.interaction = 'gray'
            ,col.high='#EDB931'
            ,col.low='steelblue4'
            ,text.rot=30
            )
    
    return(kp_mat)
  }
  
  Who_connect_whom(Co_occur)

  #Fig.5e-f: Importance of the Eukayote taxa by degree and betweenness contrality####
  # closer interactions between euka and proka
  # Calculate normalized degree and normalized betweenness centrality for eukaryotic nodes
  
  # Get network node table and adjacency matrix
  
  # Calculate importance of eukaryotic nodes
  node_table <- Co_occur$res_node_table
  adj_matrix <- Co_occur$res_adjacency_matrix
  
  # Select only eukaryotic nodes
  eukaryotic_nodes <- node_table[grepl("Euka_", rownames(node_table)), ]
  cat("Number of Eukaryotic taxa:", nrow(eukaryotic_nodes), "\n")
  
  # Calculate normalized degree centrality
  node_count <- nrow(node_table)
  eukaryotic_nodes$norm_degree <- eukaryotic_nodes$degree / (node_count - 1)
  
  rownames(eukaryotic_nodes)
  
  # Calculate normalized betweenness centrality
  betweenness <- Co_occur$res_node_table$betweenness_centrality[rownames(Co_occur$res_node_table)%in%rownames(eukaryotic_nodes)]
  eukaryotic_nodes$norm_betweenness <- betweenness/((node_count-1)*(node_count-2))

 
  # colors for each kingdom
  kingdom_colors <- c(
    "k__Amoebozoa" = "#E41A1C",        
    "k__Chloroplastida" = "#377EB8",   
    "k__Cryptophyceae" = "#4DAF4A",    
    "k__Discoba" = "#984EA3",           
    "k__Fungi" = "#FF7F00",           
    "k__norank_d__Eukaryota" = "#FFFF33" 
  )
  
  # Plot scatter plot
  par(mar = c(5, 5, 4, 8) + 0.1) 
  
  x_min <- ifelse(min(eukaryotic_nodes$norm_degree, na.rm = TRUE) * 0.8>0,min(eukaryotic_nodes$norm_degree, na.rm = TRUE) * 0.8,0)
  x_max <- max(eukaryotic_nodes$norm_degree, na.rm = TRUE) * 1.25
  y_min <- ifelse(min(eukaryotic_nodes$norm_betweenness, na.rm = TRUE) * 0.8>0,min(eukaryotic_nodes$norm_betweenness, na.rm = TRUE) * 0.8,0)
  y_max <- max(eukaryotic_nodes$norm_betweenness, na.rm = TRUE) * 1.25
  
  # Emplty canvas
  plot(eukaryotic_nodes$norm_degree, 
       eukaryotic_nodes$norm_betweenness,
       type = "n",
       xlab = "Normalized Degree Centrality",
       ylab = "Normalized Betweenness Centrality",
       main = paste(title.month, "- Eukaryotic Nodes Centrality"),
       xlim = c(x_min, x_max),
       ylim = c(y_min, y_max),
       las = 1)
  
  # Add points
  for (kingdom in names(kingdom_colors)) {
    kingdom_data <- eukaryotic_nodes[eukaryotic_nodes$Kingdom == kingdom, ]
    if (nrow(kingdom_data) > 0) {
      points(kingdom_data$norm_degree, 
             kingdom_data$norm_betweenness,
             pch = 15, 
             cex = 1.5,
             col = kingdom_colors[kingdom]
      )
    }
  }

# Legend
legend("topright", 
       inset = c(-0.25, 0),
       xpd = TRUE,
       legend = names(kingdom_colors),
       pch = 16,
       col = kingdom_colors,
       title = "Kingdom",
       cex = 0.8)
  
  eukaryotic_nodes

#Fig.6: Spearman correlation and PLS-PM ####

  OTU.Proka.new_names <- paste0('Proka_',rownames(Proka.mt$otu_table));length(rownames(Proka.mt$otu_table))
  rownames(Proka.mt$otu_table) <- OTU.Proka.new_names
  rownames(Proka.mt$tax_table) <- OTU.Proka.new_names
  rownames(Proka.mt$tax_table)

  OTU.Euka.new_names <- paste0('Euka_',rownames(Euka.mt$otu_table));length(rownames(Euka.mt$otu_table))
  rownames(Euka.mt$otu_table) <- OTU.Euka.new_names
  rownames(Euka.mt$tax_table) <- OTU.Euka.new_names
  rownames(Euka.mt$tax_table)

 
  # transform the data into numeric
  env.for.test <- apply(env.for.test, 2, as.numeric)
  env.for.test <- as.data.frame(env.for.test)
  
  Month <- Proka.mt$sample_table$Month
  
  PLS.data <- function(title.month){
    
    Proka.mt_month <- clone(Proka.mt)
    Euka.mt_month <- clone(Euka.mt)
    
    #define dataset for different month
    Proka.mt_month$sample_table <- Proka.mt_month$sample_table[Proka.mt_month$sample_table$Month==title.month,]
    Euka.mt_month$sample_table <- Euka.mt_month$sample_table[Euka.mt_month$sample_table$Month==title.month,]
    Proka.mt_month$tidy_dataset()
    Euka.mt_month$tidy_dataset()
    
    Euka.mt_month$cal_alphadiv()
    richness.euka <- Euka.mt_month$alpha_diversity$Observed
    
    #PC1, Richness
    Cal.var <- function(Proka.mt_month,Euka.mt_month){
      
      # Extract PC1
      library(vegan)
      pca_result <- prcomp(t(Proka.mt_month$otu_table), center = TRUE, scale. = TRUE)
      
      pc1_scores <- pca_result$x[, 1]  # PC1
      pc1_scores
      
      binary = Proka.mt_month$otu_table;binary[binary>0]=1
      
      
      Proka_variables <- data.frame(
        PC1 = pc1_scores
        ,richness = colSums(binary)
        ,row.names = NULL
      )
      rownames(Proka_variables)=names(pc1_scores)
      
      #euka
      Euka.mt_month$cal_alphadiv(PD=FALSE)
      richness.euka <- Euka.mt_month$alpha_diversity$Observed
      
      pca_result <- prcomp(t(Euka.mt_month$otu_table), center = TRUE, scale. = TRUE)
      pc1_scores_Euka <- pca_result$x[, 1]  # PC1
      
      Euka_variables <- data.frame(
        PC1_Euka = pc1_scores_Euka 
        ,richness = richness.euka
      )
      rownames(Euka_variables) <- names(pc1_scores_Euka)
      
      sem.data.both <- data.frame(Proka = Proka_variables,Euka=Euka_variables)
      return(sem.data.both)
    }
    Eu_Proka_variables <- Cal.var(Proka.mt_month,Euka.mt_month)
    sem.data.both <- cbind(env.for.test[Euka.mt$sample_table$Month==title.month,],Eu_Proka_variables)
    return(sem.data.both)
  }
  pls.data <- PLS.data(title.month=title.month)
  pls.data;'How many sites:';nrow(pls.data);'How many variables:';ncol(pls.data)

  # Spearman correlation ####
  library(Hmisc)
  library(corrplot)
  cor_matrix_result <- rcorr(as.matrix(pls.data), type = "spearman")
  
  # correlation and p values
  correlation_matrix <- cor_matrix_result$r
  p_value_matrix <- cor_matrix_result$P
  
  # Assume first 8 columns are environmental variables, next 4 are biological
  # Extract submatrix between environmental and biological variables (8x4)
  # Modify this based on your actual data
  env_bio_cor <- correlation_matrix[1:8, 9:12]
  env_bio_p <- p_value_matrix[1:8, 9:12]
  
  # Env-Bio correlation matrix
  corrplot(t(env_bio_cor_rm),
           method = "color",
           type = "full",
           tl.col = "black",
           tl.srt = 0,
           diag = TRUE,
           is.corr = FALSE,
           col = colorRampPalette(c("blue", "white", "red"))(200),
           mar = c(0, 0, 1, 0),
           title = title.month,
           p.mat = t(env_bio_p),  
           sig.level = c(0.001, 0.01, 0.05), 
           insig = "label_sig", 
           pch.cex = 0.8, 
           pch.col = "black",
           col.lim = c(-1, 1)
  )

  
  # Based on correlations, select environmental variables for analysis
  # Path analysis####
  library(plspm)
  #devtools::install_github("gastonstat/plspm")        
  
  path <- rbind(
     env <- c(0,0,0)
    ,Proka <- c(1,0,0)
    ,Euka <- c(1,1,0)
  )
  group <- list(env = 1,Proka = 2:3,Euka = 4:5)

  pls.data <- pls.data[,c('TP','Proka.PC1','Proka.richness','Euka.PC1_Euka','Euka.richness')]

  innerplot(path, colpos = "red", arr.width = 0.1)

  colnames(pls.data)
  
  # PLS-PM
  sat.pls <- plspm( Data = pls.data
                    ,path_matrix = path      # path of the model
                    ,blocks = group  # latent group
                    ,modes = rep('A',3)
                    ,scaled = TRUE      # whether scale the data
                    ,boot.val = TRUE    # perform bootstrap validation
                    ,br=999
  )
  print(summary(sat.pls))
  
  PLS_Jan <- PLS('Jan');#saveRDS(PLS_Jan,'PLS_Jan_TN.rds')
  PLS_Aug <- PLS('Aug');#saveRDS(PLS_Aug,'PLS_Aug_rm_TP.rds')
