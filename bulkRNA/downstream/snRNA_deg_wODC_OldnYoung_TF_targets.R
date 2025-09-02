###==================================
# Date: 01/29/2025 - 02/14/2025
# Robin JM Franklin
# Project: CNS_remyelination
# Ze
###==================================

###============###============###============###============
# NOTE Date 02/12/2025  - Update: 02/21/2025
###============###============###============###============

library(cowplot)
library(ggpubr)
library(ggrepel)
library(patchwork)
library(ggplot2)

theme_set(theme_cowplot())

library(dplyr)

###============
# PATH
outPATH <- '/dfs7/swaruplab/zechuas/Collaborations/AltosLab/RobinJMFranklin/CNS_remyelination/Analysis/'
dataPATH <- '/dfs7/swaruplab/zechuas/Collaborations/AltosLab/RobinJMFranklin/CNS_remyelination/data/processed_data/'
tfPATH <- '/dfs7/swaruplab/zechuas/Collaborations/AltosLab/RobinJMFranklin/CNS_remyelination/data/fromSam/stuff_ze/'

# ODC_Old_dir <- paste0(PATH, "/", "ODC_Old")
# ODC_Young_dir <- paste0(PATH, "/", "ODC_Young")
#
# dir.create(ODC_Old_dir)
# dir.create(ODC_Young_dir)

###============
# Load data
#------------------------------------------------------------------
# Bar plot showing the number of DEGs of each tf target type
#------------------------------------------------------------------

# NOTE ODC_Old

# Define variables
fc_cutoff <- 0
barplot_groups <- c('Mat-ODC', 'Int-ODC', 'MF-ODC', 'NFOL', 'OPC')
conditions <- c('ODC_Old')  # Add more conditions if needed, 'ODC_Young'
TFnames <- c('Bach2', 'Bhlhe41', 'Elf2', 'Foxk2', 'Nr6a1', 'Sox8', 'Stat3')  # All TFs
Days <- c('Day5', 'Day14', 'Day30', 'Naive')  # All days

# Loop through each condition, TF, and day
for (condition in conditions) {
  for (TFname in TFnames) {
    for (Day in Days) {

      # Perform your operations here
      # For example, print the combination
      cat("Condition:", condition, "\tTF:", TFname, "\tDay:", Day, "\n")

      # Read TF target data
      # tfData <- read.csv(paste0(dataPATH, condition, '_targets/', condition, '_targets_', TFname, '.csv'))
      tfData <- read.csv(paste0(tfPATH, condition, '_targets_', TFname, '.csv'))

      # Extract primary and secondary genes
      primary_genes <- tfData %>% subset(depth == "1") %>% .$target
      secondary_genes <- subset(tfData, depth == "2") %>% .$target

      # Read DEG data for the current day
      DEG_All <- read.csv(paste0(tfPATH, Day, '_DEGs.csv'))

      # Subset DEGs for the relevant groups
      DEG_sub <- subset(DEG_All, group %in% barplot_groups)

      # Initialize a data frame for plotting
      plot_df <- data.frame()

      # Loop through each group to calculate up/down-regulated genes
      for (cur_group in barplot_groups) {
        # p_val_adj < 0.05
        cur_up <- DEG_sub %>% subset(group == cur_group & p_val_adj < 0.1 & avg_log2FC >= fc_cutoff) %>% .$gene
        cur_down <- DEG_sub %>% subset(group == cur_group & p_val_adj < 0.1 & avg_log2FC <= -1 * fc_cutoff) %>% .$gene

        cur_primary <- primary_genes
        cur_secondary <- secondary_genes
        cur_secondary <- setdiff(cur_secondary, cur_primary)

        cur_primary_up <- intersect(cur_up, cur_primary)
        cur_primary_down <- intersect(cur_down, cur_primary)
        cur_secondary_up <- intersect(cur_up, cur_secondary)
        cur_secondary_down <- intersect(cur_down, cur_secondary)
        cur_other_up <- setdiff(cur_up, unique(c(cur_primary_up, cur_secondary_up)))
        cur_other_down <- setdiff(cur_down, unique(c(cur_primary_down, cur_secondary_down)))

        cur_df <- data.frame(
          group = cur_group,
          target_type = c('primary', 'primary', 'secondary', 'secondary', 'other', 'other'),
          direction = c('up', 'down', 'up', 'down', 'up', 'down'),
          n = c(
            length(cur_primary_up),
            length(cur_primary_down),
            length(cur_secondary_up),
            length(cur_secondary_down),
            length(cur_other_up),
            length(cur_other_down)
          )
        )

        plot_df <- rbind(plot_df, cur_df)
      }

      # Adjust the levels of target_type and group
      plot_df$target_type <- factor(as.character(plot_df$target_type), levels = c('primary', 'secondary', 'other'))
      plot_df$group <- factor(plot_df$group, levels = barplot_groups)

      # Adjust the sign of `n` based on direction
      plot_df$n <- ifelse(plot_df$direction == 'down', -1 * plot_df$n, plot_df$n)

      # Calculate the maximum absolute value for scaling
      plot_max <- plot_df %>%
        group_by(direction, group) %>%
        summarise(x = sum(n)) %>%
        .$x %>%
        abs %>%
        max

      # Create the bar plot
      p <- plot_df %>%
        ggplot(aes(x = n, y = group, fill = target_type)) +
        geom_bar(position = 'stack', stat = 'identity') +
        geom_vline(xintercept = 0) +
        geom_text(aes(label = abs(n)), position = position_stack(vjust = 0.5), size = 3) +
        theme(
          axis.line.y = element_blank(),
          axis.title.y = element_blank()
        ) +
        xlab(bquote("N"[genes])) +
        scale_fill_manual(
          values = c(
            'other' = 'light grey',      # Light grey for other genes
            'secondary' = '#E287E2',     # #E287E2 for secondary genes
            'primary' = 'darkorchid3'    # darkorchid3 for primary genes
          )
        )

      # Save the plot data
      write.csv(plot_df, paste0(outPATH,  'snRNA_deg_wODC_OldnYoung_TF_targets/', condition, '/', Day, '/', 'degs_scRNA_bar_', TFname, '.csv'))

      # Save the plot as a PDF
      pdf(paste0(outPATH,  'snRNA_deg_wODC_OldnYoung_TF_targets/', condition, '/', Day, '/', 'degs_scRNA_bar_', TFname, '.pdf'), width = 5, height = 2, useDingbats = FALSE)
      print(p)
      dev.off()
    }
  }
}






# NOTE ODC_Young

# Define variables
fc_cutoff <- 0
barplot_groups <- c('Mat-ODC', 'Int-ODC', 'MF-ODC', 'NFOL', 'OPC')
conditions <- c('ODC_Young')  # Add more conditions if needed, 'ODC_Young'
TFnames <- c('Bach2', 'Bhlhe41', 'Elf2', 'Foxk2', 'Nr6a1', 'Sox8', 'Stat3')  # All TFs
Days <- c('Day5', 'Day14', 'Day30', 'Naive')  # All days

# Loop through each condition, TF, and day
for (condition in conditions) {
  for (TFname in TFnames) {
    for (Day in Days) {

      # Perform your operations here
      # For example, print the combination
      cat("Condition:", condition, "\tTF:", TFname, "\tDay:", Day, "\n")

      # Read TF target data
      # tfData <- read.csv(paste0(dataPATH, condition, '_targets/', condition, '_targets_', TFname, '.csv'))
      tfData <- read.csv(paste0(tfPATH, condition, '_targets_', TFname, '.csv'))

      # Extract primary and secondary genes
      primary_genes <- tfData %>% subset(depth == "1") %>% .$target
      secondary_genes <- subset(tfData, depth == "2") %>% .$target

      # Read DEG data for the current day
      DEG_All <- read.csv(paste0(tfPATH, Day, '_DEGs.csv'))

      # Subset DEGs for the relevant groups
      DEG_sub <- subset(DEG_All, group %in% barplot_groups)

      # Initialize a data frame for plotting
      plot_df <- data.frame()

      # Loop through each group to calculate up/down-regulated genes
      for (cur_group in barplot_groups) {
        # p_val_adj < 0.05
        cur_up <- DEG_sub %>% subset(group == cur_group & p_val_adj < 0.1 & avg_log2FC >= fc_cutoff) %>% .$gene
        cur_down <- DEG_sub %>% subset(group == cur_group & p_val_adj < 0.1 & avg_log2FC <= -1 * fc_cutoff) %>% .$gene

        cur_primary <- primary_genes
        cur_secondary <- secondary_genes
        cur_secondary <- setdiff(cur_secondary, cur_primary)

        cur_primary_up <- intersect(cur_up, cur_primary)
        cur_primary_down <- intersect(cur_down, cur_primary)
        cur_secondary_up <- intersect(cur_up, cur_secondary)
        cur_secondary_down <- intersect(cur_down, cur_secondary)
        cur_other_up <- setdiff(cur_up, unique(c(cur_primary_up, cur_secondary_up)))
        cur_other_down <- setdiff(cur_down, unique(c(cur_primary_down, cur_secondary_down)))

        cur_df <- data.frame(
          group = cur_group,
          target_type = c('primary', 'primary', 'secondary', 'secondary', 'other', 'other'),
          direction = c('up', 'down', 'up', 'down', 'up', 'down'),
          n = c(
            length(cur_primary_up),
            length(cur_primary_down),
            length(cur_secondary_up),
            length(cur_secondary_down),
            length(cur_other_up),
            length(cur_other_down)
          )
        )

        plot_df <- rbind(plot_df, cur_df)
      }

      # Adjust the levels of target_type and group
      plot_df$target_type <- factor(as.character(plot_df$target_type), levels = c('primary', 'secondary', 'other'))
      plot_df$group <- factor(plot_df$group, levels = barplot_groups)

      # Adjust the sign of `n` based on direction
      plot_df$n <- ifelse(plot_df$direction == 'down', -1 * plot_df$n, plot_df$n)

      # Calculate the maximum absolute value for scaling
      plot_max <- plot_df %>%
        group_by(direction, group) %>%
        summarise(x = sum(n)) %>%
        .$x %>%
        abs %>%
        max

      # Create the bar plot
      p <- plot_df %>%
        ggplot(aes(x = n, y = group, fill = target_type)) +
        geom_bar(position = 'stack', stat = 'identity') +
        geom_vline(xintercept = 0) +
        geom_text(aes(label = abs(n)), position = position_stack(vjust = 0.5), size = 3) +
        theme(
          axis.line.y = element_blank(),
          axis.title.y = element_blank()
        ) +
        xlab(bquote("N"[genes])) +
        scale_fill_manual(
          values = c(
            'other' = 'light grey',      # Light grey for other genes
            'secondary' = '#A7E9AF',     # #E287E2 for secondary genes
            'primary' = 'seagreen'    # darkorchid3 for primary genes
          )
        )

      # Save the plot data
      write.csv(plot_df, paste0(outPATH,  'snRNA_deg_wODC_OldnYoung_TF_targets/', condition, '/', Day, '/', 'degs_scRNA_bar_', TFname, '.csv'))

      # Save the plot as a PDF
      pdf(paste0(outPATH,  'snRNA_deg_wODC_OldnYoung_TF_targets/', condition, '/', Day, '/', 'degs_scRNA_bar_', TFname, '.pdf'), width = 5, height = 2, useDingbats = FALSE)
      print(p)
      dev.off()
    }
  }
}
