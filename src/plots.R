# create folder_path vector
folder_list <- c("data/CAMP000427", "data/CAMP001339", "data/CAMP001490",
                "data/CAMP001523", "data/CAMP002274", "data/CAMP003468",
                "data/CAMP004884", "data/CAMP007136", "data/Kemp")

# readin sequence_date.csv files and join 
library(tidyverse)
library(hrbrthemes)
library(ggplot2)
library(viridis)

readin_sequence_dates <- function(folders) {
  # Initialize an empty list to store data frames from each folder
  df_list <- list()

  # Iterate through the list of folders
  for (folder in folders) {
    # Create the file path for the "sequence_dates.csv" file in the current folder
    file_path <- file.path(folder, "sequence_dates.csv")

    # Check if the file exists
    if (file.exists(file_path)) {
      # Read the "sequence_dates.csv" file into a data frame
      current_df <- read.csv(file_path, header = TRUE)

      # Extract the folder name as the Sample_Name
      current_df$Sample_Name <- basename(folder)

      # Rename the required columns
      colnames(current_df) <- c("Sequence_Name", "Date", "Px_Name")

      # Append the current data frame to the list
      df_list[[length(df_list) + 1]] <- current_df
    }
  }

  # Combine all data frames in the list into a single data frame
  final_df <- do.call(rbind, df_list)

  return(final_df)
}
date_dat <- readin_sequence_dates(folder_list)

# add patient number column
px_ns <- read.csv("data/Px_names_to_n.csv", header = TRUE)
date_dat <- full_join(date_dat, px_ns, by = "Px_Name")
date_dat$Px_n <- as.factor(date_dat$Px_n)

# dot plot of sequence dates
time_plot <- 
    date_dat %>%
    ggplot(aes(x = Date, y = Px_n, fill = Px_n)) +
        geom_line(color="grey") +
        geom_point(shape=21, color="black", size=6) +
        scale_fill_viridis(discrete = TRUE, alpha = 0.6)
        theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
        theme_bw() +
        theme(axis.ticks.x = element_text(angle = 90)
        )
        

plot(time_plot)
