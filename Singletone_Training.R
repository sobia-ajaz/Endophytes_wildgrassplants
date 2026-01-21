library(readxl)
library(rstudioapi)
library(dplyr)

setwd(dirname(getActiveDocumentContext()$path))
# Read the Excel file (replace with your file path)
data <- read_excel("feature-table.xlsx")

# Check the structure of the data to confirm column count
str(data)


# Assuming you want to count values from columns 2 to the last one
# Adjust the column range accordingly
num_columns <- ncol(data)  # Get the total number of columns

# If you want to count values greater than 0 in all columns from 2 to the last column
data$Count <- apply(data[, 2:num_columns], 1, function(x) sum(x > 0))

# Add the 'MaxValue' column: count the maximum value in columns 2 to the last one
data$MaxValue <- apply(data[, 2:num_columns], 1, max, na.rm = TRUE)

# View the updated data with the count column added
head(data)

# Assuming your data is stored in the 'data' dataframe
singletone_data <- data %>%
  filter(Count == 1, MaxValue <= 100)

# View the singletone_data
head(singletone_data)

# Write the singletone_data to a tab-delimited .txt file
write.table(singletone_data, "singletone_feature-table.txt", sep = "\t", row.names = FALSE, quote = FALSE)




# Filter the data: Count > 1 and MaxValue >= 100
filtered_data <- data %>%
  filter(Count > 1, MaxValue >= 100)

# View the filtered data
head(filtered_data)

# Optionally, write the filtered data to a tab-delimited .txt file
write.table(filtered_data, "filtered_feature-table.txt", sep = "\t", row.names = FALSE, quote = FALSE)



