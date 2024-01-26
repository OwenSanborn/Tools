library(rvest)
library(dplyr)
library(stringr)
library(ggplot2)

process_output <- function(path) {
  files <- dir(path)
  data_frames <- list()
  
  for (file in files) {
    # Read the HTML file
    html <- read_html(paste0(path, file))
  
    # Extract the XPath value
    xpath_value <- html_text(html_element(html, xpath = '//*[@id="for-use-in-r"]/pre[2]/code'))
  
    # Clean up the extracted text
    xpath_value <- gsub("\n", "", xpath_value)
    xpath_value <- gsub("\\\\", "", xpath_value)
  
    data_frame <- eval(parse(text = xpath_value))

    file <- gsub(".html", "", file)
  
    data_frame$Group <- gsub("\\..*", "", file)
    data_frame$Rep <- gsub("^.*\\.","", file )
 
    # Append the data frame to the list
    data_frames[[file]] <- data_frame
  }

  # Combine all data frames into a single data frame
  combined_df <- do.call(rbind, data_frames)

  positions <- c("3", "8", "13", "18", "23", "28", "33")
  combined_df <- subset(combined_df, guide.position %in% positions)
  
  return(combined_df.i)
}
