#install.packages("stringr")  # Use double quotes
#install.packages("e1071", dependencies = TRUE)
library(stringr)
library(e1071)

# Ex1 

# 1st method
func_reverse <- function(x) sapply(lapply(strsplit(x, ""), rev), paste, collapse = "")

func_palindrome <- function(phrase) {
  ifelse(func_reverse(phrase) == phrase, "True", "False")
}

no <- 1221
func_palindrome(as.character(no))  # Convert number to string before passing

# 2nd method 

is_palindrome <- function(num) {
  str_num <- as.character(num)  # Convert number to string
  rev_num <- paste(rev(strsplit(str_num, "")[[1]]), collapse = "")  # Reverse string
  
  return(str_num == rev_num)  # Compare original and reversed string
}

# Test cases
is_palindrome(1221)   # TRUE
is_palindrome(1234)   # FALSE
is_palindrome(0)      # TRUE
is_palindrome(101)    # TRUE


# Ex2 
text <- "seemerightnow"

# (a) Extract "see" (1st to 3rd character)
sub_a <- substr(text, 1, 3)

# (b) Extract "me" (4th to 5th character)
sub_b <- substr(text, 4, 5)

# (c) Extract "right" (7th to 11th character)
sub_c <- substr(text, 7, 11)

# Print the results
sub_a  # "see"
sub_b  # "me"
sub_c  # "right"


# Ex3
gc_fraction <- function(sequence) {
  sequence <- toupper(sequence)  # Ensure uppercase
  gc_count <- sum(strsplit(sequence, "")[[1]] %in% c("G", "C"))  # Count G and C
  return(gc_count / nchar(sequence))  # Fraction of GC bases
}

# Example usage
seq <- "ATTGCGCATAGTCCGGG"
gc_fraction(seq)


# Ex4
is_palindromic_dna <- function(seq) {
  # Define complement pairs
  complement_map <- c("A" = "T", "T" = "A", "G" = "C", "C" = "G")
  
  # Convert sequence into a character vector
  seq_vector <- strsplit(seq, "")[[1]]
  
  # Get the complementary sequence
  complement_seq <- sapply(seq_vector, function(base) complement_map[base])
  
  # Reverse the complement
  reverse_complement <- rev(complement_seq)
  
  # Join back into a string
  reverse_complement_str <- paste(reverse_complement, collapse = "")
  
  # Compare original with its reverse complement
  return(seq == reverse_complement_str)
}

# Example usage
is_palindromic_dna("TGGATCCA")  # TRUE
is_palindromic_dna("ATGC")      # FALSE


# Ex5

find_largest_words <- function(sentence) {
  # Step 1: Split the sentence into words
  words <- strsplit(sentence, " ")[[1]]  # Split by space
  
  # Step 2: Remove punctuation manually (basic method)
  words <- gsub("[,.!?;]", "", words)  # Remove common punctuation
  
  # Step 3: Compute word lengths
  word_lengths <- nchar(words)
  
  # Step 4: Find the maximum and second maximum word lengths
  max_length <- max(word_lengths)
  second_max_length <- max(word_lengths[word_lengths < max_length])
  
  # Step 5: Extract words of these lengths
  largest_words <- words[word_lengths == max_length]
  second_largest_words <- words[word_lengths == second_max_length]
  
  return(list("Largest Words" = largest_words, "Second Largest Words" = second_largest_words))
}

# Example usage
sentence <- "She sells hundreds of sea oysters on the sea shore."
find_largest_words(sentence)



# Ex6 Load the data from 'worldfloras.txt'

world_data <- read.table("/home/ibab/Downloads/worldfloras.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Check structure
str(world_data)

# Create subsets for each continent
asia_data <- subset(world_data, Continent == "Asia")
europe_data <- subset(world_data, Continent == "Europe")
africa_data <- subset(world_data, Continent == "Africa")
america_data <- subset(world_data, Continent == "America")
oceania_data <- subset(world_data, Continent == "Oceania")

# View sample data
head(asia_data)

# (6b) Boxplot for floral count distribution within each continent
par(mar = c(5, 5, 5, 5))  # Adjust margins to avoid "figure margins too large" error
boxplot(Flora ~ Continent, data = world_data, 
        main = "Floral Count Distribution by Continent", 
        xlab = "Continent", ylab = "Floral Count", 
        col = "lightgreen", border = "black")

# Statistical summary
summary_stats <- aggregate(Flora ~ Continent, data = world_data, summary)
print(summary_stats)

# Calculate mean and standard deviation
mean_floral <- mean(world_data$Flora, na.rm = TRUE)
sd_floral <- sd(world_data$Flora, na.rm = TRUE)
print(mean_floral)
print(sd_floral)

# Skewness and kurtosis
floral_skewness <- skewness(world_data$Flora, na.rm = TRUE)
floral_kurtosis <- kurtosis(world_data$Flora, na.rm = TRUE)
print(floral_skewness)  # Skewness interpretation
print(floral_kurtosis)  # Kurtosis interpretation

# (6c) Boxplot and histogram for population distribution
boxplot(Population ~ Continent, data = world_data, 
        main = "Population Distribution by Continent", 
        xlab = "Continent", ylab = "Population", 
        col = "lightgreen", border = "black")

hist(world_data$Population, 
     main = "Population Distribution", 
     xlab = "Population", 
     col = "purple", border = "black")

# Statistical summary for population
population_stats <- aggregate(Population ~ Continent, data = world_data, summary)
print(population_stats)

# Skewness and kurtosis for population
population_skewness <- skewness(world_data$Population, na.rm = TRUE)
population_kurtosis <- kurtosis(world_data$Population, na.rm = TRUE)
print(population_skewness)
print(population_kurtosis)


# Ex7 Read in the data from 'HumanBones.txt'

bones_data <- readLines("/home/ibab/Downloads/HumanBones.txt")

# Initialize vectors
categories <- c()
bone_names <- c()
bone_numbers <- c()
current_category <- NULL

# Process each line
for (line in bones_data) {
  if (!grepl("\\(", line)) {
    current_category <- trimws(line)  # Set new category
  } else {
    bone_info <- strsplit(line, "\\(")[[1]]
    bone_name <- trimws(bone_info[1])  
    bone_number <- as.numeric(str_extract(bone_info[2], "\\d+"))  # Extract only first number
    categories <- c(categories, current_category)
    bone_names <- c(bone_names, bone_name)
    bone_numbers <- c(bone_numbers, bone_number)
  }
}

# Create a dataframe
bones_info <- data.frame(category = categories, name_of_bone = bone_names, number_of_bones = bone_numbers, stringsAsFactors = FALSE)
print(head(bones_info))

# (8) Find category with the maximum number of bones
category_bones_summary <- aggregate(number_of_bones ~ category, data = bones_info, sum)
max_category <- category_bones_summary[which.max(category_bones_summary$number_of_bones), ]
print(max_category$category)

# Create a frequency table
category_frequency <- table(bones_info$category)
print(category_frequency)

# Bar plot
barplot(category_bones_summary$number_of_bones, names.arg = category_bones_summary$category, 
        main = "Number of Bones by Category", xlab = "Category", ylab = "Number of Bones", 
        col = "pink", border = "black")

# (9) Subset category "Legs" and find bones with names longer than 5 letters
legs_data <- subset(bones_info, category == "Legs")
long_bones_legs <- subset(legs_data, nchar(name_of_bone) > 5)
print(long_bones_legs$name_of_bone)

# (10) List bones starting with "M" and replace "a" with "A"
bones_starting_M <- subset(bones_info, grepl("^M", name_of_bone))
bones_starting_M$name_of_bone <- gsub("a", "A", bones_starting_M$name_of_bone)
print(bones_starting_M$name_of_bone)

# (11) List bones ending with "e" and convert them to lowercase
bones_ending_with_e <- subset(bones_info, grepl("e$", name_of_bone))
bones_ending_with_e$name_of_bone <- tolower(bones_ending_with_e$name_of_bone)
print(bones_ending_with_e$name_of_bone)

# (12) List bones with two "o"s in their names
bones_with_two_o <- subset(bones_info, grepl("o.*o", name_of_bone, ignore.case = TRUE))
print(bones_with_two_o$name_of_bone)
