library(tidyverse)

# Load complete dataset
occ <- read_tsv("./data/occurrences/all_occurrences.csv") # works now??

# Create a smaller data set
write_csv(head(occ), "data/occurrences/all_occurrences_mini.csv")

# Try to reexport
write_tsv(head(occ, 100000), "data/occurrences/all_occurrences2.csv")

# Check exported dataset
occ2 <- read_tsv("./data/occurrences/all_occurrences2.csv") # wrong number of rows??
head(occ, 100000) %>% nrow()

# Investigate problems
problem_rows <- problems(occ2) %>% pull(row)
problem_rows_buffered <- sort(c(problem_rows-1, problem_rows, problem_rows+1))
problem_occ <- slice(occ, problem_rows_buffered)
problem_occ
write_tsv(problem_occ, "./data/occurrences/all_occurrences_problems.csv")

# Re-read problems?
occ3 <- read_tsv("./data/occurrences/all_occurrences_problems.csv")
all(occ3 == problem_occ, na.rm = TRUE)
occ3 == problem_occ %>% View()
occ3[which(occ3 != problem_occ)
