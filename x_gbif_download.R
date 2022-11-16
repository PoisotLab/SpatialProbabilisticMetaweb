library(tidyverse)

# Load complete dataset
occ <- read_tsv("./data/occurrences/subset/all_occurrences.csv") # works now??

# Create a smaller data set
write_csv(head(occ), "data/occurrences/subset/all_occurrences_mini.csv")

# Try to reexport
write_tsv(head(occ, 100000), "data/occurrences/subset/all_occurrences2.csv")

# Check exported dataset
occ2 <- read_tsv("./data/occurrences/subset/all_occurrences2.csv") # wrong number of rows??
head(occ, 100000) %>% nrow()

# Investigate problems
problem_rows <- problems(occ2) %>% pull(row)
problem_rows_buffered <- sort(c(problem_rows-1, problem_rows, problem_rows+1))
problem_occ <- slice(occ, problem_rows_buffered)
problem_occ
write_tsv(problem_occ, "./data/occurrences/subset/all_occurrences_problems.csv")

# Re-read problems?
occ3 <- read_tsv("./data/occurrences/subset/all_occurrences_problems.csv")
all(occ3 == problem_occ, na.rm = TRUE)
occ3 == problem_occ %>% View()
occ3[which(occ3 != problem_occ)

## Temporary workaround
# Keep only basic columns, separate by species, export to separate files

# First use a smaller set
# occ_mini <- head(occ, 1000)

# Now select only the same columns as in your initial analyses
# occ_select <- occ_mini |>
occ_select <- occ |>
  select(species, decimalLatitude, decimalLongitude) |>
  rename(name = species,
         latitude = decimalLatitude,
         longitude = decimalLongitude)

# Separate by species
# occ_select %>%
#   group_by(name) %>%
#   # summarize(test = length(name)) %>%
#   group_split()

# Get species names
sp_list <- unique(occ_select$name)

# Export as separate files
for (sp in sp_list) {
  sp_name <- str_replace(sp, " ", "_")
  sp_path <- file.path("data", "occurrences", paste0(sp_name, ".csv"))
  occ_select |>
    filter(name == sp) |>
    write_tsv(sp_path)
}
# 67 MB
# 1,579,706 obs
