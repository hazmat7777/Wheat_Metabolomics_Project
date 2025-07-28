##### STEP 1- GET THE PATHWAYS FILE INTO A READABLE FORM

# Read the entire .dat file as one text vector
lines <- readLines("../data/metacyc_data/pathways.dat")

# Collapse into one long string and split by "//"
entries <- strsplit(paste(lines, collapse = "\n"), split = "//\\s*")[[1]]

# Function to parse an individual entry into a named list
parse_entry <- function(entry) {
  lines <- strsplit(entry, "\n")[[1]]
  lines <- trimws(lines)
  lines <- lines[lines != ""]  # remove empty lines
  
  parsed <- list()
  for (line in lines) {
    if (grepl(" - ", line)) {
      keyval <- strsplit(line, " - ", fixed = TRUE)[[1]]
      key <- keyval[1]
      val <- keyval[2]
      parsed[[key]] <- c(parsed[[key]], val)  # allows multiple values for same key
    } else {
      # Continuation lines: append to previous key
      if (length(parsed) > 0) {
        last_key <- tail(names(parsed), 1)
        parsed[[last_key]] <- paste(parsed[[last_key]], line)
      }
    }
  }
  return(parsed)
}

# Parse all entries
parsed_pathways <- lapply(entries, parse_entry)

# Example: show only entries with REACTION-LAYOUT
with_reaction_layout <- Filter(function(x) "REACTION-LAYOUT" %in% names(x), parsed_pathways)

length(parsed_pathways) # only 3646 pathways total
length(with_reaction_layout) # filtered a couple out

# View the first one
print(with_reaction_layout[[1]]) # have created a list of lists

with_reaction_layout[[1]]$"UNIQUE-ID"

# save this as a new rda file
saveRDS(parsed_pathways, "../data/metacyc_data/pathways_parsed.rds")

##### STEP 2- FILTER THE PATHWAYS DOWN TO THE ONES PREDICTED IN MY COMMUNITY

# Read in the predicted pathways file
predicted_pathways <- read.delim("../results/picrust2/farmkits_16S/picrust2_out_pipeline_split/pathways_out/path_abun_unstrat.tsv", header=TRUE)

# Extract pathways present 
pathways_present <- predicted_pathways[,1]

length(pathways_present) # 327

# Filter pathways list to only those present in samples
pathways_present_data <- Filter(function(x) x$"UNIQUE-ID" %in% pathways_present, parsed_pathways)

length(pathways_present_data) #325

## finding the ones missing from parsed_pathways
database_unique_ids <- sapply(parsed_pathways, function(x) x$"UNIQUE-ID")
missing_pathways <- setdiff(pathways_present, database_unique_ids)

# Print the missing pathways
print(missing_pathways) # dont know why these 2 are missing and idk what to do. Not that important



###### STEP 3- GET A LIST OF ALL COMPOUNDS FOUND IN THE PATHWAYS PRESENT

pathways_present_data[[length(pathways_present_data)]]

# Assuming parsed_pathways is a list of lists and each list contains information about different reactions
# We will extract the compounds from the "LEFT-PRIMARIES" and "RIGHT-PRIMARIES" sections of each reaction

# Create an empty list to store all extracted compounds
all_compounds <- list()

# Iterate through each pathway
for (pathway in pathways_present_data) {
  # Iterate through each reaction in the pathway
  for (reaction in pathway$"REACTION-LAYOUT") {
      
      # Use regular expressions to extract compounds from 'LEFT-PRIMARIES' and 'RIGHT-PRIMARIES'
      # Pattern to extract the compounds
      left_compounds <- gsub(".*:LEFT-PRIMARIES ([^\\)]+).*", "\\1", reaction)
      right_compounds <- gsub(".*:RIGHT-PRIMARIES ([^\\)]+).*", "\\1", reaction)
      
      # Split the compounds into individual compounds (assuming they are separated by commas or spaces)
      left_compounds <- strsplit(left_compounds, ",")[[1]]
      right_compounds <- strsplit(right_compounds, ",")[[1]]
      
      # Add the compounds to the list of all compounds (combining both sides)
      all_compounds <- c(all_compounds, left_compounds, right_compounds)
    }
}

# Remove duplicates from the list of compounds
all_compounds <- unique(all_compounds)

# View the first few compounds
all_compounds


# Create an empty list to store all extracted compounds
all_compounds <- list()

# Create an empty list to store all extracted compounds
all_compounds <- list()

# Iterate through each pathway
for (pathway in pathways_present_data) {
  # Iterate through each reaction in the pathway
  for (reaction in pathway$"REACTION-LAYOUT") {
      
    # Check if the reaction has valid LEFT-PRIMARIES and RIGHT-PRIMARIES
    if (grepl(":LEFT-PRIMARIES", reaction) && grepl(":RIGHT-PRIMARIES", reaction)) {
      
        # Extract compounds from LEFT-PRIMARIES and RIGHT-PRIMARIES
        left_compounds <- gsub(".*:LEFT-PRIMARIES ([^\\)]+).*", "\\1", reaction)
        right_compounds <- gsub(".*:RIGHT-PRIMARIES ([^\\)]+).*", "\\1", reaction)

        # Check if any compounds are found in either LEFT-PRIMARIES or RIGHT-PRIMARIES
        if (left_compounds == reaction || left_compounds == "") {
            left_compounds <- character(0)  # Set to empty vector if no compounds
        }

        if (right_compounds == reaction || right_compounds == "") {
            right_compounds <- character(0)  # Set to empty vector if no compounds
        }

        # Split the compounds into individual compounds only if they are not empty
        if (length(left_compounds) > 0) {
            left_compounds <- strsplit(left_compounds, ",")[[1]]
        } else {
            left_compounds <- character(0)  # Set to empty vector if no valid compounds
        }

        if (length(right_compounds) > 0) {
            right_compounds <- strsplit(right_compounds, ",")[[1]]
        } else {
            right_compounds <- character(0)  # Set to empty vector if no valid compounds
        }

        # Only add non-empty compound entries to the all_compounds list
        if (length(left_compounds) > 0) {
            all_compounds <- append(all_compounds, left_compounds)  # Use append for efficiency
        }
        if (length(right_compounds) > 0) {
            all_compounds <- append(all_compounds, right_compounds)  # Use append for efficiency
        }
    }
  }
}

# Remove duplicates from the list of compounds
all_compounds <- unique(all_compounds)

# Flatten list
all_compounds<- unlist(all_compounds)
all_compounds



###### STEP 4- Read in the metacyc compounds data, correctly parsed

# Read the entire .dat file as one text vector
lines <- readLines("../data/metacyc_data/compounds.dat", encoding = "ISO-8859-1")
lines <- iconv(lines, "UTF-8") # change encoding

# Collapse (merge) into one long string and split by "//"
entries <- strsplit(paste(lines, collapse = "\n"), split = "//\\s*")[[1]] 

length(entries) # info on 26k compounds

# Parse each entry in the dataset
all_compounds_data <- lapply(entries, parse_entry) # using fn defined for the pathways.dat

print(all_compounds_data)

# example of a formula in this list
all_compounds_data[[26175]]$`CHEMICAL-FORMULA`
class(all_compounds_data[[26175]]$`CHEMICAL-FORMULA`)

# save this as a new rda file
saveRDS(all_compounds_data, "../data/metacyc_data/compounds_parsed.rds")


##### STEP 5- read in and process the metabolomics data

# am now doing this in a separate file

mass_spec_df <- read.csv("../data/20250327_Soil_Bact_OB.csv")
    
unique(mass_spec_df$"Reference.Ion") # Calc.MW accounts for the method of ionization

colnames(mass_spec_df)

# filter out the unrelated the soybean study
soil_columns <- grep("oil", colnames(mass_spec_df), value = TRUE)

# Subset the data frame to include only soil-related columns
mass_spec_df <- mass_spec_df[, c("Tags", "Checked", "Name", "Formula", "Annot..Source..Predicted.Compositions", 
                              "Annot..Source..mzCloud.Search", "Annot..Source..mzVault.Search", "Annot..DeltaMass..ppm.", 
                              "Calc..MW", "m.z", "RT..min.", "Area..Max..", "X..mzCloud.Results", 
                              "X..mzVault.Results", "mzCloud.Best.Match", "mzCloud.Best.Match.Confidence", 
                              "mzCloud.Best.Sim..Match", "mzCloud.Best.Tree.Match", "mzVault.Best.Match", 
                              "mzCloud.Library.Match..Autoprocessed", "mzCloud.Library.Match..Reference", 
                              "mzVault.Library.Match..20220905_ALL_GNPS_v2c_with_CompoundMetadata_20240813", 
                              "mzVault.Library.Match..20240730_IMPACT_ALL_v2_with_CompoundMetadata", 
                              "MS2", "Reference.Ion", soil_columns)]

# how many entries have associated compound names?

soil_columns