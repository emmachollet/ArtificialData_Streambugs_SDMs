## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## --- Create all necessary directories
##     project of comparison of sdms on artificial data ---
##
## --- Decembre 2024 -- Emma Chollet ---
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 1. Set working directory to current file location 

# 2. Create directories

directories <- c(
    "00_utilities",
    "A1_inputs_streambugs",
    "A2_scripts_streambugs",
    "A3_outputs_streambugs",
    "B1_scripts_sdms",
    "B2_outputs_sdms",
    "C1_documentation"
)

for(directory in directories){
    dir.create(directory)
}

# Backwards: extract list of all files in all folders
list.files(path = ".", recursive = TRUE)
writeLines(list.files(".", recursive = TRUE), "folder_structure_full.txt")

