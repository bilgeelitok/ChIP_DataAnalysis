PROJECT_DIR <- "/server/Bilge_GEP" #User project location.
DATA_DIR <- "/server/GEP" #Location of raw data.

MOTIF_PATH <- file.path(DATA_DIR,"motifs/GEP.motifs.zip")

WORKSPACE_PATH <- file.path(PROJECT_DIR,"MOTIFS")

unzip(MOTIF_PATH, exdir = WORKSPACE_PATH, overwrite = TRUE)
allDIRs<- list.files(path=WORKSPACE_PATH, full.names = TRUE)

#Only known results will be used for downstream analysis.
knownResults_list <- list()

for (dir in allDIRs) {
  #Get dir name:
  dir_name <- basename(dir)
  #Get clean names
  clean_name <- gsub("_trimmed.TD.ANTIBODY.regions.txt_motifs","", dir_name)
  clean_name <- gsub("DU145", "", clean_name)
  
  knownResult <- file.path(dir, "knownResults.txt")
  
  #Check if file exists.
  if (file.exists(knownResult)) {
    #cat("\n SUCCESS: knownResults.txt found for:", clean_name)
    knownResults_list[[clean_name]] <- knownResult
    
  } else {
    cat("\n ERROR: No knownResults.txt file for:", clean_name)
  }
  
  
}

#Motif names (n=7) to look for:
GEP_motifs <- c("MOTIF1/Homer",
                "MOTIF2/Homer",
                "MOTIF3/Homer",
                "MOTIF4/Homer",
                "MOTIF5/Homer",
                "MOTIF6/Homer",
                "MOTIF/Homer"
)

#Filter knownResults and store the filtered dataframes.
knownResults_filtered <- list()
sample_moi_only <- list()

for (name in names(knownResults_list)) {
  sample <- name
  knownResult <- read.delim(knownResults_list[[name]], 
                            header = TRUE, 
                            sep = "\t", 
                            stringsAsFactors = FALSE)
  #Use only motifs with p < 0.001 and over 5% enrichment
  sample_filt <- knownResult[(knownResult$P.value < 0.0001 & knownResult$X..of.Target.Sequences.with.Motif > 5),]
  
  #Remember R is more precise and excel might overestimate.
  print(paste("Sample", name, "has", (dim(sample_filt)[1]), "motifs passed the filters."))
  
  #Store filtered results intact.
  knownResults_filtered[[name]] <- sample_filt
  
  #Continue with only the motifs of interest (moi).
  sample_moi <- sample_filt[sample_filt$Motif.Name %in% GEP_motifs, ]
  sample_moi_only[[name]] <- sample_moi
}

moi_matrix <- matrix(data = FALSE,
                     nrow = length(sample_moi_only),
                     ncol = length(GEP_motifs),
                     dimnames = list(names(sample_moi_only))
)

#Each row represents a sample:
row.names(moi_matrix) <- names(sample_moi_only)
#Each column represents a motif of interest:
colnames(moi_matrix) <- GEP_motifs

#If motif is present in one sample fill in TRUE:
for (n in names(sample_moi_only)) {
  found <- sample_moi_only[[n]]$Motif.Name
  for (m in GEP_motifs){
    moi_matrix[n, m] <- ifelse(m %in% found, TRUE, FALSE)
  }
}

#This table goes to the report.
result_path <- file.path(PROJECT_DIR,"MOTIFS_in_SAMPLES.csv")
moi_df <- as.data.frame(moi_matrix)
write.csv(moi_df, result_path, quote = FALSE)

print(moi_df)