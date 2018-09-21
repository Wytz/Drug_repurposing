# Create the feature set from the pathway data

# Load the packages
require(data.table)
require(rjson)
require(plyr)
require(doMC)
registerDoMC(cores = 20)

# Read in the files
paths = fread("/home/wytze/ODEX_usecases/Drug_repurposing_LUMC/2018_04_06/Indirect paths.csv")
setkey(paths)
paths = unique(paths)

colnames(paths) = c("Drug", "DB_id", "Drug_id", "Drug_species", "Drug_ST", "Drug_SG",
                    "AB_preds", "AB_prov", "Intermediate", "Intermediate_id", "Intermediate_species",
                    "Intermediate_ST", "Intermediate_SG", "BC_preds", "BC_prov", "Disease", "Disease_CUI",
                    "Disease_id", "Disease_species", "Disease_ST", "Disease_SG", "Direct_preds", "Direct_prov", 
                    "Status")

paths$Intermediate_ST = gsub("\'", "\"", paths$Intermediate_ST)

repoDB = read.csv("ML5/RepoDB.csv", stringsAsFactors = F)

# Select the relevant subset of repoDB
repoDB = repoDB[repoDB$status %in% c("Approved", "Terminated"),]

# Subset to the drugs and diseases we were able to extract from the platform
repoDB = repoDB[repoDB$drug_id %in% paths$DB_id & repoDB$ind_id %in% paths$Disease_CUI,]

# Create the features
st = unique(unlist(lapply(paths$Intermediate_ST, fromJSON)))
sg = unique(paths$Intermediate_SG)

# Create emtpy feature space

features = as.data.frame(matrix(nrow = 0, ncol = length(st) + length(sg) + 6))
colnames(features) = c("drug_id", "drug_id_EKP", "direct", make.names(st), make.names(sg), "disease_id", "disease_id_EKP", "Status")

# Create the feature set from the data
feature_set = foreach(i = 1:nrow(repoDB), .packages ="rjson") %dopar% {
  paths_subset = paths[paths$DB_id == repoDB$drug_id[i] & paths$Disease_CUI == repoDB$ind_id[i]]
  
  # Take only the relevant columns of this subset, and make them unique (Paths can apparently be split up sometimes)
  paths_relevant_subset = unique(paths_subset[,c("Intermediate_id", "Intermediate_ST", "Intermediate_SG"), with = F])
  
  if(nrow(paths_relevant_subset) == 0){
    out = NULL
  } else {
  subset_sts = unlist(lapply(paths_relevant_subset$Intermediate_ST, function(x){unique(fromJSON(x))}))
  subset_sts = as.data.frame(table(subset_sts))
  subset_sgs = as.data.frame(table(paths_relevant_subset$Intermediate_SG))
  add_features = as.data.frame(t(data.frame(feature = c(subset_sts$Freq, subset_sgs$Freq))))
  colnames(add_features) = make.names(c(as.character(subset_sts$subset_sts), as.character(subset_sgs$Var1)))
  out = rbind.fill(features, add_features)
  out[,c("drug_id", "drug_id_EKP", "direct", "disease_id", "disease_id_EKP", "Status")] = 
    c(repoDB$drug_id[i], toJSON(unique(paths_subset$Drug_id)), 
      unique(ifelse(paths_subset$Direct_preds == "[]", FALSE, TRUE)), repoDB$ind_id[i], 
      toJSON(unique(paths_subset$Disease_id)), repoDB$status[i])
  }
  out_line = out
}

# Remove the instances which did not have any indirect paths between them 
feature_set = feature_set[-(which(sapply(feature_set,is.null),arr.ind=TRUE))]
feature_set = do.call(rbind.data.frame, feature_set)

# Transform the NA values to 0
feature_set[is.na(feature_set)] = 0

# Check how many duplicate rows there are, and remove these
feature_set = unique(feature_set)

# Write away the data
write.csv2(feature_set, paste("ML final/Features based on semantic types and semantic groups ", Sys.Date(), ".csv", sep = ""), row.names = F)

# Remove the duplicated instances. Use the Euretos ID's instead of the DrugBank ID's, due to their deprecation. Write away the data as a second feature set
feature_set = feature_set[!duplicated(feature_set[,c("disease_id_EKP", "drug_id_EKP", "Status")]), ]
write.csv2(feature_set, paste("ML final/Features based on semantic types and semantic groups, no duplicate drugs or diseases ", Sys.Date(), ".csv", sep = ""), row.names = F)

# Clean up
paths = NULL
