################################################################################
# FilterBepiPredRes     #
#########################
# Script uses a slightly modified BepiPred-2.0 output to separate each
# epitope candidate in a line + computes its average RSA, Coil probability
# and Epitope probability + annotates the peptide.
#
# Slightly modified BepiPred-2.0 output includes a column named "End"
# that marks the end of a peptide with an x.
#
# Written by Ana Paula Vargas R.
# Bioinformatics lab - UPCH. Lima, Peru.
# Last updated: 05-jan-20.
###############################################################################

Data_MSP7C <- read.csv(file = "BepiPred_WIP_MSP7C.csv", header = TRUE, stringsAsFactors = FALSE)


# Init
df <- data.frame(
  "Protein" = "",
  "Variant" = "",
  "Peptide" = "",
  "PosI" = "",
  "PosF" = "",
  "Len" = "",
  "Exp-Bur" = "",
  "RSA_Ave" = "",
  "CoilP_Ave" = "",
  "EpiP_Ave" = ""
)

# Start

i <- 1
Amino <- c()
IPos <- c()
EB <- c()
RSA <- c()
CoilP <- c()
EP <- c()
Vari <- c()
FPos <- c()

while (i <= nrow(Data_MSP7C)){
  # Collect info
  IPos <- c(IPos, Data_MSP7C$Position[i])
  Amino <- c(Amino, Data_MSP7C$AminoAcid[i])
  EB <- c(EB, Data_MSP7C$Exposed.Buried[i])
  RSA <- c(RSA, Data_MSP7C$RelativeSurfaceAccessilibity[i])
  CoilP <- c(CoilP, Data_MSP7C$CoilProbability[i])
  EP <- c(EP, Data_MSP7C$EpitopeProbability[i])
  Vari <- c(Vari, Data_MSP7C$Entry[i])

  if (Data_MSP7C$End[i] == "x"){
    # Set last position
    FPos <- c(Data_MSP7C$Position[i])

    # Migrate data to new, transitional df: de
    de <- data.frame(
      "Protein" = "MSP7C",
      "Variant" = as.character(Vari[1]),
      "Peptide" = paste(Amino, collapse = ""),
      "PosI" = paste(IPos[1]),
      "PosF" = paste(FPos[1]),
      "Len" = paste(length(Amino)),
      "Exp-Bur" = paste(EB, collapse = ""),
      "RSA_Ave" = paste(mean(RSA)),
      "CoilP_Ave" = paste(mean(CoilP)),
      "EpiP_Ave" = paste(mean(EP))
    )
    # Concat de + df
    df <- rbind(df, de)

    # Reset vars
    Amino <- c()
    IPos <- c()
    EB <- c()
    RSA <- c()
    CoilP <- c()
    EP <- c()
    Vari <- c()
    FPos <- c()
  }

  i <- i + 1
}

# Remove first empty row
df = df[-1,]

# Export as CSV
write.csv(df, "BP_Sorted.csv", row.names = FALSE)

###########################
