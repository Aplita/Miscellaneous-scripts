
################################################################################
# RNAseq change names     #
###########################
# This script annonates gene names after an RNA-seq analysis.
# Uses as input a .csv result table from edgeR.
#
# Note: working, but can be optimized.
#
# Written by Ana Paula Vargas R.
# Bioinformatics lab - UPCH. Lima, Peru.
# Last updated: 13-feb-20.
###############################################################################


# Read csv file generated with edgeR
File <- read.csv("pH7-PZA100vs0.csv", stringsAsFactors = FALSE, header = TRUE, fill = TRUE)
ANOT <- read.table("NEW_h37rv.gff3", stringsAsFactors = FALSE, fill = TRUE, sep = "\t", header = FALSE, quote = "",
                   col.names = c("Chrom", "Ena", "Cod", "StartPos", "EndPos", "S1", "S2", "S3", "Info"))
ANOT <- ANOT[grep("gene", ANOT$Cod),]

# Init variables for loop
C1 <- 1
Start_VNames <- c()
End_VNames <- c()

while (C1 <= length(File$GeneID)) {
  # Parse start/end positions by transcript
  SPos <- unlist(strsplit(File$Start[C1], ";"))
  EPos <- unlist(strsplit(File$End[C1], ";"))
  StartIndex <- c()
  EndIndex <- c()

  C2 <- 1

  # Parse row number of each start/end position in ANOT table
  while (C2 <= length(SPos)) {
    # Get value for C2 index
    StartIndex1 <- which(ANOT$StartPos == as.integer(SPos[C2]))
    EndIndex1 <- which(ANOT$EndPos == as.integer(EPos[C2]))

    # Check whether it exists; if it does, add value; if not, add flag
    if (length(StartIndex1) == 0) {
      StartIndex <- c(StartIndex, "NotFound")
    } else {
      StartIndex <- c(StartIndex, StartIndex1)
    }

    if (length(EndIndex1) == 0) {
      EndIndex <- c(EndIndex, "NotFound")
    } else {
      EndIndex <- c(EndIndex, EndIndex1)
    }

    C2 <- C2 + 1
  }

  # Are there undetermined values (Start Index)?
  if ("NotFound" %in% StartIndex) {
    C3 <- 1
    Start_VNames1 <- c()
    while(C3 <= length(StartIndex)){
      if (StartIndex[C3] == "NotFound") {
        Start_VNames1 <- c(Start_VNames1, "NotFound")
      } else {
        x <- as.integer(StartIndex[C3])
        a <- strsplit(ANOT$Info[x], ";")
        a <- strsplit(unlist(lapply(a, "[[", 1)), ":")
        a <- unlist(lapply(a, "[[", 2))
        a <- paste(a, collapse = ";")
      }
      C3 <- C3 + 1
    }
    a <- paste(c(Start_VNames1, b), collapse = ";")
    Start_VNames <- c(Start_VNames, a)
  } else {
    # Get gene names for each pos
    a <- strsplit(ANOT$Info[StartIndex], ";")
    a <- strsplit(unlist(lapply(a, "[[", 1)), ":")
    a <- unlist(lapply(a, "[[", 2))
    # Merge names in one vector, separated by ; for easy parsing
    a <- paste(a, collapse = ";")
    # Add vector of names to a master vector
    Start_VNames <- c(Start_VNames, a)
  }

  # Are there undetermined values (End Index)?
  if ("NotFound" %in% EndIndex) {
    C3 <- 1
    End_VNames1 <- c()
    while(C3 <= length(EndIndex)){
      if (EndIndex[C3] == "NotFound") {
        End_VNames1 <- c(End_VNames1, "NotFound")
      } else {
        x <- as.integer(EndIndex[C3])
        b <- strsplit(ANOT$Info[x], ";")
        b <- strsplit(unlist(lapply(b, "[[", 1)), ":")
        b <- unlist(lapply(b, "[[", 2))
        b <- paste(b, collapse = ";")
      }
      C3 <- C3 + 1
    }
    b <- paste(c(End_VNames1, b), collapse = ";")
    End_VNames <- c(End_VNames, b)
  } else {
    b <- strsplit(ANOT$Info[EndIndex], ";")
    b <- strsplit(unlist(lapply(b, "[[", 1)), ":")
    b <- unlist(lapply(b, "[[", 2))
    b <- paste(b, collapse = ";")
    End_VNames <- c(End_VNames, b)
  }
  C1 <- C1 + 1
}


# Bind master vector to df and export as csv.
File <- cbind(File, data.frame(StartPosNames = Start_VNames, stringsAsFactors = FALSE))
File <- cbind(File, data.frame(EndPosNames = End_VNames, stringsAsFactors = FALSE))
write.csv(File, "pH7-PZA100vs0_NAMES.csv", row.names = FALSE)
