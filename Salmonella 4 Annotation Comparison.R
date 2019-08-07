library(magrittr)
library(tibble)
library(tools)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(readxl)
library(MASS)
library(Matrix)
library(reshape2)
library(gplots)
library(devtools)
library(plyr)
#source("https://bioconductor.org/biocLite.R")
# biocLite("muscle")
# biocLite("Biostrings")
library(Biostrings)
library(muscle)
library(tidyverse)

StartTime <- Sys.time()

###############
## Functions ##
###############

'%nin%' <- Negate('%in%')

uniquelength <- function(x) {x %>% unique() %>% length()}

# Rotate the elements of a vector
rot <- function(x, l) (x[c(2:l,1)])

# Filter a blast hit table for the best unique entries for each query and subject - qseqid and ssseqid can only appear once
# Blast hit table from comparison of 2 sets of CDS
FindBestBlastHits <- function(Hits) {
    # Assumes Hits is a df of blast hits

    # Remove any functions that do not match
    BestHits <- Hits %>% filter(qfunction == sfunction)

    # Pick best hit by qfraction * sfraction and by mismatch - ties are unlikely but possible
    BestHits %<>% mutate (qs = qfraction * sfraction)
    # BestHits %<>% group_by(qseqid) %>% filter(qs == max(qs))
    # BestHits %<>% group_by(sseqid) %>% filter(qs == max(qs))
    BestHits %<>% group_by(qseqid) %>% filter(qs == max(qs), mismatch == min(mismatch))
    BestHits %<>% group_by(sseqid) %>% filter(qs == max(qs), mismatch == min(mismatch))

    # Remove matches with low qfraction/sfraction
    BestHits %<>% filter(qs > 0.25)

    # Check for multiple qseqid entries and choose one at random
    DupRowNrs <- BestHits$qseqid %>% duplicated() %>% which()
    DupIDs <- BestHits[DupRowNrs,] %>% .$qseqid
    AllDupRowNrs <- (BestHits$qseqid %in% DupIDs) %>% which()
    if (length(DupIDs) > 0) {
        Keep <- vector()
        for (ID in DupIDs) {
            DupRowNrs2 <- (BestHits$qseqid == ID) %>% which() # now includes the first duplicate
            Keep <- c(Keep, sample(DupRowNrs2,1)) # choose one to keep
        }
        Lose <- AllDupRowNrs[-Keep] # from the duplicate row nrs in the original df, remove the one we want to keep
        BestHits %<>% .[-Lose,] # remove the others
        #BestHits$qseqid %>% unique() %>% length()
    }

    # Check for multiple sseqid entries and choose one at random
    DupRowNrs <- BestHits$sseqid %>% duplicated() %>% which()
    DupIDs <- BestHits[DupRowNrs,] %>% .$sseqid
    AllDupRowNrs <- (BestHits$sseqid %in% DupIDs) %>% which()
    if (length(DupIDs) > 0) {
        Keep <- vector()
        for (ID in DupIDs) {
            DupRowNrs2 <- (BestHits$sseqid == ID) %>% which() # now includes the first duplicate
            Keep <- c(Keep, sample(DupRowNrs2,1)) # choose one to keep
        }
        Lose <- AllDupRowNrs[-Keep] # from the duplicate row nrs in the original df, remove the one we want to keep
        BestHits %<>% .[-Lose,] # remove the others
        #BestHits$qseqid %>% unique() %>% length()
    }

    BestHits
}

# Filter a blast hit table for the best unique entries for each query and subject - qseqid and ssseqid can only appear once
# Blast hit table from comparison of a sets of CDS to a set of contigs
#Hits <- DNABlasts2[["KC_109_vs_T-Chr"]] %>% filter(cdsid == "cdsid0338")
FindBestBlastHits2 <- function(Hits) {
    # Assumes Hits is a df of blast hits
    BestHits <- Hits

    # Pick best hit by qfraction * "sfraction" - ties are unlikely but possible
    # Cannot use the lengths of subject contigs, so use length of subject match / length of query sequence
    #BestHits %<>% mutate (qs = qfraction * (slength/qfulllength))
    #BestHits %<>% mutate (qs = qfraction * min(slength/qfulllength, 1))
    #BestHits %<>% group_by(qseqid) %>% filter(qs == max(qs))
    BestHits %<>% group_by(qseqid) %>% filter(qfraction == max(qfraction))
    BestHits %<>% group_by(qseqid) %>% filter(qfraction == max(qfraction), pident == max(pident))

    # Remove matches with low qfraction/sfraction
    BestHits %<>% filter(qfraction > 0.25)

    # Check for multiple qseqid entries and choose one at random
    DupRowNrs <- BestHits$qseqid %>% duplicated() %>% which()
    DupIDs <- BestHits[DupRowNrs,] %>% .$qseqid
    AllDupRowNrs <- (BestHits$qseqid %in% DupIDs) %>% which()
    if (length(DupIDs) > 0) {
        Keep <- vector()
        for (ID in DupIDs) {
            DupRowNrs2 <- (BestHits$qseqid == ID) %>% which() # now includes the first duplicate
            Keep <- c(Keep, sample(DupRowNrs2,1)) # choose one to keep
        }
        Lose <- AllDupRowNrs[-Keep] # from the duplicate row nrs in the original df, remove the one we want to keep
        BestHits %<>% .[-Lose,] # remove the others
        #BestHits$qseqid %>% unique() %>% length()
    }

    BestHits
}

# Align 2 sequences, either amino acids or nucleotides
doMuscleComp <- function (Strings, Names, Quiet, Type) {
    if (Type == "AA") {
        Strings <- Biostrings::AAStringSet(Strings)
    } else if (Type == "DNA") {
        Strings <- Biostrings::DNAStringSet(Strings)
    } else {
        stop(sprintf("Expected type of 'AA' or 'DNA', got '%s'", Type))
    }
    names(Strings) <- Names
    Align <- muscle::muscle(Strings, quiet = Quiet)
}

# Load a subsequence from a contig(s) crr file by name of contig and position
GetSeqFromDecodeCRR <- function(Ref, Contig, Start, Stop) {
    # Takes a crr file, contig within the file, start pos and stop pos (1 based)
    # Extracts the sequence from the crr file using cgatools decodecrr
    # If start > stop, gets the complement
    # If there is an error, more than one line in the vector is returned
    # At the moment, function dies if this happens

    Range <- paste(Contig, min(Start, Stop)-1, max(Start, Stop), sep = ",")
    Snt <- system2(command = paste0("cgatools"), args = c( "decodecrr", "--reference", Ref, "--range", Range), stdout = TRUE, stderr = TRUE)

    if(length(Snt) == 1) {
        if (Start > Stop) {
            Snt <- chartr("ATGC", "TACG", Snt)
            Snt <- stringi::stri_reverse(Snt)
        }
    } else {
        print("Error retrieving sequence")
        stop()
    }

    return(Snt)
}

# Extract and format protein variants by comparing 2 protein sequences
FormatProteinVariant <- function(type, currseq, currseq2, pos) {
    type <- type %>% unique()
    currseq <- paste0(currseq, collapse = "")
    currseq2 <- paste0(currseq2, collapse = "")

    if (length(type) == 1) {
        if (type == "Del") {
            var <- sprintf("del:%d:-/%s", pos-nchar(currseq2)+1, currseq2)
        }
        else if (type == "Ins") {
            var <- sprintf("ins:%d:%s/-", pos-nchar(currseq)+1, currseq)
        }
        else if (type == "Sub") {
            var <- sprintf("sub:%d:%s/%s", pos-nchar(currseq)+1, currseq, currseq2)
        }
    } else {
        var <- sprintf("delins:%d:%s/%s", pos-nchar(currseq)+1, gsub("-", "", currseq), gsub("-", "", currseq2))
    }
}

# Extract and format DNA variants by comparing 2 DNA sequences
FormatDNAVariant <- function(type, currseq, currseq2, pos) {
    type <- type %>% unique()
    currseq <- paste0(currseq, collapse = "")
    currseq2 <- paste0(currseq2, collapse = "")

    if (length(type) == 1) {
        if (type == "Del") {
            var <- sprintf("del:%d:-/%s", pos-nchar(currseq2)+1, currseq2)
        }
        else if (type == "Ins") {
            var <- sprintf("ins:%d:%s/-", pos-nchar(currseq)+1, currseq)
        }
        else if (type == "Sub") {
            var <- sprintf("snp:%d:%s/%s", pos-nchar(currseq)+1, currseq, currseq2)
        }
    } else { # multiple types, must be a delins
        var <- sprintf("delins:%d:%s/%s", pos-nchar(currseq)+1, gsub("-", "", currseq), gsub("-", "", currseq2))
    }
}

####################
## Set components, genomes, paths ##
###################
# Component usually Chr then Plasmid_1 etc. Starts with cap, use tolower() to change
Components <- c("Chr", "Plasmid_1", "Plasmid_2", "Plasmid_3", "Plasmid_4")
NrComponents <- Components %>% length()

# Temp store of genome names
x <- c("KC_109", "T", "TM", "TMu")
NrGenomes <- length(x)

# Create a df of rotated genomes
Numbers <- c("One","Two","Three","Four","Five","Six","Seven","Eight","Nine","Ten") # if you need more make it longer

Genomes = matrix(NA_character_, nrow = NrGenomes, ncol = NrGenomes)
colnames(Genomes) <- Numbers[1:NrGenomes]

Genomes[,1] <- x
for (n in 2:NrGenomes) { Genomes[,n] <- rot(Genomes[,n-1], NrGenomes) }

Genomes %<>% as_data_frame()
Genomes

# Set file prefix
FilePrefix <- "Salmonella_4"

# Set paths
PATRICDirIn <- "/Users/rtearle/Documents/Roseworthy_Projects/6.Students/Nitish_Salmonella/Apr_2018/Out/PATRIC" # PATRIC annotations
DNABlastDirIn <- "/Users/rtearle/Documents/Roseworthy_Projects/6.Students/Nitish_Salmonella/Apr_2018/Out/BLAST_DNA_CDS_vs_CDS" # results of NCBI website blasts of feature DNAs against each other
DNABlastDirIn2 <- "/Users/rtearle/Documents/Roseworthy_Projects/6.Students/Nitish_Salmonella/Apr_2018/Out/BLAST_DNA_CDS_vs_Genome" # # results of NCBI website blasts of feature DNAs against other genomes
CrrDirIn <- "/Users/rtearle/Documents/Roseworthy_Projects/6.Students/Nitish_Salmonella/Holt_Pipeline_Feb_2018/Out/CG_CRR_Out" # CG crr files for each genome component file found in the bandage output
GGDirIn <- "/Users/rtearle/Documents/Roseworthy_Projects/6.Students/Nitish_Salmonella/May_2018/Out/Genome-vs-Genome" # results of genome-genome comparisons with minimap2

DNAAnnotDirOut <- "/Users/rtearle/Documents/Roseworthy_Projects/6.Students/Nitish_Salmonella/Apr_2018/Out/R_DNA_ANNOT" # DNA analysis contained in this script
SummaryDirOut <- "/Users/rtearle/Documents/Roseworthy_Projects/6.Students/Nitish_Salmonella/May_2018/Out/R_Summaries"
#dir.create(SummaryDirOut)

# NB crr files should already exist.
# If they don't, run "Make CRR Files from Fasta Files"

rm(x, Numbers, n)

##########################
## Get annotation files ##
##########################
FileSuffix <- "Annotation.txt"
ColNames <- c("ContigID", "FeatureID", "Type", "Location", "Start", "Stop", "Strand", "Function", "Aliases", "FigFam", "EvidenceCodes", "Nucleotide", "Protein")

AnnotationFilesIn <- vector("character", NrGenomes * NrComponents)
Annotations <- list() # converted to a df later

# Get annotation files
c <- 1
for (Component in Components) {
    for (n in 1:NrGenomes) {
        # Get file names
        AnnotationFilesIn[c] <- paste0("Salmonella-", Genomes$One[n], "-", Component, "-", FileSuffix)

        # Open files - note use of double brackets to extract name for use
        Name <- paste(Genomes$One[n], Component, sep = "-")
        Annotations[[Name]] = read_tsv(file.path(PATRICDirIn, Genomes$One[n], Component, AnnotationFilesIn[c]), col_names = TRUE, skip = 0) %>% as_data_frame()
        colnames(Annotations[[Name]]) <- ColNames

        # Remove repeats
        Annotations[[Name]] %<>% filter(Type == "CDS")

        # Add length, component, genome
        Annotations[[Name]] %<>% mutate(Length = abs(Stop - Start)+1)
        Annotations[[Name]] %<>% mutate(Component = Component)
        Annotations[[Name]] %<>% mutate(Genome = Genomes$One[n])

        # Fix Start and Stop so Start < Stop
        Annotations[[Name]] %<>% mutate(Start2 = ifelse(Start < Stop, Start, Stop))
        Annotations[[Name]] %<>% mutate(Stop2 = ifelse(Start > Stop, Start, Stop))
        Annotations[[Name]] %<>% mutate(Start = Start2, Stop = Stop2)
        Annotations[[Name]] %<>% select(-c(Start2, Stop2))

        # Sort
        Annotations[[Name]] %<>% arrange(Genome, Component, ContigID, Start)

        # Mark CDS at begin and end of contig - will use to check if it explains truncated CDS
        Annotations[[Name]]$EndOfContig = NA
        Annotations[[Name]]$EndOfContig %<>% as.integer()
        Contigs <- Annotations[[Name]]$ContigID %>% unique()

        for (Contig in Contigs) {
            Subset <- Annotations[[Name]] %>% filter(ContigID == Contig)
            # First
            FID <- Subset[1,]$FeatureID # first CDS
            Distance <- Subset[1,]$Start - 1 # distance to begin of contig
            Pos <- (Annotations[[Name]]$FeatureID == FID) %>% which() # row in original df
            Annotations[[Name]][Pos, "EndOfContig"] <- Distance # assign distance
            #y <- Annotations[[Name]] %>% head()
            # Last
            FID <- Subset[nrow(Subset),]$FeatureID # last CDS
            Length <- gsub(".+length_(.+)_cov.+", "\\1", Contig) %>% as.integer # length of contig from title
            Distance = Length - Subset[nrow(Subset),]$Stop # distance to end of contig
            Pos <- (Annotations[[Name]]$FeatureID == FID) %>% which() # row in original df
            Annotations[[Name]][Pos, "EndOfContig"] <- Distance # assign distance
            #y <- Annotations[[Name]] %>% tail()
        }

        Annotations[[Name]] %<>% select(Genome, Component, everything())

        c <- c+1
    }
}

AnnotationFilesIn %>% length()
Annotations %>% length()
rm(FileSuffix, ColNames, AnnotationFilesIn, Subset, FID, Length, Distance, Pos, c)
rm(Contig, Contigs)

# Make one df from annotations list
df <- matrix(nrow = 0, ncol = Annotations[[Name]] %>% ncol()) %>% as_data_frame()
colnames(df) <- colnames(Annotations[[Name]])
for (n in 1:length(Annotations)) {
    Name <- names(Annotations[n])
    df <- bind_rows(df, Annotations[[Name]])
}
Annotations <- df

rm(df)

Annotations %>% dim() # 20566  17

colnames(Annotations)
#  [1] "Genome"        "Component"     "ContigID"      "FeatureID"     "Type"
#  [6] "Location"      "Start"         "Stop"          "Strand"        "Function"
# [11] "Aliases"       "FigFam"        "EvidenceCodes" "Nucleotide"    "Protein"
# [16] "Length"        "EndOfContig"

# FeatureID is specific to each genome, cannot use
# Function is not unique because of repeated use of 'hypothetical protein' as a label
# Have no unique identifier

# Quick look at size of datasets function ids across genomes
DataSizes <- tibble(Genome = character(0), Component = character(0), Rows = integer(0), Columns = integer(0))
for (C in Components) {
    for (G in Genomes$One) {
        Sub <- Annotations %>% filter(Genome == G & Component == C)
        DataSizes <- add_row(DataSizes, Genome = G, Component = C, Rows = Sub %>% nrow(), Columns = Sub %>% ncol())
    }
}

DataSizes
#    Genome Component  Rows Columns
#    <chr>  <chr>     <int>   <int>
#  1 KC_109 Chr        4843      17
#  2 T      Chr        4839      17
#  3 TM     Chr        4839      17
#  4 TMu    Chr        4840      17
#  5 KC_109 Plasmid_1   145      17
#  6 T      Plasmid_1   142      17
#  7 TM     Plasmid_1   148      17
#  8 TMu    Plasmid_1   142      17
#  9 KC_109 Plasmid_2   155      17
# 10 T      Plasmid_2   147      17
# 11 TM     Plasmid_2   146      17
# 12 TMu    Plasmid_2   144      17
# 13 KC_109 Plasmid_3     3      17
# 14 T      Plasmid_3     3      17
# 15 TM     Plasmid_3     3      17
# 16 TMu    Plasmid_3     6      17
# 17 KC_109 Plasmid_4     6      17
# 18 T      Plasmid_4     6      17
# 19 TM     Plasmid_4     4      17
# 20 TMu    Plasmid_4     5      17

DataSizes %<>% select(-c(Columns)) %>% rename(CDSCount = Rows)
FileOut <- paste0(FilePrefix,"-CDS_Counts_Genome_Component.tsv")
write_tsv(DataSizes, file.path(DNAAnnotDirOut, FileOut), col_names = TRUE, na = "")

Annotations %>% filter(Genome == "KC_109", Component == "Chr") %>% nrow()

Annotations$Function %>% length() # 20566
Annotations$Function %>% unique() %>% length() # 4064
# Clearly not usable to compare genones, will need to blast

rm(n, Sub, DataSizes, G, C, Name, FileOut, Component)

##################################
## Get CDS vs CDS DNA blast data ##
##################################
# These are the results of blasting Plasmid_1 CDSs (DNA) for each genome against each other
# Will load, extract the best hits and then delete

# Set file names
DNABlastFiles <- tibble(Comparison = character(0), Query = character(0), Subject = character(0), Component = character(0), File = character(0))
for (Component in Components) {
    for (n in 1:NrGenomes) {
        for (m in 1:NrGenomes) {

            if (Genomes$One[n] == Genomes$One[m]) {next} # no comparisons of self to self

            # Get file name eg Chr-CDS-TM-vs-KC_109.tsv
            File <- paste0(Component, "-CDS_DNA-", Genomes$One[n], "-vs-", Genomes$One[m], ".tsv")
            DNABlastFiles %<>% add_row(Comparison = paste0(Genomes$One[n], "_vs_", Genomes$One[m]), Query = Genomes$One[n],
                                       Subject = Genomes$One[m], Component = Component, File = File)
        }
    }
}

DNABlastFiles %>% dim() # 60  5
DNABlastFiles %>% head()
rm(Component, n, m, File)

# Get blast results - note use of double brackets to extract name for use
DNABlasts <- list() # leaving as a list
BlastColNames <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sstrand")

# Load blast data
for (n in 1:nrow(DNABlastFiles)) {
    Name <- paste(DNABlastFiles$Comparison[n], DNABlastFiles$Component[n], sep = "-")

    print(Name)
    DNABlasts[[Name]] <- read_tsv(file.path(DNABlastDirIn, DNABlastFiles$Component[n], DNABlastFiles$File[n]), col_names = FALSE, skip = 0) %>% as_data_frame()
    colnames(DNABlasts[[Name]]) <- BlastColNames
}
DNABlasts %>% length() # 60

# Modify blast data, adding extra data
for (n in 1:nrow(DNABlastFiles)) {
    #Name <- DNABlastFiles$Comparison[n]
    Name <- DNABlasts[n] %>% names()
    Query <- DNABlastFiles$Query[n]
    Subject <- DNABlastFiles$Subject[n]
    Component <- DNABlastFiles$Component[n]
    QAnnot <- Annotations %>% filter(Genome == Query & Component == Component)
    SAnnot <- Annotations %>% filter(Genome == Subject & Component == Component)

    # Remove repeats
    DNABlasts[[Name]] %<>% filter(!grepl("repeat", qseqid) & !grepl("repeat", sseqid))

    # Remove "gnl\\|" from Blast sseqid columns, so ids match those in annotations
    DNABlasts[[Name]]$sseqid %<>% gsub("gnl\\|", "", .)

    # Add number of qseqid and sseqid for sorting
    DNABlasts[[Name]]$qseqidnr <- DNABlasts[[Name]]$qseqid %>% gsub("fig.+peg.", "", .) %>% as.integer()
    DNABlasts[[Name]]$sseqidnr <- DNABlasts[[Name]]$sseqid %>% gsub("fig.+peg.", "", .) %>% as.integer()

    # add function for query, subject
    DNABlasts[[Name]] %<>% left_join(QAnnot %>% select(FeatureID, Function), by = c("qseqid" = "FeatureID"))
    DNABlasts[[Name]] %<>% dplyr::rename(qfunction = Function)
    DNABlasts[[Name]] %<>% left_join(SAnnot %>% select(FeatureID, Function), by = c("sseqid" = "FeatureID"))
    DNABlasts[[Name]] %<>% dplyr::rename(sfunction = Function)

    # Add query and subject match lengths
    DNABlasts[[Name]] %<>% mutate(qlength = abs(qend - qstart)+1)
    DNABlasts[[Name]] %<>% mutate(slength = abs(send - sstart)+1)

    # Add query and subject full lengths
    DNABlasts[[Name]] %<>% left_join(QAnnot %>% select(FeatureID, Length), by = c("qseqid" = "FeatureID"))
    DNABlasts[[Name]] %<>% dplyr::rename(qfulllength = Length)
    DNABlasts[[Name]] %<>% left_join(SAnnot %>% select(FeatureID, Length), by = c("sseqid" = "FeatureID"))
    DNABlasts[[Name]] %<>% dplyr::rename(sfulllength = Length)

    # Add query and subject fractions length/fulllength
    DNABlasts[[Name]] %<>% mutate(qfraction = qlength/qfulllength)
    DNABlasts[[Name]] %<>% mutate(sfraction = slength/sfulllength)

    # Add query and subject genomes, component, move to front
    DNABlasts[[Name]] %<>% mutate(qgenome = Query)
    DNABlasts[[Name]] %<>% mutate(sgenome = Subject)
    DNABlasts[[Name]] %<>% mutate(component = Component)

    NrCols <- DNABlasts[[Name]] %>% ncol()
    DNABlasts[[Name]] %<>% select(qgenome, sgenome, component, everything())
}
DNABlasts %>% length() # 60

# Sizes of datasets
DataSizes <- tibble(Comparison = character(0), Rows = integer(0), Columns = integer(0), UniqueSeqIds = integer(0))
for (n in 1:length(DNABlasts)) {
    x <- DNABlasts[[names(DNABlasts[n])]]
    UniqueLength <- x$qseqid %>% unique() %>% length()
    DataSizes <- add_row(DataSizes, Comparison = names(DNABlasts[n]), Rows = x %>% nrow(), Columns = x %>% ncol(), UniqueSeqIds = UniqueLength)
}

DataSizes
#    Comparison         Rows Columns UniqueSeqIds
#    <chr>              <int>   <int>        <int>
#  1 KC_109_vs_T-Chr   140786      26         4955
#  2 KC_109_vs_TM-Chr  140834      26         4955
#  3 KC_109_vs_TMu-Chr 140842      26         4955
#  4 T_vs_KC_109-Chr   140798      26         4951
#  5 T_vs_TM-Chr       140820      26         4951
#  6 T_vs_TMu-Chr      140828      26         4951
#  7 TM_vs_KC_109-Chr  140822      26         4951
#  8 TM_vs_T-Chr       140795      26         4951
#  9 TM_vs_TMu-Chr     140851      26         4951
# 10 TMu_vs_KC_109-Chr 140802      26         4952
# # ... with 50 more rows

rm(n, DataSizes, UniqueLength, Name, Query, Subject, Component, x, QAnnot, SAnnot, NrCols)

########################################
## Get CDS cs CDS DNA blast best hits ##
########################################
# Select the best hit for each query
BestDNABlasts <- list() # will be converted to df later
for (n in 1:length(DNABlasts)) {
    Name <- DNABlasts[n] %>% names() %>% .[1] # get name of comparison

    # Get best unique matches for each query and subject id
    BestDNABlasts[[Name]] <- DNABlasts[[Name]] %>% FindBestBlastHits()
}

# Sizes of datasets
DataSizes <- tibble(Comparison = character(0), Rows = integer(0), Columns = integer(0), UniqueQseqIds = integer(0), UniqueSseqIds = integer(0))
for (n in 1:length(BestDNABlasts)) {
    Name <- BestDNABlasts[n] %>% names() %>% .[1]
    qseqidUniqueLength <- BestDNABlasts[[Name]]$qseqid %>% uniquelength()
    sseqidUniqueLength <- BestDNABlasts[[Name]]$sseqid %>% uniquelength()
    DataSizes <- add_row(DataSizes, Comparison = Name, Rows = BestDNABlasts[[Name]] %>% nrow(), Columns = BestDNABlasts[[Name]] %>% ncol(),
                         UniqueQseqIds = qseqidUniqueLength, UniqueSseqIds = sseqidUniqueLength)
}

DataSizes # if nr rows, unique qseqids and unique sseqids should have the same value. If not, check code
#    Comparison         Rows Columns UniqueQseqIds UniqueSseqIds
#    <chr>             <int>   <int>         <int>         <int>
#  1 KC_109_vs_T-Chr    4825      27          4825          4825
#  2 KC_109_vs_TM-Chr   4815      27          4815          4815
#  3 KC_109_vs_TMu-Chr  4829      27          4829          4829
#  4 T_vs_KC_109-Chr    4827      27          4827          4827
#  5 T_vs_TM-Chr        4825      27          4825          4825
#  6 T_vs_TMu-Chr       4835      27          4835          4835
#  7 TM_vs_KC_109-Chr   4817      27          4817          4817
#  8 TM_vs_T-Chr        4825      27          4825          4825
#  9 TM_vs_TMu-Chr      4827      27          4827          4827
# 10 TMu_vs_KC_109-Chr  4829      27          4829          4829
# # ... with 50 more rows

(DataSizes$Rows == DataSizes$UniqueQseqIds) %>% summary()
# Mode    TRUE    NA's
# logical      60       0
(DataSizes$Rows == DataSizes$UniqueSseqIds) %>% summary()
# Mode    TRUE    NA's
# logical      60       0

rm(n,Name, qseqidUniqueLength, sseqidUniqueLength, DataSizes)

# Make one df from best hits list
Name <- names(BestDNABlasts[1])
df <- matrix(nrow = 0, ncol = BestDNABlasts[[Name]] %>% ncol()) %>% as_data_frame()
colnames(df) <- colnames(BestDNABlasts[[Name]])

for (n in 1:length(BestDNABlasts)) {
    Name <- BestDNABlasts[n] %>% names()

    # Add to df
    df <- bind_rows(df, BestDNABlasts[[Name]])
}

BestDNABlasts <- df
BestDNABlasts %>% dim() # 61296    27
Annotations %>% dim() # 20566  17

rm(df, Name, n)

###########################################
## Compare DNA blasts across genomes ##
###########################################
# Compare number of ids in Annotations and in BestBlast, and collect missing genes
# Store counts of genes in/not in Annotations and BestBlast
IDCounts <-  tibble(Comparison = character(0), Component = character(0),
                    Query = character(0), QueryCount = integer(0), CompQueryCount = integer(0),
                    AnnotInQuery = integer(0), AnnotNinQuery = integer(0),
                    Subject = character(0), SubjectCount = integer(0), CompSubjectCount = integer(0),
                    AnnotInSubject = integer(0), AnnotNinSubject = integer(0))

# Lists of query genes missing from blasts
Names <- c("Comparison", "Query", "Subject", names(Annotations))
MissingCDS <-  matrix(nrow = 0, ncol = length(Names)) %>% as_data_frame() # list of dfs for each comparison
names(MissingCDS) <- Names

# Getting data using the Blast Files
for (n in 1:nrow(DNABlastFiles)) {
    Name <- DNABlastFiles$Comparison[n]
    Query <- DNABlastFiles$Query[n]
    Subject <- DNABlastFiles$Subject[n]
    Component <- DNABlastFiles$Component[n]

    QueryFeatureIDs <- Annotations %>% filter(Genome == DNABlastFiles$Query[n] & Component == DNABlastFiles$Component[n]) %>% .$FeatureID
    SubjectFeatureIDs <- Annotations %>% filter(Genome == DNABlastFiles$Subject[n]  & Component == DNABlastFiles$Component[n]) %>% .$FeatureID

    Comps <- BestDNABlasts %>% filter(qgenome == Query, sgenome == Subject, component == Component)

    CompQueryIDs <- Comps$qseqid %>% unique()
    CompSubjectIDs <- Comps$sseqid %>% unique()

    AnnotInQuery <- (QueryFeatureIDs %in% CompQueryIDs) %>% sum(na.rm=TRUE)
    AnnotNinQuery <- (QueryFeatureIDs %nin% CompQueryIDs) %>% sum(na.rm=TRUE)
    AnnotInSubject <- (SubjectFeatureIDs %in% CompSubjectIDs) %>% sum(na.rm=TRUE)
    AnnotNinSubject <- (SubjectFeatureIDs %nin% CompSubjectIDs) %>% sum(na.rm=TRUE)

    IDCounts %<>% add_row(Comparison = Name, Component = Component,
                          Query = Query, QueryCount = uniquelength(QueryFeatureIDs), CompQueryCount = length(CompQueryIDs),
                          Subject = Subject, SubjectCount = uniquelength(SubjectFeatureIDs), CompSubjectCount = length(CompSubjectIDs),
                          AnnotInQuery = AnnotInQuery, AnnotNinQuery = AnnotNinQuery,
                          AnnotInSubject = AnnotInSubject, AnnotNinSubject = AnnotNinSubject)

    if(AnnotNinQuery > 0) {
        t <- (QueryFeatureIDs %nin% CompQueryIDs) %>% which()
        x <- Annotations %>% filter(Genome == Query & Component == Component) %>% .[t,]
        y <- tibble(Comparison = c(rep(Name,nrow(x))), Query = c(rep(Query,nrow(x))), Subject = c(rep(Subject,nrow(x))))
        x <- cbind(y, x) %>% as_data_frame()
        MissingCDS %<>% bind_rows(x)
    }
}

IDCounts %<>% arrange(Component, Subject)
IDCounts %>% dim() # 60 12

# Save CDS counts
FileOut <- paste0(FilePrefix, "-CDS_vs_CDS-Counts.tsv")
write_tsv(IDCounts, file.path(DNAAnnotDirOut, FileOut), col_names = TRUE, na = "")

# Save missing CDS
MissingCDS %<>% select(-Genome)
FileOut <- paste0(FilePrefix,"-CDS_vs_CDS-Missing.tsv")
write_tsv(MissingCDS, file.path(DNAAnnotDirOut, FileOut), col_names = TRUE, na = "")

# Make a feature table linking ids across genomes
FeatureIDTableList <- list()
x<- FeatureIDTableList$T
for (g in 1:NrGenomes){

    Query <- Genomes[g,1] %>% .[[1]]

    FeatureIDTable <- tibble(Component = Annotations %>% filter (Genome == Query) %>% .$Component,
                             Name = Annotations %>% filter (Genome == Query) %>% .$FeatureID)
    colnames(FeatureIDTable) <- c("Component", Query)
    #FeatureIDTable$Component %>% as.factor() %>% summary()

    for (g2 in 2:ncol(Genomes)){
        Subject <- Genomes[g,g2] %>% .[[1]]

        if (Query == Subject) {next} # do not do self to self

        # Join Query/Subject seqids from BestDNABlasts
        Comps <- BestDNABlasts %>% filter(qgenome == Query, sgenome == Subject) %>% dplyr::select(component, qseqid, sseqid)
        if (nrow(Comps) > 0) {
            FeatureIDTable %<>% left_join(Comps, by = c("Component" = "component", setNames("qseqid", Query)))
            colnames(FeatureIDTable)[ncol(FeatureIDTable)] <- Subject
        }
    }
    FeatureIDTableList[[Query]] <- FeatureIDTable
    #FeatureIDTable$Component %>% as.factor() %>% summary()

}

# Merge feature tables
FeatureIDTable <- FeatureIDTableList[[Genomes$One[1]]] # get first table
for (n in 2:NrGenomes) {
    FeatureIDTable %<>% union(FeatureIDTableList[[Genomes$One[n]]] %>% select(Component, Genomes$One))
    #print(FeatureIDTable %>% dim())
}

FeatureIDTable %<>% select(Component, eval(Genomes$One))
FeatureIDTable %<>% arrange(Component, eval(as.name(Genomes$One[1])))
FeatureIDTable %>% dim() # 5209

# Sort feature table by component and then each genome feature if in turn
for (G in Genomes$One) {
    Name <- paste0(G,"_Nr")
    #FeatureIDTable %<>% mutate(eval(as.name(Name)) = gsub("fig.+peg.", "", .) as.integer)
    FeatureIDTable %<>% mutate(!!(Name) := gsub("fig.+peg.", "", eval(as.name(G))) %>% as.integer())
}
SortNames <- c("Component", paste0(Genomes$One,"_Nr"))
FeatureIDTable %<>% .[do.call(order, .[, SortNames]), ]
FeatureIDTable <- FeatureIDTable[,-((ncol(FeatureIDTable)-3):ncol(FeatureIDTable))]

# Resolve dups - NOT USING FOR NOW
# Unresolved <- FeatureIDTable %>% filter(Component == "XXX") # generate empty df
# for (g in 1:NrGenomes){
#     Genome <- Genomes$One[g]
#
#     # Find duplicate entries in Feature table for this genome
#     Dups <- FeatureIDTable %>% group_by(eval(as.name(Genome))) %>% # creates a new col, dont know why
#         filter(!is.na(eval(as.name(Genome)))) %>%
#         filter(n()>1) %>%
#         arrange(eval(as.name(Genome))) %>%
#         ungroup() %>%
#         select(-(ncol(.))) # remove the last col, created by group_by command
#     print(Genome)
#     print(Dups)
#
#     # Remove clear duplicated records - this may remove some true dups?
#     Remove <- FeatureIDTable %>% filter(Component == "XXX") # generate empty df
#     if(nrow(Dups) > 0) {
#         FeatureIDs <- Dups[,Genome] %>% unique() %>% unlist() %>% as.vector()
#
#         for (FeatureID in FeatureIDs) {
#             Dups2 <- Dups %>% filter(eval(as.name(Genome)) == FeatureID)
#             # Test 1: if Feature IDs the same, keep row with most Feature IDs
#             for (n in 1:(nrow(Dups2)-1)) {
#                 for (m in (n+1):nrow(Dups2)) {
#                     Match <- (Dups2[n,] == Dups2[m,]) %>% .[1,2:ncol(Dups)] %>% summary()
#                     if ("FALSE" %nin% names(Match)) { # a value does not match, have to keep both
#                         NotNAs <- rowSums(!is.na(Dups2[c(n,m),]))
#                         if (NotNAs[1] > NotNAs[2]) {
#                             Remove %<>% rbind(Dups2[m,])
#                         } else {
#                             Remove %<>% rbind(Dups2[n,])
#                         }
#                     }
#                 }
#             }
#             # Test 2: what to do with dups where other genomes have different Feature IDs?
#         }
#         if(nrow(Remove) > 0) {
#             print(paste0("Removing ", nrow(Remove), " rows from data based on duplicate entries for genome ", Genome))
#             FeatureIDTable %<>% setdiff(Remove)
#             Dups %<>% setdiff(Remove) # not sure this works
#         }
#         Unresolved %<>% rbind(Dups)
#     }
# }
#
FeatureIDTable %>% dim()
FeatureIDTable %>% unique() %>% dim()
# Unresolved %<>% unique()
# Unresolved %>% dim()


rm(n, g, g2, t, x, y, AnnotInQuery, AnnotInSubject,AnnotNinQuery, AnnotNinSubject, Comps, CompQueryIDs, CompSubjectIDs,
   FileOut, Name, Names, Query, QueryFeatureIDs, Subject, SubjectFeatureIDs, FeatureIDTableList, Component)
rm(Genome, SortNames, G, Name)

##################################
## Create and add unique CDS id ##
##################################
# Create unique id for each gene
t <- nrow(FeatureIDTable)
c <- t %>% as.character() %>% nchar()
CDSIDs <- vector(mode = "character", t)
for (n in 1:t) {
    #print(n)
    m <- c - (n %>% as.character() %>% nchar())
    if (m > 0) {
        Pad <- paste0(rep("0", m), collapse = "")
    } else {
        Pad <- ""
    }
    CDSIDs[n] <- paste0("cdsid", Pad, n)
}
CDSIDs %<>% as_data_frame()
colnames(CDSIDs) <- c("cdsid")
CDSIDs %>% head()
CDSIDs %>% tail()

# Add CDS IDs to feature table
FeatureIDTable <- cbind(CDSIDs,FeatureIDTable) %>% as_data_frame()

# Add CDS IDs to annotations
x <- melt(FeatureIDTable, id = c("cdsid")) %>% filter(!is.na(value)) %>% as_data_frame() # easier to join
colnames(x) <- c("CDSID", "Genome", "FeatureID")
x %<>% dplyr::select(-Genome)

Annotations %<>% left_join(x, by = "FeatureID")

# Add CDS IDs to blast data
colnames(x) <- c("cdsid", "qseqid") # redo col names so they match
for (n in 1:length(DNABlasts)) {
    Name <- DNABlasts[n] %>% names()
    DNABlasts[[Name]] <- left_join(DNABlasts[[Name]], x, by = "qseqid")
    DNABlasts[[Name]] %<>% select(cdsid, everything())
}

BestDNABlasts <- left_join(BestDNABlasts, x, by =  "qseqid")
BestDNABlasts %<>% select(cdsid, everything())

rm(n, x, Name)
rm(c, m, Pad, t)

######################
## Presence/Absence ##
######################
# Get patterns of presence/absence across genomes
FeaturePresenceTable <- FeatureIDTable %>% select(-c(cdsid, Component)) %>% as.matrix(nrow = nrow(.), ncol = ncol(.))
FeaturePresenceTable %>% dim() # 5207    4
FeaturePresenceTable[!is.na(FeaturePresenceTable)] <- TRUE # presence/absence of an id to T/F
FeaturePresenceTable[is.na(FeaturePresenceTable)] <- FALSE
FeaturePresenceTable %>% dim()
FeaturePresenceTable %>% summary()
#    KC_109         T            TM          TMu
#  FALSE:  51   FALSE:  67   FALSE:  65   FALSE:  66
#  TRUE :5156   TRUE :5140   TRUE :5142   TRUE :5141
class(FeaturePresenceTable) <- "logical" # change to logical
FeaturePresenceTable %<>% as_data_frame()

# Create a string by joining each row of truth values, so can use to summarise
String <- unite(FeaturePresenceTable, sep = " ", remove = TRUE)
colnames(String) <- c("String")
FeaturePresenceTable <- cbind(FeatureIDTable[,1], FeaturePresenceTable, String) %>% as_data_frame() # join cdids, truth table and truth string

x <- FeaturePresenceTable$String %>% as.factor() %>% summary() # get a summary
x <- cbind(names(x), x %>% as.vector()) %>% as_data_frame() # convert names and values to df
colnames(x) <- c("Present", "Counts") # rename the cols
x$Counts %<>% as.integer()
x %<>% arrange(desc(Present)) # sort by the truth string
#    Present                Counts
#    <chr>                   <int>
#  1 TRUE TRUE TRUE TRUE      5086
#  2 TRUE TRUE TRUE FALSE        3
#  3 TRUE TRUE FALSE TRUE       20
#  4 TRUE TRUE FALSE FALSE       5
#  5 TRUE FALSE TRUE TRUE        4
#  6 TRUE FALSE TRUE FALSE       9
#  7 TRUE FALSE FALSE TRUE       2
#  8 TRUE FALSE FALSE FALSE     27
#  9 FALSE TRUE TRUE TRUE       14
# 10 FALSE TRUE TRUE FALSE       5
# 11 FALSE TRUE FALSE TRUE       6
# 12 FALSE TRUE FALSE FALSE      1
# 13 FALSE FALSE TRUE TRUE       5
# 14 FALSE FALSE TRUE FALSE     16
# 15 FALSE FALSE FALSE TRUE      4

# Save feature presence table summary
x <- str_split_fixed(x$Present, " ", 4) %>% cbind(x$Counts) %>% as_data_frame()
colnames(x) <- c(Genomes$One, "Counts")

AnnotDirOut <- "/Users/rtearle/Documents/Roseworthy_Projects/6.Students/Nitish_Salmonella/Apr_2018/Out/R_DNA_ANNOT"
FileOut <- paste0(FilePrefix,"-Matching_CDS_IDs-Truth_Table_Summary.tsv")
write_tsv(x, file.path(AnnotDirOut, FileOut), col_names = TRUE, na = "")

rm(x, String, FileOut)

##############################################
## Compare matching CDS sizes and sequences ##
##############################################
# Create vectors of types char, int, real
VChar <- rep(NA_character_, nrow(FeatureIDTable))
VInt <- rep(NA_integer_, nrow(FeatureIDTable))
VReal <- rep(NA_real_, nrow(FeatureIDTable))

# Create empty table for comparisons
CDSComparisons <- tibble(CDSID = VChar, Component = VChar, Function = VChar)
for (T in c("_ID", "_AA_Length", "_NT_Length", "_AA", "_NT")) {
    for (G in Genomes$One) {
        x <- paste0(G, T)
        if (T == "_ID" | T == "_DNA_Length" | T == "_AA_Length") {
            CDSComparisons %<>% add_column(!!(x) := VInt)
        } else if (T == "_AA") {
            CDSComparisons %<>% add_column(!!(x) := VChar)
        } else if (T == "_NT") {
            CDSComparisons %<>% add_column(!!(x) := VChar)
        }
    }
}

CDSComparisons %<>% add_column(Mean_AA := VReal)
CDSComparisons %<>% add_column(SD_AA := VReal)
CDSComparisons %<>% add_column(Mean_NT := VReal)
CDSComparisons %<>% add_column(SD_NT := VReal)

# Populate the comparisons table
for (n in 1:nrow(FeatureIDTable)) {
    #FeatureIDTable[n,]
    #CDSComparisons[n,]

    # Add CDS ID
    CDSComparisons$CDSID[n] <- CDS <- FeatureIDTable$cdsid[n]
    if (is.na(CDS)) {next}

    # Add component
    CDSComparisons$Component[n] <- Component <- FeatureIDTable$Component[n]

   # Loop through the genomes to extract data from annotations and feature table
    for (m in 2:(NrGenomes+1)) {
        Name <- colnames(FeatureIDTable)[m+1] # eg KC_109
        Name2 <- paste0(Name, "_ID") # eg KC_109_ID

        # Add feature ID from feature table
        CDSComparisons[n, Name2] <- FID <- FeatureIDTable[n,m+1] %>% .[[1]]

        if (is.na(FID)) {next}

        # Get annotation record for this feature id
        Annot <- Annotations %>% dplyr::filter(Genome == Name, FeatureID == FID)

        # Add protein from comparisons
        Protein <- paste0(Name, "_AA") # eg KC_109_AA
        CDSComparisons[n, Protein] <- Annot$Protein[[1]]

        # Add DNA from comparisons
        DNA <- paste0(Name, "_NT") # eg KC_109_AA
        CDSComparisons[n, DNA] <- Annot$Nucleotide[[1]]

        # Add length from comparisons
        Length <- paste0(Name, "_AA_Length") # eg KC_109_Length
        CDSComparisons[n, Length] <- CDSComparisons[n, Protein] %>% nchar() #%>% as.integer()
        Length <- paste0(Name, "_NT_Length") # eg KC_109_Length
        CDSComparisons[n, Length] <- CDSComparisons[n, DNA] %>% nchar() #%>% as.integer()

    }
    # Add function from annotations
    CDSComparisons$Function[n] <- Annot$Function[[1]]
}

# Add means and sds
CDSComparisons$Mean_AA <- CDSComparisons %>% select(ends_with("_AA_Length")) %>% rowMeans(na.rm = TRUE)
CDSComparisons$SD_AA <- apply(CDSComparisons %>% select(ends_with("_AA_Length")), 1, sd, (na.rm = TRUE))
CDSComparisons$Mean_NT <- CDSComparisons %>% select(ends_with("_NT_Length")) %>% rowMeans(na.rm = TRUE)
CDSComparisons$SD_NT <- apply(CDSComparisons %>% select(ends_with("_NT_Length")), 1, sd, (na.rm = TRUE))

CDSComparisons %>% dim()

rm(VChar, VInt, VReal, T, G, x, n, m, CDS, Name, Name2, Annot, Length, FID, DNA)
rm(na.rm, Protein)
rm(Component)

######################################
## Get CDS vs Genome DNA blast data ##
######################################
# These are the results of blasting CDSs (DNA) for each genome against the other genomes
# Will load, extract the best hits and then delete

# Set file names
DNABlastFiles2 <- tibble(Comparison = character(0), Query = character(0), Subject = character(0), Component = character(0), File = character(0))
for (n in 1:length(Genomes$One)) {
    for (m in 1:length(Genomes$One)) {
        for (Component in Components) {
            if (Genomes$One[n] == Genomes$One[m]) {next} # no comparisons of self to self

            # Get file name eg Plasmid_1-CDS-T-vs-KC_109-Genome.tsv
            File <- paste0(Component, "-CDS_DNA-", Genomes$One[n], "-vs-", Genomes$One[m], "-Genome.tsv")
            DNABlastFiles2 %<>% add_row(Comparison = paste0(Genomes$One[n], "_vs_", Genomes$One[m], "-", Component), Query = Genomes$One[n], Subject = Genomes$One[m], Component = Component, File = File)
        }
    }
}

DNABlastFiles2 %>% dim() # 60  5
DNABlastFiles2$File

rm(n, m, File)

# Get blast results - note use of double brackets to extract name for use
DNABlasts2 <- list() # leaving as a list
BlastColNames <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sstrand")

# Load blast data
for (n in 1:nrow(DNABlastFiles2)) {
    Name <- DNABlastFiles2$Comparison[n]
    print(Name)
    DNABlasts2[[Name]] <- read_tsv(file.path(DNABlastDirIn2, DNABlastFiles2$Component[n], DNABlastFiles2$File[n]), col_names = FALSE, skip = 0) %>% as_data_frame()
    colnames(DNABlasts2[[Name]]) <- BlastColNames
}

# Melt feature table, makes it easier to join CDS IDs to blast results
x <- melt(FeatureIDTable, id = c("cdsid")) %>% filter(!is.na(value)) %>% as_data_frame()
colnames(x) <- c("cdsid", "Genome", "qseqid")
x %<>% dplyr::select(-Genome)

# Modify blast data, adding extra data
for (n in 1:nrow(DNABlastFiles2)) {

    Name <- DNABlastFiles2$Comparison[n]
    Query <- DNABlastFiles2$Query[n]
    Subject <- DNABlastFiles2$Subject[n]
    Component <- DNABlastFiles2$Component[n]
    QAnnot <- Annotations %>% filter(Genome == Query)
    SAnnot <- Annotations %>% filter(Genome == Subject)

    # Remove repeats and rna
    DNABlasts2[[Name]] %<>% filter(!grepl("repeat", qseqid))
    DNABlasts2[[Name]] %<>% filter(!grepl("rna", qseqid))
    DNABlasts2[[Name]] %<>% filter(!grepl("crispr", qseqid))


    # Remove "gnl\\|" from Blast sseqid columns, so ids match those in annotations
    DNABlasts2[[Name]]$sseqid %<>% gsub("gnl\\|", "", .)

    # Add number of qseqid and sseqid for sorting
    DNABlasts2[[Name]]$qseqidnr <- DNABlasts2[[Name]]$qseqid %>% gsub("fig.+peg.", "", .) %>% as.integer()
    DNABlasts2[[Name]]$sseqidnr <- DNABlasts2[[Name]]$sseqid %>% gsub("fig.+peg.", "", .) %>% as.integer()

    # add function for query and subject
    DNABlasts2[[Name]] %<>% left_join(QAnnot %>% select(FeatureID, Function), by = c("qseqid" = "FeatureID"))
    DNABlasts2[[Name]] %<>% dplyr::rename(qfunction = Function)
    DNABlasts2[[Name]] %<>% left_join(SAnnot %>% select(FeatureID, Function), by = c("sseqid" = "FeatureID"))
    DNABlasts2[[Name]] %<>% dplyr::rename(sfunction = Function)

    # Add query and subject match lengths
    DNABlasts2[[Name]] %<>% mutate(qlength = abs(qend - qstart)+1)
    DNABlasts2[[Name]] %<>% mutate(slength = abs(send - sstart)+1)

    # Add query and subject full lengths
    DNABlasts2[[Name]] %<>% left_join(QAnnot %>% select(FeatureID, Length), by = c("qseqid" = "FeatureID"))
    DNABlasts2[[Name]] %<>% dplyr::rename(qfulllength = Length)
    DNABlasts2[[Name]] %<>% left_join(SAnnot %>% select(FeatureID, Length), by = c("sseqid" = "FeatureID"))
    DNABlasts2[[Name]] %<>% dplyr::rename(sfulllength = Length)

    # Add query and subject fractions length/fulllength
    DNABlasts2[[Name]] %<>% mutate(qfraction = qlength/qfulllength)
    DNABlasts2[[Name]] %<>% mutate(sfraction = slength/sfulllength)

    # Add cdsid, query and subject genomes, component, move to front
    DNABlasts2[[Name]] <- left_join(DNABlasts2[[Name]], x, by = "qseqid")
    DNABlasts2[[Name]] %<>% mutate(qgenome = Query)
    DNABlasts2[[Name]] %<>% mutate(sgenome = Subject)
    DNABlasts2[[Name]] %<>% mutate(component = Component)

    DNABlasts2[[Name]] %<>% select(cdsid, qgenome, sgenome, component, everything())

}
DNABlasts2 %>% length() # 60

# Sizes of datasets
DataSizes <- tibble(Genome = character(0), Rows = integer(0), Columns = integer(0), UniqueSeqIds = integer(0))
for (n in 1:length(DNABlasts2)) {
    x <- DNABlasts2[[names(DNABlasts2[n])]]
    UniqueLength <- x$qseqid %>% unique() %>% length()
    DataSizes <- add_row(DataSizes, Genome = names(DNABlasts2[n]), Rows = x %>% nrow(), Columns = x %>% ncol(), UniqueSeqIds = UniqueLength)
}
DataSizes
#    Genome                   Rows Columns UniqueSeqIds
#    <chr>                   <int>   <int>        <int>
#  1 KC_109_vs_T-Chr        134354      27         4843
#  2 KC_109_vs_T-Plasmid_1    1888      27          145
#  3 KC_109_vs_T-Plasmid_2    1808      27          155
#  4 KC_109_vs_T-Plasmid_3      16      27            3
#  5 KC_109_vs_T-Plasmid_4      16      27            6
#  6 KC_109_vs_TM-Chr       134360      27         4843
#  7 KC_109_vs_TM-Plasmid_1   1885      27          145
#  8 KC_109_vs_TM-Plasmid_2   1712      27          155
#  9 KC_109_vs_TM-Plasmid_3     17      27            3
# 10 KC_109_vs_TM-Plasmid_4     16      27            6
# # ... with 50 more rows

# Annotations %>% filter(Genome == Genomes$One[1], Component == "Chr") %>% .$FeatureID %>% length() # 4844
# Annotations %>% filter(Genome == Genomes$One[1], Component == "Chr") %>% .$FeatureID %>% uniquelength() # 4843
# x <- DNABlasts2[[names(DNABlasts2[1])]]
# y1 <- Annotations %>% filter(Genome == Genomes$One[1], Component == "Chr") %>% .$FeatureID %>% unique()
# y2 <- DNABlasts2[[names(DNABlasts2[1])]] %>% .$qseqid %>% unique()
#
# y1 %in% y2 %>% summary()
# #    Mode    TRUE    NA's
# # logical    4843       0
#
# y2 %in% y1 %>% summary()
# #    Mode   TRUE    NA's
# # logical   4843       0

rm(n, DataSizes, UniqueLength, Name, Query, Subject, x, QAnnot, SAnnot)
rm(Component)

###########################################
## Get CDS vs genome DNA blast best hits ##
###########################################
# Select the best hit for each query
BestDNABlasts2 <- list() # will be converted to df later
for (n in 1:length(DNABlasts2)) {
    Name <- DNABlasts2[n] %>% names() %>% .[1] # get name of comparison

    # Get best unique matches for each query and subject id
    BestDNABlasts2[[Name]] <- DNABlasts2[[Name]] %>% FindBestBlastHits2()
}

# Sizes of datasets
DataSizes <- tibble(Comparison = character(0), Rows = integer(0), Columns = integer(0), UniqueQseqIds = integer(0), UniqueSseqIds = integer(0))
for (n in 1:length(BestDNABlasts2)) {
    Name <- BestDNABlasts2[n] %>% names() %>% .[1]
    qseqidUniqueLength <- BestDNABlasts2[[Name]]$qseqid %>% uniquelength()
    sseqidUniqueLength <- BestDNABlasts2[[Name]]$sseqid %>% uniquelength()
    DataSizes <- add_row(DataSizes, Comparison = Name, Rows = BestDNABlasts2[[Name]] %>% nrow(), Columns = BestDNABlasts2[[Name]] %>% ncol(),
                         UniqueQseqIds = qseqidUniqueLength, UniqueSseqIds = sseqidUniqueLength)
}
DataSizes # nr rows, unique qseqids should have the same value
#    Comparison              Rows Columns UniqueQseqIds UniqueSseqIds
#    <chr>                  <int>   <int>         <int>         <int>
#  1 KC_109_vs_T-Chr         4842      27          4842            79
#  2 KC_109_vs_T-Plasmid_1    145      27           145             3
#  3 KC_109_vs_T-Plasmid_2    148      27           148             7
#  4 KC_109_vs_T-Plasmid_3      3      27             3             1
#  5 KC_109_vs_T-Plasmid_4      6      27             6             1
#  6 KC_109_vs_TM-Chr        4842      27          4842            79
#  7 KC_109_vs_TM-Plasmid_1   145      27           145             4
#  8 KC_109_vs_TM-Plasmid_2   145      27           145             7
#  9 KC_109_vs_TM-Plasmid_3     3      27             3             1
# 10 KC_109_vs_TM-Plasmid_4     6      27             6             1
# # ... with 50 more rows

(DataSizes$Rows - DataSizes$UniqueQseqIds == 0) %>% summary()
#    Mode    TRUE    NA's
# logical      60       0

rm(n, Name, qseqidUniqueLength, sseqidUniqueLength, DataSizes)

# Make one df from best hits list
Name <- names(BestDNABlasts2[1])
df <- matrix(nrow = 0, ncol = length(colnames(BestDNABlasts2[[Name]]))+1) %>% as_data_frame()
colnames(df) <- c("comp", colnames(BestDNABlasts2[[Name]]))

for (n in 1:length(BestDNABlasts2)) {
    Name <- BestDNABlasts2[n] %>% names()

    # create a df with Name as only col
    Name2 <- matrix(Name, ncol = 1, nrow = nrow(BestDNABlasts2[[Name]])) %>% as_data_frame()
    colnames(Name2) <- c("comp")

    # Add to front
    x <- bind_cols(Name2, BestDNABlasts2[[Name]])
    df <- bind_rows(df, x)
}

BestDNABlasts2 <- df
BestDNABlasts2 %>% dim() # 61647   29

# y1 <- Annotations %>% filter(Genome == Genomes$One[1], Component == "Chr") %>% .$FeatureID %>% unique()
# y1 %>% length() # 4843
# y2 <- BestDNABlasts2 %>% filter(comp == "KC_109_vs_T-Chr") %>% .$qseqid %>% unique()
# y2 %>% length() # 4842
#
# y1 %in% y2 %>% summary()
# #    Mode   FALSE    TRUE    NA's
# # logical       1    4842       0
#
# y2 %in% y1 %>% summary()
# #    Mode    TRUE    NA's
# # logical    4842       0
#
# y1 %nin% y2 %>% which() %>% head()
# y1[4677]

rm(n, Name, Name2, df, x)

###########################################
## Compare DNA blasts across genomes ##
###########################################
# Compare number of ids in Annotations and in BestBlast, and collect missing genes
# Store counts of genes in/not in Annotations and BestBlast
IDCounts2 <-  tibble(Comparison = character(0), Component = character(0),
                     Query = character(0), QueryCount = integer(0), CompQueryCount = integer(0),
                     AnnotInQuery = integer(0), AnnotNinQuery = integer(0),
                     Subject = character(0))

# Lists of query genes missing from blasts
Names <- c("Comparison", "Query", "Subject", names(Annotations))
MissingCDS2 <-  matrix(nrow = 0, ncol = length(Names)) %>% as_data_frame() # list of dfs for each comparison
names(MissingCDS2) <- Names

for (n in 1:nrow(DNABlastFiles2)) {

    Name <- DNABlastFiles2$Comparison[n]
    Query <- DNABlastFiles2$Query[n]
    Subject <- DNABlastFiles2$Subject[n]
    Component <- DNABlastFiles2$Component[n]

    Annot <- Annotations %>% filter(Genome == Query & Component == Component)
    QueryFeatureIDs <- Annot$FeatureID
    QueryFeatureIDs %>% length()
    QueryFeatureIDs %>% unique() %>% length()
    Annot %>% tail()

    Comps <- BestDNABlasts2 %>% filter(qgenome == Query, sgenome == Subject, component == Component)
    CompQueryIDs <- Comps$qseqid %>% unique()
    CompQueryIDs %>% length()

    AnnotInQuery <- (QueryFeatureIDs %in% CompQueryIDs) %>% sum(na.rm=TRUE)
    AnnotNinQuery <- (QueryFeatureIDs %nin% CompQueryIDs) %>% sum(na.rm=TRUE)

    IDCounts2 %<>% add_row(Comparison = Name, Component = Component,
                           Query = Query, QueryCount = uniquelength(QueryFeatureIDs), CompQueryCount = length(CompQueryIDs),
                           AnnotInQuery = AnnotInQuery, AnnotNinQuery = AnnotNinQuery,
                           Subject = Subject)

    if(AnnotNinQuery > 0) {
        t <- (QueryFeatureIDs %nin% CompQueryIDs) %>% which()
        x <- Annot[t,]
        y <- tibble(Comparison = c(rep(Name,nrow(x))), Query = c(rep(Query,nrow(x))), Subject = c(rep(Subject,nrow(x))))
        x <- cbind(y, x) %>% as_data_frame()
        MissingCDS2 %<>% bind_rows(x)
    }
}
IDCounts2 %<>% arrange(Component, Subject)

# Save CDS counts
FileOut <- paste0(FilePrefix,"-CDS_vs_Genomes-Counts.tsv")
write_tsv(IDCounts2, file.path(DNAAnnotDirOut, FileOut), col_names = TRUE, na = "")

# Save missing CDS
FileOut <- paste0(FilePrefix,"-CDS_vs_Genomes-Missing.tsv")
write_tsv(MissingCDS2, file.path(DNAAnnotDirOut, FileOut), col_names = TRUE, na = "")

# Extend feature table to add genome info
FeatureIDTable2 <- FeatureIDTable
for (G in Genomes$One) {
    Name <- paste0(G, "_vs_Genomes")
    FeatureIDTable2[,Name] <- ""
}

# Add presence/absence of CDS when blasted against other whole genomes
for (n in 1:NrGenomes){

    q <- Genomes[n,1] %>% .[[1]] # CDS query genome name
    q2 <- paste0(q, "_vs_Genomes")
    f <- colnames(FeatureIDTable2) %in% q %>% which() # col pos to get CDS ID
    f2 <- colnames(FeatureIDTable2) %in% q2 %>% which() # col pos to put result

    for (m in 2:NrGenomes) {

        s <- Genomes[n,m] %>% .[[1]] # CDS subject genome name

        for (r in 1:nrow(FeatureIDTable2)) {

            # Special case, missing id
            if (is.na(FeatureIDTable2[r,f])) {
                FeatureIDTable2[r,f2] = "NA NA NA"
                next
            }

            # Have id
            x <- BestDNABlasts2 %>% filter(qseqid == FeatureIDTable2[r,f] %>% .[[1]], qgenome == q, sgenome == s)

            if (nrow(x) == 1) {
                if (FeatureIDTable2[r,f2] == "") {
                    FeatureIDTable2[r,f2] = "TRUE"
                } else {
                    FeatureIDTable2[r,f2] = paste(FeatureIDTable2[r,f2], "TRUE")
                }
            } else if (nrow(x) == 0) {
                if (FeatureIDTable2[r,f2] == "") {
                    FeatureIDTable2[r,f2] = "FALSE"
                } else {
                    FeatureIDTable2[r,f2] = paste(FeatureIDTable2[r,f2], "FALSE")
                }
            } else {
                print("Error: multiple results when matching query id to best blast results")
                print(q, q2, f, f2)
                print(s)
                print(r)
                print(x)
                stop()
            }
        }
    }
}

Cols <- colnames(FeatureIDTable2) %>% grepl("_Genome", .) %>% which()

x <- FeatureIDTable2[,Cols] %>% as.matrix() %>% as.vector() %>% as.factor() %>% table() %>% as_data_frame()
colnames(x) <- c("Presence", "Count")
x %>% arrange(desc(Presence))
#   Presence          Count
#   <chr>             <int>
# 1 TRUE TRUE TRUE    20543
# 2 TRUE TRUE FALSE       3
# 3 TRUE FALSE TRUE       3
# 4 NA NA NA            249
# 5 FALSE TRUE TRUE       3
# 6 FALSE FALSE FALSE    27  # Vast majority of CDS present in other genomes, a few not. C/f CDS more missing

rm(n, t, x, AnnotInQuery, AnnotNinQuery, Comps, CompQueryIDs, Annot,
   FileOut, Name, Names, Query, QueryFeatureIDs, Subject, f, f2)
rm(G, m, q, q2, r, s, Component, Cols)
rm(y)

#######################################################
## Find changes between AA and DNA sequences ##
#######################################################
# Compare AA sequences across genomes
# When they differ, compare DNA differences

# Df for variant data
Variants <- tibble(CDSID = character(0), Comp = character(0), Query = character(0), QueryID = character(0),
                   Subject  = character(0), SubjectID = character(0), SubjectContig = character(0),
                   EndOfContig = logical(0), ProteinVars = character(0), DNAVars = character(0))

# This takes a long time!
for (n in 1:nrow(CDSComparisons)) { # row of the feature table
    #CDSComparisons[n,]
    CDSID <- CDSComparisons$CDSID[n]
    Component <- CDSComparisons$Component[n]

    for (g in 1:NrGenomes) { # query genome
        Query <- Genomes[g,1][[1]] # first genome in col
        Qid <- CDSComparisons[n, paste0(Query, "_ID")][[1]]

        if(is.na(Qid)) {next} # do not process if no query id

        # Labels, lengths, sequences for Query
        Qlabel <- paste(CDSID, Query, sep = " ")
        Qaa <- CDSComparisons[n, paste0(Query,"_AA")][[1]]
        Qaalen <- CDSComparisons[n, paste0(Query,"_AA_Length")][[1]]
        Qnt <- CDSComparisons[n, paste0(Query,"_NT")][[1]]
        Qntlen <- CDSComparisons[n, paste0(Query,"_NT_Length")][[1]]

        for (g2 in 2:NrGenomes){ # subject genomes in turn
            Subject <- Genomes[g,g2][[1]]
            Sid <- CDSComparisons[n, paste0(Subject, "_ID")][[1]]

            # Compare protein sequences for CDS vs CDS, extracting differences
            ProteinVars <- vector()

            if(!is.na(Sid)) { # only generate protein variants if there is a subject id
                # Labels, lengths, sequences for Subject
                Slabel <- paste(CDSID, Subject, sep = " ")
                Saa <- CDSComparisons[n, paste0(Subject,"_AA")][[1]]
                Saalen <- CDSComparisons[n, paste0(Subject,"_AA_Length")][[1]]

                if( Qaa != Saa ) { # aa seqs different
                    # Find differences between query and subject protein sequences
                    Mat <- doMuscleComp(c(Qaa, Saa), c(Qlabel, Slabel), Quiet = TRUE, Type = "AA")
                    Mat2 <- as.matrix(Mat)

                    CurrType <- vector()
                    CurrSeq <- vector()
                    CurrSeq2 <- vector()

                    # Convert differences to variant strings
                    for (c in 1:ncol(Mat2)) {
                        base <- Mat2[1,c] %>% .[[1]]
                        base2 <- Mat2[2,c] %>% .[[1]]
                        if (base != Mat2[2,c]) {
                            if (base == "-") {
                                CurrType <- c(CurrType, "Del")
                            }  else if (base2 == "-") {
                                CurrType <- c(CurrType, "Ins")
                            } else {
                                CurrType <- c(CurrType, "Sub")
                            }
                            # Add to current seqs
                            CurrSeq <- c(CurrSeq, base)
                            CurrSeq2 <- c(CurrSeq2, base2)
                        } else {
                            # Process current, if any
                            if (length(CurrType) > 0) {
                                ProteinVars %<>% append(FormatProteinVariant(CurrType, CurrSeq, CurrSeq2, c-1))
                                CurrType <- vector()
                                CurrSeq <- vector()
                                CurrSeq2 <- vector()
                            } else {
                            }
                        }
                    }
                }
            }

            # Process last, if any
            if (length(ProteinVars) > 0) {
                ProteinVars %<>% paste(collapse = "|")
            } else {
                ProteinVars = ""
            }

            # DNA variants extracted by comparing CDS seq to matching seq from other genome
            Comp <- paste0(Query,"_vs_",Subject, "-", Component)
            Blast <- BestDNABlasts2 %>% filter(comp == Comp & qseqid == Qid)

            DNAVars <- vector()
            if(nrow(Blast) > 0) { # only generate protein variants if there is a blast result
                # Labels, lengths, sequences for Subject
                Sctg <- Blast$sseqid[[1]]
                Path <- file.path(CrrDirIn, Subject, paste0(tolower(Component),".crr"))
                Snt <- GetSeqFromDecodeCRR(Path, Blast$sseqid[1], Blast$sstart[1], Blast$send[1])
                Sntlen <- Snt %>% nchar()
                if (Sntlen == 0) {
                    sprintf("Cannot find sequence by position in subject genome %s matching %s query %s from genome %s", Subject, Component, CDSID, Query)
                    stop()
                }

                if( Qnt != Snt ) { # nt seqs different
                    # Find differences between query and subject nt sequences
                    Mat <- doMuscleComp(c(Qnt, Snt), c(Qlabel, Slabel), Quiet = TRUE, Type = "DNA")
                    Mat2 <- as.matrix(Mat)

                    CurrType <- vector()
                    CurrSeq <- vector()
                    CurrSeq2 <- vector()

                    # Convert differences to variant strings
                    for (c in 1:ncol(Mat2)) {
                        base <- Mat2[1,c] %>% .[[1]]
                        base2 <- Mat2[2,c] %>% .[[1]]
                        if (base != base2) {
                            if (base == "-") {
                                CurrType <- c(CurrType, "Del")
                            }  else if (base2 == "-") {
                                CurrType <- c(CurrType, "Ins")
                            } else {
                                CurrType <- c(CurrType, "Sub")
                            }
                            # Add to current seqs
                            CurrSeq <- c(CurrSeq, base)
                            CurrSeq2 <- c(CurrSeq2, base2)
                        } else {
                            # Process current, if any
                            if (length(CurrType) > 0) {
                                DNAVars %<>% append(FormatDNAVariant(CurrType, CurrSeq, CurrSeq2, c-1))
                                CurrType <- vector()
                                CurrSeq <- vector()
                                CurrSeq2 <- vector()
                            } else {
                            }
                        }
                    }
                }
            }

            # Process last, if any
            if (length(DNAVars) > 0) {
                DNAVars %<>% paste(collapse = "|")
            } else {
                DNAVars = ""
            }

            # Add end of contig
            EndOfContig <- Annotations %>% filter(FeatureID == Qid) %>% .$EndOfContig %>% .[1]

            Variants %<>% add_row(CDSID = CDSID, Comp = Comp, Query = Query, QueryID = Qid,
                                  Subject = Subject, SubjectID = Sid, SubjectContig = Sctg, EndOfContig = EndOfContig,
                                  ProteinVars = ProteinVars, DNAVars = DNAVars)
        }
    }
}

rm(n, Query, Qid, CDSID, Subject, Sid, Comp, ProteinVars, Qlabel, Slabel, Qaalen, Saalen, Qaa, Saa, Qnt, Snt, Qntlen, Sntlen,
    Mat, Mat2, CurrType, CurrSeq, CurrSeq2, DNAVars, Sctg, Path, c)
rm(EndOfContig, g, g2, base, base2, Component, first_time)

####################################
## Consolidate and save data sets ##
####################################
# Consolidating
# 1. Annotations
# 2. CDSComparisons = CDSComparisons + FeatureTable2
# 3. CDSBlastsVariants = BestDNABlasts + BestDNABlasts2 + Variants

# Save annotations
Annotations %<>% dplyr::select(CDSID, everything()) # put CDS ID first

FileOut <- paste0(FilePrefix,"-CDS_Annotations.tsv")
write_tsv(Annotations, file.path(SummaryDirOut, FileOut), col_names = TRUE, na = "")

# Save feature table plus comparisons
CDSComparisons %>% dim() # 5152   27
FeatureIDTable2 %>% dim() # 5152   10
CDSComparisons2 <- left_join(CDSComparisons, FeatureIDTable2 %>% dplyr::select(-c(one_of(Genomes$One), Component)), by = c("CDSID" = "cdsid"))
CDSComparisons2 %>% dim() # 5152   31
CDSComparisons2 %>% colnames()

FileOut <- paste0(FilePrefix,"-CDS-Multigenome_Annotations.tsv")
write_tsv(CDSComparisons2, file.path(SummaryDirOut, FileOut), col_names = TRUE, na = "")


# Subset of BestDNABlasts with col name changes
BDB <- BestDNABlasts %>% dplyr::select(-c(qseqidnr, sseqidnr, qfunction, sfunction))
colnames(BDB)
Names <- c("CDSID", "QueryGenome", "SubjectGenome", "Component", "QueryID", "SubjectID", "CDSIdentity", "CDSCompLength", "CDSMismatch", "CDSGaps",
           "CDSQueryStart", "CDSQueryEnd", "CDSSubjectStart", "CDSSubjectEnd", "CDSEValue", "CDSBitScore", "CDSStrand",
           "CDSQueryLength", "CDSSubjectLength", "CDSQueryFullLength", "CDSSubjectFullLength",
           "CDSQueryFraction", "CDSSubjectFraction", "CDSFractionProduct") # for BestDNABlasts
colnames(BDB) <- Names
BDB %<>% mutate(Comparison = paste0(QueryGenome, "_vs_", SubjectGenome, "-", Component))
BDB %<>% dplyr::select(Comparison, everything())
BDB %>% dim() # 61296    25

# Subset and split of BestDNABlasts2 with col name changes
BDB2 <- BestDNABlasts2 %>% dplyr::select(-c(qseqidnr, sseqidnr, sfunction, sfraction))
colnames(BDB2)
Names <- c("Comparison", "CDSID", "QueryGenome", "SubjectGenome", "Component", "QueryID", "SubjectContig", "GenIdentity", "GenCompLength", "GenMismatch", "GenGaps",
           "GenQueryStart", "GenQueryEnd", "GenSubjectStart", "GenSubjectEnd", "GenEValue", "GenBitScore", "GenStrand",
           "Function", "GenQueryLength", "GenSubjectLength", "GenQueryFullLength", "GenSubjectFullLength",
           "GenQueryFraction") # for BestDNABlasts2
colnames(BDB2) <- Names
BDB2 %>% dim() # 61674    24

# Variants with col name changes
colnames(Variants) <- c("CDSID", "Comparison", "QueryGenome", "QueryID", "SubjectGenome", "SubjectID", "SubjectContig", "EndOfContig", "ProteinVars", "DNAVars")

# New df with key cols
#CDSBlastsVariants <- BDB2 %>% select(Comparison, CDSID)
CDSBlastsVariants <- BDB2 %>% dplyr::select(Comparison, CDSID, QueryGenome, SubjectGenome, Component, QueryID)
CDSBlastsVariants %>% dim() # 61674    6

#BDB %>% select(Comparison, CDSID, QueryID, Component) %>% arrange(CDSID) %>% filter()

CDSBlastsVariants %<>% left_join(BDB %>% dplyr::select(-c(QueryGenome, SubjectGenome, Component)), by = c("Comparison", "CDSID", "QueryID"))
CDSBlastsVariants %>% dim() # 61674    25
colnames(CDSBlastsVariants)

CDSBlastsVariants %<>% left_join(BDB2 %>% dplyr::select(-c(QueryGenome, SubjectGenome, Component)), by = c("Comparison", "CDSID", "QueryID"))
CDSBlastsVariants %>% dim() # 61674    43
colnames(CDSBlastsVariants)

CDSBlastsVariants %<>% left_join(Variants %>% dplyr::select(-c(QueryGenome, SubjectGenome, SubjectID, SubjectContig)), by = c("Comparison", "CDSID", "QueryID"))
CDSBlastsVariants %>% dim() # 61674    46
colnames(CDSBlastsVariants)

FileOut <- paste0(FilePrefix,"-CDS_Blasts_Variants.tsv")
write_tsv(CDSBlastsVariants, file.path(SummaryDirOut, FileOut), col_names = TRUE, na = "")

rm(BDB, BDB2, FileOut, Names)
#rm(BestDNABlasts, BestDNABlasts2, Variants, FeatureIDTable)




StopTime <- Sys.time()

StartTime
StopTime
StopTime - StartTime

#######################################
## Consolidate Genome vs Genome data ##
#######################################
# Genomes were compared with minimap2
# Two datasets analysed: comparisons as paf datasets, stats

# Data stores
PAFs <- tibble(Query = character(0), QueryLength = character(0), QueryStart = character(0), QueryEnd = character(0),
               Strand = character(0), Target = character(0), TargetLength = character(0), TargetStart = character(0),
               TargetEnd = character(0), Matches = character(0), AlignLength = character(0), MapQ = character(0), Other = character(0))

PAFStats <- tibble(Types = c("Bases", "Substitutions", "Del_1bp", "Ins_1bp", "Del_2bp", "Ins_2bp", "Del_3_50bp", "Ins_3_50bp", "Del_>50bp", "Ins_>50bp"))

# Process genome-genome comparisons
for (n in 1:NrGenomes){
    Subject <- Genomes[n,1] %>% .[[1]] # CDS query genome name

    for (m in 2:NrGenomes) {
        Query <- Genomes[n,m] %>% .[[1]] # CDS subject genome name
        Comp <- paste0(Query, "_Genome-vs-", Subject, "_Genome")
        ShortComp <- paste0(Query, "-vs-", Subject)

        # T_Genome-vs-KC_109_Genome.paf
        File <- paste0(Comp, ".paf")

        PAF <- read_tsv(file.path(GGDirIn, Subject, File), col_names = FALSE)
        PAF %>% dim()
        PAF %<>% unite("X13", colnames(PAF)[13:ncol(PAF)], sep = "")
        colnames(PAF) <- c("QueryID", "QueryLength", "QueryStart", "QueryEnd", "Strand", "SubjectID", "SubjectLength", "SubjectStart", "SubjectEnd", "Matches", "AlignLength", "MapQ", "Other")

        PAF$Comp <- Comp
        PAF$Query <- Query
        PAF$Subject <- Subject

        PAF %<>% select(Comp, Query, Subject, everything())

        PAFs %<>% rbind(PAF)

        # T_Genome-vs-KC_109_Genome-stats.tsv
        File <- paste0(Comp, "-stats.tsv")

        Stats <- read_tsv(file.path(GGDirIn, Subject, File), col_names = FALSE)
        Stats$X1 %<>% gsub("^([0-9]+) .+", "\\1", .)
        colnames(Stats) = c(ShortComp)

        PAFStats %<>% cbind(Stats)
    }
}

PAFs %>% dim() # 1223   16

# Save
FileOut <- paste0(FilePrefix,"-Genome_vs_Genome.paf")
write_tsv(PAFs, file.path(GGDirIn, FileOut), col_names = TRUE, na = "")

FileOut <- paste0(FilePrefix,"-Genome_vs_Genome_PAFStats.tsv")
write_tsv(PAFStats, file.path(GGDirIn, FileOut), col_names = TRUE, na = "")



















##########################################
## Create bed files from CDS vs genomes ##
##########################################
# Create bed files for each comparison
rgb = "0,0,255"
HalveIntensity <- function(rgb) {
    rgb <- strsplit(rgb, ",") %>% unlist()
    for (n in length(rgb)) {
        rgb[n] <- (as.integer(rgb[n])/2.0) %>% as.integer()
    }
    rgb <- paste(rgb, collapse = ",")
    print(rgb)
    rgb
}

BedFileDirOut <- sprintf("/Users/rtearle/Documents/Roseworthy_Projects/6.Students/Nitish_Salmonella/Apr_2018/Out/R_DNA_ANNOT/%s/BED", Component)
dir.create(BedFileDirOut)
n = 1; m = 2

for (n in 1:NrGenomes){

    Query <- Genomes[n,1] %>% .[[1]] # CDS query genome name

    for (m in 2:NrGenomes) {

        Subject <- Genomes[n,m] %>% .[[1]] # CDS subject genome name
        Comp <- paste0(Query, "_vs_", Subject)

        BedLabel <- sprintf('track name="%s" description="%s" itemRgb="On"', Comp, Comp)
        Bed <- CDSBlastsVariants %>% filter(Comp == Comp, !is.na(SubjectContig)) %>%
            mutate(Start = ifelse(GenSubjectStart < GenSubjectEnd, GenSubjectStart, GenSubjectEnd),
                   End = ifelse(GenSubjectStart > GenSubjectEnd, GenSubjectStart, GenSubjectEnd),
                   Start2 = Start,
                   End2 = End,
                   RGB = ifelse(GenStrand == "plus", "255,0,0", "0,0,255"),
                   RGB2 = HalveIntensity(RGB),
                   X = ifelse(startsWith(Function, "hypothetical"), RGB2, RGB))  %>%
    select(SubjectContig, Start, End, CDSID, CDSIdentity, GenStrand, Start2, End2, RGB2, QueryID) %>%
            unite(CDSID, CDSID, QueryID, sep = "|", remove = TRUE) %>%
            arrange(factor(SubjectContig), Start)

        # Bed <- CDSBlastsVariants %>% filter(Comp == Comp, !is.na(SubjectContig)) %>%
        #     mutate(Start = ifelse(GenSubjectStart < GenSubjectEnd, GenSubjectStart, GenSubjectEnd),
        #            End = ifelse(GenSubjectStart > GenSubjectEnd, GenSubjectStart, GenSubjectEnd),
        #            Start2 = Start, End2 = End, RGB = ifelse(GenStrand == "plus", "255,0,0", "0,0,255")) %>%
        #     select(SubjectContig, Start, End, CDSID, CDSIdentity, GenStrand, Start2, End2, RGB, QueryID, Function) %>%
        #     unite(CDSID, CDSID, QueryID, Function, sep = "|", remove = TRUE) %>%
        #     mutate(CDSID = paste0('\"', CDSID, '\"'))

        FileOut <- sprintf("%s_CDS-vs-%s_Genome.bed", Query, Subject)
        writeLines(BedLabel, file.path(BedFileDirOut, FileOut))
        write_tsv(Bed, file.path(BedFileDirOut, FileOut), col_names = FALSE, na = "", append = TRUE)

    }
}






z <- BestDNABlasts2 %>% filter(cdsid == "cdsid0338")
