library(magrittr)
library(tibble)
library(tools)
library(stringr)
library(RColorBrewer)
library(readxl)
library(MASS)
library(Matrix)
library(reshape2)
library(gplots)
library(devtools)
library(plyr)
library(Biostrings)
library(muscle)
library(tibble)
library(tidyverse)

StartTime <- Sys.time()

###############
## Functions ##
###############

'%nin%' <- Negate('%in%')

uniquelength <- function(x) {x %>% unique() %>% length()}

# Rotate the elements of a vector
rot <- function(x, l) (x[c(2:l,1)])

# Filter a blast hit table for the best unique entries for each query and subject
# Blast hit table from comparison of 2 sets of CDS
#Hits <- DNABlasts[[1]]
FindBestBlastHits <- function(Hits) {
    # Assumes Hits is a df of blast hits
    BestHits <- Hits

    # Group by gseqid and keep top hit, based on qs
    BestHits %<>% group_by(qseqid) %>% filter(qs == max(qs, na.rm = TRUE))
    BestHits %<>% group_by(qseqid) %>% filter(pident == max(pident, na.rm = TRUE))

    # Remove matches qfraction and sfraction one less than 0.95, other less than 0.7
    BestHits %<>% filter((qfraction > 0.95 & sfraction > 0.7) | (qfraction > 0.7 & sfraction > 0.95)) # arbitrary threshold for matching lengths

    BestHits
}

# Filter a blast hit table for the best unique entries for each query and subject
# Blast hit table from comparison of a sets of CDS to a set of contigs
#Hits <- DNABlasts2[["KC_109_vs_T-Chr"]] %>% filter(cdsid == "cdsid0338")
FindBestBlastHits2 <- function(Hits) {
    # Assumes Hits is a df of blast hits
    BestHits <- Hits

    # Pick best hit by highest qfraction and idenitity
    # BestHits %<>% group_by(qseqid) %>% filter(qfraction == max(qfraction, na.rm = TRUE),
    #                                          )
    BestHits %<>% group_by(qseqid) %>% filter(qfraction == max(qfraction, na.rm = TRUE))
    BestHits %<>% group_by(qseqid) %>% filter(pident == max(pident, na.rm = TRUE))

    # Remove matches with low qfraction/sfraction
    #BestHits %<>% filter((qfraction > 0.95 & slength/qfulllength > 0.7 ) | (qfraction > 0.7 & slength/qfulllength > 0.95 )) # arbitrary threshold for matching lengths
    #BestHits %<>% filter(qfraction > 0.5) # arbitrary threshold for matching lengths
    BestHits %<>% filter(qfraction > 0.25) # arbitrary threshold for matching lengths

    # Check for multiple qseqid entries and choose one at random
    # DupRowNrs <- BestHits$qseqid %>% duplicated() %>% which()
    # DupIDs <- BestHits[DupRowNrs,] %>% .$qseqid
    # AllDupRowNrs <- (BestHits$qseqid %in% DupIDs) %>% which()
    # if (length(DupIDs) > 0) {
    #     Keep <- vector()
    #     for (ID in DupIDs) {
    #         DupRowNrs2 <- (BestHits$qseqid == ID) %>% which() # now includes the first duplicate
    #         Keep <- c(Keep, sample(DupRowNrs2,1)) # choose one to keep
    #     }
    #     Lose <- AllDupRowNrs[-Keep] # from the duplicate row nrs in the original df, remove the one we want to keep
    #     BestHits %<>% .[-Lose,] # remove the others
        #BestHits$qseqid %>% unique() %>% length()
    #}

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

AnnotDirOut <- "/Users/rtearle/Documents/Roseworthy_Projects/6.Students/Nitish_Salmonella/May_2018/Out/R_Annotation" # DNA analysis contained in this script
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
        Annotations[[Name]] %<>% dplyr::select(-c(Start2, Stop2))

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

        Annotations[[Name]] %<>% dplyr::select(Genome, Component, everything())

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

DataSizes %<>% dplyr::select(-c(Columns)) %>% dplyr::rename(CDSCount = Rows)

FileOut <- paste0(FilePrefix,"-CDS_Counts_Genome_Component.tsv")
write_tsv(DataSizes, file.path(AnnotDirOut, FileOut), col_names = TRUE, na = "")

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
DNABlastsList <- list() # leaving as a list
BlastColNames <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sstrand")

# Load blast data
for (n in 1:nrow(DNABlastFiles)) {
    Name <- paste(DNABlastFiles$Comparison[n], DNABlastFiles$Component[n], sep = "-")

    print(Name)
    DNABlastsList[[Name]] <- read_tsv(file.path(DNABlastDirIn, DNABlastFiles$Component[n], DNABlastFiles$File[n]), col_names = FALSE, skip = 0) %>% as_data_frame()
    colnames(DNABlastsList[[Name]]) <- BlastColNames
    DNABlastsList[[Name]] %<>% filter(grepl("peg", qseqid))
}
DNABlastsList %>% length() # 60

# Modify blast data, adding extra data
for (n in 1:nrow(DNABlastFiles)) {
    #Name <- DNABlastFiles$Comparison[n]
    Name <- DNABlastsList[n] %>% names()
    Query <- DNABlastFiles$Query[n]
    Subject <- DNABlastFiles$Subject[n]
    Component <- DNABlastFiles$Component[n]
    QAnnot <- Annotations %>% filter(Genome == Query & Component == Component)
    SAnnot <- Annotations %>% filter(Genome == Subject & Component == Component)

    # Remove repeats
    DNABlastsList[[Name]] %<>% filter(!grepl("repeat", qseqid) & !grepl("repeat", sseqid))

    # Remove "gnl\\|" from Blast sseqid columns, so ids match those in annotations
    DNABlastsList[[Name]]$sseqid %<>% gsub("gnl\\|", "", .)

    # Add number of qseqid and sseqid for sorting
    DNABlastsList[[Name]]$qseqidnr <- DNABlastsList[[Name]]$qseqid %>% gsub("fig.+peg.", "", .) %>% as.integer()
    DNABlastsList[[Name]]$sseqidnr <- DNABlastsList[[Name]]$sseqid %>% gsub("fig.+peg.", "", .) %>% as.integer()

    # add function for query, subject
    DNABlastsList[[Name]] %<>% left_join(QAnnot %>% dplyr::select(FeatureID, Function), by = c("qseqid" = "FeatureID"))
    DNABlastsList[[Name]] %<>% dplyr::rename(qfunction = Function)
    DNABlastsList[[Name]] %<>% left_join(SAnnot %>% dplyr::select(FeatureID, Function), by = c("sseqid" = "FeatureID"))
    DNABlastsList[[Name]] %<>% dplyr::rename(sfunction = Function)

    # Add query and subject match lengths
    DNABlastsList[[Name]] %<>% mutate(qlength = abs(qend - qstart)+1)
    DNABlastsList[[Name]] %<>% mutate(slength = abs(send - sstart)+1)

    # Add query and subject full lengths
    DNABlastsList[[Name]] %<>% left_join(QAnnot %>% dplyr::select(FeatureID, Length), by = c("qseqid" = "FeatureID"))
    DNABlastsList[[Name]] %<>% dplyr::rename(qfulllength = Length)
    DNABlastsList[[Name]] %<>% left_join(SAnnot %>% dplyr::select(FeatureID, Length), by = c("sseqid" = "FeatureID"))
    DNABlastsList[[Name]] %<>% dplyr::rename(sfulllength = Length)

    # Add query and subject fractions length/fulllength
    DNABlastsList[[Name]] %<>% mutate(qfraction = qlength/qfulllength)
    DNABlastsList[[Name]] %<>% mutate(sfraction = slength/sfulllength)

    # Add query and subject genomes, component, move to front
    DNABlastsList[[Name]] %<>% mutate(qgenome = Query)
    DNABlastsList[[Name]] %<>% mutate(sgenome = Subject)
    DNABlastsList[[Name]] %<>% mutate(component = Component)

    NrCols <- DNABlastsList[[Name]] %>% ncol()
    DNABlastsList[[Name]] %<>% dplyr::select(qgenome, sgenome, component, everything())
}
DNABlastsList %>% length() # 60

# Sizes of datasets
DataSizes <- tibble(Comparison = character(0), Rows = integer(0), Columns = integer(0), UniqueSeqIds = integer(0))
for (n in 1:length(DNABlastsList)) {
    x <- DNABlastsList[[names(DNABlastsList[n])]]
    UniqueLength <- x$qseqid %>% unique() %>% length()
    DataSizes <- add_row(DataSizes, Comparison = names(DNABlastsList[n]), Rows = x %>% nrow(), Columns = x %>% ncol(), UniqueSeqIds = UniqueLength)
}

DataSizes
#    Comparison          Rows Columns UniqueSeqIds
#    <chr>              <int>   <int>        <int>
#  1 KC_109_vs_T-Chr   138514      26         4843
#  2 KC_109_vs_TM-Chr  138546      26         4843
#  3 KC_109_vs_TMu-Chr 138568      26         4843
#  4 T_vs_KC_109-Chr   138526      26         4839
#  5 T_vs_TM-Chr       138531      26         4839
#  6 T_vs_TMu-Chr      138553      26         4839
#  7 TM_vs_KC_109-Chr  138550      26         4839
#  8 TM_vs_T-Chr       138522      26         4839
#  9 TM_vs_TMu-Chr     138576      26         4839
# 10 TMu_vs_KC_109-Chr 138533      26         4840
# # ... with 50 more rows

# Save Sizes
FileOut <- paste0(FilePrefix, "-CDS_Blasts-Counts.tsv")
write_tsv(DataSizes, file.path(AnnotDirOut, FileOut), col_names = TRUE, na = "")

rm(n, DataSizes, UniqueLength, Name, Query, Subject, Component, x, QAnnot, SAnnot, NrCols)

############################################
## Extract CDS vs CDS DNA blast best hits ##
############################################
# Select the best hit for each query
BestDNABlastsList <- list() # will be converted to df later
for (n in 1:length(DNABlastsList)) {
    Name <- names(DNABlastsList[n])

    # add qs to blast data
    DNABlastsList[[Name]] %<>% mutate(qs = qfraction * sfraction)
    # print(DNABlasts[[Name]])

    # Get best unique matches for each query and subject id
    BestDNABlastsList[[Name]] <- DNABlastsList[[Name]] %>% FindBestBlastHits()
    BestDNABlastsList[[Name]] %<>% unique()
}

# Sizes of datasets
DataSizes <- tibble(Comparison = character(0), Rows = integer(0), Columns = integer(0), UniqueQseqIds = integer(0), UniqueSseqIds = integer(0))
for (n in 1:length(BestDNABlastsList)) {
    Name <- BestDNABlastsList[n] %>% names() %>% .[1]
    qseqidUniqueLength <- BestDNABlastsList[[Name]]$qseqid %>% uniquelength()
    sseqidUniqueLength <- BestDNABlastsList[[Name]]$sseqid %>% uniquelength()
    DataSizes <- add_row(DataSizes, Comparison = Name, Rows = BestDNABlastsList[[Name]] %>% nrow(), Columns = BestDNABlastsList[[Name]] %>% ncol(),
                         UniqueQseqIds = qseqidUniqueLength, UniqueSseqIds = sseqidUniqueLength)
}

DataSizes
#    Comparison         Rows Columns UniqueQseqIds UniqueSseqIds
#    <chr>             <int>   <int>         <int>         <int>
#  1 KC_109_vs_T-Chr    4826      27          4826          4825
#  2 KC_109_vs_TM-Chr   4813      27          4813          4812
#  3 KC_109_vs_TMu-Chr  4829      27          4829          4829
#  4 T_vs_KC_109-Chr    4827      27          4827          4827
#  5 T_vs_TM-Chr        4822      27          4822          4822
#  6 T_vs_TMu-Chr       4835      27          4835          4835
#  7 TM_vs_KC_109-Chr   4814      27          4814          4814
#  8 TM_vs_T-Chr        4822      27          4822          4822
#  9 TM_vs_TMu-Chr      4824      27          4824          4824
# 10 TMu_vs_KC_109-Chr  4829      27          4829          4829
# # ... with 50 more rows

# Save Sizes
FileOut <- paste0(FilePrefix, "-CDS_Best_Blasts-Counts.tsv")
write_tsv(DataSizes, file.path(AnnotDirOut, FileOut), col_names = TRUE, na = "")

(DataSizes$Rows == DataSizes$UniqueQseqIds) %>% summary()
#    Mode    TRUE    NA's
# logical      60       0
(DataSizes$Rows == DataSizes$UniqueSseqIds) %>% summary()
#    Mode   FALSE    TRUE    NA's
# logical       7      53       0

DataSizes$Comparison[(DataSizes$Rows != DataSizes$UniqueSseqIds) %>% which()]

rm(n,Name, qseqidUniqueLength, sseqidUniqueLength, DataSizes)

# Make one df from best hits list
BestDNABlasts <- BestDNABlastsList[[1]]
for (n in 2:length(BestDNABlastsList)) {
    BestDNABlasts %<>% rbind(BestDNABlastsList[[n]])
}
BestDNABlasts %>% dim() # 61215    27
BestDNABlasts %<>% ungroup()

BestDNABlasts %>% dplyr::select(qgenome, component, qseqid, sseqid) %>% dim() # 61215     4
BestDNABlasts %>% dplyr::select(qgenome, component, qseqid, sseqid) %>% unique() %>% dim() # 61215     4

BestDNABlasts$qs %>% summary()
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.7212  1.0000  1.0000  0.9991  1.0000  1.0000

(BestDNABlasts$qs < 1) %>% summary()
#    Mode   FALSE    TRUE    NA's
# logical   60571     644       0


# x <- BestDNABlasts %>% filter(qs < 1)
# x <- DNABlasts[[1]] %>% filter(qseqid == "fig|590.17911.peg.72")
# x <- BestDNABlasts %>% filter(qseqid == "fig|590.17911.peg.72")
# DNABlasts[[1]] %>% filter(qseqid == "fig|590.17911.peg.72") %>%
#     group_by(qseqid) %>% dplyr::filter(qs == max(qs, na.rm = TRUE))
# DNABlasts[[1]] %>% filter(qseqid == "fig|590.17911.peg.72") %>% .$qs
#     max()
#
# DNABlasts[[1]] %>% filter(qseqid == "fig|590.17911.peg.73")
# BestDNABlasts[[1]] %>% filter(qseqid == "fig|590.17911.peg.73")

(BestDNABlasts$qfunction == BestDNABlasts$sfunction) %>% summary()
#    Mode   FALSE    TRUE    NA's
# logical      18   61197       0

#x <- BestDNABlasts %>% filter(qfunction != sfunction)

rm(df, n)

########################################
## Group by DNA blasts across genomes ##
########################################
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

# Generating counts
for (n in 1:nrow(DNABlastFiles)) {
    Name <- DNABlastFiles$Comparison[n]
    Query <- DNABlastFiles$Query[n]
    Subject <- DNABlastFiles$Subject[n]
    Component <- DNABlastFiles$Component[n]

    # Feature ids for query genome and component, or subject genome/component
    QueryFeatureIDs <- Annotations %>% filter(Genome == DNABlastFiles$Query[n] & Component == DNABlastFiles$Component[n]) %>% .$FeatureID
    SubjectFeatureIDs <- Annotations %>% filter(Genome == DNABlastFiles$Subject[n]  & Component == DNABlastFiles$Component[n]) %>% .$FeatureID

    # Find blasts for this comparison
    Comps <- BestDNABlasts %>% filter(qgenome == Query, sgenome == Subject, component == Component)

    # query and subject comparison feature ids
    CompQueryIDs <- Comps$qseqid %>% unique()
    CompSubjectIDs <- Comps$sseqid %>% unique()

    # Count how many found not found for each
    AnnotInQuery <- (QueryFeatureIDs %in% CompQueryIDs) %>% sum(na.rm=TRUE)
    AnnotNinQuery <- (QueryFeatureIDs %nin% CompQueryIDs) %>% sum(na.rm=TRUE)
    AnnotInSubject <- (SubjectFeatureIDs %in% CompSubjectIDs) %>% sum(na.rm=TRUE)
    AnnotNinSubject <- (SubjectFeatureIDs %nin% CompSubjectIDs) %>% sum(na.rm=TRUE)

    # Add to counts
    IDCounts %<>% add_row(Comparison = Name, Component = Component,
                          Query = Query, QueryCount = uniquelength(QueryFeatureIDs), CompQueryCount = length(CompQueryIDs),
                          Subject = Subject, SubjectCount = uniquelength(SubjectFeatureIDs), CompSubjectCount = length(CompSubjectIDs),
                          AnnotInQuery = AnnotInQuery, AnnotNinQuery = AnnotNinQuery,
                          AnnotInSubject = AnnotInSubject, AnnotNinSubject = AnnotNinSubject)

        # Collect feature ids not found in other genomes
        if(AnnotNinQuery > 0) {
        t <- (QueryFeatureIDs %nin% CompQueryIDs) %>% which()
        x <- Annotations %>% filter(Genome == Query & Component == Component) %>% .[t,]
        y <- tibble(Comparison = c(rep(Name,nrow(x))), Query = c(rep(Query,nrow(x))), Subject = c(rep(Subject,nrow(x))))
        x <- cbind(y, x) %>% as_data_frame()
        MissingCDS %<>% bind_rows(x)
    }
}
rm(Names, n, Name, Query, Subject, Component, QueryFeatureIDs, SubjectFeatureIDs, Comps, CompQueryIDs, CompSubjectIDs, AnnotInQuery,
   AnnotNinQuery, AnnotInSubject, AnnotNinSubject, t, x, y)

IDCounts %<>% arrange(Component, Subject)
IDCounts %>% dim() # 60 12

# Save CDS counts
FileOut <- paste0(FilePrefix, "-CDS_vs_CDS-Counts.tsv")
write_tsv(IDCounts, file.path(AnnotDirOut, FileOut), col_names = TRUE, na = "")

# Save missing CDS
MissingCDS %<>% dplyr::select(-Genome)
FileOut <- paste0(FilePrefix,"-CDS_vs_CDS-Missing.tsv")
write_tsv(MissingCDS, file.path(AnnotDirOut, FileOut), col_names = TRUE, na = "")

# Make a feature table linking ids across genomes
FeatureIDTableList <- list()
for (g in 1:NrGenomes){

    Query <- Genomes[g,1] %>% .[[1]]

    FeatureIDTable <- tibble(Component = Annotations %>% filter (Genome == Query) %>% .$Component,
                             Name = Annotations %>% filter (Genome == Query) %>% .$FeatureID)
    colnames(FeatureIDTable) <- c("Component", Query)

    for (g2 in 2:ncol(Genomes)){

        Subject <- Genomes[g,g2] %>% .[[1]]

        #if (Query == Subject) {next} # do not do self to self

        # Join Query/Subject seqids from BestDNABlasts
        Comps <- BestDNABlasts %>% filter(qgenome == Query, sgenome == Subject) %>% dplyr::select(component, qseqid, sseqid)

        # y <- DNABlasts[[1]] %>% filter(qseqid == "fig|590.17911.peg.72")
        # BestDNABlasts %>% filter(qseqid == "fig|590.17911.peg.72")
        # Comps %>% filter(qseqid == "fig|590.17911.peg.72")
        if (nrow(Comps) > 0) {
            FeatureIDTable %<>% left_join(Comps, by = c("Component" = "component", setNames("qseqid", Query)))
            colnames(FeatureIDTable)[ncol(FeatureIDTable)] <- Subject
        }
    }
    FeatureIDTableList[[Query]] <- FeatureIDTable
}
rm(g, Query, g2, Subject, Comps)

# Merge feature tables
FeatureIDTable <- FeatureIDTableList[[Genomes$One[1]]] # get first table
for (n in 2:NrGenomes) {
    FeatureIDTable %<>% union(FeatureIDTableList[[Genomes$One[n]]] %>% dplyr::select(Component, Genomes$One))
}
rm(n)

FeatureIDTable %<>% dplyr::select(Component, eval(Genomes$One)) # order cols
FeatureIDTable %>% dim() # 5226    5

# Sort feature table by component and then each genome feature if in turn
# New cols of feature ids reduced to just the last element (integer) so we can sort
for (G in Genomes$One) {
    Name <- paste0(G,"_Nr")
    FeatureIDTable %<>% mutate(!!(Name) := gsub("fig.+peg.", "", eval(as.name(G))) %>% as.integer())
}
SortNames <- c("Component", paste0(Genomes$One,"_Nr")) # create names to sort with
FeatureIDTable %<>% .[do.call(order, .[, SortNames]), ] # sort
FeatureIDTable <- FeatureIDTable[,-((ncol(FeatureIDTable)-3):ncol(FeatureIDTable))] # remove the feature id integer cols
rm(G, Name, SortNames)

# Count all and unique feature ids
Counts <- tibble(Genome = character(0), Component = character(0), NrFeatureIDs = integer(0), NrUniqueFeatureIDs = integer(0))
for (G in Genomes$One) {
    for (C in Components) {
        x <- FeatureIDTable %>% filter(Component == C) %>% .[,G] %>% unlist %>% .[!is.na(.)]
        c <- x %>% length()
        c2 <- x %>% unique() %>% length()
        Counts %<>% add_row(Genome = G, Component = C, NrFeatureIDs = c, NrUniqueFeatureIDs = c2)

    }
}
Counts
#    Genome Component NrFeatureIDs NrUniqueFeatureIDs
#    <chr>  <chr>            <int>              <int>
#  1 KC_109 Chr               4845               4843
#  2 KC_109 Plasmid_1          145                145
#  3 KC_109 Plasmid_2          157                155
#  4 KC_109 Plasmid_3            3                  3
#  5 KC_109 Plasmid_4            8                  6
#  6 T      Chr               4840               4839
#  7 T      Plasmid_1          142                142
#  8 T      Plasmid_2          149                147
#  9 T      Plasmid_3            3                  3
# 10 T      Plasmid_4            7                  6
# 11 TM     Chr               4840               4839
# 12 TM     Plasmid_1          148                148
# 13 TM     Plasmid_2          148                146
# 14 TM     Plasmid_3            3                  3
# 15 TM     Plasmid_4            5                  4
# 16 TMu    Chr               4842               4840
# 17 TMu    Plasmid_1          142                142
# 18 TMu    Plasmid_2          145                144
# 19 TMu    Plasmid_3            6                  6
# 20 TMu    Plasmid_4            6                  5
rm(G, C, x, c, c2, Counts)

# Resolve dups
KeepDups <- FeatureIDTable %>% filter(Component == "XXX") # generate empty df
RemoveDups <- FeatureIDTable %>% filter(Component == "XXX") # generate empty df
for (g in 1:NrGenomes){
    G <- Genomes$One[g]

    # Find duplicate entries in Feature table for this genome
    Dups <- FeatureIDTable %>% group_by(eval(as.name(G))) %>% # creates a new col, dont know why
        filter(!is.na(eval(as.name(G)))) %>%
        filter(n()>1) %>%
        arrange(eval(as.name(G))) %>%
        ungroup() %>%
        dplyr::select(-(ncol(.))) # remove the last col, created by group_by command

    # Remove clear duplicated records - this may remove some true dups?
    if(nrow(Dups) < 2) {next}

    FeatureIDs <- Dups[,G] %>% unique() %>% unlist() %>% as.vector()
    #FID <- "fig|590.17913.peg.152"

    for (FID in FeatureIDs) {
        Dups2 <- Dups %>% filter(eval(as.name(G)) == FID)

        # If doing head to head comps the same, keep row with most Feature IDs
        r <- nrow(Dups2)
        c <- ncol(Dups2)
        for (n in 1:(r-1)) {
            for (m in (n+1):r) {
                # Check to see if either is already on the remove list, if so, do nothing
                # print(c(n,m))
                # print(Dups2[c(n,m),])
                # print(nrow(anti_join(Dups2[n,], RemoveDups)) == 0 | nrow(anti_join(Dups2[m,], RemoveDups)) == 0)
                if(nrow(anti_join(Dups2[n,], RemoveDups)) == 0 | nrow(anti_join(Dups2[m,], RemoveDups)) == 0) {
                    # at least one already in remove list, do nothing
                } else {
                    # neither in remove list so choose one to remove
                    Match <- (Dups2[n,] == Dups2[m,]) %>% .[1,2:c] %>% summary() # do genome FIDs match?
                    if ("FALSE" %in% names(Match)) { # at least FID across genomes does not match, keep both
                        KeepDups %<>% rbind(Dups2[c(n,m),])
                    } else { # FIDs match apart from NAs
                        NotNAs <- rowSums(!is.na(Dups2[c(n,m),]))
                        if (NotNAs[1] > NotNAs[2]) {
                            RemoveDups %<>% rbind(Dups2[m,])
                        } else {
                            RemoveDups %<>% rbind(Dups2[n,])
                        }
                    }
                }
            }
        }
    }
}
rm(g, G, Dups, FeatureIDs, FID, Dups2, r, c, n, m, Match)

KeepDups %<>% unique() %>% arrange(Component, eval(as.name(Genomes$One)))
RemoveDups %<>% unique() %>% arrange(Component, eval(as.name(Genomes$One)))
KeepDups %>% dim() # 5 5
RemoveDups %>% dim() # 5 5
FeatureIDTable %>% dim() # 5226    5

# Remove dups marked to be removed
if(nrow(RemoveDups) > 0) {
    print(paste("Removing", nrow(RemoveDups), "rows from data as duplicate entries", sep = " "))
    FeatureIDTable %<>% setdiff(RemoveDups)
}
FeatureIDTable %>% dim() # 5221    5

Counts <- tibble(Genome = character(0), Component = character(0), NrFeatureIDs = integer(0), NrUniqueFeatureIDs = integer(0))
for (G in Genomes$One) {
    for (C in Components) {
        x <- FeatureIDTable %>% filter(Component == C) %>% .[,G] %>% unlist %>% .[!is.na(.)]
        c <- x %>% length()
        c2 <- x %>% unique() %>% length()
        Counts %<>% add_row(Genome = G, Component = C, NrFeatureIDs = c, NrUniqueFeatureIDs = c2)
    }
}
Counts
#    Genome Component NrFeatureIDs NrUniqueFeatureIDs
#    <chr>  <chr>            <int>              <int>
#  1 KC_109 Chr               4844               4843
#  2 KC_109 Plasmid_1          145                145
#  3 KC_109 Plasmid_2          156                155
#  4 KC_109 Plasmid_3            3                  3
#  5 KC_109 Plasmid_4            6                  6
#  6 T      Chr               4840               4839
#  7 T      Plasmid_1          142                142
#  8 T      Plasmid_2          147                147
#  9 T      Plasmid_3            3                  3
# 10 T      Plasmid_4            6                  6
# 11 TM     Chr               4840               4839
# 12 TM     Plasmid_1          148                148
# 13 TM     Plasmid_2          146                146
# 14 TM     Plasmid_3            3                  3
# 15 TM     Plasmid_4            4                  4
# 16 TMu    Chr               4841               4840
# 17 TMu    Plasmid_1          142                142
# 18 TMu    Plasmid_2          144                144
# 19 TMu    Plasmid_3            6                  6
# 20 TMu    Plasmid_4            5                  5


FileOut <- paste0(FilePrefix,"-CDS-FeatureIDTable-Counts.tsv")
write_tsv(Counts, file.path(AnnotDirOut, FileOut), col_names = TRUE, na = "")

FileOut <- paste0(FilePrefix,"-FeatureIDTable-Duplicate_Entries.tsv")
write_tsv(KeepDups, file.path(AnnotDirOut, FileOut), col_names = TRUE, na = "")

rm(G, c, c2, Counts)
rm(FileOut, KeepDups, RemoveDups)

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
# CDSIDs %>% head()
# CDSIDs %>% tail()

# Add CDS IDs to feature table
FeatureIDTable <- cbind(CDSIDs,FeatureIDTable) %>% as_data_frame()

# Add CDS IDs to annotations
CDSIDFeaturesTable <- melt(FeatureIDTable %>% dplyr::select(-Component), id = c("cdsid")) %>% filter(!is.na(value)) %>% as_data_frame() # easier to join
colnames(CDSIDFeaturesTable) <- c("CDSID", "Genome", "FeatureID")
CDSIDFeaturesTable %<>% dplyr::select(-Genome)

Annotations %>% dim() # 20566    17
Annotations %<>% left_join(CDSIDFeaturesTable, by = "FeatureID")
Annotations %>% dim() # 20571    18 # now have dups

BestDNABlasts %>% dim() # 61215    27

BestDNABlasts$cdsid <- NA_character_
for (n in 1:nrow(BestDNABlasts)) {
    BestDNABlasts$cdsid[n] <- FeatureIDTable %>% filter(eval(as.name(BestDNABlasts$qgenome[n])) == BestDNABlasts$qseqid[n],
                                                        eval(as.name(BestDNABlasts$sgenome[n])) == BestDNABlasts$sseqid[n]) %>% .$cdsid %>% .[[1]]
}

BestDNABlasts %<>% dplyr::select(cdsid, everything())
BestDNABlasts %>% dim() # 61215    28

grepl("cdsid", BestDNABlasts$cdsid) %>% summary() # should not be any
#    Mode    TRUE    NA's
# logical   61215       0

rm(n, CDSIDFeaturesTable, m, Pad, t)

######################
## Presence/Absence ##
######################
# Get patterns of presence/absence across genomes
FeaturePresenceTable <- FeatureIDTable %>% dplyr::select(-c(cdsid, Component)) %>% as.matrix(nrow = nrow(.), ncol = ncol(.))
FeaturePresenceTable %>% dim() # 5221    4
FeaturePresenceTable[!is.na(FeaturePresenceTable)] <- TRUE # presence/absence of an id to T/F
FeaturePresenceTable[is.na(FeaturePresenceTable)] <- FALSE
FeaturePresenceTable %>% dim() # 5221    4
FeaturePresenceTable %>% summary()
#    KC_109         T            TM          TMu
#  FALSE:  67   FALSE:  83   FALSE:  80   FALSE:  83
#  TRUE :5154   TRUE :5138   TRUE :5141   TRUE :5138

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
x
#    Present                Counts
#    <chr>                   <int>
#  1 TRUE TRUE TRUE TRUE      5075
#  2 TRUE TRUE TRUE FALSE        3
#  3 TRUE TRUE FALSE TRUE       22
#  4 TRUE TRUE FALSE FALSE       5
#  5 TRUE FALSE TRUE TRUE        4
#  6 TRUE FALSE TRUE FALSE      11
#  7 TRUE FALSE FALSE TRUE       2
#  8 TRUE FALSE FALSE FALSE     32
#  9 FALSE TRUE TRUE TRUE       16
# 10 FALSE TRUE TRUE FALSE       5
# 11 FALSE TRUE FALSE TRUE       8
# 12 FALSE TRUE FALSE FALSE      4
# 13 FALSE FALSE TRUE TRUE       4
# 14 FALSE FALSE TRUE FALSE     23
# 15 FALSE FALSE FALSE TRUE      7

# Save feature presence table summary
x <- str_split_fixed(x$Present, " ", 4) %>% cbind(x$Counts) %>% as_data_frame()
colnames(x) <- c(Genomes$One, "Counts")

FileOut <- paste0(FilePrefix,"-Matching_CDS_IDs-Truth_Table.tsv")
write_tsv(FeaturePresenceTable, file.path(AnnotDirOut, FileOut), col_names = TRUE, na = "")

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
CDSComparisons$Mean_AA <- CDSComparisons %>% dplyr::select(ends_with("_AA_Length")) %>% rowMeans(na.rm = TRUE)
CDSComparisons$SD_AA <- apply(CDSComparisons %>% dplyr::select(ends_with("_AA_Length")), 1, sd, (na.rm = TRUE))
CDSComparisons$Mean_NT <- CDSComparisons %>% dplyr::select(ends_with("_NT_Length")) %>% rowMeans(na.rm = TRUE)
CDSComparisons$SD_NT <- apply(CDSComparisons %>% dplyr::select(ends_with("_NT_Length")), 1, sd, (na.rm = TRUE))

CDSComparisons %>% dim() # 5221   27

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
rm(n, m, File)

DNABlastFiles2 %>% dim() # 60  5
#DNABlastFiles2$File

# Get blast results - note use of double brackets to extract name for use
DNABlasts2List <- list() # leaving as a list
BlastColNames <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sstrand")

# Load blast data
for (n in 1:nrow(DNABlastFiles2)) {
    Name <- DNABlastFiles2$Comparison[n]
    print(Name)
    DNABlasts2List[[Name]] <- read_tsv(file.path(DNABlastDirIn2, DNABlastFiles2$Component[n], DNABlastFiles2$File[n]), col_names = FALSE, skip = 0) %>% as_data_frame()
    colnames(DNABlasts2List[[Name]]) <- BlastColNames
    DNABlasts2List[[Name]] %<>% filter(grepl("peg", qseqid))
}
rm(Name)

# Melt feature table, makes it easier to join CDS IDs to blast results
CDSIDFeaturesTable <- melt(FeatureIDTable, id = c("cdsid")) %>% filter(!is.na(value)) %>% as_data_frame()
colnames(CDSIDFeaturesTable) <- c("cdsid", "Genome", "qseqid")
CDSIDFeaturesTable %<>% dplyr::select(-Genome)

# Modify blast data, adding extra data
for (n in 1:nrow(DNABlastFiles2)) {

    Name <- DNABlastFiles2$Comparison[n]
    Query <- DNABlastFiles2$Query[n]
    Subject <- DNABlastFiles2$Subject[n]
    Component <- DNABlastFiles2$Component[n]
    QAnnot <- Annotations %>% filter(Genome == Query)
    SAnnot <- Annotations %>% filter(Genome == Subject)

    # Remove repeats and rna
    DNABlasts2List[[Name]] %<>% filter(!grepl("repeat", qseqid))
    DNABlasts2List[[Name]] %<>% filter(!grepl("rna", qseqid))
    DNABlasts2List[[Name]] %<>% filter(!grepl("crispr", qseqid))


    # Remove "gnl\\|" from Blast sseqid columns, so ids match those in annotations
    DNABlasts2List[[Name]]$sseqid %<>% gsub("gnl\\|", "", .)

    # Add number of qseqid and sseqid for sorting
    DNABlasts2List[[Name]]$qseqidnr <- DNABlasts2List[[Name]]$qseqid %>% gsub("fig.+peg.", "", .) %>% as.integer()
    DNABlasts2List[[Name]]$sseqidnr <- DNABlasts2List[[Name]]$sseqid %>% gsub("fig.+peg.", "", .) %>% as.integer()

    # add function for query and subject
    DNABlasts2List[[Name]] %<>% left_join(QAnnot %>% dplyr::select(FeatureID, Function), by = c("qseqid" = "FeatureID"))
    DNABlasts2List[[Name]] %<>% dplyr::rename(qfunction = Function)
    DNABlasts2List[[Name]] %<>% left_join(SAnnot %>% dplyr::select(FeatureID, Function), by = c("sseqid" = "FeatureID"))
    DNABlasts2List[[Name]] %<>% dplyr::rename(sfunction = Function)

    # Add query and subject match lengths
    DNABlasts2List[[Name]] %<>% mutate(qlength = abs(qend - qstart)+1)
    DNABlasts2List[[Name]] %<>% mutate(slength = abs(send - sstart)+1)

    # Add query and subject full lengths
    DNABlasts2List[[Name]] %<>% left_join(QAnnot %>% dplyr::select(FeatureID, Length), by = c("qseqid" = "FeatureID"))
    DNABlasts2List[[Name]] %<>% dplyr::rename(qfulllength = Length)
    DNABlasts2List[[Name]] %<>% left_join(SAnnot %>% dplyr::select(FeatureID, Length), by = c("sseqid" = "FeatureID"))
    DNABlasts2List[[Name]] %<>% dplyr::rename(sfulllength = Length)

    # Add query and subject fractions length/fulllength
    DNABlasts2List[[Name]] %<>% mutate(qfraction = qlength/qfulllength)
    DNABlasts2List[[Name]] %<>% mutate(sfraction = slength/sfulllength)

    # Add cdsid, query and subject genomes, component, move to front
    DNABlasts2List[[Name]] <- left_join(DNABlasts2List[[Name]], CDSIDFeaturesTable, by = "qseqid")
    DNABlasts2List[[Name]] %<>% mutate(qgenome = Query)
    DNABlasts2List[[Name]] %<>% mutate(sgenome = Subject)
    DNABlasts2List[[Name]] %<>% mutate(component = Component)

    DNABlasts2List[[Name]] %<>% dplyr::select(cdsid, qgenome, sgenome, component, everything())

}
DNABlasts2List %>% length() # 60

# Sizes of datasets
DataSizes <- tibble(Genome = character(0), Rows = integer(0), Columns = integer(0), UniqueSeqIds = integer(0))
for (n in 1:length(DNABlasts2List)) {
    x <- DNABlasts2List[[names(DNABlasts2List[n])]]
    UniqueLength <- x$qseqid %>% unique() %>% length()
    DataSizes <- add_row(DataSizes, Genome = names(DNABlasts2List[n]), Rows = x %>% nrow(), Columns = x %>% ncol(), UniqueSeqIds = UniqueLength)
}
DataSizes
#    Genome                   Rows Columns UniqueSeqIds
#    <chr>                   <int>   <int>        <int>
#  1 KC_109_vs_T-Chr        134354      27         4843
#  2 KC_109_vs_T-Plasmid_1    1888      27          145
#  3 KC_109_vs_T-Plasmid_2    1607      27          155
#  4 KC_109_vs_T-Plasmid_3      16      27            3
#  5 KC_109_vs_T-Plasmid_4      16      27            6
#  6 KC_109_vs_TM-Chr       134360      27         4843
#  7 KC_109_vs_TM-Plasmid_1   1885      27          145
#  8 KC_109_vs_TM-Plasmid_2   1549      27          155
#  9 KC_109_vs_TM-Plasmid_3     17      27            3
# 10 KC_109_vs_TM-Plasmid_4     16      27            6
# # ... with 50 more rows

# Save Sizes
FileOut <- paste0(FilePrefix, "-CDS_vs_Genomes-Blasts-Counts.tsv")
write_tsv(DataSizes, file.path(AnnotDirOut, FileOut), col_names = TRUE, na = "")


rm(n, DataSizes, UniqueLength, Name, Query, Subject, x, QAnnot, SAnnot)
rm(Component)

###########################################
## Get CDS vs genome DNA blast best hits ##
###########################################
# Select the best hit for each query
BestDNABlasts2List <- list() # will be converted to df later
for (n in 1:length(DNABlasts2List)) {
    # Get best match(es) for each query and subject id
    Name <- DNABlasts2List[n] %>% names() %>% .[1] # get name of comparison
    BestDNABlasts2List[[Name]] <- DNABlasts2List[[Name]] %>% FindBestBlastHits2()
    BestDNABlasts2List[[Name]] %<>% unique()
}

# Sizes of datasets
DataSizes <- tibble(Comparison = character(0), Rows = integer(0), Columns = integer(0), UniqueQseqIds = integer(0), UniqueSseqIds = integer(0))
for (n in 1:length(BestDNABlasts2List)) {
    Name <- BestDNABlasts2List[n] %>% names() %>% .[1]
    qseqidUniqueLength <- BestDNABlasts2List[[Name]]$qseqid %>% uniquelength()
    sseqidUniqueLength <- BestDNABlasts2List[[Name]]$sseqid %>% uniquelength()
    DataSizes <- add_row(DataSizes, Comparison = Name, Rows = BestDNABlasts2List[[Name]] %>% nrow(), Columns = BestDNABlasts2List[[Name]] %>% ncol(),
                         UniqueQseqIds = qseqidUniqueLength, UniqueSseqIds = sseqidUniqueLength)
}
DataSizes
#    Comparison              Rows Columns UniqueQseqIds UniqueSseqIds
#    <chr>                  <int>   <int>         <int>         <int>
#  1 KC_109_vs_T-Chr         4844      27          4843            79
#  2 KC_109_vs_T-Plasmid_1    145      27           145             3
#  3 KC_109_vs_T-Plasmid_2    151      27           150             7
#  4 KC_109_vs_T-Plasmid_3      3      27             3             1
#  5 KC_109_vs_T-Plasmid_4      6      27             6             1
#  6 KC_109_vs_TM-Chr        4844      27          4843            79
#  7 KC_109_vs_TM-Plasmid_1   145      27           145             4
#  8 KC_109_vs_TM-Plasmid_2   148      27           147             7
#  9 KC_109_vs_TM-Plasmid_3     3      27             3             1
# 10 KC_109_vs_TM-Plasmid_4     6      27             6             1
# # ... with 50 more rows

# Save Sizes
FileOut <- paste0(FilePrefix, "-CDS_vs_Genomes-Best_Blasts-Counts.tsv")
write_tsv(DataSizes, file.path(AnnotDirOut, FileOut), col_names = TRUE, na = "")

(DataSizes$Rows == DataSizes$UniqueQseqIds) %>% summary()
#    Mode    TRUE    NA's
# logical      15      45

rm(n, Name, qseqidUniqueLength, sseqidUniqueLength, DataSizes)

# Make one df from best hits list

BestDNABlasts2 <- BestDNABlasts2List[[1]]
for (n in 2:length(BestDNABlasts2List)) {
    BestDNABlasts2 %<>% rbind(BestDNABlasts2List[[n]])
}

BestDNABlasts2 %>% dim() # 61689    27
BestDNABlasts2 %>% unique() %>% dim() # 61689     27

BestDNABlasts2 %>% dplyr::select(qseqid, sseqid) %>% dim() # 61689     2
BestDNABlasts2 %>% dplyr::select(qseqid, sseqid) %>% unique() %>% dim() # 59389     2

BestDNABlasts2 %>% dplyr::select(qseqid, sgenome) %>% dim() # 61689     2
BestDNABlasts2 %>% dplyr::select(qseqid, sgenome) %>% unique() %>% dim() # 61674     2
# very surprised there are dups for qseqid/sseqid but not qseqid/sgenome - can there be duplicated contig names?

# Sort this out later
##
ContigNames <- vector("character")
# for (n in 1:length(BestDNABlasts2List)) {
#     G <- BestDNABlasts2List[[n]]
#     ContigNames <- c(ContigNames, BestDNABlasts2List[[n]]$sseqid %>% unique())
# }
#
# ContigNames %>% length() # 1086
# ContigNames %>% unique() %>% length() # 361
#
# t <- ContigNames %>% duplicated() %>% which()
# x <- ContigNames[t] %>% unique() %>% sort()
#
# BestDNABlasts2 %<>% ungroup()
# x <- BestDNABlasts2 %>% dplyr::select(qgenome, sgenome, component, sseqid) %>% unique()
#
# x$sseqid %>% length() # 1086
# x$sseqid %>% unique() %>% length() # 361
# t <- x$sseqid %>% duplicated %>% which()
# x$sseqid
# BestDNABlasts2$sseqid %>% duplicated() %>% summary()
#
# x <- BestDNABlasts2 %>% filter(qgenome == "KC_109", sgenome == "T", component == "Chr")
# x$qseqid %>% duplicated() %>% summary()
# t <- x$qseqid %>% duplicated() %>% which()
# x$qseqid[t]
# y <- x %>% filter(qseqid == "fig|590.17911.peg.259")
# y <- x %>% filter(qseqid == "fig|590.17911.peg.4677")
##

BestDNABlasts2$qfraction %>% summary()
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.4240  1.0000  1.0000  0.9998  1.0000  1.0000

(BestDNABlasts2$qfraction < 1) %>% summary()
#    Mode   FALSE    TRUE    NA's
# logical   61631      58       0

x <- BestDNABlasts2 %>% filter(qfraction < 1)
rm(n, Name2, df, x)

###############################################
## Group DNA vs genome blasts across genomes ##
###############################################
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

BestDNABlasts2 %<>% ungroup()
BestDNABlasts2 %<>% group_by(qgenome, sgenome, component, add = TRUE)

for (n in 1:nrow(DNABlastFiles2)) {

    Name <- DNABlastFiles2$Comparison[n]
    print(Name)
    Q <- DNABlastFiles2$Query[n]
    S <- DNABlastFiles2$Subject[n]
    C <- DNABlastFiles2$Component[n]

    Annot <- Annotations %>% filter(Genome == Q & Component == C)
    QueryFeatureIDs <- Annot$FeatureID %>% unique()
    QueryFeatureIDs %>% length()
    QueryFeatureIDs %>% tail()
    Annot %>% tail()

    Comps <- BestDNABlasts2 %>% filter(qgenome == Q, sgenome == S, component == C)
    Comps %>% tail()
    CompQueryIDs <- Comps$qseqid %>% unique()
    CompQueryIDs %>% length()

    AnnotInQuery <- (QueryFeatureIDs %in% CompQueryIDs) %>% sum(na.rm=TRUE)
    AnnotNinQuery <- (QueryFeatureIDs %nin% CompQueryIDs) %>% sum(na.rm=TRUE)

    IDCounts2 %<>% add_row(Comparison = Name, Component = C,
                           Query = Q, QueryCount = length(QueryFeatureIDs), CompQueryCount = length(CompQueryIDs),
                           AnnotInQuery = AnnotInQuery, AnnotNinQuery = AnnotNinQuery,
                           Subject = S)

    if(AnnotNinQuery > 0) {
        t <- (QueryFeatureIDs %nin% CompQueryIDs) %>% which()
        x <- Annot[t,]
        x$Comparison <- Name
        x$Query <- Q
        x$Subject <- paste0(S)
        x %<>% dplyr::select(Comparison, Query, Subject, everything())
        # y <- tibble(Comparison = c(rep(Name,nrow(x))), Query = c(rep(Q,nrow(x))), S = c(rep(Subject,nrow(x))))
        # x <- cbind(y, x) %>% as_data_frame()
        MissingCDS2 %<>% bind_rows(x)
    }
}
IDCounts2 %<>% dplyr::select(Comparison, Component, Query, Subject, everything()) %>% arrange(Component, Subject)

# Save CDS counts
FileOut <- paste0(FilePrefix,"-CDS_vs_Genomes-Counts.tsv")
write_tsv(IDCounts2, file.path(AnnotDirOut, FileOut), col_names = TRUE, na = "")

# Save missing CDS
FileOut <- paste0(FilePrefix,"-CDS_vs_Genomes-Missing.tsv")
write_tsv(MissingCDS2, file.path(AnnotDirOut, FileOut), col_names = TRUE, na = "")

# Extend feature table to add genome info
FeatureIDTable2 <- FeatureIDTable
for (G in Genomes$One) {
    Name <- paste0(G, "_vs_Genomes")
    FeatureIDTable2[,Name] <- ""
}

BestDNABlasts2 %<>% ungroup()
BestDNABlasts2 %<>% group_by(cdsid, qseqid, add = TRUE)
for (n in 1:NrGenomes){

    QueryGenome <- Genomes[n,1] %>% .[[1]] # CDS query genome name
    Col1 <- colnames(FeatureIDTable2) %in% QueryGenome %>% which() # col pos to get CDS ID
    Col2 <- colnames(FeatureIDTable2) %in% paste0(QueryGenome, "_vs_Genomes") %>% which() # col pos to put result

    for (m in 2:NrGenomes) {

        SubjectGenome <- Genomes[n,m] %>% .[[1]] # CDS subject genome name
        print(c(QueryGenome, SubjectGenome))
        Blasts <- BestDNABlasts2 %>% filter(qgenome == QueryGenome, sgenome == SubjectGenome)

        for (r in 1:nrow(FeatureIDTable2)) {

            CDSID <- FeatureIDTable2$cdsid[r] %>% .[[1]]

            # Missing feature ID, dont process
            if (is.na(FeatureIDTable2[r,Col1])) {
                FeatureIDTable2[r,Col2] = NA
                next
            }

            # Have id, so see if there is a blast result
            Q <- FeatureIDTable2[r,Col1] %>% .[[1]]
            #c <- Blasts %>% dplyr::filter(cdsid == CDSID, qseqid == Q) %>% nrow()
            c <- ((Blasts$cdsid == CDSID & Blasts$qseqid == Q) == TRUE) %>% which() %>% length()
            #c <- Blasts[,Blasts$cdsid == CDSID & Blasts$qseqid == Q] %>% nrow()
            #c <- Blasts %>% dplyr::filter(cdsid == CDSID) %>% nrow()

            if (c == 1) {
                if (m == 2) {
                    FeatureIDTable2[r,Col2] = "TRUE"
                } else {
                    FeatureIDTable2[r,Col2] = paste(FeatureIDTable2[r,Col2], "TRUE")
                }
            } else if (c == 0) {
                if (m == 2) {
                    FeatureIDTable2[r,Col2] = "FALSE"
                } else {
                    FeatureIDTable2[r,Col2] = paste(FeatureIDTable2[r,Col2], "FALSE")
                }
            } else {
                print("Error: multiple results when matching query id to best blast results")
                print(c(QueryGenome, SubjectGenome))
                print(FeatureIDTable2[r,])
                print(Blasts %>% dplyr::filter(cdsid == CDSID, qseqid == Q))
                stop()
            }
        }
    }
}

# Redo, making a single truth table for all genomes
All <- FeatureIDTable %>% dplyr::select(cdsid)
for (G in Genomes$One) {
    All[,G] <- NA_character_
} # create df with cds id, empty cols for each genome

for (r in 1:nrow(FeatureIDTable2)) {
    CDSID <- FeatureIDTable2$cdsid[r]

    for (n in 1:NrGenomes) {
        Name <- paste0(Genomes$One[n], "_vs_Genomes")
        if (!is.na(FeatureIDTable2[r,Name])) {
            x <- str_split(FeatureIDTable2[r,Name], " ") %>% unlist() %>% as.vector() # genomes matched as vector
            pos <- c(1:4) %>% .[-(n)] # genomes reported in x

            x2 <- vector("logical", NrGenomes) # store
            c <- 1
            for (m in pos) # put truth in to new vector, in genome numerical order
            {
                x2[m] <- x[c]
                c <- c + 1
            }
            x2[n] <- TRUE # add current genome
            All[All$cdsid == CDSID,][2:(NrGenomes+1)] <- x2 # add to All
            break # do not need to do more cols from FeatureIDTable2
        }
    }
}

All %<>% unite("All", Genomes$One, sep = " ")
FeatureIDTable2 %<>% left_join(All, by = "cdsid")
FeatureIDTable2$cdsid %>% length() # 5221
FeatureIDTable2$cdsid %>% unique() %>% length() # 5221

x <- FeatureIDTable2$All %>% as.matrix() %>% as.vector() %>% as.factor() %>% table() %>% as_data_frame()
colnames(x) <- c("Presence", "Count")
x %<>% arrange(desc(Presence))
#   Presence               Count
#   <chr>                  <int>
# 1 TRUE TRUE TRUE TRUE     5213
# 2 TRUE TRUE FALSE TRUE       3
# 3 TRUE FALSE FALSE FALSE     5 # Most CDS present in other genomes

# y <- FeatureIDTable2 %>% filter(All != "TRUE TRUE TRUE TRUE")
# y <- FeatureIDTable2 %>% filter(All == "TRUE TRUE TRUE TRUE")
# y$KC_109 %>% length()
# y <- FeatureIDTable2 %>% filter(is.na(TMu))

# Save CDS vs genomes truth table
FileOut <- paste0(FilePrefix,"-CDS_vs_Genomes-Truth_Table.tsv")
write_tsv(FeatureIDTable2, file.path(AnnotDirOut, FileOut), col_names = TRUE, na = "")

# Save CDS vs genomes truth table summary
FileOut <- paste0(FilePrefix,"-CDS_vs_Genomes-Truth_Table_Summary.tsv")
write_tsv(x, file.path(AnnotDirOut, FileOut), col_names = TRUE, na = "")

rm(n, t, x, AnnotInQuery, AnnotNinQuery, C, CompQueryIDs, Annot,
   FileOut, Name, Names, Q, QueryFeatureIDs, S, c, NotNAs)
rm(G, m, r, Component, pos)
rm(y, CDSID, Col1, Col2, QueryGenome, SubjectGenome, s, x2)

#######################################################
## Find changes between AA and DNA sequences ##
#######################################################
# Compare AA sequences across genomes
# When they differ, compare DNA differences

#Load contigs
FastaDirIn <- "/Users/rtearle/Documents/Roseworthy_Projects/6.Students/Nitish_Salmonella/Holt_Pipeline_Feb_2018/Out/Bandage_Out"

FastAs <- tibble(Genome = character(0), Component = character(0), ContigID = character(0), Seq = character(0), Length = integer(0))
for (G in Genomes$One){
    for (C in Components){
        # Input file
        In <- file.path(FastaDirIn, G, paste0(tolower(C), ".fasta"))
        if (!file.exists(In)) {
            print(paste0("Cannot find input file \'",basename(In), "\'"))
            stop()}

        s = readDNAStringSet(In)

        for (n in 1:length(s)) {
            FastAs %<>% add_row(Genome = G, Component = C, ContigID = names(s[n]), Seq = toString(s[n]), Length = width(s[n]))
        }
    }
}

# Modifications to BestDNABlasts2
BestDNABlasts2 %<>% mutate(comparison = paste0(qgenome,"_vs_", sgenome, "-", component))
BestDNABlasts2 %<>% dplyr::select(cdsid, comparison, qgenome, sgenome, component, everything())
BestDNABlasts2 %<>% ungroup()
BestDNABlasts2 %<>% group_by(comparison, qseqid, add = TRUE)

# Fix order of subject start and end
BestDNABlasts2 %<>% mutate(newsstart = ifelse(sstart < send, sstart, send))
BestDNABlasts2 %<>% mutate(newssend = ifelse(sstart < send, send, sstart))
BestDNABlasts2 %<>% mutate(sstart = newsstart)
BestDNABlasts2 %<>% mutate(send = newssend)
BestDNABlasts2 %<>% dplyr::select(-c(newsstart, newssend))

#Use <- BestDNABlasts2 %>% select(cdsid, comparison, qseqid, sseqid, sstart, send, sstrand, pident) %>% filter(pident < 100)
Use <- BestDNABlasts2 %>% select(cdsid, comparison, qseqid, sseqid, sstart, send, sstrand, pident)
Use %<>% group_by(cdsid, comparison, qseqid, add = TRUE)

# Df for variant data
Variants <- tibble(CDSID = character(0), Comp = character(0), Query = character(0), QueryID = character(0),
                   Subject  = character(0), SubjectID = character(0), SubjectContig = character(0),
                   EndOfContig = logical(0), SubjectDNA = character(0), ProteinVars = character(0), DNAVars = character(0))

Start <- Sys.time()

# New
for (n in 1:nrow(CDSComparisons)) { # row of the feature table

    if (n %% 1000 == 0) {print(c(n, Sys.time()-Start))}

    CDSID <- CDSComparisons$CDSID[n]
    C <- CDSComparisons$Component[n]

    for (g in 1:NrGenomes) { # query genome
        Query <- Genomes[g,1][[1]]

        Qid <- CDSComparisons[n, paste0(Query, "_ID")][[1]]
        if(is.na(Qid)) {next} # do not process if no query id

        Qaa <- CDSComparisons[n, paste0(Query,"_AA")][[1]]

        for (g2 in 2:NrGenomes){ # subject genomes in turn

            Subject <- Genomes[g,g2][[1]]

            Sid <- CDSComparisons[n, paste0(Subject, "_ID")][[1]]
            if(is.na(Sid)) {next} # do not process if no subject id

            Saa <- CDSComparisons[n, paste0(Subject,"_AA")][[1]]

            Comp <- paste0(Query,"_vs_",Subject, "-", C)
            Sctg <- NA
            Snt <- NA

            ProteinVars <- vector()

            # Compare protein sequences for CDS vs CDS, extracting differences
            if( Qaa != Saa ) { # aa seqs different - rare, so set vars after knowing this
                # Labels, lengths, sequences for Subject

                Qlabel <- paste(CDSID, Query, sep = " ")
                Qaalen <- CDSComparisons[n, paste0(Query,"_AA_Length")][[1]]
                Slabel <- paste(CDSID, Subject, sep = " ")
                Saalen <- CDSComparisons[n, paste0(Subject,"_AA_Length")][[1]]

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

            # Concatenate
            if (length(ProteinVars) > 0) {
                ProteinVars %<>% paste(collapse = "|")
            } else {
                ProteinVars = ""
            }

            # DNA variants extracted by comparing CDS seq to matching seq from other genome
            Blast <- Use %>% filter(cdsid == CDSID, comparison == Comp, qseqid == Qid)

            DNAVars <- vector()

            if(nrow(Blast) > 0) { # only generate protein variants if there is a blast result

                Qnt <- CDSComparisons[n, paste0(Query,"_NT")][[1]]

                # Get substring of contig sequence
                Sctg <- Blast$sseqid[[1]]
                Snt <- substr(FastAs %>% filter(Genome == Subject, Component == C, ContigID == Sctg|NA) %>% .$Seq, Blast$sstart[[1]], Blast$send[[1]])
                Sntlen <- nchar(Snt)
                if (Sntlen == 0) {
                    sprintf("Cannot find sequence by position in subject genome %s matching %s query %s from genome %s", Subject, C, CDSID, Query)
                    stop()
                }

                # Complement if necessary
                if (Blast$sstrand[1] == "minus") {
                    Snt <- chartr("ATGC", "TACG", Snt)
                    Snt <- stringi::stri_reverse(Snt)
                }

                # Compare nt seqs
                if( Qnt != Snt ) {
                    Qntlen <- CDSComparisons[n, paste0(Query,"_NT_Length")][[1]]
                    Sntlen <- Snt %>% nchar()

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
            EndOfContig <- Annotations %>% filter(FeatureID == Qid) %>% .$EndOfContig %>% .[[1]]

            Variants %<>% add_row(CDSID = CDSID, Comp = Comp, Query = Query, QueryID = Qid,
                                  Subject = Subject, SubjectID = Sid, SubjectContig = Sctg, EndOfContig = EndOfContig, SubjectDNA = Snt,
                                  ProteinVars = ProteinVars, DNAVars = DNAVars)
        }
    }
}

rm(n, Query, Qid, CDSID, Subject, Sid, Comp, ProteinVars, Qlabel, Slabel, Qaalen, Saalen, Qaa, Saa, Qnt, Snt, Qntlen, Sntlen,
    Mat, Mat2, CurrType, CurrSeq, CurrSeq2, DNAVars, Sctg, Path, c)
rm(EndOfContig, g, g2, base, base2, C, first_time)

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
write_tsv(Annotations, file.path(AnnotDirOut, FileOut), col_names = TRUE, na = "")

# Save feature table plus comparisons
CDSComparisons %>% dim() # 5221   27
FeatureIDTable2 %>% dim() # 5221   11
CDSComparisons2 <- left_join(CDSComparisons, FeatureIDTable2 %>% dplyr::select(-c(one_of(Genomes$One), Component)), by = c("CDSID" = "cdsid"))
CDSComparisons2 %>% dim() # 5221   32
CDSComparisons2 %>% colnames()

FileOut <- paste0(FilePrefix,"-CDS-Multigenome_Annotations.tsv")
write_tsv(CDSComparisons2, file.path(AnnotDirOut, FileOut), col_names = TRUE, na = "")


# Subset of BestDNABlasts with col name changes
BestDNABlasts %<>% ungroup()
BDB <- BestDNABlasts %>% dplyr::select(-c(qseqidnr, sseqidnr, qfunction, sfunction))
BDB %<>% dplyr::select(cdsid, everything())
colnames(BDB)
Names <- c("CDSID", "QueryGenome", "SubjectGenome", "Component", "QueryID", "SubjectID", "CDSIdentity", "CDSCompLength", "CDSMismatch", "CDSGaps",
           "CDSQueryStart", "CDSQueryEnd", "CDSSubjectStart", "CDSSubjectEnd", "CDSEValue", "CDSBitScore", "CDSStrand",
           "CDSQueryLength", "CDSSubjectLength", "CDSQueryFullLength", "CDSSubjectFullLength",
           "CDSQueryFraction", "CDSSubjectFraction", "CDSFractionProduct") # for BestDNABlasts
colnames(BDB) <- Names
BDB %<>% mutate(Comparison = paste0(QueryGenome, "_vs_", SubjectGenome, "-", Component))
BDB %<>% dplyr::select(Comparison, everything())
BDB %>% dim() # 61215    25

# Subset and split of BestDNABlasts2 with col name changes
BestDNABlasts2 %<>% ungroup()
BDB2 <- BestDNABlasts2 %>% dplyr::select(-c(qseqidnr, sseqidnr, sfunction, sfraction))
colnames(BDB2)
Names <- c("CDSID", "Comparison", "QueryGenome", "SubjectGenome", "Component", "QueryID", "SubjectContig", "GenIdentity", "GenCompLength", "GenMismatch", "GenGaps",
           "GenQueryStart", "GenQueryEnd", "GenSubjectStart", "GenSubjectEnd", "GenEValue", "GenBitScore", "GenStrand",
           "Function", "GenQueryLength", "GenSubjectLength", "GenQueryFullLength", "GenSubjectFullLength",
           "GenQueryFraction") # for BestDNABlasts2
colnames(BDB2) <- Names
#BDB2 %<>% mutate(Comparison = paste0(QueryGenome, "_vs_", SubjectGenome, "-", Component))
BDB2 %<>% dplyr::select(Comparison, CDSID, everything())
BDB2 %>% dim() # 61689    24

# Variants with col name changes
colnames(Variants) <- c("CDSID", "Comparison", "QueryGenome", "QueryID", "SubjectGenome", "SubjectID", "SubjectContig", "EndOfContig", "SubjectDNA", "ProteinVars", "DNAVars")

# New df with key cols
CDSBlastsVariants <- BDB2 %>% dplyr::select(Comparison, CDSID, QueryGenome, SubjectGenome, Component, QueryID)
CDSBlastsVariants %>% dim() # 61689     6

CDSBlastsVariants %<>% left_join(BDB %>% dplyr::select(-c(QueryGenome, SubjectGenome, Component)), by = c("Comparison", "CDSID", "QueryID"))
CDSBlastsVariants %>% dim() # 61689    25
colnames(CDSBlastsVariants)

CDSBlastsVariants %<>% left_join(BDB2 %>% dplyr::select(-c(QueryGenome, SubjectGenome, Component)), by = c("Comparison", "CDSID", "QueryID"))
CDSBlastsVariants %>% dim() # 61689    43
colnames(CDSBlastsVariants)

CDSBlastsVariants %<>% left_join(Variants %>% dplyr::select(-c(QueryGenome, SubjectGenome, SubjectID, SubjectContig)), by = c("Comparison", "CDSID", "QueryID"))
CDSBlastsVariants %>% dim() # 61689    47
colnames(CDSBlastsVariants)

CDSBlastsVariants %<>% dplyr::select(c(Comparison, CDSID, QueryGenome, SubjectGenome, Component, Function, QueryID, SubjectID), everything())

FileOut <- paste0(FilePrefix,"-CDS_Blasts_Variants.tsv")
write_tsv(CDSBlastsVariants, file.path(AnnotDirOut, FileOut), col_names = TRUE, na = "")

rm(BDB, BDB2, FileOut, Names)
#rm(BestDNABlasts, BestDNABlasts2, Variants, FeatureIDTable)

# Save DNABlasts
Name <- names(DNABlastsList[1])
DNABlasts <- matrix(nrow = 0, ncol = DNABlastsList[[Name]] %>% ncol()) %>% as_data_frame()
colnames(DNABlasts) <- colnames(DNABlastsList[[Name]])

for (n in 1:length(DNABlastsList)) {
    Name <- DNABlastsList[n] %>% names()
    DNABlasts <- bind_rows(DNABlasts, DNABlastsList[[Name]])
}

DNABlasts %>% dim() # 1703581      27
x <- DNABlasts %>% filter(qs > 0.9)
x %>% dplyr::select(qseqid, sseqid) %>% dim() # 62116
x %>% dplyr::select(qseqid, sseqid) %>% unique() %>% dim() # 62116

FileOut <- paste0(FilePrefix,"-CDS_Blasts.tsv")
write_tsv(DNABlasts, file.path(AnnotDirOut, FileOut), col_names = TRUE, na = "")

FileOut <- paste0(FilePrefix,"-CDS_Best_Blasts.tsv")
write_tsv(BestDNABlasts, file.path(AnnotDirOut, FileOut), col_names = TRUE, na = "")

# Save DNABlasts2
Name <- names(DNABlasts2List[1])
DNABlasts2 <- matrix(nrow = 0, ncol = DNABlasts2List[[Name]] %>% ncol()) %>% as_data_frame()
colnames(DNABlasts2) <- colnames(DNABlasts2List[[Name]])

for (n in 1:length(DNABlasts2List)) {
    Name <- DNABlasts2List[n] %>% names()
    DNABlasts2 <- bind_rows(DNABlasts2, DNABlasts2List[[Name]])
}

DNABlasts2 %>% dim() # 1650315      27

FileOut <- paste0(FilePrefix,"-CDS_vs_Genome_Blasts.tsv")
write_tsv(DNABlasts2, file.path(AnnotDirOut, FileOut), col_names = TRUE, na = "")

FileOut <- paste0(FilePrefix,"-CDS_vs_Genome_Best_Blasts.tsv")
write_tsv(BestDNABlasts2, file.path(AnnotDirOut, FileOut), col_names = TRUE, na = "")

rm(n, Name)










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
            dplyr::select(SubjectContig, Start, End, CDSID, CDSIdentity, GenStrand, Start2, End2, RGB2, QueryID) %>%
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
