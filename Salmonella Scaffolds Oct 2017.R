library(ggplot2)
library(iterators)
library(lubridate)
library(magrittr)
library(RColorBrewer)
library(readr)
library(readxl)
library(reshape2)
library(scales)
library(statmod)
library(stringr)
library(tibble)
library(tools)
library(VennDiagram)
library(car)
library(corrplot)
library(dplyr)

###############
## Functions ##
###############

'%nin%' <- Negate('%in%')

####################
## Load scaffolds ##
####################
Scaffolds <- list()

# KC_109
FileIn <- "/Users/rtearle/Documents/Roseworthy_Projects/Students/Talia_Salmonella/Salmonella_Oct_2017/Out/TMu-Scaffold_Summary.tsv"
Scaffolds$KC_109 <- read_tsv(FileIn, col_names = FALSE) %>% as_data_frame()
colnames(Scaffolds$KC_109) = c("Name", "Node", "Length", "Coverage")
Scaffolds$KC_109 %>% dim() # 345 x 4

FileIn <- "/Users/rtearle/Documents/Roseworthy_Projects/Students/Talia_Salmonella/Salmonella_Oct_2017/Out/T-Scaffold_Summary.tsv"
Scaffolds$T <- read_tsv(FileIn, col_names = FALSE) %>% as_data_frame()
colnames(Scaffolds$T) = c("Name", "Node", "Length", "Coverage")
Scaffolds$T %>% dim() # 177   4

FileIn <- "/Users/rtearle/Documents/Roseworthy_Projects/Students/Talia_Salmonella/Salmonella_Oct_2017/Out/TM-Scaffold_Summary.tsv"
Scaffolds$TM <- read_tsv(FileIn, col_names = FALSE) %>% as_data_frame()
colnames(Scaffolds$TM) = c("Name", "Node", "Length", "Coverage")
Scaffolds$TM %>% dim() # 251   4

FileIn <- "/Users/rtearle/Documents/Roseworthy_Projects/Students/Talia_Salmonella/Salmonella_Oct_2017/Out/KC_109-Scaffold_Summary.tsv"
Scaffolds$TMu <- read_tsv(FileIn, col_names = FALSE) %>% as_data_frame()
colnames(Scaffolds$TMu) = c("Name", "Node", "Length", "Coverage")
Scaffolds$TMu %>% dim() # 340   4
for (n in 1:length(Scaffolds)) {
    #print(Scaffolds[[n]])
    Scaffolds[[n]] %>% filter(Coverage < 10) %>% nrow() %T>% print()


}

##############################
## Load Blast vs Salmonella ##
##############################
# Data is from bandage, copied from Blast table
SalmonellaBlast <- list()
# KC_109
FileIn <- "/Users/rtearle/Documents/Roseworthy_Projects/Students/Talia_Salmonella/Salmonella_Oct_2017/KC_109-blast_vs_LT2_SJTUF10584.tsv"
SalmonellaBlast$KC_109 <- read_tsv(FileIn, col_names = FALSE) %>% as_data_frame()
colnames(SalmonellaBlast$KC_109) = c("Query", "Node", "Pct_Identity", "Length", "Pct_Cover", "Mismatches", "Gap_Open", "Query_Start", "Query_End", "Node_Start", "Node_End", "E_Value", "Bit_Score")

SalmonellaBlast$KC_109 %<>%
    mutate(Pct_Cover = gsub("%", "", Pct_Cover)) %>% # remove % from Cover
    mutate(Strand = gsub("[0-9]+", "", Node)) %>% # new col with strand
    mutate(Node = gsub("[+-]", "", Node)) %>% # remove strand from Node
    select(Query, Node, Strand, everything()) # reorder

SalmonellaBlast$KC_109$Node %>% unique()
SalmonellaBlast$KC_109$Node %<>% as.integer()
SalmonellaBlast$KC_109$Pct_Cover %<>% as.numeric()
SalmonellaBlast$KC_109$Length %<>% as.integer()
SalmonellaBlast$KC_109$Query %<>% as.factor()
SalmonellaBlast$KC_109

SalmonellaBlast$KC_109 %>% dim() # 2654   14

# Add a column with the gap to the next entry for each query and node
SalmonellaBlast$KC_109$QueryGap <- NA_real_
SalmonellaBlast$KC_109$NodeGap <- NA_real_

SalmonellaBlast$KC_109 %<>% arrange(Query, Query_Start)
for(n in 1:(nrow(SalmonellaBlast$KC_109)-1)){
    # Query gap
    if(SalmonellaBlast$KC_109$Query[n] == SalmonellaBlast$KC_109$Query[n+1]) {
        SalmonellaBlast$KC_109$QueryGap[n] = SalmonellaBlast$KC_109$Query_Start[n+1] - SalmonellaBlast$KC_109$Query_End[n]
    }
}

SalmonellaBlast$KC_109 %<>% arrange(Node, Node_Start)
for(n in 1:(nrow(SalmonellaBlast$KC_109)-1)){
    # Node gap
    if(!is.na(SalmonellaBlast$KC_109$Node[n]) & !is.na(SalmonellaBlast$KC_109$Node[n+1]) & SalmonellaBlast$KC_109$Node[n] == SalmonellaBlast$KC_109$Node[n+1]) {
        SalmonellaBlast$KC_109$NodeGap[n] = SalmonellaBlast$KC_109$Node_Start[n+1] - SalmonellaBlast$KC_109$Node_End[n]
    }
}

# Generate a bed file
DirOut <- "/Users/rtearle/Documents/Roseworthy_Projects/Students/Talia_Salmonella/Salmonella_Oct_2017/Out"
FileOut <- "KC109_Nodes_vs_LT2_SJTUF10584.bed"
colnames(SalmonellaBlast$KC_109)
BedFile <- SalmonellaBlast$KC_109 %>% select(Query, Query_Start, Query_End, Node, Length, Strand)
BedFile$Query_Start <- BedFile$Query_Start - 1
Header <- 'track name="KC109 Nodes vs LT2_SJTUF10584" description="KC109 Nodes vs LT2_SJTUF10584"'
writeLines(Header, file.path(DirOut, FileOut))
write_tsv(BedFile, file.path(DirOut, FileOut), append = TRUE)

##############################
## Load Blast vs Salmonella ##
##############################
# Data is from NCBI web blast
NTBlastTable <- list()

# KC_109
FileIn <- "/Users/rtearle/Documents/Roseworthy_Projects/Students/Talia_Salmonella/Salmonella_Oct_2017/Out/Blast/YUM2UUGH015-Alignment-2.txt"
NTBlastTable$KC_109 <- read_tsv(FileIn, col_names = FALSE, skip = 7) %>% as_data_frame()
NTBlastTable$KC_109 %>% dim()

colnames(NTBlastTable$KC_109) = c("QueryID", "SubjectID", "QueryAccessn", "SubjectAccessn", "Pct_Identity", "Length", "Mismatches", "Gap_Open", "QueryStart", "QueryEnd", "SubjectStart", "SubjectEnd", "EValue", "BitScore")
NTBlastTable$KC_109 %<>% filter(!grepl("^#", QueryID))
NTBlastTable$KC_109 %>% dim()

NTBlastText <- list()

#NTBlastText$KC_109 <- readLines(file.path(DirIn, File), skipNul = TRUE) %>%
FileIn <- "/Users/rtearle/Documents/Roseworthy_Projects/Students/Talia_Salmonella/Salmonella_Oct_2017/Out/Blast/YUM2UUGH015-Alignment.txt"
DF <- readLines(FileIn, skipNul = TRUE, n = 100000) %>%
    grep("(Query  )|(\\|{3,})|(Sbjct  )", ., invert = TRUE, value = TRUE) %>%
    gsub("^[0-9]+", "", .) %>%
    stringr::str_trim(side = "both") %>%
    as_data_frame()
DF %<>% filter(value != "")

DF %>% nrow()
DF2 <- DF[8:nrow(DF),]
DF2 %>% nrow()



