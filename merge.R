library(stringr)
library(magrittr)

#### import ABIDE data: area and phenotype ####
# downloaded from https://www.nitrc.org/frs/download.php/7721/abide.cc.measurements.20150118.zip
AB_area <- read.csv("data/abide/abide.cc.area.csv")
AB_phenotype <- read.csv("data/abide/Phenotypic_V1_0b.csv")

# in AB_area create SUB_ID variable with patient's number extracted from ID
AB_area$SUB_ID <- as.numeric(str_extract(AB_area$ID, "\\d{7}"))
# merge AB_area and AB_phenotype by matching SUB_ID values
ABIDE <- merge(AB_area, AB_phenotype, by = "SUB_ID")
rm(AB_area, AB_phenotype)

# export/import
write.csv(ABIDE, "data/ABIDE.csv", row.names = FALSE)
# ABIDE <- read.csv("data/ABIDE.csv")


#### get PCP FILE_ID list for downloading ####
# saves FILE_ID values to /data/pcp/FILE_ID.txt
# this file then can be used as input for download.sh script to download PCP data
read.csv("data/pcp/summary.csv") %$%
    FILE_ID %>%
    # ignore patients with no file
    setdiff("no_filename") %>%
    as.character %>%
    write("./data/pcp/FILE_ID.txt")


#### convert PCP FS aseg.stats to aseg.csv ####
aseg.stats_colnames <- c("Index", "SegId", "NVoxels", "Volume_mm3",
                         "StructName", "normMean", "normStdDev", "normMin",
                         "normMax", "normRange")

aseg.stats_path <- "./data/pcp/aseg.stats"
files <- list.files(path = aseg.stats_path)
for(file in files) {
    file_noext <- str_extract(file, ".+\\.aseg")
    write.csv(
        read.table(paste0(aseg.stats_path, "/", file),
                   col.names = aseg.stats_colnames),
        file = paste0("./data/pcp/aseg.csv/", file_noext, ".csv"),
        row.names = FALSE
    )
}

# clean up
rm(aseg.stats_path, file_noext, file, files, aseg.stats_colnames)




#### convert PCP FS aparc.stats to aparc.csv ####
aparc.stats_colnames <- c("StructName","NumVert","SurfArea","GrayVol",
                          "ThickAvg","ThickStd","MeanCurv","GausCurv",
                          "FoldInd","CurvInd")
for(hemi in c("rh", "lh")) {
    aparc.stats_path <- paste0("./data/pcp/", hemi, ".aparc.stats")
    files <- list.files(path = aparc.stats_path)
    for(file in files) {
        file_noext <- str_extract(file, ".+\\.aparc")
        write.csv(read.table(paste0(aparc.stats_path, "/", file),
                             col.names = aparc.stats_colnames),
                  file = paste0("./data/pcp/", hemi, ".aparc.csv/", file_noext, ".csv"),
                  row.names = FALSE)
    }
}
# clean up
rm(aparc.stats_path, file_noext, hemi, file, files)





#### import PCP FS ASEG data - Volume_mm3, normMean, normStdDev ####
aseg.csv_path <- "./data/pcp/aseg.csv"
files <- list.files(path = aseg.csv_path)

# columns in aseg.csv to be imported
# 4, 6, 7 correspond to Volume_mm3, normMean, normStdDev
ASEG_cols <- c(4, 6, 7)

# read first file to create varaible names for ASEG
# !!! assuming first file has all structures
file1_table <- read.csv(paste0(aseg.csv_path, "/", files[1]))

ASEG_StructName <- file1_table$StructName

ASEG_StructName_measure <- NULL
for (meas in colnames(file1_table)[ASEG_cols]) {
    for (struct in ASEG_StructName) {
        ASEG_StructName_measure <- c(ASEG_StructName_measure,
                                     paste0(struct, "_", meas))
    }
}

# initialize empty ASEG df
ASEG <- data.frame(matrix(NA,
                          nrow = length(files),
                          ncol = 1 + length(ASEG_StructName_measure)))
names(ASEG) <- c("SUB_ID", ASEG_StructName_measure)

# add data to ASEG df
i <- 1
for (file in files) {
    file_table <- read.csv(paste0(aseg.csv_path, "/", file))
    file_id <- as.numeric(str_extract(file, "\\d{7}"))
    ASEG$SUB_ID[i] <- file_id
    for (row in 1:nrow(file_table)) {
        for (col in ASEG_cols) {
            ASEG[[match(paste0(file_table$StructName[row], "_",
                               colnames(file_table)[col]),
                        names(ASEG))]][i] <- 
                file_table[row, col]
        }
    }
    i <- i + 1
}
# clean up
rm (file1_table, file, files, file_id, file_table, i, row, col, aseg.csv_path,
    ASEG_StructName_measure, ASEG_cols, meas, struct, ASEG_StructName)


#### append global volumes to ASEG table - using asegstats2table output ####
# load global values and store in ASEG_GLOBALS keeping colnames
read.table("./data/pcp/asegstats2table/abide_as2t_vol.txt",
           header = TRUE)[ ,47:56] -> ASEG_GLOBALS

# load SUB_ID values from list.txt and append to ASEG_GLOBALS
read.table("./data/pcp/asegstats2table/list.txt")[[1]] %>%
    str_extract("\\d{7}") %>%
    as.numeric -> ASEG_GLOBALS$SUB_ID

# clean up
ASEG <- merge(ASEG, ASEG_GLOBALS, by = "SUB_ID")
ABIDE.ASEG <- merge(ABIDE, ASEG, by = "SUB_ID")
rm(ASEG_GLOBALS)

# export/import
write.csv(ASEG, "data/ASEG.csv", row.names = FALSE)
write.csv(ABIDE.ASEG, "data/ABIDE.ASEG.csv", row.names = FALSE)
# ASEG <- read.csv("data/ASEG.csv")
# ABIDE.ASEG <- read.csv("data/ABIDE.ASEG.csv")





#### import PCP FS APARC data - all measures ####
# create variable names based on first file
# !!! assuming first file has all structures
file1_table <- read.csv("./data/pcp/rh.aparc.csv/Caltech_0051456.rh.aparc.csv")

APARC_StructName <- file1_table$StructName

APARC_StructName_hemi_meas <- NULL
for (meas in colnames(file1_table)[-1]) {
    for (struct in APARC_StructName) {
        for (hemi in c("rh", "lh")) {
            APARC_StructName_hemi_meas <- c(
                APARC_StructName_hemi_meas,
                paste0(struct, "_", hemi, "_", meas)
            )
        }
    }
}

# initialize empty APARC df
APARC <- data.frame(matrix(NA,
                           nrow = 1035,
                           ncol = 1 + length(APARC_StructName_hemi_meas)))
names(APARC) <- c("SUB_ID", APARC_StructName_hemi_meas)

# add data to APARC df
# suspiciously slow              
for (hemi in c("rh", "lh")) {
    aparc.csv_path <- paste0("./data/pcp/", hemi, ".aparc.csv")
    files <- list.files(path = aparc.csv_path)
    
    i <- 1
    for (file in files) {
        file_table <- read.csv(paste0(aparc.csv_path, "/", file))
        file_id <- as.numeric(str_extract(file, "\\d{7}"))
        
        # only set SUB_ID when doing RH
        # confirm that order of subjects is the same in RH and LH
        if (hemi == "rh") {
            APARC$SUB_ID[i] <- file_id
        } else if (APARC$SUB_ID[i] != file_id) {
            print("Error: different order of subjects in RH and LH!")
        }
        
        for (row in 1:nrow(file_table)) {
            for (col in 2:ncol(file_table)) {
                APARC[[match(paste0(file_table$StructName[row],"_",
                                  hemi, "_",
                                  colnames(file_table)[col]),
                           names(APARC))]][i] <- file_table[row, col]
            }
        }
        i <- i + 1
    }
}

# clean up
rm(file1_table, meas, struct, hemi, aparc.csv_path, files, file,
   file_table, file_id, i, col, row, APARC_StructName_hemi_meas,
   APARC_StructName)


#### append global vars to APARC table - using extract_aparc_globals.sh output ####
# rh
read.table("data/pcp/rh.aparc.SUB_ID.txt")[[1]] %>%
    as.numeric %>%
    data.frame -> APARC_GLOBALS_RH

names(APARC_GLOBALS_RH)[1] <- "SUB_ID"

# VARIABLE RENAME
read.table("data/pcp/rh.aparc.NumVert.txt")[[1]] %>%
    as.numeric -> APARC_GLOBALS_RH$NumVert_rh

# VARIABLE RENAME
read.table("data/pcp/rh.aparc.WhiteSurfArea.txt")[[1]] %>%
    as.numeric -> APARC_GLOBALS_RH$WhiteSurfArea_rh

# VARIABLE RENAME
read.table("data/pcp/rh.aparc.MeanThickness.txt")[[1]] %>%
    as.numeric -> APARC_GLOBALS_RH$MeanThickness_rh

# lh
read.table("data/pcp/lh.aparc.SUB_ID.txt")[[1]] %>%
    as.numeric %>%
    data.frame -> APARC_GLOBALS_LH

names(APARC_GLOBALS_LH)[1] <- "SUB_ID"

# VARIABLE RENAME
read.table("data/pcp/lh.aparc.NumVert.txt")[[1]] %>%
    as.numeric -> APARC_GLOBALS_LH$NumVert_lh

# VARIABLE RENAME
read.table("data/pcp/lh.aparc.WhiteSurfArea.txt")[[1]] %>%
    as.numeric -> APARC_GLOBALS_LH$WhiteSurfArea_lh

# VARIABLE RENAME
read.table("data/pcp/lh.aparc.MeanThickness.txt")[[1]] %>%
    as.numeric -> APARC_GLOBALS_LH$MeanThickness_lh

# clean up
APARC_GLOBALS <- merge(APARC_GLOBALS_RH, APARC_GLOBALS_LH, by = "SUB_ID")
APARC <- merge(APARC, APARC_GLOBALS, by = "SUB_ID")
ABIDE.ASEG.APARC <- merge(ABIDE.ASEG, APARC, by = "SUB_ID")

rm(APARC_GLOBALS_RH, APARC_GLOBALS_LH, APARC_GLOBALS)

# export/import
write.csv(APARC, "data/APARC.csv", row.names = FALSE)
write.csv(ABIDE.ASEG.APARC, "data/ABIDE.ASEG.APARC.csv", row.names = FALSE)
# APARC <- read.csv("data/APARC.csv")
# ABIDE.ASEG.APARC <- read.csv("data/ABIDE.ASEG.APARC.csv")



#### case-control matching ####
A <- ABIDE.ASEG.APARC

# more understandable variables
# VARIABLE RENAME
pairClass <- factor(ifelse(A$DX_GROUP == 1, "case", "control"))
A$SEX <- factor(ifelse(A$SEX == 1, "M", "F"))
names(A)[match("AGE_AT_SCAN", names(A))] <- "AGE"

# A$SEX %<>% '==' 1 %>% ifelse("M", "F") %>% factor
# pairClass <- A$DX_GROUP %>% '==' 1 %>% ifelse("case", "control") %>% factor

# new variables for matching. Unmatched pairNumber will be set to NA later
pairNumber <- rep(0, nrow(A))
matchType <- rep(NA, nrow(A))
matched <- rep(FALSE, nrow(A))

# using data.frame() rather than cbind() changes Aseg var names ("-" to "\.", "^\d" to "^X\d")
A <- cbind(A, pairNumber, pairClass, matched, matchType)
rm(pairNumber, pairClass, matched, matchType)

# maximum allowed age difference for Close Match and Distant Match
age_CM <- 0.5
age_DM <- 1

# optionally merge different samples from same sites as same sites
# VARIABLE RENAME
if (TRUE) {
    SITE <- as.character(A$SITE_ID)
    SITE <- ifelse(xor(SITE == "LEUVEN_1", SITE == "LEUVEN_2"),
                   "LEUVEN",
                   SITE)
    SITE <- ifelse(xor(SITE == "UCLA_1", SITE == "UCLA_2"),
                   "UCLA",
                   SITE)
    SITE <- ifelse(xor(SITE == "UM_1", SITE == "UM_2"),
                   "UM",
                   SITE)
    SITE <- factor(SITE)
    A <- cbind(A, SITE)
    rm(SITE)
} else {
    SITE <- A$SITE_ID
    A <- cbind(A, SITE)
    rm(SITE)
}

# 1. matching pass: Gender + Close Match age + Site
for (i in 1:nrow(A)) {
    if (A$pairClass[i] == "case" & !A$matched[i]) {
        # list of SUB_ID: at same site & of same sex & not matched & controls
        control_id <-
            A$SUB_ID[A$SITE == A$SITE[i] & A$SEX == A$SEX[i] &
                     !A$matched & A$pairClass == "control"]
        
        # TODO rest of the loop could be a function repeated in 2. and 3. pass
        
        # only continue if any potential matches found
        if (length(control_id) > 0) {
            control_agediff <-
                abs(A$AGE[i] - A$AGE[A$SUB_ID %in% control_id])
            
            # only ontinue if the best match is a Close Match
            if (min(control_agediff) <= age_CM) {
                control_id <- control_id[which.min(control_agediff)]
                
                # make a pair
                A$matched[i] <- TRUE
                A$matched[A$SUB_ID == control_id] <- TRUE
                A$pairNumber[i] <- max(A$pairNumber) + 1
                A$pairNumber[A$SUB_ID == control_id] <- A$pairNumber[i]
                A$matchType[i] <- "GCMS"
                A$matchType[A$SUB_ID == control_id] <- "GCMS"
            }
        }
    }
}

# 2. matching pass: Gender + Close Match age
for (i in 1:nrow(A)) {
    if (A$pairClass[i] == "case" & !A$matched[i]) {
        # list of SUB_ID: of same sex & not matched & controls
        control_id <- 
            A$SUB_ID[A$SEX == A$SEX[i] & !A$matched &
                         A$pairClass == "control"]
        
        # only continue if any potential matches found
        if (length(control_id) > 0) {
            control_agediff <-
                abs(A$AGE[i] - A$AGE[A$SUB_ID %in% control_id])
            
            # only ontinue if the best match is a Close Match
            if (min(control_agediff) <= age_CM) {
                control_id <- control_id[which.min(control_agediff)]
                
                # make a pair
                A$matched[i] <- TRUE
                A$matched[A$SUB_ID == control_id] <- TRUE
                A$pairNumber[i] <- max(A$pairNumber) + 1
                A$pairNumber[A$SUB_ID == control_id] <- A$pairNumber[i]
                A$matchType[i] <- "GCM"
                A$matchType[A$SUB_ID == control_id] <- "GCM"
            }
        }
    }
}

# 3. matching pass: Gender + Distant Match age
for (i in 1:nrow(A)) {
    if (A$pairClass[i] == "case" & !A$matched[i]) {
        # list of SUB_ID: of same sex & not matched & controls
        control_id <- 
            A$SUB_ID[A$SEX == A$SEX[i] & !A$matched &
                         A$pairClass == "control"]
        
        # only continue if any potential matches found
        if (length(control_id) > 0) {
            control_agediff <-
                abs(A$AGE[i] - A$AGE[A$SUB_ID %in% control_id])
            
            # only ontinue if the best match is Distant Match
            if (min(control_agediff) <= age_DM) {
                control_id <- control_id[which.min(control_agediff)]
                
                # make a pair
                A$matched[i] <- TRUE
                A$matched[A$SUB_ID == control_id] <- TRUE
                A$pairNumber[i] <- max(A$pairNumber) + 1
                A$pairNumber[A$SUB_ID == control_id] <- A$pairNumber[i]
                A$matchType[i] <- "GDM"
                A$matchType[A$SUB_ID == control_id] <- "GDM"
            }
        }
    }
}

# clean up
rm(control_id, control_agediff, i)
A$pairNumber <- ifelse(A$pairNumber == 0, NA, A$pairNumber)
A$matchType <- factor(A$matchType,
                      levels = c("GCMS", "GCM", "GDM"),
                      ordered = TRUE)

# check how many matched, how many pairs, summary of matchType
sum(A$matched)
max(A$pairNumber, na.rm = TRUE)
summary(A$matchType)

# export/import
write.csv(A, "data/A.csv", row.names = FALSE)
# A <- read.csv("data/A.csv")

# check mean age difference in pairs
pairAgeDifference <- rep(NA, max(A$pairNumber, na.rm = TRUE))
for (i in 1:length(pairAgeDifference)) {
    pairAgeDifference[i] <-
        max(A$AGE[A$pairClass == "case" & A$pairNumber == i],
            na.rm = TRUE) -
        max(A$AGE[A$pairClass == "control" & A$pairNumber == i],
            na.rm = TRUE)
}
rm(i)
mean(pairAgeDifference, na.rm = TRUE)
sd(pairAgeDifference, na.rm = TRUE)
rm(pairAgeDifference)

#### miscellaneous ####
# functions to retrieve Aseg and Aparc data
GetAsegValue <- function(sub_id, struct, meas){
    return(A[match(sub_id, A$SUB_ID),
             match(paste0(struct, "_", meas), names(A))])
}

GetAparcValue <- function(sub_id, hemi, struct, meas){
    return(A[match(sub_id, A$SUB_ID),
             match(paste0(struct, "_", hemi, "_", meas), names(A))])
}


#### analysis ####
# sort by pairNumber, exclude not matched
B <- A[order(A$pairNumber, A$pairClass), ]
B <- B[B$matched, ]

# verify proper sorting
all(B$pairNumber[B$pairClass == "control"] ==
    B$pairNumber[B$pairClass == "case"])

# case-vs-control wrappers for wilcox.test and t.test
wilcox.cvc <- function(x) {
    return(wilcox.test(x[B$pairClass == "case"],
                       x[B$pairClass == "control"],
                       paired = TRUE))
}
t.cvc <- function(x) {
    return(t.test(x[B$pairClass == "case"],
                  x[B$pairClass == "control"],
                  paired = TRUE))
}

# compare WhiteSurfArea
WhiteSurfArea_total <- B$WhiteSurfArea_rh + B$WhiteSurfArea_lh
wilcox.cvc(WhiteSurfArea_total)
wilcox.cvc(B$WhiteSurfArea_rh)
wilcox.cvc(B$WhiteSurfArea_lh)

# compare WhiteSurfArea realtive to CC_area
wilcox.cvc(WhiteSurfArea_total / B$CC_area)
wilcox.cvc(B$WhiteSurfArea_rh / B$CC_area)
wilcox.cvc(B$WhiteSurfArea_lh / B$CC_area)

# compare MeanThickness
wilcox.cvc(B$MeanThickness_rh)
wilcox.cvc(B$MeanThickness_lh)

# compare TotalGrayVol
wilcox.cvc(B$TotalGrayVol)
wilcox.cvc(B$TotalGrayVol / B$CC_area)

