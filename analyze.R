library(magrittr)
library(effsize)

#### import merged data ####
A <- read.csv("data/A.csv", check.names = FALSE)

ASEG_structures <- unlist(read.table("data/ASEG_structures.txt"),
                          use.names = FALSE)
ASEG_measures <- unlist(read.table("data/ASEG_measures.txt"),
                          use.names = FALSE)
ASEG_globals <- unlist(read.table("data/ASEG_globals.txt"),
                          use.names = FALSE)

APARC_structures <- unlist(read.table("data/APARC_structures.txt"),
                          use.names = FALSE)
APARC_measures <- unlist(read.table("data/APARC_measures.txt"),
                           use.names = FALSE)
APARC_globals <- unlist(read.table("data/APARC_globals.txt"),
                           use.names = FALSE)

# add summed structures as variables to A
# BrainSegNotVent - same as BrainSeg without ventricles (lateral, inferior lateral, 3rd, 4th, 5th), CSF, and choroid plexus.
# Ventricles
A$Ventricles_Volume_mm3 <-
  A$`Left-Lateral-Ventricle_Volume_mm3` +
  A$`Left-Inf-Lat-Vent_Volume_mm3` +
  A$`Right-Lateral-Ventricle_Volume_mm3` +
  A$`Right-Inf-Lat-Vent_Volume_mm3` +
  A$`3rd-Ventricle_Volume_mm3` +
  A$`4th-Ventricle_Volume_mm3` +
  A$`5th-Ventricle_Volume_mm3`

# Fluid
A$Fluid_Volume_mm3 <- 
  A$Ventricles_Volume_mm3 +
  A$`Left-choroid-plexus_Volume_mm3` +
  A$`Right-choroid-plexus_Volume_mm3` +
  A$CSF_Volume_mm3

# Pallidum
A$Pallidum_Volume_mm3 <-
  A$`Left-Pallidum_Volume_mm3` +
  A$`Right-Pallidum_Volume_mm3`


# sort by pairNumber, exclude not matched
B <- A[order(A$pairNumber, A$pairClass), ]
B <- B[B$matched, ]


# verify proper sorting
all(B$pairNumber[B$pairClass == "control"] ==
        B$pairNumber[B$pairClass == "case"])


#### functions ####
# retrieve Aseg and Aparc data
GetAsegValue <- function(sub_id, struct, meas){
    return(A[match(sub_id, A$SUB_ID),
             match(paste0(struct, "_", meas), names(A))])
}

GetAparcValue <- function(sub_id, hemi, struct, meas){
    return(A[match(sub_id, A$SUB_ID),
             match(paste0(struct, "_", hemi, "_", meas), names(A))])
}


# case-vs-control wrapper for wilcox.test
wilcox.cvc <- function(x) {
    return(wilcox.test(x[B$pairClass == "case"],
                       x[B$pairClass == "control"],
                       paired = TRUE))
}

# print case-vs-control differences: wilcox.test, cohen.d, mean ± sd
PrintDifference <- function(x, paired = FALSE) {
    w <- wilcox.test(x ~ B$pairClass, paired = paired)$p.value %>% round(4)
    d <- cohen.d(x ~ B$pairClass, paired = paired)$estimate %>% round(4)
    m0 <- mean(x[B$pairClass == "control"]) %>% round()
    m1 <- mean(x[B$pairClass == "case"]) %>% round()
    s0 <- sd(x[B$pairClass == "control"]) %>% round()
    s1 <- sd(x[B$pairClass == "case"]) %>% round()
    
    cat("Wilcox:", w, "|",
        "Cohen's d:", d, "|",
        "control:", m0, "±", s0, "|",
        "case:", m1, "±", s1, "\n")
}

#### interesting things ####
# Ventricles/CSF
PrintDifference(B$Ventricles_Volume_mm3)
PrintDifference(B$Fluid_Volume_mm3)

# Pallidum
PrintDifference(B$Pallidum_Volume_mm3)
PrintDifference(B$`Right-Pallidum_Volume_mm3`)
PrintDifference(B$`Left-Pallidum_Volume_mm3`)


#### CC area ####
wilcox.cvc(B$CC_area)

# difference realive to control
B$CC_area[B$pairClass == "case"] %>%
    `-`(B$CC_area[B$pairClass == "control"]) %>%
    `/`(B$CC_area[B$pairClass == "control"]) -> CC_area_RelDiff
mean(CC_area_RelDiff)
sd(CC_area_RelDiff)
# plot(B$AGE[B$pairClass == "control"], CC_area_RelDiff)



#### CC area ~ volume of brain structures ####
## volume only - individual structures
for (structure in ASEG_structures) {
    control_all <- GetAsegValue(B$SUB_ID[B$pairClass == "control"],
                                struct = structure,
                                meas = "Volume_mm3")
    case_all <- GetAsegValue(B$SUB_ID[B$pairClass == "case"],
                             struct = structure,
                             meas = "Volume_mm3")
    
    # remove pairs where one value is NA
    control <- control_all[!is.na(control_all) & !is.na(case_all)]
    case <- case_all[!is.na(control_all) & !is.na(case_all)]
    
    # TODO cannot compute exact p-value with zeroes
    test <- wilcox.test(control, case, paired = TRUE)$p.value
    
    if (test < 0.05) {
        cat("control:\t", mean(control), "\t",
            "case:\t", mean(case), "\t",
            test, "\t", structure, "\n")
    } else {
        cat("nonsignificant\t", test, "\t", structure, "\n")
    }
}

## volume only - globals
for (global in ASEG_globals) {
    control <- B[ , match(global, names(B))][B$pairClass == "control"]
    case <- B[ , match(global, names(B))][B$pairClass == "case"]
    
    test <- wilcox.test(control, case, paired = TRUE)$p.value
    
    if (test < 0.05) {
        cat("control:\t", mean(control), "\t",
            "case:\t", mean(case), "\t",
            test, "\t", global, "\n")
    } else {
        cat("nonsignificant\t", test, "\t", global, "\n")
    }
}

## relative to CC area - individual structures
for (structure in ASEG_structures) {
    control_all <- GetAsegValue(B$SUB_ID[B$pairClass == "control"],
                                struct = structure,
                                meas = "Volume_mm3")
    case_all <- GetAsegValue(B$SUB_ID[B$pairClass == "case"],
                             struct = structure,
                             meas = "Volume_mm3")
    
    # relative to CC area
    control_all <- control_all / B$CC_area[B$pairClass == "control"]
    case_all <- case_all / B$CC_area[B$pairClass == "case"]
    
    # remove pairs where one value is NA
    control_CC <- control_all[!is.na(control_all) & !is.na(case_all)]
    case_CC <- case_all[!is.na(control_all) & !is.na(case_all)]
    
    # TODO cannot compute exact p-value with zeroes
    test <- wilcox.test(control_CC, case_CC, paired = TRUE)$p.value
    
    if (test < 0.05) {
        cat("control:\t", mean(control_CC), "\t",
            "case:\t", mean(case_CC), "\t",
            test, "\t", structure, "\n")
    } else {
        cat("nonsignificant\t", test, "\t", structure, "\n")
    }
}

## relative to CC area - globals
for (global in ASEG_globals) {
    control <- B[ , match(global, names(B))][B$pairClass == "control"]
    case <- B[ , match(global, names(B))][B$pairClass == "case"]
    
    # relative to CC area
    control_CC <- control / B$CC_area[B$pairClass == "control"]
    case_CC <- case / B$CC_area[B$pairClass == "case"]
    
    test <- wilcox.test(control_CC, case_CC, paired = TRUE)$p.value
    
    if (test < 0.05) {
        cat("control:\t", mean(control_CC), "\t",
            "case:\t", mean(case_CC), "\t",
            test, "\t", global, "\n")
    } else {
        cat("nonsignificant\t", test, "\t", global, "\n")
    }
}




#### CC area ~ surface area of cortical structures ####
## surface only - individual structures
for (structure in APARC_structures) {
    for (hemisphere in c("rh", "lh")) {
        control_all <- GetAparcValue(B$SUB_ID[B$pairClass == "control"],
                                     hemi = hemisphere,
                                     struct = structure,
                                     meas = "SurfArea")
        case_all <- GetAparcValue(B$SUB_ID[B$pairClass == "case"],
                                  hemi = hemisphere,
                                  struct = structure,
                                  meas = "SurfArea")
        
        # remove pairs where one value is NA
        control <- control_all[!is.na(control_all) & !is.na(case_all)]
        case <- case_all[!is.na(control_all) & !is.na(case_all)]
        
        test <- wilcox.test(control, case, paired = TRUE)$p.value
        
        if (test < 0.05) {
            cat("control:\t", mean(control), "\t",
                "case:\t", mean(case), "\t",
                test, "\t", structure, "_", hemisphere, "\n",
                sep = "")
        } else {
            cat("nonsignificant\t", test, "\t", structure, "_", hemisphere,"\n",
                sep = "")
        }
    }
}

## surface only - globals
WhiteSurfArea_total <- B$WhiteSurfArea_rh + B$WhiteSurfArea_lh
wilcox.cvc(WhiteSurfArea_total)$p.value
wilcox.cvc(B$WhiteSurfArea_rh)$p.value
wilcox.cvc(B$WhiteSurfArea_lh)$p.value

## relative to CC area - individual structures
for (structure in APARC_structures) {
    for (hemisphere in c("rh", "lh")) {
        control_all <- GetAparcValue(B$SUB_ID[B$pairClass == "control"],
                                     hemi = hemisphere,
                                     struct = structure,
                                     meas = "SurfArea")
        case_all <- GetAparcValue(B$SUB_ID[B$pairClass == "case"],
                                  hemi = hemisphere,
                                  struct = structure,
                                  meas = "SurfArea")
        
        # relative to CC area
        control_all <- control_all / B$CC_area[B$pairClass == "control"]
        case_all <- case_all / B$CC_area[B$pairClass == "case"]
        
        # remove pairs where one value is NA
        control_CC <- control_all[!is.na(control_all) & !is.na(case_all)]
        case_CC <- case_all[!is.na(control_all) & !is.na(case_all)]
        
        test <- wilcox.test(control_CC, case_CC, paired = TRUE)$p.value
        
        if (test < 0.05) {
            cat("control:\t", mean(control_CC), "\t",
                "case:\t", mean(case_CC), "\t",
                test, "\t", structure, "_", hemisphere, "\n",
                sep = "")
        } else {
            cat("nonsignificant\t", test, "\t", structure, "_", hemisphere,"\n",
                sep = "")
        }
    }
}

## relative to CC area - globals
wilcox.cvc(WhiteSurfArea_total / B$CC_area)$p.value
wilcox.cvc(B$WhiteSurfArea_rh / B$CC_area)$p.value
wilcox.cvc(B$WhiteSurfArea_lh / B$CC_area)$p.value
