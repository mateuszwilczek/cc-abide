library(magrittr)
library(effsize)
library(ggplot2)
# library(nortest)


# import merged data ------------------------------------------------------

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


# create variables: summed structures -------------------------------------

# TODO simpler but less consistent names without "_Volume_mm3"?
# Lateral Ventricles
A$LatVentricles_Volume_mm3 <- 
    A$`Left-Lateral-Ventricle_Volume_mm3` +
    A$`Left-Inf-Lat-Vent_Volume_mm3` +
    A$`Right-Lateral-Ventricle_Volume_mm3` +
    A$`Right-Inf-Lat-Vent_Volume_mm3`

# Ventricles
A$Ventricles_Volume_mm3 <-
    A$LatVentricles_Volume_mm3 +
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

# Amygdala
A$Amygdala_Volume_mm3 <-
    A$`Left-Amygdala_Volume_mm3` +
    A$`Right-Amygdala_Volume_mm3`

# precise Total Brain Volume
# attempt to construct absent BrainSegNotVent from other values
A$pTBV <- 
    A$SupraTentorialVol +
    A$`Left-Cerebellum-White-Matter_Volume_mm3` +
    A$`Left-Cerebellum-Cortex_Volume_mm3` +
    A$`Right-Cerebellum-White-Matter_Volume_mm3` +
    A$`Right-Cerebellum-Cortex_Volume_mm3` +
    A$`Brain-Stem_Volume_mm3` -
    A$LatVentricles_Volume_mm3 -
    A$`Left-choroid-plexus_Volume_mm3` -
    A$`Right-choroid-plexus_Volume_mm3`
    
# for (s in ASEG_structures) {
#     cat("A$`", s, "_Volume_mm3` +\n", sep = "")
# }
# for (g in ASEG_globals) {
#     cat("A$", g, " +\n", sep = "")
# }

# estimated Total Brain Volume
# based on eTIV
A$eTBV <- A$IntraCranialVol - A$Fluid_Volume_mm3

# estimated Total Intracranial Volume
# new name for consistency with pTBV and eTBV
A$eTIV <- A$IntraCranialVol

# sort and exclude unmatched ----------------------------------------------

B <- A[order(A$pairNumber, A$pairClass), ]
B <- B[B$matched, ]

# verify proper sorting
all(B$pairNumber[B$pairClass == "control"] ==
        B$pairNumber[B$pairClass == "case"])



# functions ---------------------------------------------------------------

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

# print case-vs-control differences: wilcox.test, cohen.d, median [1Q - 3Q]
PrintDifference <- function(x, round = 0, wp = FALSE, dp = FALSE) {
    w <- wilcox.test(x ~ B$pairClass, paired = wp)$p.value %>% round(4)
    d <- cohen.d(x ~ B$pairClass, paired = dp)
    
    q_control <- quantile(x[B$pairClass == "control"]) %>% round(round)
    q_case <- quantile(x[B$pairClass == "case"]) %>% round(round)
    
    cat("Wilcox: \t", w, "\n",
        
        "d:      \t", round(d$estimate, 4),
        "\t[", round(d$conf.int[1], 4), " : ", round(d$conf.int[2], 4), "]\n",
        
        "control:\t", q_control[3],
        "\t[", q_control[2], " : ", q_control[4], "]\n",
        
        "case:   \t", q_case[3],
        "\t[", q_case[2], " : ", q_case[4], "]\n",
        
        sep = "")
}

# generate a table with results
GenerateReport <- function(variable,
                           wilcox_paired = FALSE,
                           d_paired = FALSE,
                           p = 0.05,
                           digits = 4) {
    
    x <- B[ , match(variable, names(B))]
    
    w <- wilcox.test(x ~ B$pairClass, paired = wilcox_paired)
    d <- cohen.d(x ~ B$pairClass, paired = d_paired)
    
    q_control <- quantile(x[B$pairClass == "control"])
    q_case <- quantile(x[B$pairClass == "case"])
    
    c(
        "Variable" = variable,
        
        "Wilcox_p.value" = w$p.value %>% round(digits),
        "Wilcox_significant" = w$p.value < p,
        
        "Cohens_d" = d$estimate[[1]] %>% round(digits),
        "Cohens_d_CI_inf" = d$conf.int[[1]] %>% round(digits),
        "Cohens_d_CI_sup" = d$conf.int[[2]] %>% round(digits),
        
        "control_median" = q_control[[3]],
        "control_1Q" = q_control[[2]],
        "control_3Q" = q_control[[4]],
        
        "case_median" = q_case[[3]],
        "case_1Q" = q_case[[2]],
        "case_3Q" = q_case[[4]]
    ) %>% return
}


# interesting things ------------------------------------------------------    

# wilcox.test and cohen.d paired
wp = TRUE
dp = FALSE

# generate a report table
c(
    "LatVentricles_Volume_mm3",
    "Ventricles_Volume_mm3",
    "Fluid_Volume_mm3",
    "pTBV",
    "eTBV",
    "eTIV",
    "Pallidum_Volume_mm3",
    "Right-Pallidum_Volume_mm3",
    "Left-Pallidum_Volume_mm3",
    "Amygdala_Volume_mm3"
) %>% 
    sapply(GenerateReport, wilcox_paired = wp, d_paired = dp) %>%
    t -> report
View(report)
write.csv(report, file = "results/report.csv", row.names = FALSE)
# read.csv("results/report.csv") %>% View("report")

# print report to console
# Ventricles/CSF
PrintDifference(B$LatVentricles_Volume_mm3, wp = wp, dp = dp)
PrintDifference(B$Ventricles_Volume_mm3, wp = wp, dp = dp)
PrintDifference(B$Fluid_Volume_mm3, wp = wp, dp = dp)

# Pallidum
PrintDifference(B$Pallidum_Volume_mm3, wp = wp, dp = dp)
PrintDifference(B$`Right-Pallidum_Volume_mm3`, wp = wp, dp = dp)
PrintDifference(B$`Left-Pallidum_Volume_mm3`, wp = wp, dp = dp)

# Amygdala
PrintDifference(B$Amygdala_Volume_mm3,
                wp = wp, dp = dp)


# plots -------------------------------------------------------------------

p1 <- ggplot(B, aes(x = pairClass, y = Fluid_Volume_mm3)) +
    geom_boxplot()
p1

# p2 <- ggplot(report, aes(x = Variable))

# fishing expedition below ------------------------------------------------

# CC area -----------------------------------------------------------------

wilcox.cvc(B$CC_area)

# difference realive to control
B$CC_area[B$pairClass == "case"] %>%
    `-`(B$CC_area[B$pairClass == "control"]) %>%
    `/`(B$CC_area[B$pairClass == "control"]) -> CC_area_RelDiff
mean(CC_area_RelDiff)
sd(CC_area_RelDiff)
# plot(B$AGE[B$pairClass == "control"], CC_area_RelDiff)



# CC area ~ volume of brain structures ------------------------------------

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




# CC area ~ surface area of cortical structures ---------------------------

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
