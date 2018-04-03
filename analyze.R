library(magrittr)
library(effsize)
library(ggplot2)
library(tidyr)
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

# # print all existing ASEG variable names
# for (s in ASEG_structures) {
#     cat("A$`", s, "_Volume_mm3` +\n", sep = "")
# }
# for (g in ASEG_globals) {
#     cat("A$", g, " +\n", sep = "")
# }

# Lateral Ventricles
A$LatVentricles <- 
    A$`Left-Lateral-Ventricle_Volume_mm3` +
    A$`Left-Inf-Lat-Vent_Volume_mm3` +
    A$`Right-Lateral-Ventricle_Volume_mm3` +
    A$`Right-Inf-Lat-Vent_Volume_mm3`

# Ventricles
A$Ventricles <-
    A$LatVentricles +
    A$`3rd-Ventricle_Volume_mm3` +
    A$`4th-Ventricle_Volume_mm3` +
    A$`5th-Ventricle_Volume_mm3`

# Fluid
A$Fluid <- 
    A$Ventricles +
    A$`Left-choroid-plexus_Volume_mm3` +
    A$`Right-choroid-plexus_Volume_mm3` +
    A$CSF_Volume_mm3

# Pallidum
A$Pallidum <-
    A$`Left-Pallidum_Volume_mm3` +
    A$`Right-Pallidum_Volume_mm3`

# Amygdala
A$Amygdala <-
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
    A$LatVentricles -
    A$`Left-choroid-plexus_Volume_mm3` -
    A$`Right-choroid-plexus_Volume_mm3`

# estimated Total Brain Volume
# based on eTIV
A$eTBV <- A$IntraCranialVol - A$Fluid

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


# generate a table with results
Report <- function(variable,
                   wilcox_paired = FALSE,
                   d_paired = FALSE,
                   alpha = 0.05,
                   digits = 4) {
    
    x <- B[ , match(variable, names(B))]
    
    w <- wilcox.test(x ~ B$pairClass, paired = wilcox_paired)
    d <- cohen.d(x ~ B$pairClass, paired = d_paired)
    
    q_control <- quantile(x[B$pairClass == "control"])
    q_case <- quantile(x[B$pairClass == "case"])
    
    c(
        "Variable" = variable,
        
        "Wilcox_p.value" = w$p.value %>% round(digits),
        "Wilcox_significant" = w$p.value < alpha,
        
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

# difference from control
RelDifference <- function(x) {
    x[B$pairClass == "case"] %>%
        `-`(x[B$pairClass == "control"]) %>%
        `/`(x[B$pairClass == "case"] + x[B$pairClass == "control"]) %>%
        `*`(200) %>%
        return()
}

AbsDifference <- function(x) {
    x[B$pairClass == "case"] %>%
        `-`(x[B$pairClass == "control"]) %>%
        return()
}


# create df with deltas ---------------------------------------------------

D <- B[B$pairClass == "control", match("AGE", names(B))] %>% data.frame()
names(D)[1] <- "AGE_control"

# verify
all(D$AGE_control == B$AGE[B$pairClass == "control"])

# add more columns to D
# basic metrics
D$AGE_case <- B$AGE[B$pairClass == "case"]
D$SEX <- B$SEX[B$pairClass == "control"]

# Fluid
D$LatVentricles_Rel <- RelDifference(B$LatVentricles)
D$LatVentricles_Abs <- AbsDifference(B$LatVentricles)

D$Ventricles_Rel <- RelDifference(B$Ventricles)
D$Ventricles_Abs <- AbsDifference(B$Ventricles)

D$Fluid_Rel <- RelDifference(B$Fluid)
D$Fluid_Abs <- AbsDifference(B$Fluid)

# Total Volulmes
D$pTBV_Rel <- RelDifference(B$pTBV)
D$pTBV_Abs <- AbsDifference(B$pTBV)

D$eTBV_Rel <- RelDifference(B$eTBV)
D$eTBV_Abs <- AbsDifference(B$eTBV)

D$eTIV_Rel <- RelDifference(B$eTIV)
D$eTIV_Abs <- AbsDifference(B$eTIV)

# Pallidum
D$Pallidum_Rel <- RelDifference(B$Pallidum)
D$Pallidum_Abs <- AbsDifference(B$Pallidum)

D$R_Pallidum_Rel <- RelDifference(B$`Right-Pallidum_Volume_mm3`)
D$R_Pallidum_Abs <- AbsDifference(B$`Right-Pallidum_Volume_mm3`)

D$L_Pallidum_Rel <- RelDifference(B$`Left-Pallidum_Volume_mm3`)
D$L_Pallidum_Abs <- AbsDifference(B$`Left-Pallidum_Volume_mm3`)

# Amygdala
D$Amygdala_Rel <- RelDifference(B$Amygdala)
D$Amygdala_Abs <- AbsDifference(B$Amygdala)


# report ------------------------------------------------------------------    

# wilcox.test and cohen.d paired
wp = TRUE
dp = FALSE

# generate a report table
c(
    "LatVentricles",
    "Ventricles",
    "Fluid",
    "pTBV",
    "eTBV",
    "eTIV",
    "Pallidum",
    "Right-Pallidum_Volume_mm3",
    "Left-Pallidum_Volume_mm3",
    "Amygdala"
) %>% 
    sapply(Report, wilcox_paired = wp, d_paired = dp, USE.NAMES = FALSE) %>%
    t -> report
View(report)
write.csv(report, file = "results/report.csv", row.names = FALSE)
# read.csv("results/report.csv") %>% View("report.csv")


# plots -------------------------------------------------------------------

# boxplot with volumes
# select which variables to plot
C <- gather(B,
            value = "Volume_mm3",
            key = "Structure",
            LatVentricles,
            Ventricles,
            Fluid,
            pTBV,
            eTBV,
            eTIV,
            Pallidum,
            `Right-Pallidum_Volume_mm3`,
            `Left-Pallidum_Volume_mm3`,
            Amygdala)

# define order of Structures on plot
C$Structure %<>% factor(levels = c(
    "LatVentricles",
    "Ventricles",
    "Fluid",
    "pTBV",
    "eTBV",
    "eTIV",
    "Pallidum",
    "Right-Pallidum_Volume_mm3",
    "Left-Pallidum_Volume_mm3",
    "Amygdala"
))

boxplot_volumes <- ggplot(C,
                          aes(x = pairClass,
                              y = Volume_mm3)) +
    geom_boxplot() +
    facet_wrap( ~ Structure, scales = "free", ncol = 3) +
    scale_x_discrete(limits = rev(levels(B$pairClass)))
boxplot_volumes

ggsave("results/boxplot_volumes.pdf", width = 6, height = 12, scale = 2)

# sex by color
boxplot_volumes_sex <- boxplot_volumes + aes(color = SEX)
boxplot_volumes_sex
ggsave("results/boxplot_volumes_sex.pdf", width = 6, height = 12, scale = 2)

# boxplot with relative deltas
# select variables to plot
Dr <- gather(D,
             value = "Relative_Difference",
             key = "Structure",
             LatVentricles_Rel,
             Ventricles_Rel,
             Fluid_Rel,
             pTBV_Rel,
             eTBV_Rel,
             eTIV_Rel,
             Pallidum_Rel,
             R_Pallidum_Rel,
             L_Pallidum_Rel,
             Amygdala_Rel)

# define order of Structures
Dr$Structure %<>% factor(levels = c(
    "LatVentricles_Rel",
    "Ventricles_Rel",
    "Fluid_Rel",
    "pTBV_Rel",
    "eTBV_Rel",
    "eTIV_Rel",
    "Pallidum_Rel",
    "R_Pallidum_Rel",
    "L_Pallidum_Rel",
    "Amygdala_Rel"
))


boxplot_reldiff <- ggplot(Dr,
                          aes(x = Structure,
                              y = Relative_Difference)) +
    geom_boxplot() +
    scale_y_continuous(name = "Relative Difference as % of mean\n(case - control) / (case + control) * 200")
boxplot_reldiff

ggsave("results/boxplot_reldiff.pdf", width = 6, height = 3, scale = 3)

# sex by color
boxplot_reldiff_sex <- boxplot_reldiff + aes(color = SEX)
boxplot_reldiff_sex
ggsave("results/boxplot_reldiff_sex.pdf", width = 6, height = 3, scale = 3)

# Fluid_Rel vs AGE
plot_fluidrel_age <- ggplot(D,
                            aes(x = AGE_control,
                                y = Fluid_Rel)) +
    geom_point() +
    geom_smooth(span = 1) +
    scale_x_continuous(name = "Age") +
    scale_y_continuous(name = "Fluid Volume: Relative Difference as % of mean\n(case - control) / (case + control) * 200")
plot_fluidrel_age

ggsave("results/plot_fluidrel_age.pdf", width = 6, height = 3, scale = 3)


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
