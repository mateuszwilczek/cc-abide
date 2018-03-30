library(stringr)
library(magrittr)

#### import merged data ####
# unlist(use.names = TRUE)

# sort by pairNumber, exclude not matched
B <- A[order(A$pairNumber, A$pairClass), ]
B <- B[B$matched, ]

# verify proper sorting
all(B$pairNumber[B$pairClass == "control"] ==
        B$pairNumber[B$pairClass == "case"])


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


#### CC area ~ volume of structures ####
#### CC area ~ curvature, surface area of cortical structures ####
#### random ####
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
