# read in corpus callosum areas
cc <- read.csv("abide.cc.area.csv")

# read in brain volume
bvol <- read.csv("abide.brain.vol.csv",h=F)
names(bvol)[2] <- "bvol.mm"

# read in phenotypic data
demo <- read.csv("Phenotypic_V1_0b.csv")

# read in QA assessments
qa <- read.csv("abide_scan_rating.csv")

# make common id column across groups so we can merge cc area, brain vol and phenotypic data
for (i in 1:dim(cc)[1]) cc$id[i] <- strsplit(as.character(cc$ID[i]),split = "/")[[1]][3]
for (i in 1:dim(bvol)[1]) bvol$id[i] <- strsplit(as.character(bvol$V1[i]),split = "/")[[1]][3]
for (i in 1:dim(demo)[1]) demo$id[i] <- formatC(x = demo$SUB_ID[i], width = 7, format = "d", flag = "0")
for (i in 1:dim(qa)[1]) qa$id[i] <- formatC(x = qa$sub.id[i], width = 7, format = "d", flag = "0")

# merge these three
cc.bvol <- merge(cc,bvol, by = "id")
cc.bvol.demo <- merge(cc.bvol,demo, by = "id")
cc.bvol.demo.qa <- merge(cc.bvol.demo, qa, by = "id")

# make diagnosis and dsm category columns
cc.bvol.demo.qa$diagnosis <- factor(cc.bvol.demo.qa$DX_GROUP, labels = c("autism","control"))

# convert -9999 entries to NA
for (i in 1:dim(cc.bvol.demo.qa)[1]) if (cc.bvol.demo.qa$DSM_IV_TR[i] == -9999) cc.bvol.demo.qa$dsm.iv[i] <- NA else cc.bvol.demo.qa$dsm.iv[i] <- cc.bvol.demo.qa$DSM_IV_TR[i]
cc.bvol.demo.qa$dsm.category <- factor(cc.bvol.demo.qa$dsm.iv, labels = c("control","autism","aspergers","pdd.nos","aspergers.or.nos"))

# make data frame that is only scans rated 1 or 2
cc.bvol.demo.qa.no.poor.quality <- cc.bvol.demo.qa[ !((cc.bvol.demo.qa$scan.quality == 1) | (cc.bvol.demo.qa$scan.quality == 2)),]

# make a simpler name
abide <- cc.bvol.demo.qa.no.poor.quality
rm(cc.bvol.demo.qa.no.poor.quality)

# total CC ASD vs controls
quad.fit <- lm( CC_area ~ diagnosis + poly(AGE_AT_SCAN,2) + diagnosis:poly(AGE_AT_SCAN,2) + as.factor(SEX) + bvol.mm + diagnosis:bvol.mm + SITE_ID, data = abide)

# brain volume in ASD
intracranial.vol.asd <- lm(bvol.mm ~ diagnosis + poly(AGE_AT_SCAN,2) + diagnosis:poly(AGE_AT_SCAN,2) + as.factor(SEX) + as.factor(SITE_ID), data = abide)

# brain volume in subcategories

intracranial.vol.asd.subcat <- lm(bvol.mm ~ dsm.category + poly(AGE_AT_SCAN,2) + dsm.category:poly(AGE_AT_SCAN,2) + as.factor(SEX) + as.factor(SITE_ID), data = abide)

# Corpus callosum cohen's d
cohens.d <- (mean(abide$CC_area[ abide$diagnosis == "autism" ])-mean(abide$CC_area[ abide$diagnosis == "control" ]))/sd(abide$CC_area)

# total CC ASD subtypes vs controls
quad.fit.subcat <- lm( CC_area ~ dsm.category + poly(AGE_AT_SCAN,2) + dsm.category:poly(AGE_AT_SCAN,2) + SEX + bvol.mm + SITE_ID, data = abide)

# Witelson subregions
w1.fit.asd <- lm( W1 ~ diagnosis + poly(AGE_AT_SCAN,2) + diagnosis:poly(AGE_AT_SCAN,2) + SEX + bvol.mm + diagnosis:bvol.mm + SITE_ID, data = abide)
w1.fit.subtypes <- lm( W1 ~ dsm.category + poly(AGE_AT_SCAN,2) + dsm.category:poly(AGE_AT_SCAN,2) + SEX + bvol.mm + SITE_ID, data = abide)

w2.fit.asd <- lm( W2 ~ diagnosis + poly(AGE_AT_SCAN,2) + diagnosis:poly(AGE_AT_SCAN,2) + SEX + bvol.mm + diagnosis:bvol.mm + SITE_ID, data = abide)
w2.fit.subtypes <- lm( W2 ~ dsm.category + poly(AGE_AT_SCAN,2) + dsm.category:poly(AGE_AT_SCAN,2) + SEX + bvol.mm + SITE_ID, data = abide)

w3.fit.asd <- lm( W3 ~ diagnosis + poly(AGE_AT_SCAN,2) + diagnosis:poly(AGE_AT_SCAN,2) + SEX + bvol.mm + diagnosis:bvol.mm + SITE_ID, data = abide)
w3.fit.subtypes <- lm( W3 ~ dsm.category + poly(AGE_AT_SCAN,2) + dsm.category:poly(AGE_AT_SCAN,2) + SEX + bvol.mm + SITE_ID, data = abide)

w4.fit.asd <- lm( W4 ~ diagnosis + poly(AGE_AT_SCAN,2) + diagnosis:poly(AGE_AT_SCAN,2) + SEX + bvol.mm + diagnosis:bvol.mm + SITE_ID, data = abide)
w4.fit.subtypes <- lm( W4 ~ dsm.category + poly(AGE_AT_SCAN,2) + dsm.category:poly(AGE_AT_SCAN,2) + SEX + bvol.mm + SITE_ID, data = abide)

w5.fit.asd <- lm( W5 ~ diagnosis + poly(AGE_AT_SCAN,2) + diagnosis:poly(AGE_AT_SCAN,2) + SEX + bvol.mm + diagnosis:bvol.mm + SITE_ID, data = abide)
w5.fit.subtypes <- lm( W5 ~ dsm.category + poly(AGE_AT_SCAN,2) + dsm.category:poly(AGE_AT_SCAN,2) + SEX + bvol.mm + SITE_ID, data = abide)

w6.fit.asd <- lm( W6 ~ diagnosis + poly(AGE_AT_SCAN,2) + diagnosis:poly(AGE_AT_SCAN,2) + SEX + bvol.mm + diagnosis:bvol.mm + SITE_ID, data = abide)
w6.fit.subtypes <- lm( W6 ~ dsm.category + poly(AGE_AT_SCAN,2) + dsm.category:poly(AGE_AT_SCAN,2) + SEX + bvol.mm + SITE_ID, data = abide)

w7.fit.asd <- lm( W7 ~ diagnosis + poly(AGE_AT_SCAN,2) + diagnosis:poly(AGE_AT_SCAN,2) + SEX + bvol.mm + diagnosis:bvol.mm + SITE_ID, data = abide)
w7.fit.subtypes <- lm( W7 ~ dsm.category + poly(AGE_AT_SCAN,2) + dsm.category:poly(AGE_AT_SCAN,2) + SEX + bvol.mm + SITE_ID, data = abide)

# make table summarizing main p vals
# with layout like
# total CC W1 W2... W7
# ASD
# Autism
# Aspergers
# PDD
# PDD- NOS

cc.overall.p <- c(summary(quad.fit)$coefficients[2,4],summary(quad.fit.subcat)$coefficients[2,4],summary(quad.fit.subcat)$coefficients[3,4],summary(quad.fit.subcat)$coefficients[4,4],summary(quad.fit.subcat)$coefficients[5,4])
w1.p <- c(summary(w1.fit.asd)$coefficients[2,4],summary(w1.fit.subtypes)$coefficients[2,4],summary(w1.fit.subtypes)$coefficients[3,4],summary(w1.fit.subtypes)$coefficients[4,4],summary(w1.fit.subtypes)$coefficients[5,4])
w2.p <- c(summary(w2.fit.asd)$coefficients[2,4],summary(w2.fit.subtypes)$coefficients[2,4],summary(w2.fit.subtypes)$coefficients[3,4],summary(w2.fit.subtypes)$coefficients[4,4],summary(w2.fit.subtypes)$coefficients[5,4])
w3.p <- c(summary(w3.fit.asd)$coefficients[2,4],summary(w3.fit.subtypes)$coefficients[2,4],summary(w3.fit.subtypes)$coefficients[3,4],summary(w3.fit.subtypes)$coefficients[4,4],summary(w3.fit.subtypes)$coefficients[5,4])
w4.p <- c(summary(w4.fit.asd)$coefficients[2,4],summary(w4.fit.subtypes)$coefficients[2,4],summary(w4.fit.subtypes)$coefficients[3,4],summary(w4.fit.subtypes)$coefficients[4,4],summary(w4.fit.subtypes)$coefficients[5,4])
w5.p <- c(summary(w5.fit.asd)$coefficients[2,4],summary(w5.fit.subtypes)$coefficients[2,4],summary(w5.fit.subtypes)$coefficients[3,4],summary(w5.fit.subtypes)$coefficients[4,4],summary(w5.fit.subtypes)$coefficients[5,4])
w6.p <- c(summary(w6.fit.asd)$coefficients[2,4],summary(w6.fit.subtypes)$coefficients[2,4],summary(w6.fit.subtypes)$coefficients[3,4],summary(w6.fit.subtypes)$coefficients[4,4],summary(w6.fit.subtypes)$coefficients[5,4])
w7.p <- c(summary(w7.fit.asd)$coefficients[2,4],summary(w7.fit.subtypes)$coefficients[2,4],summary(w7.fit.subtypes)$coefficients[3,4],summary(w7.fit.subtypes)$coefficients[4,4],summary(w7.fit.subtypes)$coefficients[5,4])

p.val.summary <- data.frame(cc.overall.p,w1.p,w2.p,w3.p,w4.p,w5.p,w6.p,w7.p)
row.names(p.val.summary) <- c("asd","autism","aspergers","pdd.nos","aspergers.or.nos")
