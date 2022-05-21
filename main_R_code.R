# destdir= "./data"
dir.create("data")

load(file.path(getwd(), "../data/OISCohorts_processed_20220401.rda"))
load(file.path(getwd(), "../data/OIS_sigMatrix_processed_20220422.rda"))
load(file.path(getwd(), "../data/TCGA.pancan_processed_20211101.rda"))
load("HRASKRASNRASMAP2K1BRAF.all.maf.collapse.rda")  # excluded silent mutation

require(ADAPTS)
require(ADAPTSdata2)
require(xCell)
library(preprocessCore)
require(ggplot2)
# require(maftools)


# Build a new signature matrix containing the OIS phenotype
onlyFibroblast = F
onlyMelanocyte = F
onlyRAS = F
onlyRASRAF = F
filter.sd.quantile = NULL
postNorm = T

s <- OISCohorts_proccessed
s1.cbind = do.call("cbind", lapply(s[c(1)], function(x)
  x$x))
s1.rbind.y = do.call("rbind", lapply(s[c(1)], function(x)
  x$y))
colnames(s1.cbind) <-
  ifelse(s1.rbind.y[, 1] == "C", "melanocytes", "OIS")
rownames(s1.cbind) <- mapg(rownames(s1.cbind))

if (onlyRAS && !onlyRASRAF) {
  sel = c(2:7, 9)
  s2.cbind = do.call("cbind", lapply(s[sel], function(x)
    x$x))
  s2.rbind.y = do.call("rbind", lapply(s[sel], function(x)
    x$y))
} else if (onlyRASRAF && !onlyRAS) {
  s2.cbind = do.call("cbind", lapply(s[c(2:7, 9:10)], function(x)
    x$x))
  s2.rbind.y = do.call("rbind", lapply(s[c(2:7, 9:10)], function(x)
    x$y))
} else {
  s2.cbind = do.call("cbind", lapply(s[c(2:10)], function(x)
    x$x))
  s2.rbind.y = do.call("rbind", lapply(s[c(2:10)], function(x)
    x$y))
}
colnames(s2.cbind) <-
  ifelse(s2.rbind.y[, 1] == "C", "fibroblast", "OIS")
rownames(s2.cbind) <- mapg(rownames(s2.cbind))
dim(s2.cbind)

if (onlyFibroblast &&
    !onlyMelanocyte) {
  s.cbind = s2.cbind
} else if (!onlyFibroblast &&
           onlyMelanocyte) {
  s.cbind = s1.cbind
} else {
  s.cbind = cbind(s1.cbind, s2.cbind)
}
LM22 <- ADAPTS::LM22
LM22.source <- ADAPTSdata2::fullLM22
tar <- s.cbind
inx <- intersect(rownames(LM22.source), rownames(tar))
allData <- cbind(LM22.source[inx, ], tar[inx, ])

if (!is.null(filter.sd.quantile)) {
  allData.qn <- normalize.q(allData, filter.sd.quantile)
} else {
  allData.qn = data.frame(preprocessCore::normalize.quantiles(as.matrix(allData)))
  rownames(allData.qn) = rownames(allData)
  colnames(allData.qn) = colnames(allData)
  dim(allData)
  dim(allData.qn)
  peep(allData.qn)
  peep(allData)
}
  seedSigMatrix <- LM22

colnames(seedSigMatrix) <-
  sub('\\.[0-9]+$', '', colnames(seedSigMatrix))
seedSigMatrix.pre <-
  lapply(unique(colnames(seedSigMatrix)), function(x) {
    rowMeans(seedSigMatrix[, colnames(seedSigMatrix) == x, drop = F], na.rm =
               TRUE)
  })
seedSigMatrix.pre <- do.call(cbind, seedSigMatrix.pre)
colnames(seedSigMatrix.pre) <- unique(colnames(seedSigMatrix))
seedSigMatrix <- seedSigMatrix.pre

# Rank genes for each cell type: Use a t-statistic to rank to features for each cell type
gList <- rankByT(allData.qn)
ref <- allData.qn[, colnames(LM22.source)]
if (onlyFibroblast &&
    !onlyMelanocyte) {
  addCellTypeColInx = grep("fibroblast|OIS", colnames(allData.qn))
} else if (!onlyFibroblast && onlyMelanocyte) {
  addCellTypeColInx = grep("melanocytes|OIS", colnames(allData.qn))
} else {
  addCellTypeColInx = grep("melanocytes|fibroblast|OIS", colnames(allData.qn))
}
addCellType <- allData.qn[, addCellTypeColInx, drop = F]

colnames(ref) <- sub('\\.[0-9]+$', '', colnames(ref))
colnames(addCellType)  <-
  sub('\\.[0-9]+$', '', colnames(addCellType))

inxSeedTar = intersect(rownames(seedSigMatrix), rownames(tar))
seedSigMatrix = seedSigMatrix[inxSeedTar,]

newSigMatrix <-
  AugmentSigMatrix(
    seedSigMatrix,
    ref,
    addCellType,
    gList,
    plotToPDF = T,
    imputeMissing = F,
    nGenes = 1:100,
    postNorm = postNorm,
    minSumToRem = NA,
    addTitle = NULL,
    calcSpillOver = T,
    pdfDir = "./"
  )



# Apply the newSigmatrix to 4 cancer types with the highest frequency of Ras/Raf/MEK mutations to explore the gene expression of OIS phenotype. 
onlyRASRAFMEK = F
excludeRASRAFMEK = F
onlypaired = F
onlyPrimaryTumor = F
excludeNormalAdjacentTissue = T
CDKN2Aorder = F
minInx = 12

results <- list()
tarCAlist = unique(TCGA.pancan@pheno$cancerType)

for (tarCA in tarCAlist) {
  RASRAFMEKInx <-
    grep(paste0(substr(
      unique(HRASKRASNRASMAP2K1BRAF[[tarCA]]$maf$Tumor_Sample_Barcode),
      1,
      12
    ), collapse = "|"), colnames(TCGA.pancan@x))
  
  if (onlyRASRAFMEK) {
    inx2 = colnames(TCGA.pancan@x)[RASRAFMEKInx]
    inx2
  } else {
    inx2 <-
      excludeInx <-
      rownames(TCGA.pancan@pheno)[which(TCGA.pancan@pheno$cancerType %in% tarCA)]
  }
  
  if (excludeRASRAFMEK) {
    inx2 = setdiff(excludeInx, inx2)
  }
  
  if (onlypaired) {
    inx3 = inx2[duplicated(substr(unique(inx2), 1, 12))]
    grepInx = paste0(substr(unique(inx3), 1, 12), collapse = "|")
    if (grepInx != "") {
      inx2 = inx2[grep(grepInx, inx2)]
    } else {
      inx2 = NULL
    }
    
    if (onlyPrimaryTumor)
      inx2 = inx2[which(substr(unique(inx2), 14, 15) == "01")]
  }
  
  if (onlyPrimaryTumor)
    inx2 = inx2[which(substr(unique(inx2), 14, 15) == "01")]
  if (excludeNormalAdjacentTissue &&
      length(which(substr(unique(inx2), 14, 15) == "11")) != 0)
    inx2 = inx2[-which(substr(unique(inx2), 14, 15) == "11")]
  
  if (is.null(inx2))
    next   
  expr <- TCGA.pancan@x[, inx2]
  if (dim(expr)[2] < minInx)
    next  
  anno <-
    data.frame(
      SampleID = colnames(expr),
      ExperimentID = "",
      SampleType = TCGA.pancan@pheno[inx2, "sample_type"]
    )
  anno
  rownames(anno) <- anno$SampleID
  dim(anno)
  tmp <-
    list(expr.list = list(tmp.expr = expr),
         anno.list = list(tmp.anno = anno))
  rownames(expr) <- mapg(rownames(expr))
  
  estCellPercentDCQ <-
    estCellPercent.DCQ(refExpr = newSigMatrix$sigMatrix, geneExpr = expr)
  
  inxOIShigh <-
    colnames(estCellPercentDCQ)[which(estCellPercentDCQ["OIS", ] > mean(estCellPercentDCQ["OIS", ]))]
  inxOISlow <-
    colnames(estCellPercentDCQ)[-which(estCellPercentDCQ["OIS", ] > mean(estCellPercentDCQ["OIS", ]))]
  inxOrder <-
    colnames(estCellPercentDCQ)[order(estCellPercentDCQ["OIS", ])]
  results[[tarCA]]$inxOrder.list <- inxOrder
  estCellPercentDCQ.pre <-  estCellPercentDCQ

  if (xcell.val) {
    require(xCell)
    estCellPercentDCQ <- data.frame(xCellAnalysis(expr))
    colnames(estCellPercentDCQ) <-
      gsub("\\.", "-", colnames(estCellPercentDCQ))
  }
  
  bdf = anno[inxOrder,]
  bdf$'OIS score' <-
    t(scale01(t(estCellPercentDCQ.pre["OIS", inxOrder, drop = F])))[1, ]
  
  SampleTypes = unique(TCGA.pancan@pheno$sample_type)[c(2, 1, 4, 7, 5, 6, 3)]
  SampleTypesCol = c("#1f77b4","#FF9900","#FF0000","#d62728","#d62728","#ff7f0e","#2ca02c")
  names(SampleTypesCol) = SampleTypes 
  
  bdf$SampleType = factor(bdf$SampleType)
  bdf$SampleType = factor(bdf$SampleType, levels = intersect(SampleTypes, bdf$SampleType))
  SampleTypeColIn = SampleTypesCol[levels(bdf$SampleType)]
  bdf$'CDKN2A' <- expr["CDKN2A", rownames(bdf)]
  
  
  load("HRASKRASNRASMAP2K1BRAF.all.maf.collapse.rda")  # excluded silent mutation
  bdf$SampleType
  rownames(bdf)
  bdf$Variant_Classification <-
    ifelse(
      bdf$SampleType %in% c("Primary Tumor", "Metastatic"),
      HRASKRASNRASMAP2K1BRAF.all.maf.collapse[rownames(bdf), c("Variant_Classification")],
      NA
    )
  bdf$Oncogene <-
    ifelse(
      bdf$SampleType  %in% c("Primary Tumor", "Metastatic"),
      HRASKRASNRASMAP2K1BRAF.all.maf.collapse[rownames(bdf), "Hugo_Symbol"],
      NA
    )
  oncogenes[[tarCA]] <- bdf$Oncogene
  bdf$Oncogene <- factor(bdf$Oncogene)
  bdf$Mutation <-
    factor(bdf$Oncogene, levels = intersect(names(laslafmakCol2), bdf$Oncogene))
  MutationColIn <- laslafmakCol2[levels(bdf$Oncogene)]
  
  if (!xcell.val) {
    ICBs = c(
      "PDCD1",
      "CD274",
      "CTLA4",
      "PDCD1LG2",
      "HAVCR2",
      "LAG3",
      "CD244"
    ) #
    for (i in ICBs) {
      dat = data.frame(ois = bdf$'OIS score', tmp = expr[i, rownames(bdf)])
      OISICBcorES[[i]] <- cor.test(dat$tmp, dat$ois)$estimate
      OISICBcorPVal[[i]] <- cor.test(dat$tmp, dat$ois)$p.value
      OISICBres[[i]] <-
        gplot(
          tmp ~ ois,
          dat,
          col = cancerTypeCol[tarCA],
          fill = "black",
          shape = 1,
          shape.alpha = .4,
          add.params = list(
            col = cancerTypeCol[tarCA],
            size = .4,
            alpha = .7
          )
        ) + scale_x_continuous(n.breaks = 3) + labs(
          subtitle = tarCA,
          x = expression(bold("OIS score")),
          y = substitute(paste(bolditalic(i)), list(i = i))
        )
    }
    results[[tarCA]]$resOISICBcorES <- OISICBcorES
    results[[tarCA]]$resOISICBcorPVal <- OISICBcorPVal
    results[[tarCA]]$resOISICBScatterPlot <- OISICBres
  }
  
  if (CDKN2Aorder)
    inxOrder = rownames(bdf)[order(bdf$'CDKN2A')]
  
  colnames(bdf)
  bdf = bdf[inxOrder, 3:9, drop = F]
  bdf <-
    add.gs(
      estCellPercentDCQ.pre[, inxOrder],
      bdf,
      list(tar1 = "OIS"),
      add.gs.quantile = list(
        median = c(0, 1 / 2, 1),
        tertiles = c(0, 1 / 3, 2 / 3, 1),
        quartiles = c(0, 1 / 4, 2 / 4, 3 / 4, 1),
        quintiles = c(0, 1 / 5, 2 / 5, 3 / 5, 4 / 5, 1)
      )
    )$y
  
  results[[tarCA]]$bdf.di <- bdf
  results[[tarCA]]$icb.di <-
    add.gs(expr[, inxOrder], bdf, list(tar2 = "PDCD1"), add.gs.quantile = list(median =
                                                                                 c(0, 1 / 2, 1)))$y
}

rsp <-
  do.call("c", lapply(results, function(x)
    x$resOISICBScatterPlot[1:4]))
  do.call(ggpubr::ggarrange, c(rsp, list(
  ncol = 4, nrow = 4, align = "h"
)))

daES = sapply(results, function(x)
  unlist(x$resOISICBcorES))[4:10,]
daPVal = sapply(results, function(x)
  unlist(x$resOISICBcorPVal))[4:10,]
da = melt(daES, 1)
colnames(da)[1:2] <- c("IC Genes", "Cancer Type")
da$`IC Genes` <- gsub("\\.cor", "", da$`IC Genes`)
da$sig <- "***"
da$`Cancer Type` = factor(da$`Cancer Type`, levels = rev(tarCAlist))

p <-
  gbarplot(
    value ~ `IC Genes`,
    da[dim(da):1,],
    2,
    sd = F,
    beside = T,
    rotate = T,
    border = F,
    col = cancerTypeCol
  )
p + geom_text(
  aes(label = sig, group = `Cancer Type`),
  position = position_dodge(width = .8),
  size = 3,
  vjust = 1,
  hjust = 1.2
)
tcgaCl.tmp <- data.frame(TCGA.pancan@meta$PancanClinical2)
rownames(tcgaCl.tmp) <-
  TCGA.pancan@meta$PancanClinical2$bcr_patient_barcode
sample_type = unique(TCGA.pancan@pheno$sample_type)
inxOrder.list = lapply(results, function(x)
  x$inxOrder.list)
canInx = names(inxOrder.list)



#  Examine whether the stratification of patients by their OIS index allows disease outcome prediction through survival analysis
monthScale = T
onlyEarlyStage = F
onlyLateStage = F
RFS = F
require(survival)
require(survminer)
res <- ggRes <- ggResPval <- bdeRes <- coxfit <- msr <- list()
for (i in canInx) {
  if (!RFS) {
    OSorRFS = 3:4
  } else {
    OSorRFS = 13:14
  }
  cl <-
    cbind(TCGA.pancan@pheno[inx2 <-
                              rownames(TCGA.pancan@pheno)[intersect(
                                which(TCGA.pancan@pheno$sample_type %in% sample_type),
                                which(TCGA.pancan@pheno$X_primary_disease %in% i)
                              )], names(TCGA.pancan@pheno)[c(OSorRFS, 19, 20)]] , data.frame(tcgaCl.tmp)[substr(inx2, 1, 12), c("pathologic_stage"), drop =
                                                                                                           F])
  
  colnames(cl) <- c("time", "status", "age", "sex", "stage")
  cl <- cl[intersect(inxOrder.list[[i]], rownames(cl)),]
  dim(cl)
  if (monthScale)
    cl$time <- cl$time / 30.5
  
  if (onlyEarlyStage &&
      !onlyLateStage) {
    cl = cl[grep(
      "Stage I$|Stage IA$|Stage IB$|Stage II$|Stage IIA$|Stage IIB$|Stage IIC$",
      cl$stage
    ), , drop = F]
  } else if (!onlyEarlyStage && onlyLateStage) {
    cl = cl[grep(
      "Stage III$|Stage IIIA$|Stage IIIB$|Stage IIIC$|Stage IV$|Stage IVA$|Stage IVB$|Stage IVC$",
      cl$stage
    ), , drop = F]
  } else if (!onlyEarlyStage &&
             !onlyLateStage) {
    cl = cl
  } else {
    cat.box(, "CHECK!!")
  }
  
  exp = TCGA.pancan@x[, rownames(cl)]
  rownames(exp) <- mapg(rownames(exp))

  cl$sex <- ifelse(cl$sex == "", NA, cl$sex)
  cl$stage <- factor(cl$stage)
  
  bdf.di.list  = lapply(results, function(x)
    x$bdf.di)
  colnames(bdf)
  bdf = bdf.di.list[[i]]
  rownames(bdf) == colnames(exp)
  
  cl$Subtype <- bdf[rownames(cl), "tar1.median"]
  # cl$Subtype <- bdf[rownames(cl),"tar1.quartiles" ]
  
  bde <-
    cbind(cl, t(exp))
  bde
  bde$sex <-
    as.factor(bde$sex)
  bde$stage <- as.factor(bde$stage)
  dim(bde)
  
  bdeRes[[i]] <- bde
  
  fit <-
    tryCatch(
      survfit(Surv(time, status) ~ Subtype, data = bde),
      error = function(e) {
        NULL
      },
      warning = function(w) {
        cat("undefined columns selected")
      }
    )
  coxfit[[i]] <-
    tryCatch(
      coxph(Surv(time, status) ~ ., bde[, c(1:6)]),
      error = function(e) {
        NULL
      },
      warning = function(w) {
        cat("No (non-missing) observations")
      }
    )
  
  surv_diff <-
    tryCatch(
      survdiff(Surv(time, status) ~ Subtype, bde),
      error = function(e) {
        NULL
      },
      warning = function(w) {
        cat("undefined columns selected")
      }
    )
  ggResPval[[i]] <-
    tryCatch(
      1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1),
      error = function(e) {
        NULL
      },
      warning = function(w) {
        cat("undefined columns selected")
      }
    )
  
  risk.table = "nrisk_cumevents"
  risk.table = F
  surv.median.line = "none"
  ggRes[[i]] <-
    tryCatch(
      ggsurvplot(
        fit,
        data = bde,
        pval = T,
        palette = c("#467EB3","#D4352B"),
        risk.table = risk.table,
        surv.median.line = surv.median.line,
        tables.height = .3,
        conf.int = F,
        size = 0.5,
        pval.size = 3.2,
        ggtheme = theme_minimal(),
        add.all = F,
        risk.table.fontsize = 3.5,
        font.main = c(11, "plain", "black"),
        font.x = c(11),
        font.y = c(11),
        font.tickslab = c(9),
        fontsize = 4.5,
        tables.y.text = T,
        risk.table.title = ""
      ) + labs(subtitle = paste(i, ifelse(RFS, "RFS", "OS"))),
      error = function(e) {
        NULL
      },
      warning = function(w) {
        cat("check ggsurvplot")
      }
    )
  
}
