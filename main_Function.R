options(stringsAsFactors = FALSE) 

loadUrl <-
  function(url,
           downloadPath = NA,
           sep = c("RData", " ", "," , "\t", ";", "xls", "gsheet"),
           onedrive = T,
           ...) {
    cat(
      'onedrive: copy link\n googlesheet: share-> Anyone with the link\n sep: "RData", ..."xls", "gsheet"\n'
    )
    if (!is.na(downloadPath))  {
      tmpFile <- downloadPath
      
    } else {
      tmpFile <-
        file.path(getwd(), paste0(substr(Sys.time(), 1, 10), '.rda'))
    }
    if (onedrive) {
      url2 <- gsub("e=.*", "download=1", url)
    } else {
      url2 <- url
    }
    
    download.file(url2, tmpFile, mode = "wb") 
    sep <- match.arg(sep)
    if (sep == "RData") {
      print(tmpFile)
      tmpFile <-  gsub("\\\\", "/", tmpFile)
      justLoaded <- try(load(tmpFile), silent = T)
      ;
      try(assign(justLoaded, eval(as.symbol(justLoaded)), .GlobalEnv), silent = T)
      ;
      if (class(justLoaded) == "try-error") {
        justLoaded <-
          try(read.delim(tmpFile, ...), silent = T)
        ; message("Need 'sep' argument, is it txt file?")
      }
    } else if (sep == "xls") {
      is.installed('readxl')
      justLoaded <- try(readxl::read_excel(tmpFile, ...), silent = T)
      
    } else if (sep == "gsheet") {
      is.installed("gsheet")
      cat('gsheet should be public, click share-> Anyone with the link\n')
      justLoaded <- gsheet::gsheet2tbl(url, ...)
    } else {
      justLoaded <- try(read.delim(tmpFile, sep = sep, ...), silent = T)
    }
    justLoaded
  }


Capitalize <-
  function(v,
           captalOption = c("toupper",
                            "tolower",
                            "firstCapital",
                            "multiCapital",
                            "titleCase"),
           str.width = NULL) {
    captalOption = match.arg(captalOption)
    
    if (captalOption == "toupper") {
      bb <-
        lapply(as.list(v), function(x)
          toupper(gsub("^ |^  |^   | $|  $|   $|", "", x)))
      cc <- bb
    }
    if (captalOption == "tolower") {
      bb <-
        lapply(as.list(v), function(x)
          tolower(gsub("^ |^  |^   | $|  $|   $|", "", x)))
      cc <- bb
    }
    
    if (captalOption == "firstCapital") {
      bb <-
        lapply(as.list(v), function(x)
          gsub("^ |^  |^   | $|  $|   $|", "", x))
      cc <-
        lapply(bb, function(x) {
          capped <- grep("^[A-Z]", x, invert = T)
          substr(x[capped], 1, 1) <- toupper(substr(x[capped],
                                                    1, 1))
          paste(x, collapse = " ")
        })
    }
    
    if (captalOption == "multiCapital") {
      bb <-
        lapply(strsplit(v, " "), function(x)
          gsub("^ |^  |^   | $|  $|   $|", "", x))
      
      cc <-
        lapply(bb, function(x) {
          capped <- grep("^[A-Z]", x, invert = T)
          substr(x[capped], 1, 1) <- toupper(substr(x[capped],
                                                    1, 1))
          paste(x, collapse = " ")
        })
    }
    
    if (captalOption == "titleCase") {
      bb <-
        lapply(as.list(v), function(x)
          tolower(gsub("^ |^  |^   | $|  $|   $|", "", x)))
      cc <- lapply(bb, function(x)
        tools::toTitleCase(x))
    }
    if (!is.null(str.width)) {
      cc <- do.call("c", .str(cc, str.width))
    } else {
      cc <- do.call("c", cc)
    }
    cc
  }    





mapg <-
  function(keys,
           probeSeq = F,
           probePlatform = c("hgu133a", "hgu133plus2", "GPL10558")) {
    require(org.Hs.eg.db) 
    cat.box(
      "mapG pipeline" ,
      " 1. EntrezID \n 2. HGNC official symbols were mapped to entrez ID (cf. in case of asNA => show alias to HGNC)"
    )
    
    if (length(keys) == 0)
      return(NA)
    keys <-  gsub("///.*", "", keys)
    temp <-
      tryCatch(
        as.numeric(keys),
        warning = function(x)
          print("Assumes that terms are all HGNC approved symbols")
      )
    probePlatform = match.arg(probePlatform)
    
    if (probeSeq) {
      if (any(probePlatform == "hgu133a")) {
        if (!is.numeric(temp)) {
          require(hgu133a.db)
          columns(hgu133a.db)
          a = mapIds(hgu133a.db, keys, c("PROBEID"), "SYMBOL", multiVals =
                       "list")
          print(a)
        } else {
          a = mapIds(hgu133a.db, keys, c("PROBEID"), "ENTREZID", multiVals = "list")
          print(a)
        }
        require(hgu133aprobe) 
        data(hgu133aprobe)
        return(hgu133aprobe[which(hgu133aprobe[, "Probe.Set.Name"] == a[[1]][1]), ])
        
      }
      if (any (probePlatform == "hgu133aplus2")) {
        if (!is.numeric(temp)) {
          require(hgu133plus2.db)
          a = mapIds(hgu133plus2.db,
                     keys,
                     c("PROBEID"),
                     "SYMBOL",
                     multiVals = "list")
          print(a)
        } else {
          a = mapIds(hgu133a.db, keys, c("PROBEID"), "ENTREZID", multiVals = "list")
          print(a)
        }
        require(hgu133aprobe)
        data(hgu133aprobe)
        return(hgu133aprobe[which(hgu133aprobe[, "Probe.Set.Name"] == a[[1]][1]), ])
      }
      
      if (any (probePlatform == "GPL10558")) {
        # illumina GPL10558 Illumina HumanHT-12 V4.0 expression beadchip
        load("~/rawTemp/GPL10558-50081.RData")
        colnames(illu.seq)
        print(keys)
        return(illu.seq[which(illu.seq[, "Symbol"] == keys),])
      }
    }
    if (is.numeric(temp)) {
      hgnc <-
        mapIds(
          org.Hs.eg.db,
          keys = keys,
          column = "SYMBOL",
          keytype = "ENTREZID",
          multiVals = "asNA"
        )
      duplicated.inx <<- which(duplicated(hgnc))
      cat.box(
        "Handling multiVals = asNA",
        paste0(
          paste(keys[duplicated.inx], collapse = ","),
          "\n\nIndex of above ID wewe now assigned as global variable 'duplicated.inx'\n\n"
        )
      )
      return(hgnc) 
    } else {
      entrez <-
        mapIds(
          org.Hs.eg.db,
          keys = keys,
          column = "ENTREZID",
          keytype = "SYMBOL",
          multiVals = "asNA"
        )
      
      multi.hgnc <-
        mapIds(
          org.Hs.eg.db,
          keys = keys,
          column = "SYMBOL",
          keytype = "ALIAS",
          multiVals = "list"
        )
      a = as.list(keys)
      b = multi.hgnc
      c = c()
      for (i in 1:length(a)) {
        c[i] <- isTRUE(a[[i]] == b[[i]])
      }
      diff.aliasToHGNC <- multi.hgnc[which(c == F)]
      attributes(entrez)$keytype <- "SYMBOL"
      attributes(entrez)$multiVals <- "asNA"
      attributes(entrez)$aliasToMulti.hgnc <- multi.hgnc
      attributes(entrez)$diff.aliasToHGNC <- diff.aliasToHGNC
      
      if (any(is.na(entrez))) {
        first.hgnc <-
          mapIds(
            org.Hs.eg.db,
            keys = keys,
            column = "SYMBOL",
            keytype = "ALIAS"
          )
        aliasToFirstHGNCToEntrez <-
          mapIds(
            org.Hs.eg.db,
            keys = first.hgnc,
            column = "ENTREZID",
            keytype = "SYMBOL",
            multiVals = "asNA"
          )
        entrez[which(is.na(entrez))] <-
          aliasToFirstHGNCToEntrez[which(is.na(entrez))]
        return(entrez)
      }   else  {
        return(entrez)
      }
    }
  }



toTCGAAcronym <- function(vec) {
  tcgaStudyAbbreviations <<-
    data.frame(
      Abbreviation = c(
        "LAML",
        "ACC",
        "BLCA",
        "LGG",
        "BRCA",
        "CESC",
        "CHOL",
        "LCML",
        "COAD",
        "CNTL",
        "ESCA",
        "FPPP",
        "GBM",
        "HNSC",
        "KICH",
        "KIRC",
        "KIRP",
        "LIHC",
        "LUAD",
        "LUSC",
        "DLBC",
        "MESO",
        "MISC",
        "OV",
        "PAAD",
        "PCPG",
        "PRAD",
        "READ",
        "SARC",
        "SKCM",
        "STAD",
        "TGCT",
        "THYM",
        "THCA",
        "UCS",
        "UCEC",
        "UVM"
      ),
      Study = c(
        "Acute Myeloid Leukemia",
        "Adrenocortical carcinoma",
        "Bladder Urothelial Carcinoma",
        "Brain Lower Grade Glioma",
        "Breast invasive carcinoma",
        "Cervical squamous cell carcinoma and endocervical adenocarcinoma",
        "Cholangiocarcinoma",
        "Chronic Myelogenous Leukemia",
        "Colon adenocarcinoma",
        "Controls",
        "Esophageal carcinoma",
        "FFPE Pilot Phase II",
        "Glioblastoma multiforme",
        "Head and Neck squamous cell carcinoma",
        "Kidney Chromophobe",
        "Kidney renal clear cell carcinoma",
        "Kidney renal papillary cell carcinoma",
        "Liver hepatocellular carcinoma",
        "Lung adenocarcinoma",
        "Lung squamous cell carcinoma",
        "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma",
        "Mesothelioma",
        "Miscellaneous",
        "Ovarian serous cystadenocarcinoma",
        "Pancreatic adenocarcinoma",
        "Pheochromocytoma and Paraganglioma",
        "Prostate adenocarcinoma",
        "Rectum adenocarcinoma",
        "Sarcoma",
        "Skin Cutaneous Melanoma",
        "Stomach adenocarcinoma",
        "Testicular Germ Cell Tumors",
        "Thymoma",
        "Thyroid carcinoma",
        "Uterine Carcinosarcoma",
        "Uterine Corpus Endometrial Carcinoma",
        "Uveal Melanoma"
      ),
      Surrogate = c(
        "acute",
        "adrenocortical",
        "bladder",
        "brain",
        "breast",
        "cervical",
        "cholangiocarcinoma",
        "chronic",
        "colon",
        "controls",
        "esophageal",
        "ffpe",
        "glioblastoma",
        "head",
        "chromophobe",
        "clear",
        "papillary",
        "liver",
        "lung.*adenocarcinoma",
        "lung.*squamous",
        "lymphoma",
        "mesothelioma",
        "miscellaneous",
        "ovarian",
        "pancreatic",
        "pheochromocytoma",
        "prostate",
        "rectum",
        "sarcoma",
        "skin",
        "stomach",
        "testicular",
        "thymoma",
        "thyroid",
        "carcinosarcoma",
        "corpus",
        "uveal"
      )
    )
  tcgaStudyAbbreviations.pre <<- tcgaStudyAbbreviations
  print(tcgaStudyAbbreviations)
  
  cat("\ntcgaStudyAbbreviations(df) is now assinged as global variable\n")
  inx <- c()
  if (any(vec %in% tcgaStudyAbbreviations$Abbreviation)) {
    for (i in 1:dim(tcgaStudyAbbreviations)[1]) {
      vec[k <-
            grep(tcgaStudyAbbreviations$Abbreviation[i], vec)] <-
        tcgaStudyAbbreviations$Study[i]
      if (length(k) != 0)
        inx <- c(inx, i)
    }
    vec = Capitalize(Capitalize(vec, "tolower"), "firstCapital")
  } else {
    vec <- tolower(vec)
    for (i in 1:dim(tcgaStudyAbbreviations)[1]) {
      vec[k<-grep(tcgaStudyAbbreviations$Surrogate[i], vec)] <-
        tcgaStudyAbbreviations$Abbreviation[i]
      if (length(k) != 0)
        inx <- c(inx, i)
    }
  }
  
  cat("\nAbbreviations: ")
  cat(sprintf(
    "%s, %s;",
    tcgaStudyAbbreviations$Abbreviation[rev(inx)],
    tolower(tcgaStudyAbbreviations$Study[rev(inx)])
  ))
  cat("\n\n")
  print(vec)
}


factorizeDf <- function(df, lengthOfUniqueCut = 10) {
  if (!is.null(df)) {
    df.tmp <-
      data.frame(lapply(df, function(x) {
        if (length(unique(x)) < lengthOfUniqueCut) {
          as.factor(x)
        } else {
          x
        }
      }))
    colnames(df.tmp) <-
      colnames(df)
    rownames(df.tmp) <- rownames(df)
    df.tmp
  }
}


cat.box = function(title = "CAUTION", text, type = 2) {
  require(crayon)
  title = paste(unlist(strsplit(toupper(title) , "")), collapse = " ")
  if (type == 1) {
    cat(
      bold('\n\n\\') %+% white(
        '----------------------------------------------------//\n'
      ) %+% bold(paste0(" ", title)) %+% silver(
        '\n ----------------------------------------------------\n'
      ) %+% yellow$italic(text) %+% white(
        '\n//----------------------------------------------------'
      ) %+% bold('\\\n\n\n')
    )
  } else {
    cat(
      red$bold('\n\n\\___') %+% blue$bold('___') %+% silver$bold('______________________________________________\n\n') %+% bold(paste0(" ", title)) %+% silver(
        '\n ----------------------------------------------------\n'
      ) %+% yellow$italic(text) %+% white$bold('\n ___') %+% white$bold(
        '_________________________________________________\n\n\n'
      )
    )
  }
}






options(mySeedV=76543217)
biPDSmlv2 <- function (d.merged,
                       devCohorts,
                       pw.in = NULL,
                       path = "",
                       PDS.min_exp = -10,
                       PDS.min_std = 0.2,
                       className = 'class',
                       classesTrain = c("C", "E"),
                       familyName = "binomial",
                       Screening = F,
                       screen.name = "",
                       CVal = c("LOOCV", "LOSOCV", "nfold"),
                       nfold = 10,
                       GlobalOp = c("epsgo", "GA",  "none"),
                       alphas.seq = seq(0.03, 1, length = 10),
                       type.measure = "deviance",
                       type.min = c("lambda.min", "lambda.1se"),
                       geneSymbol = TRUE,
                       populations = 30,
                       generation = 5,
                       penalty.factor = rep(1, dim(d.merged[[1]]$x)[1]),
                       authors = paste(names(d.merged), "authors"),
                       seed = 29365917,
                       prior.pds.results = NULL,
                       prior.predsList = NULL,
                       browseTable = F,
                       selectedStat = T,
                       preCalpredsList = NULL,
                       valBatchCorrectMethod = "combat",
                       ...) {
  require(foreach)
  require(sva)
  require(glmnet)
  require(c060)
  require(pROC)
  require(Epi)
  require(plot3D)
  require(ggplot2)
  require(precrec)
  require(reshape2)
  require(RColorBrewer)
  require(dplyr)
  require(beeswarm)
  require(ROCR)
  require(pamr)
  require(pheatmap)
  require(gridExtra)
  require(GSAR)
  require(org.Hs.eg.db)
  require(AnnotationDbi)
  
  try(dev.off(dev.list()["RStudioGD"]), silent = TRUE)
  
  op <- par(no.readonly = T)
  set.seed(seed)
  
  GlobalOp = match.arg(GlobalOp)
  GlobalOp.pre <- GlobalOp
  type.min = match.arg(type.min)
  CVal = match.arg(CVal)
  p0 <- p1 <- p2 <- p3 <- p4 <- p5 <- p6 <- p7 <- p8 <- NULL
  pp1 <- pp2 <- pp3 <- pp4 <- pp5 <- pp6 <- pp7 <- pp8 <- NULL
  out.st2 <- out.t2 <-  NULL
  
  valCohorts = 1 - devCohorts
  pdf.options(family = 'Helvetica', useDingbats = F)
  
  ifelse(!dir.exists("./results"),
         dir.create("./results"),
         "Folder exists already")
  
  if (!is.null(pw.in)) {
    pw.in.pre <- NULL
    pw.in.pre$pathway <- as.list(names(pw.in))
    pw.in.pre$entrez_gene_ids <- pw.in
    pw.in <- pw.in.pre
    PDS = T
  } else {
    PDS = F
  }
  
  for (i in 1:length(d.merged))
    d.merged[[i]][[2]] <-  d.merged[[i]][[2]][, 1]
  study <-
    foreach(a = 1:length(d.merged), .combine = "c") %do% {
      rep(names(d.merged)[a], length(d.merged[[a]][[2]]))
    }
  sample <-
    foreach(a = 1:length(d.merged), .combine = "c") %do% {
      colnames(d.merged[[a]][[1]])
    }
  class <-
    foreach(a = 1:length(d.merged), .combine = "c") %do% {
      d.merged[[a]][[2]]
    }
  ematList <-
    foreach(a = 1:length(d.merged)) %do% {
      d.merged[[a]][[1]]
    }
  names(ematList) <- names(d.merged)
  if (any(rownames(ematList[[a]]) == "")) {
    for (a in 1:length(ematList)) {
      ematList[[a]] = ematList[[a]][-which(rownames(ematList[[a]]) == ""), ]
    }
  }
  meta2Tmp <-
    data.frame(
      study = as.character(study),
      sample = as.character(sample),
      class = as.factor(class),
      stringsAsFactors = F
    )
  str(ematList)
  
  
  meta <-
    data.frame(
      study = as.character(names(ematList)),
      devCohorts = devCohorts,
      valCohorts = valCohorts,
      stringsAsFactors = F
    )
  intercept = TRUE
  rownames(meta) = meta[, 'study']
  meta[, 'devCohorts'] = meta[, 'devCohorts'] == 1
  meta[, 'valCohorts'] = meta[, 'valCohorts'] == 1
  rownames(meta2Tmp) = meta2Tmp[, 'sample']
  meta2 = meta2Tmp[meta2Tmp[, 'study'] %in% meta[, 'study'], ]
  mergedClass = mergeStudyData(ematList[meta[meta[, 'devCohorts'], 'study']], meta2, covariateName =
                                 NA)
  dim(mergedClass)
  
  mergedClass.pre <<- mergedClass
  
  devNames = meta[meta[, 'devCohorts'], 'study']
  idxDev = (meta2[, 'study'] %in% devNames) &
    (meta2[, className] %in% classesTrain)
  devSampleNames = rownames(meta2)[idxDev]
  meta2Dev = meta2[idxDev, ]
  meta2Dev[, className] = factor(meta2Dev[, className], levels = classesTrain)
  devMerged = mergedClass[, meta2Dev[, 'sample']]
  s.Title <- sprintf("%sAND%s", devNames[1], length(devNames) - 1)
  glmnetArgs = makeGlmnetArgs(meta2Dev)
  print(glmnetArgs)
  print(table(meta2$class))
  print(table(meta2$study))
  
  if (!is.null(prior.pds.results) || PDS) {
    if (!is.null(prior.pds.results)) {
      exp <- devMerged
      PDS.res = prior.pds.results
      insel <- meta2[which(meta2[, 'study'] %in% devNames), ]
      PDSmatrix <- do.call(rbind.data.frame, PDS.res$scores)
      colnames(PDSmatrix) <- colnames(exp)
      PDSmatrix <- as.matrix(PDSmatrix)
      PDS = T
      gene_ids <- PDS.res$gene_ids
      
    } else {
      exp <- devMerged
      pw.in <- pw.in
      gene_ids <- unique(do.call("c", pw.in$entrez_gene_ids))
      nor <- rep(F, dim(exp)[2])
      nor[which(colnames(exp) %in% meta2$sample[which(meta2$class == "C")])] <-
        T
      normal.cv.merged <- nor
      print(dim(exp))
      print(pw.in$entrez_gene_ids[1:2])
      print(pw.in$pathway[1:2])
      print(normal.cv.merged)
      exp.cv.merged <- exp
      
      if (file.exists(file.path(
        "./results",
        sprintf("%s%s_PDS.res.RData", s.Title, path)
      ))) {
        load(file.path(
          "./results",
          sprintf("%s%s_PDS.res.RData", s.Title, path)
        ))
      } else {
        PDS.res <-
          wrap.PDS(
            exp.cv.merged,
            rownames(exp.cv.merged),
            pw.in$entrez_gene_ids,
            pw.in$pathway,
            normal.cv.merged,
            attempts = 2,
            min_exp = PDS.min_exp ,
            min_std = PDS.min_std
          )
      }
      PDSmatrix <- do.call(rbind.data.frame, PDS.res$scores)
      colnames(PDSmatrix) <- colnames(exp)
      PDSmatrix <- as.matrix(PDSmatrix)
      insel <- meta2[which(meta2[, 'study'] %in% devNames), ]
      PDS.res$PDSmatrix <- PDSmatrix
      PDS.res$insel <- insel
      PDS.res$gene_ids <- gene_ids
      if (is.null(prior.pds.results) &&
          !file.exists(file.path(
            "./results",
            sprintf("%s%s_PDS.res.RData", s.Title, path)
          )))
        save(PDS.res, file = file.path(
          "./results",
          sprintf("%s%s_PDS.res.RData", s.Title, path)
        ))
      cat(sprintf("\n\nPDSmatrix dim: %s\n\n", dim(PDSmatrix)))
    }
  }
  
  if (PDS) {
    Input.D <- "PDS"
  } else {
    Input.D <- "expr"
  }
  
  if (Input.D == "PDS") {
    Input.D.matrix <- PDSmatrix
    geneSymbol = F
  } else {
    Input.D.matrix <- devMerged
  }
  
  if (CVal == "LOOCV") {
    nRepeats = 1
    nFolds = ncol(Input.D.matrix)
    foldid = NA
  } else if (CVal == "LOSOCV") {
    nRepeats = 2
    nFolds = NA
    foldid = glmnetArgs$foldid
  } else {
    nRepeats = 1
    nFolds = nfold
    foldid = NA
  }
  
  if (GlobalOp == "none") {
    alphas = alphas.seq
    GlobalOp = NA
  }
  
  CVconfusion <-
    function(cvFit,
             lambda,
             mergedMatrix,
             meta2,
             className = 'class',
             classLevels = NA) {
      if (is.na(classLevels[1])) {
        classLevels = names(cvFit$glmnet.fit$beta)
      }
      if (class(cvFit$glmnet.fit)[1] == "lognet") {
        cvProbs = cvFit$fit.preval[, which.min(abs(cvFit$lambda - lambda)), drop =
                                     F]
        cvProbs <- cbind(1 - cvProbs, cvProbs)
      } else {
        cvProbs <-
          cvFit$fit.preval[, , which.min(abs(cvFit$lambda - lambda))]
      }
      rownames(cvProbs) <- colnames(mergedMatrix)
      colnames(cvProbs) <- cvFit$glmnet.fit$classnames
      cvProb <- cvProbs
      preds = colnames(cvProbs)[apply(cvProbs, MARGIN = 1, function(x)
        which.max(x))]
      predsFactor = factor(preds, levels = classLevels)
      trueClasses =  factor(meta2[colnames(mergedMatrix), className], levels =
                              classLevels)
      cv.preds <- preds
      cv.predsFactor <- predsFactor
      cv.trueClasses <- trueClasses
      cv.trueClasses.pre <<-
        list(trueClasses,
             predsFactor,
             meta2,
             colnames(mergedMatrix),
             className,
             classLevels)
      confus.cv <- table(trueClasses, predsFactor)
      list(
        cvProb = cvProb,
        cv.preds = cv.preds,
        cv.predsFactor = cv.predsFactor,
        cv.trueClasses = cv.trueClasses,
        confus.cv = confus.cv
      )
    }
  
  if (Screening) {
    GlobalOp = NA
    alphas.default <-
      list(0, 1, seq(0.03, 1, length = 10), 0.993)
    screen.names <-
      c("Ridge", "Lasso", "Elastic-Grid", "Elastic-EPSGO")
    compare.penalizedRegression <-
      function(alphas.default,  screen.names, seed) {
        sel.cvFitList.screen <- NULL
        alpha.screen <- NULL
        cvFit.screen <- NULL
        fitResult.screen <- NULL
        lambda.screen <- NULL
        for (i in 1:length(alphas.default)) {
          alphas = alphas.default[[i]]
          sel.cvFitList.screen[[i]] <-
            cv.merged(
              Input.D.matrix,
              meta2,
              yName = className,
              weights = glmnetArgs$weights,
              nRepeats = nRepeats,
              nFolds = nFolds,
              alphas = alphas,
              family = familyName,
              intercept = intercept,
              keep = TRUE,
              GlobalOp = GlobalOp,
              foldid = foldid,
              type.measure = type.measure ,
              type.min = type.min,
              standardize = F,
              popSize = populations,
              maxiter = generation,
              seed = seed
            )
          if (length(alphas) == 1) {
            if (is.na(foldid[1]))  {
              cvFit.screen[[i]] = sel.cvFitList.screen[[i]][[1]][[1]]
            } else {
              cvFit.screen[[i]] = sel.cvFitList.screen[[i]][[1]]
            }
            alpha.screen[[i]] = alphas
            fitResult.screen[[i]] <- cvFit.screen[[i]]$glmnet.fit
            lambda.screen[[i]] = cvFit.screen[[i]][type.min][[1]]
          } else {
            if (is.na(foldid[1]))  {
              sel.cvFitList <-  sel.cvFitList.screen[[i]]
            } else {
              sel.cvFitList <-  sel.cvFitList.screen
            }
            Crude3Dmatrix = NULL
            for (j in 1:length(sel.cvFitList[[1]])) {
              Crude3Dmatrix = rbind(Crude3Dmatrix,
                                    c(
                                      alphas[j],
                                      sel.cvFitList[[1]][[j]]$lambda.min,
                                      min(sel.cvFitList[[1]][[j]]$cvm)
                                    ))
            }
            D = as.data.frame(Crude3Dmatrix)
            scatter3D(
              D$V1,
              D$V2,
              D$V3,
              phi = 10,
              theta = 33,
              bty = "g",
              type = "p",
              cex.lab = 1.5,
              xlab = "alpha.grid",
              ylab = "lambda.min",
              zlab = "cvm.min",
              pch = 20,
              cex = c(1.2, 1.2, 1.2),
              xlim = c(0, 1.2),
              ylim = c(-2, 30),
              zlim = c(min(D$V3) * .99, max(D$V3) * 1.02)
            )
            D.min = D[which(D[, 3] == min(D[, 3])), ]
            alpha.min = which(D[, 3] == min(D[, 3]))
            alpha.sol = D[alpha.min, 1]
            alpha.screen[[i]] = alpha.sol
            cvFit.screen[[i]] = sel.cvFitList[[1]][[alpha.min]]
            fitResult.screen[[i]] <- cvFit.screen[[i]]$glmnet.fit
            lambda.screen[[i]] = cvFit.screen[[i]][type.min][[1]]
          }
        }
        roc.obj.screen <- NULL
        cv.statistics.screen <- NULL
        cvPrList.screen <- NULL
        cvClaLust.screen <- NULL
        for (k in 1:length(alphas.default)) {
          CVconf.screen <-
            CVconfusion(
              cvFit.screen[[k]],
              lambda = lambda.screen[[k]],
              Input.D.matrix,
              meta2,
              className = className,
              classesTrain
            )
          cvProb <- CVconf.screen$cvProb
          cv.trueClasses <- CVconf.screen$cv.trueClasses
          cv.predsFactor <- CVconf.screen$cv.predsFactor
          cv.preds <- CVconf.screen$cv.preds
          confus.cv <- CVconf.screen$confus.cv
          
          roc.obj.screen[[k]] <-
            roc(cv.trueClasses, cvProb[, 2], direction = "<")
          cvPrList.screen[[k]] <- cvProb[, 2]
          cvClaLust.screen[[k]] <- cv.trueClasses
          cv <- new(
            "binary",
            trueclass  = as.numeric(cv.trueClasses) - 1,
            predclass =  as.numeric(cv.predsFactor) - 1,
            predprob = cvProb[, 2]
          )
          cv.statistics.screen[[k]] <-
            round(as.data.frame(c(
              APRFMscore(cv),
              AUROC = AUC(cv),
              BRIER = brier(cv)
            ), drop = F), 3)
        }
        names(cvClaLust.screen) <- screen.names
        names(cv.statistics.screen) <- names(cvClaLust.screen)
        cv.statistics.screen.table <-
          do.call("rbind", cv.statistics.screen)
        
        save(cv.statistics.screen.table,
             file = file.path(
               "./results",
               sprintf(
                 "%s_%s_cv.statistics.screen.table.RData",
                 s.Title,
                 CVal
               )
             ))
        mx <-  cv.statistics.screen.table
        ylim <- c(0.5, 1)
        barplot(
          as.matrix(mx),
          beside = T,
          col = col.func(dim(mx)[1], "blue"),
          ylim = ylim,
          xpd = F,
          legend = dimnames(mx)[[1]]
        )
        pdf(file.path(
          "./results",
          sprintf("%s_%s_cv.screen.pdf", s.Title, CVal)
        ),
        width = 12,
        height = 4)
        par(mfrow = c(1, 2))
        mx <-  cv.statistics.screen.table[,-7]
        ylim <- c(0.5, 1)
        barplot(
          as.matrix(mx),
          beside = T,
          col = col.func(dim(mx)[1], "blue"),
          ylim = ylim,
          xpd = F
        )
        mx <-  cv.statistics.screen.table[, 7, drop = F]
        ylim <- c(0, 0.15)
        barplot(
          as.matrix(mx),
          beside = T,
          col = col.func(dim(mx)[1], "blue"),
          ylim = c(0, 1.3 * max(unlist(mx))),
          xpd = F,
          legend = dimnames(mx)[[1]]
        )
        dev.off()
        
        require(ggplot2)
        require(precrec)
        msmdat <-
          mmdata(
            cvPrList.screen,
            cvClaLust.screen,
            modnames = names(cvClaLust.screen) ,
            dsids = 1:length(cvClaLust.screen)
          )
        sscurves <- precrec::evalmod(msmdat)
        plot(sscurves)
        autoplot(sscurves)
        
        pdf(file.path(
          "./results",
          sprintf("%s_%s_cv.roc.screen.pdf", s.Title, CVal)
        ),
        width = 7,
        height = 4)
        autoplot(sscurves)
        dev.off()
      }
    compare.penalizedRegression(
      alphas.default = alphas.default,
      screen.names = screen.names,
      seed = seed
    )
    GlobalOp <- GlobalOp.pre
  }
  res <- list()
  n.fold = sprintf("nfold%s", nFolds)
  Exp. <-
    sprintf("%s.%s.%s.%s.%s.%s.%s",
            Input.D,
            s.Title,
            path,
            CVal,
            n.fold,
            GlobalOp,
            screen.name)
  savePath = file.path(getwd(), "results", Exp.)
  if (!file.exists(sprintf("%s_mergedMeta.RData", savePath))) {
    mergedMeta = list(
      x = Input.D.matrix,
      geneMatrix = devMerged,
      y = meta2[colnames(Input.D.matrix), "class"],
      meta = meta2
    )
    save(mergedMeta, file = sprintf("%s_mergedMeta.RData", savePath))
  }
  
  set.seed(seed)
  
  sel.cvFitList <-
    cv.merged(
      Input.D.matrix,
      meta2,
      yName = className,
      weights = glmnetArgs$weights,
      nRepeats = nRepeats,
      nFolds = nFolds,
      alphas = alphas,
      family = familyName,
      intercept = intercept,
      keep = TRUE,
      GlobalOp = GlobalOp,
      type.measure = type.measure,
      type.min = type.min,
      standardize = F,
      popSize = populations,
      maxiter = generation,
      seed = seed,
      foldid = foldid,
      penalty.factor = penalty.factor,
      ...
    )
  save(sel.cvFitList, file = sprintf("%s_sel.cvFitList.RData", savePath))
  
  
  GlobalOp <- GlobalOp.pre
  require(plot3D)
  if (GlobalOp == "none") {
    Crude3Dmatrix = NULL
    if (CVal == "LOSOCV")
      sel.cvFitList = list(sel.cvFitList)
    for (i in 1:length(sel.cvFitList[[1]])) {
      Crude3Dmatrix = rbind(Crude3Dmatrix,
                            c(
                              alphas[i],
                              sel.cvFitList[[1]][[i]]$lambda.min,
                              min(sel.cvFitList[[1]][[i]]$cvm)
                            ))
    }
    D = as.data.frame(Crude3Dmatrix)
    scatter3D(
      D$V1,
      D$V2,
      D$V3,
      phi = 10,
      theta = 33,
      bty = "g",
      type = "p",
      cex.lab = 1.5,
      xlab = "alpha.grid",
      ylab = "lambda.min",
      zlab = "cvm.min",
      pch = 20,
      cex = c(1.2, 1.2, 1.2),
      xlim = c(0, 1.2),
      ylim = c(-2, 30),
      zlim = c(min(D$V3) * .99, max(D$V3) * 1.02)
    )
    p1 <- recordPlot()
    D.min = D[which(D[, 3] == min(D[, 3])), ]
    alpha.min = which(D[, 3] == min(D[, 3]))
    alpha.sol = D[alpha.min, 1]
    alpha = alpha.sol
    cvFit = sel.cvFitList[[1]][[alpha.min]]
    cvFit.pre <<-  cvFit
    cvFit.pre$fit.preval
    fitResult = cvFit$glmnet.fit
    lambda = cvFit$lambda.min
  }
  if (GlobalOp == "GA") {
    DD = as.data.frame(sel.cvFitList[[3]])
    DD.min = DD[which(DD[, 3] == min(DD[, 3])), ]
    DD.min
    scatter3D(
      DD$V1,
      DD$V2,
      DD$V3,
      phi = 10,
      theta = 35,
      bty = "g",
      type = "p",
      cex.lab = 1.5,
      xlab = "alpha",
      ylab = "lambda.min",
      zlab = "cvm.min",
      pch = 20,
      cex = c(1.2, 1.2, 1.2),
      xlim = c(0, 1.2),
      ylim = c(-2, 30),
      zlim = c(min(DD$V3) * .99, max(DD$V3) * 1.02),
      main = ""
    )
    
    p1 <- recordPlot()
    alpha = round(DD[, 1], 3)
    GA.sol = round(sel.cvFitList[[2]], 3)
    cvFit = sel.cvFitList[[1]][[which(alpha == GA.sol[[1]])[1]]]
    alpha = GA.sol
    fitResult = cvFit$glmnet.fit
    lambda = cvFit$lambda.min
  }
  if (GlobalOp == "epsgo") {
    sumfit <- summary(sel.cvFitList, verbose = TRUE)
    plot(sumfit)
    p0 <- recordPlot()
    sumfit.pre <<- sumfit
    alpha = sumfit$opt.alpha
    cvFit = sumfit$opt.models$model$cvreg
    fitResult = cvFit$glmnet.fit
    lambda = sumfit$opt.lambda
    epsgo.cvFitList <- sel.cvFitList
    Crude3Dmatrix = NULL
    for (i in 1:length(epsgo.cvFitList$model.list)) {
      Crude3Dmatrix = rbind(
        Crude3Dmatrix,
        c(
          epsgo.cvFitList$model.list[[i]]$model$alpha,
          epsgo.cvFitList$model.list[[i]]$model$cvreg$lambda.min,
          min(epsgo.cvFitList$model.list[[i]]$model$cvreg$cvm)
        )
      )
    }
    DDD = as.data.frame(Crude3Dmatrix)
    scatter3D(
      DDD$V1,
      DDD$V2,
      DDD$V3,
      phi = 7,
      theta = 35,
      bty = "g",
      type = "p",
      cex.lab = 1.5,
      xlab = "alpha.grid",
      ylab = "lambda.min",
      zlab = "cvm.min",
      pch = 20,
      cex = c(1.2, 1.2, 1.2),
      xlim = c(0, 1.2),
      ylim = c(-2, 30),
      zlim = c(min(DDD$V3) * .99, max(DDD$V3) * 1.02)
    )
    p1 <- recordPlot()
    DDD.min = DDD[which(DDD[, 3] == min(DDD[, 3])), ]
  }
  
  
  alpha.lambda.sol.temp <- round(c(alpha, lambda), 4)
  cat(sprintf(
    "\n\nalpha.lambda.sol: \n %s\n",
    paste0(alpha.lambda.sol.temp, collapse = ", ")
  ))
  
  coefDf <-
    data.frame(as.matrix(coef(fitResult, s = lambda)[coef(fitResult, s = lambda)[, 1] != 0, ]))
  
  save(coefDf, file = sprintf("%s_coefDf.rda", savePath))
  coefDf <-
    data.frame('Nonzero.Features' = rownames(coefDf)[-1], Coef = coefDf[-1 , 1])
  coefDf <- coefDf[order(coefDf[, "Coef"], decreasing = T),]
  print(data.frame(
    'Nonzero.Features' = substr(coefDf[, 1], 1, 45),
    Coef = coefDf[, 2]
  ))
  attributes(coefDf)$meta <-
    paste0(savePath,
           "_alphaAndlamda_",
           paste(alpha.lambda.sol.temp, collapse = "_"))
  nonzero <- coefDf
  CVconf <-
    CVconfusion(cvFit,
                lambda,
                Input.D.matrix,
                meta2,
                className = className,
                classesTrain)
  
  print(CVconf$confus.cv)
  cvProb <- CVconf$cvProb
  cv.trueClasses <- CVconf$cv.trueClasses
  cv.predsFactor <- CVconf$cv.predsFactor
  cv.preds <- CVconf$cv.preds
  confus.cv <- CVconf$confus.cv
  
  cv <- new(
    "binary",
    trueclass  = as.numeric(cv.trueClasses) - 1,
    predclass =  as.numeric(cv.predsFactor) - 1,
    predprob = cvProb[, 2]
  )
  cv.statistics <-
    round(as.data.frame(c(
      AUROC = AUC(cv), BRIER = brier(cv), APRFMscore(cv)
    ), drop = F), 3)
  rownames(cv.statistics) <- "CV_statistics"
  cv.statistics
  cvPrList <- list()
  cvClList <- list()
  for (i in 1:length(devNames)) {
    inx10 <- meta2$sample[which(meta2$study %in% devNames[i])]
    inx11 <- which(devSampleNames %in% inx10)
    cvPrList[[i]] <- cvProb[, 2][inx10]
    cvClList[[i]] <- cv.trueClasses[inx11]
  }
  cvPrList <- c(cvPrList, list(cvProb[, 2]))
  cvClList <- c(cvClList, list(cv.trueClasses))
  names(cvPrList) <- c(devNames, "Overall")
  names(cvClList) <- c(devNames, "Overall")
  
  require(precrec)
  msmdat <-
    mmdata(
      cvPrList,
      cvClList,
      modnames =  gsub("^s.*_", "", names(cvClList)) ,
      dsids = 1:length(cvClList)
    )
  sscurves.cv <- precrec::evalmod(msmdat)
  sscurves.cv.pre <<- precrec::evalmod(msmdat)
  print(autoplot(sscurves.cv))
  
  p2 <- recordPlot()
  eval.md <- precrec::evalmod(msmdat)
  aucs <- precrec::auc(sscurves.cv)
  prauc.cv <- subset(aucs, curvetypes == "PRC")
  print(knitr::kable(aucs))
  AUPRC.cv = c(prauc.cv$aucs[1:length(prauc.cv$aucs)])
  
  cv.each.val <- list()
  cutoff = .5
  for (i in names(cvClList)) {
    cv.each.val[[i]] <-
      data.frame(
        predsFactor = as.factor(ifelse(cvPrList[[i]] > cutoff, "E", "C")),
        trueClasses = cvClList[[i]],
        predsProb = cvPrList[[i]]
      )
  }
  cv.stat <- extract.stat(cv.each.val)
  
  cv.stat.all <- cv.stat
  print(cv.stat)
  cv.plot.list <-
    list(
      p0 = p0,
      p1 = p1,
      p2 = p2,
      p3 = p3,
      p4 = p4,
      p5 = p5,
      p6 = p6,
      p7 = p7
    )
  
  if (length(which(valCohorts == 1)) == 0) {
    if (Input.D == "PDS") {
      result.list <-
        list(
          Exp. = Exp.,
          savePath = savePath,
          pw.in = pw.in,
          PDS.res = PDS.res,
          insel.discovery = insel,
          coefDf = coefDf,
          alpha = alpha,
          lambda = lambda,
          meta2 = meta2,
          cv.stat.all = cv.stat.all,
          CVconf = CVconf,
          cv.plot.list = cv.plot.list,
          nonzero = nonzero,
          table.cv = out.st2
        )
    } else {
      result.list  <-
        list(
          Exp. = Exp.,
          savePath = savePath,
          coefDf = coefDf,
          alpha = alpha,
          lambda = lambda,
          meta2 = meta2,
          cv.stat.all = cv.stat.all,
          CVconf = CVconf,
          cv.plot.list = cv.plot.list,
          nonzero = nonzero,
          table.cv = out.st2
        )
    }
    
    res.n = sprintf('%s_alpha_%s_lambda_%s_res.RData',
                    Exp.,
                    round(alpha, 3),
                    round(lambda, 3))
    res.path = file.path(getwd(), "results", res.n)
    save(result.list, file = res.path)
    return(result.list)
  }
  
  predValCohorts <-
    function(ematList,
             meta,
             meta2,
             devSampleNames,
             classesTrain,
             alpha,
             lambda,
             weights,
             batchColname = 'study',
             covariateName = NA,
             className = 'class',
             familyName = 'binomial',
             predType = 'response',
             intercept = TRUE,
             val.list = NA,
             type.m = c("m", "PDS"),
             pw.in = NULL,
             gene_ids = NULL,
             min_std = 0.2,
             attempts = 2,
             min_exp = -3,
             preCalpredsList = NULL,
             valBatchCorrectMethod = "combat",
             ...) {
      argss <- list(...)
      penalty.factor <-  argss[['penalty.factor']]
      type.m <- match.arg(type.m)
      if (!is.na(val.list)) {
        dinx = meta[meta[, 'devCohorts'], 'study']
        vinx <- val.list
      } else {
        dinx = meta[meta[, 'devCohorts'], 'study']
        vinx = meta[meta[, 'valCohorts'], 'study']
      }
      if (type.m == "m") {
        predsList = foreach(valinx = vinx) %dopar% {
          if (!is.na(val.list))
            valinx <- unlist(valinx)
          idxValidation = meta2[, 'study'] %in% valinx &
            (meta2[, className] %in% classesTrain)
          valsample = meta2[idxValidation, 'sample']
          
          ematListNow = ematList[c(dinx, valinx)]
          d = .w.mergeStudyData(
            ematListNow,
            meta2,
            batchColname = batchColname,
            covariateName = covariateName,
            merge.mathod = valBatchCorrectMethod,
            ...
          )
          fitResult = glmnet(
            t(d[, devSampleNames]),
            meta2[devSampleNames, className],
            alpha = alpha,
            lambda = lambda,
            weights = weights[devSampleNames],
            family = familyName,
            standardize = FALSE,
            intercept = intercept,
            penalty.factor = penalty.factor
          )
          preds = predict(fitResult,
                          newx = t(d[, valsample]),
                          s = lambda,
                          type = predType)
          list(preds = preds, dev = d[, devSampleNames], d[, valsample])
        }
        names(predsList) = vinx
        return(predsList)
      } else {
        predsList = foreach(valinx = vinx) %dopar% {
          if (!is.na(val.list))
            valinx <- unlist(valinx)
          idxValidation = meta2[, 'study'] %in% valinx &
            (meta2[, className] %in% classesTrain)
          valsample = meta2[idxValidation, 'sample']
          
          ematListNow = ematList[c(dinx, valinx)]
          ematMergedDevVal = .w.mergeStudyData(
            ematListNow,
            meta2,
            batchColname = batchColname,
            covariateName = covariateName,
            merge.mathod = valBatchCorrectMethod,
            ...
          )
          ematMergedDev = ematMergedDevVal[, devSampleNames]
          exp <- ematMergedDevVal
          exp.merged = exp[which(rownames(exp) %in% gene_ids), ]
          dim(exp)
          nor <- rep(F, dim(ematMergedDevVal)[2])
          nor[which(colnames(ematMergedDev) %in% meta2$sample[which(meta2$class == "C")])] <-
            T
          normal.merged <- nor
          
          if (is.null(preCalpredsList)) {
            PDS.result <-
              wrap.PDS(
                exp.merged,
                rownames(exp.merged),
                pw.in$entrez_gene_ids,
                pw.in$pathway,
                normal.merged,
                attempts = attempts,
                min_exp = min_exp,
                min_std = min_std
              )
            PDSmatrix <- do.call(rbind.data.frame, PDS.result$scores)
            print(dim(PDSmatrix))
            d <- as.matrix(PDSmatrix)
            colnames(d) <- colnames(ematMergedDevVal)
          } else {
            d = preCalpredsList[[valinx]][["dev"]]
            d.pree <<- d
          }
          
          fitResult = glmnet(
            t(d[, devSampleNames]),
            meta2[devSampleNames, className],
            alpha = alpha,
            lambda = lambda,
            weights = weights[devSampleNames],
            family = familyName,
            standardize = FALSE,
            intercept = intercept,
            penalty.factor = penalty.factor
          )
          preds = predict(fitResult,
                          newx = t(d[, valsample]),
                          s = lambda,
                          type = predType)
          list(preds = preds, dev = d[, devSampleNames], d[, valsample])
        }
        names(predsList) = vinx
        return(predsList)
      }
    }
}




resamplingOfPamSpearman <-
  function (d = NULL,
            maxK = 7,
            reps = 10,
            pItem = 0.8,
            title = "untitled_consensus_cluster",
            seed = 1234,
            verbose = T) {
    require(cluster)
    # ConsensusClusterPlus(d, seed=1234, clusterAlg = "pam",  maxK = 7,  distance = "spearman", innerLinkage = "complete")  # same results
    clusterAlg = "pam"
    distance = "spearman"
    innerLinkage = "complete"
    finalLinkage = "average"
    pFeature = 1
    corUse = "everything"
    weightsItem = NULL
    weightsFeature = NULL
    ml = NULL
    
    if (!class(d) %in% c("matrix"))
      stop("d must be a matrix")
    set.seed(seed)
    connectivityMatrix <-
      function (clusterAssignments, m, sampleKey) {
        names(clusterAssignments) <- sampleKey
        cls <-
          lapply(unique(clusterAssignments), function(i)
            as.numeric(names(clusterAssignments[clusterAssignments %in% i])))
        for (i in 1:length(cls)) {
          nelts <- 1:ncol(m)
          cl <- as.numeric(nelts %in% cls[[i]])
          updt <- outer(cl, cl)
          m <- m + updt
        }
        return(m)
      }
    sampleCols <-
      function (d,
                pSamp = NULL,
                pRow = NULL,
                weightsItem = NULL,
                weightsFeature = NULL)  {
        space <- ifelse(inherits(d, "dist"), ncol(as.matrix(d)),
                        ncol(d))
        sampleN <- floor(space * pSamp)
        sampCols <- sort(sample(space, sampleN, replace = FALSE,
                                prob = weightsItem))
        this_sample <- sampRows <- NA
        if (inherits(d, "matrix")) {
          if ((!is.null(pRow)) &&
              ((pRow < 1) || (!is.null(weightsFeature)))) {
            space = nrow(d)
            sampleN = floor(space * pRow)
            sampRows = sort(sample(space, sampleN, replace = FALSE,
                                   prob = weightsFeature))
            this_sample <- d[sampRows, sampCols]
            dimnames(this_sample) <- NULL
          } else {
            
          }
        }
        return(list(
          submat = this_sample,
          subrows = sampRows,
          subcols = sampCols
        ))
      }
    
    myPal <-  function (n = 10) {
      seq = rev(seq(0, 255, by = 255 / (n)))
      palRGB = cbind(seq, seq, 255)
      rgb(palRGB, maxColorValue = 255)
    }
    setClusterColors <- function (past_ct, ct, colorU, colorList) {
      newColors = c()
      if (length(colorList) == 0) {
        newColors = colorU[ct]
        colori = 2
      }
      else {
        newColors = rep(NULL, length(ct))
        colori = colorList[[2]]
        mo = table(past_ct, ct)
        m = mo / apply(mo, 1, sum)
        for (tci in 1:ncol(m)) {
          maxC = max(m[, tci])
          pci = which(m[, tci] == maxC)
          if (sum(m[, tci] == maxC) == 1 & max(m[pci,]) ==
              maxC & sum(m[pci,] == maxC) == 1) {
            newColors[which(ct == tci)] = unique(colorList[[1]][which(past_ct ==
                                                                        pci)])
          }
          else {
            colori = colori + 1
            newColors[which(ct == tci)] = colorU[colori]
          }
        }
      }
      return(list(newColors, colori, unique(newColors)))
    }
    
    CDF <- function (ml, breaks = 100) {
      areaK = c()
      cdfDf = list()
      for (i in 2:length(ml)) {
        v = triangle(ml[[i]], mode = 1)
        h = hist(v, plot = FALSE, breaks = seq(0, 1, by = 1 / breaks))
        h$counts = cumsum(h$counts) / sum(h$counts)
        thisArea = 0
        for (bi in 1:(length(h$breaks) - 1)) {
          thisArea = thisArea + h$counts[bi] * (h$breaks[bi + 1] - h$breaks[bi])
          bi = bi + 1
        }
        areaK = c(areaK, thisArea)
        cdfDf[[i]] = list(
          cdf = data.frame(
            K = i,
            x = h$mids,
            y = h$counts
          ),
          hist = list(
            v = v,
            h = h,
            breaks = breaks
          )
        )
      }
      deltaK = areaK[1]
      for (i in 2:(length(areaK))) {
        deltaK = c(deltaK, (areaK[i] - areaK[i - 1]) / areaK[i - 1])
      }
      for (i in 2:(length(deltaK) + 1)) {
        cdfDf[[i]] = c(cdfDf[[i]], list(deltaK = deltaK[i - 1]))
      }
      cdfDf
    }
    triangle <- function (m, mode = 1) {
      n = dim(m)[1]
      nm = matrix(0, ncol = n, nrow = n)
      fm = m
      nm[upper.tri(nm)] = m[upper.tri(m)]
      fm = t(nm) + nm
      diag(fm) = diag(m)
      nm = fm
      nm[upper.tri(nm)] = NA
      diag(nm) = NA
      vm = m[lower.tri(nm)]
      if (mode == 1) {
        return(vm)
      } else if (mode == 3) {
        return(fm)
      } else if (mode == 2) {
        return(nm)
      }
    }
    ccRun <-
      function (d = d,
                maxK = NULL,
                repCount = NULL,
                diss = inherits(d, "dist"),
                pItem = NULL,
                pFeature = NULL,
                innerLinkage = NULL,
                distance = NULL,
                clusterAlg = NULL,
                weightsItem = NULL,
                weightsFeature = NULL,
                verbose = NULL,
                corUse = NULL) {
        m = vector(mode = "list", repCount)
        ml = vector(mode = "list", maxK)
        n <- ncol(d)
        mCount = mConsist = matrix(c(0), ncol = n, nrow = n)
        ml[[1]] = c(0)
        if (is.null(distance))
          distance <- "euclidean"
        acceptable.distance <- c(
          "euclidean",
          "maximum",
          "manhattan",
          "canberra",
          "binary",
          "minkowski",
          "pearson",
          "spearman"
        )
        main.dist.obj <-
          as.dist(1 - cor(d, method = distance, use = corUse))
        attr(main.dist.obj, "method") <- distance
        
        for (i in 1:repCount) {
          if (verbose) {
            message(paste("random subsample", i))
          }
          sample_x = sampleCols(d, pItem, pFeature, weightsItem,
                                weightsFeature)
          
          this_dist = NA
          boot.cols <- sample_x$subcols
          this_dist <- as.matrix(main.dist.obj)[boot.cols,
                                                boot.cols]
          this_dist <- as.dist(this_dist)
          attr(this_dist, "method") <- attr(main.dist.obj,
                                            "method")
          
          this_cluster = NA
          mCount <- connectivityMatrix(rep(1, length(sample_x[[3]])),
                                       mCount, sample_x[[3]])
          for (k in 2:maxK) {
            if (verbose) {
              message(paste("  k =", k))
            }
            if (i == 1) {
              ml[[k]] = mConsist
            }
            this_assignment = NA
            if (clusterAlg == "pam") {
              this_assignment <- pam(
                x = this_dist,
                k,
                diss = TRUE,
                metric = distance,
                cluster.only = TRUE
              )
            } else {
              this_assignment <- get(clusterAlg)(this_dist, k)
            }
            ml[[k]] <-
              connectivityMatrix(this_assignment, ml[[k]], sample_x[[3]])
          }
        }
        res = vector(mode = "list", maxK)
        for (k in 2:maxK) {
          tmp = triangle(ml[[k]], mode = 3)
          tmpCount = triangle(mCount, mode = 3)
          res[[k]] = tmp / tmpCount
          res[[k]][which(tmpCount == 0)] = 0
        }
        message("end fraction")
        return(res)
      }
    
    if (is.null(ml) == TRUE) {
      ml <-
        ccRun(
          d = d,
          maxK = maxK,
          repCount = reps,
          diss = inherits(d, "dist"),
          pItem = pItem,
          pFeature = pFeature,
          innerLinkage = innerLinkage,
          clusterAlg = clusterAlg,
          weightsFeature = weightsFeature,
          weightsItem = weightsItem,
          distance = distance,
          verbose = verbose,
          corUse = corUse
        )
    }
    res = list()
    res[[1]] = 1:2
    for (tk in 2:maxK) {
      if (verbose) {
        message(paste("consensus ", tk))
      }
      fm = ml[[tk]]
      hc = hclust(as.dist(1 - fm), method = finalLinkage)
      message("clustered")
      ct = cutree(hc, tk)
      names(ct) = colnames(d)
      if (class(d) == "dist") {
        names(ct) = colnames(as.matrix(d))
      }
      res[[tk]] = list(
        consensusMatrix = fm,
        consensusTree = hc,
        consensusClass = ct,
        ml = ml[[tk]]
      )
    }
    cdf = CDF(ml)
    res[["cdf"]] = cdf
    return(res)
  }






consensusP <-
  function(d,
           maxK = 7,
           dendColFunc = NULL,
           seed = 1234,
           reps = 10,
           verbose = T,
           heatmap = T,
           fixedK = NULL,
           ...) {
    try(dev.off(), silent = T)
    consensus.pam = resamplingOfPamSpearman(
      d,
      maxK = maxK,
      reps = reps,
      verbose = verbose,
      seed = seed,
      ...
    )
    cdfDf2 = do.call("rbind", lapply(consensus.pam[["cdf"]], function(x)
      x$cdf))
    cdfDf2 = factorizeDf(cdfDf2)
    cdfFnLower <- cdfFnUpper <- PAC <- list()
    for (i in 1:length(levels(cdfDf2$K)) + 1) {
      cdfFnLower[[i]] <-
        cdfDf2[which(cdfDf2$K == i & cdfDf2$x == 0.095), , drop = F]$y
      cdfFnUpper[[i]] <-
        cdfDf2[which(cdfDf2$K == i & cdfDf2$x == 0.895), , drop = F]$y
      PAC[[i]] = round(cdfFnUpper[[i]] - cdfFnLower[[i]], 2)
    }
    orderInx <- order(do.call("c", PAC))
    for (k in 2:length(PAC))
      if (PAC[[k]] == min(do.call("c", PAC)))
        optK = k
    if (!is.null(fixedK))
      optK = fixedK
    list(res = consensus.pam, optK = optK)
  }


