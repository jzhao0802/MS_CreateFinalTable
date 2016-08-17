rm(list=ls())
library(dplyr)
library(plyr)

options(scipen=999)
ENetModelDir <- 
  "F:\\Jie\\MS\\03_Result\\2016-08-08\\2016-08-08 08.19.05\\"
subVarsENetModelDir <- 
  "F:\\Jie\\MS\\03_Result\\2016-08-08\\2016-08-08 09.24.44\\"
GlmModelDir <- 
  "F:\\Jie\\MS\\03_Result\\2016-08-09\\2016-08-09 06.11.50\\"

varDescDir <- "F:/Jie/MS/01_Data/ModelData/data4Model/"
varDescFile <- "lookup_20160816.csv"

outcomeList <- c("relapse_fu_any_01", "edssprog", "edssconf3",
                           "relapse_or_prog", "relapse_and_prog", "relapse_or_conf")

cohList <- c('Cmp')

n.repeat <- 1
BpartOutput <- F

timeStamp <- as.character(Sys.time())
timeStamp <- gsub(":", ".", timeStamp)  # replace ":" by "."
resultDir <- paste("./Results/", timeStamp, "/", sep = '')
dir.create(resultDir, showWarnings = TRUE, recursive = TRUE, mode = "0777")


myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}

# calculate the acu confidence interval using formula
# SE (AUC)= √((AUC ( 1-AUC)+(N_1-1)  (Q_1- (AUC)^2 )+(N_2-1)(Q_2-(AUC)^2) )/(N_1 N_2 ))
getCohortNumAndRespPosNum <- function(outcomeName, coh, dir){
  resp <- read.table(paste0(dir, outcomeName, '\\', coh, '_data_for_model.csv')
                     , sep=','
                     , header = T)[, "y"]
  n <- length(resp)
  resp.pos.n <- sum(resp)
  num_summary <- c(n, resp.pos.n)
  return(num_summary)
}
getAucCi <- function(coh, outcomeName, dir, cohDir, EnetTestAuc, subVarsCohDir, subVarsEnetTestAuc)
{
  resp <- read.table(paste0(dir, outcomeName, '\\', coh, '_data_for_model.csv')
                     , sep=','
                     , header = T)[, "y"]
  n.pos <- sum(resp)
  n.neg <- sum(1-resp)
  if(dir==cohDir){
    auc <- EnetTestAuc[EnetTestAuc$Outcome==outcomeName, 'AUC on test']
  
  }else if(dir == subVarsCohDir){
    auc <- subVarsEnetTestAuc[subVarsEnetTestAuc$Outcome==outcomeName, 'AUC on test']
  }
  q1 <- auc/(2-auc)
  q2 <- (2*auc^2)/(1+auc)
  
  SE <- sqrt((auc*(1-auc)+(n.pos-1)*(q1-auc^2)+(n.neg-1)*(q2-auc^2))/(n.pos*n.neg))
  
  aucCiDiff <- 1.96*SE
  return(aucCiDiff)
}

getQuintile <- function(outcomeName, dir, coh, inFile, n.bkt){
  subDir <- paste0(dir, outcomeName, '\\')
  predLab <- read.table(paste0(subDir, coh, inFile)
                        , sep=','
                        , header = T
                        , stringsAsFactors = F)[, -1] 
  # breaks <- seq(0, 1, 1/n.bkt)
  breaks <- quantile(predLab[, 1], probs=seq(0,1, by=1/n.bkt), na.rm=TRUE)
  bucket <- cut(predLab[, 1], breaks=breaks, include.lowest=F,right=T)
  predLabBkt <-cbind(bucket, predLab)
  names(predLabBkt) <- c("Bucket", 'Prediction', 'Label')
  quintile <- predLabBkt %>%
    group_by(Bucket) %>%
    {
      detach("package:plyr", character.only = T)
      library('dplyr')
      .
    } %>%
    summarize_each(., funs(getPosRate(.)), one_of("Label"))
  bktLevels <- sort(levels(bucket))
  quintileAllBkt <- data.frame(bktLevels
                          , format(
                            round(quintile[match(bktLevels, quintile$Bucket), -1], 2)
                            , nsmall = 2))
  library(plyr)
  # return(as.vector(as.data.frame(quintileAllBkt)[, 2]))
  return(quintileAllBkt)
}

getPosRate <- function(vct){
  return(sum(vct, na.rm = T)/length(vct))  
}

generateTables <- function(coh, iRepeat){
  #   1.	Table type 1: AUC by model type 
  #   a.	Columns:
  #     i.	Elastic-net (EN) LR model with ~80 variables
  #     ii.	EN LR model with ~10 variables
  #     iii.	Unconstrained LR model with ~10 variables
  #     iv.	[AUC and 95% CI for each model]. 
  #   b.	Rows: One row per outcome 
  #   c.	Total: 3 columns * 6 outcomes. Each cell (i.e. row / column combination has AUC and 95%ci).
  #   d.	Variable selection: 
  #     i.	Initial ~80 variables to be based on domain knowledge
  #     ii.	Top ~10 variables to be based on variable importance 
  # 
  iRepeatDir <- paste0(ENetModelDir, iRepeat, '\\')
  cohDir <- paste0(iRepeatDir, coh, '\\')
  
  subVarsiRepeatDir <- paste0(subVarsENetModelDir, iRepeat, '\\')
  subVarsCohDir <- paste0(subVarsiRepeatDir, coh, '\\')
  
  resultCohDir <- paste0(resultDir, iRepeat, '\\', coh, '\\')
  if(!dir.exists(resultCohDir))
    dir.create(resultCohDir, recursive = T)
  
  EnetTestAuc <- read.table(paste0(iRepeatDir, 'AvAUCs_test.csv')
                        , sep=','
                        , header = F
                        , stringsAsFactors = F)
  subVarsEnetTestAuc <- read.table(paste0(subVarsiRepeatDir, 'AvAUCs_test.csv')
                                   , sep=','
                                   , header = F
                                   , stringsAsFactors = F)
  names(EnetTestAuc) <- c('Cohort', 'Outcome', 'AUC on test')
  names(subVarsEnetTestAuc) <- c('Cohort', 'Outcome', 'AUC on test')
  
  
  EnetTestAuc <- filter(EnetTestAuc, Cohort==coh)
  subVarsEnetTestAuc <- filter(subVarsEnetTestAuc, Cohort==coh)
  
  testAuc_glm <- ldply(lapply(outcomeList, function(iOutcome){
    outcomeDir <- paste0(GlmModelDir, coh, '\\', iOutcome, '\\')
    aucCi <- read.table(paste0(outcomeDir, 'ci_auc_glm.csv')
                        , sep=','
                        , header = F)
    return(aucCi)
  }), rbind)
  
  rownames(testAuc_glm) <- outcomeList
  testAuc_glm <- testAuc_glm[match(EnetTestAuc$Outcome, outcomeList),]
  
  EnetAucCi <- unlist(lapply(outcomeList, getAucCi, coh=coh, 
                             dir=cohDir, cohDir=cohDir, 
                             EnetTestAuc=EnetTestAuc, subVarsCohDir=subVarsCohDir, 
                             subVarsEnetTestAuc=subVarsEnetTestAuc))
  subVarsEnetAucCi <- unlist(lapply(outcomeList, getAucCi, coh=coh,
                                    dir=subVarsCohDir, cohDir=cohDir, 
                                    EnetTestAuc=EnetTestAuc, subVarsCohDir=subVarsCohDir,
                                    subVarsEnetTestAuc=subVarsEnetTestAuc))
  GlmAucCi <- (testAuc_glm[3]-testAuc_glm[1])/2
  
  num_summary <- ldply(lapply(outcomeList, getCohortNumAndRespPosNum, dir=cohDir, coh=coh), quickdf)
  num_summary <- num_summary[match(EnetTestAuc$Outcome, outcomeList), ]
  names(num_summary) <- c("n.cohort","resp.pos.n")
  tb1 <- data.frame(Outcome=EnetTestAuc$Outcome
                    , n=num_summary[, 1]
                    , resp.pos.n=num_summary[, 2]
                    , Elastic_net=paste0(format(round(EnetTestAuc$`AUC on test`, 2), nsmall = 2)
                                         , '+/-'
                                         , format(round(EnetAucCi, 2), nsmall = 2))
                    , Elastic_net_10_variables=paste0(format(round(subVarsEnetTestAuc$`AUC on test`, 2), nsmall = 2)
                                                      , '+/-'
                                                      , format(round(subVarsEnetAucCi, 2), nsmall = 2))
                    , Unconstrained_LR=paste0(format(round(testAuc_glm[, 2], 2), nsmall = 2)
                                              , '+/-'
                                              , format(round(GlmAucCi[, 1], 2), nsmall = 2))
                    )
  tb1[tb1$resp.pos.n < 100, -match(c('Outcome', 'n', 'resp.pos.n'), names(tb1))] <- NA
  # change the names of outcomes
  lookupTabel <- data.frame(oldNm = c("relapse_fu_any_01"
                                      , "edssprog"
                                      , "edssconf3"
                                      , "relapse_or_prog"
                                      , "relapse_and_prog"
                                      , "relapse_or_conf")
                            , newNm = c("Relapse"
                                        , "EDSS Progression"
                                        , "Confirmed EDSS Progression"
                                        , "Relapse or EDSS Progression"
                                        , "Relapse and EDSS Progression"
                                        , "Relapse or confirmed EDSS Progression")
                            , flag = c('a', 'b', 'c', 'd', 'e', 'f')
                            )
  tb1$Outcome <- lookupTabel$newNm[match(tb1$Outcome, lookupTabel$oldNm)]
  
  # change the colnames of tb1
  names(tb1) <- c("Outcome"
                  , "Number of Patients in This Cohort"
                  , "Number of Patients with Positive Outcome"
                  , "Logistic regression with Elastic-Net, extended variable list (+/- 95% CI)"
                  , "Logistic regression with Elastic-Net, most important ten variables (+/- 95% CI)"
                  , "Standard logistic regression, most important ten variables (+/- 95% CI)")
  
  
  write.table(tb1
              , paste0(resultCohDir, 'Table1_old.csv')
              , sep=','
              , row.names=F
               )
  varsWithoutDescLst <- list()
  i=1
  for(iOutcome in outcomeList){
    cat(iOutcome, '\n\n')
    outcomeDir <- paste0(cohDir, iOutcome, '\\')
    
    # 2.	Table type 2: Variable importance for full EN model based on ~48 variables (sorted by descending variable importance)
    #   a.	Note there are separate tables for each outcome i.e. Tables 2A-2F.
    #   b.	Columns:
    #     i.	Variable rank [start at 1 and go to ~48]
    #     ii.	Variable description 
    #     iii.	Score 
    #     iv.	Number of times retained 
    #     v.	Average coefficient 
    #     vi.	Average odds ratio [footnote describing meaning!] 
    
    avgCoef <- read.table(paste0(outcomeDir, "av_coefs_", coh, ".csv")
                          , sep=','
                          , header = T
                          , stringsAsFactors = F)
    
    avgRank <- read.table(paste0(outcomeDir, "av_ranking_", coh, ".csv")
                          , sep=','
                          , header = T
                          , stringsAsFactors = F)
    
    varDesc <- read.csv(paste0(varDescDir, varDescFile)
                          , sep=','
                          , header = T
                          , stringsAsFactors = F)
    
    coefAllFoldAlpha <- read.table(paste0(outcomeDir, "coefs_", coh, ".csv")
                                   , sep=','
                                   , header = T
                                   , stringsAsFactors = F)
    
    # avgRank <- cbind(vars=rownames(avgRank), avgRank)
    avgCoefRank <- cbind(score=avgRank, Coef=avgCoef[match(rownames(avgRank), rownames(avgCoef)),])
    avgCoefRank <- format(round(avgCoefRank[order(avgCoefRank$x),], 3), nsmall=3)
    rank <- 1:nrow(avgCoefRank)
    OR <- format(round(exp(as.numeric(avgCoefRank$Coef)),3), nsmall = 3)
    
    varsWithoutDesc <- rownames(avgRank)[is.na(match(rownames(avgRank), varDesc[, 1]))] 
    varsWithoutDescLst[[i]] <- varsWithoutDesc
    i <- i+1
    desc <- varDesc[match(rownames(avgRank), varDesc[, 1]), 2]
    coefAllFoldAlpha <- coefAllFoldAlpha[match(rownames(avgRank), rownames(coefAllFoldAlpha)), ]
    n.retained <-apply(apply(coefAllFoldAlpha, 2, function(x)x!=0), 1, sum)
    n.retained <- paste0(n.retained, '/', ncol(coefAllFoldAlpha))
    
    tb2 <- cbind(avgCoefRank, Rank=rank, Desc=desc, OR=OR, n_retained=n.retained)
    names(tb2)[1] <- c("Score")
    
    
    
    tb2 <- tb2[, c('Rank', 'Desc', 'n_retained', 'Coef', 'OR')]
    colnames(tb2) <- c('Variable Rank (based on average rank across all models)'
                       , 'Variable Description' 
                       , 'Number of Times Variable Retained'
                       , 'Average Coefficient'
                       , 'Average Odds Ratio*')
    if(num_summary$resp.pos.n[match(iOutcome, outcomeList)] >= 100 | BpartOutput==F){
      write.table(tb2
                  , paste0(resultCohDir
                           , 'Table3'
                           , lookupTabel$flag[match(iOutcome, lookupTabel$oldNm)]
                           , '_'
                           , lookupTabel$newNm[match(iOutcome, lookupTabel$oldNm)]
                           , '_old.csv')
                  , sep=','
                  # , col.names=NA
                  , row.names = F
      )
      
    }
    
    #   3.	Table type 3: Odds ratio for unconstrained LR based on most important 10 variables
    #     a.	One table for each outcome e.g. Table 3A-3F.
    #     b.	Columns:
    #       i.	Variable description 
    #     ii.	Odds ratio 
    #     iii.	95% CI for odds ratio 
    #     iv.	P-value 
    cohDir_glm <- paste0(GlmModelDir, coh, '\\', iOutcome, '\\')
    coefInf_GLM_lst <- myTryCatch(read.table(paste0(cohDir_glm, 'coef_info.csv')
                              , sep=','
                              , header = T
                              , stringsAsFactors = F
                              , check.names =F)[-1, ])
    
    if(!is.null(coefInf_GLM_lst$value)){
      coefInf_GLM <- coefInf_GLM_lst$value
      desc <- varDesc[match(coefInf_GLM[, 1], varDesc[, 1]), 2]
      descOrder <- varDesc[match(rownames(avgRank), varDesc[, 1]), 2]
      # orDiff <- coefInf_GLM[, grepl("odds_97.5", names(coefInf_GLM))]-coefInf_GLM[, grepl("odds_2.5", names(coefInf_GLM))]
      
      # orInf <- paste0(coefInf_GLM[, 'odds'], "+/-", orDiff/2)
      #     tb3 <- data.frame(Desc=desc
      #                       , OR=coefInf_GLM$odds
      #                       , OR_2.5=coefInf_GLM$`odds_2.5%`
      #                       , OR_97.5=coefInf_GLM$`odds_97.5%`
      #                       , Pvalue=coefInf_GLM$`Pr(>|z|)`)
      
      OR_CI <- paste0("[", format(round(coefInf_GLM$`odds_2.5%`, 3), nsmall = 3)
                      , ",", format(round(coefInf_GLM$`odds_97.5%`, 3), nsmall = 3)
                      , "]")
      tb3 <- data.frame(desc, format(round(coefInf_GLM$odds, 3), nsmall = 3), OR_CI
                        , Pvalue=format(round(coefInf_GLM$`Pr(>|z|)`, 3), nsmall=3))
      
      tb3$Pvalue <- as.character(tb3$Pvalue)
      tb3$Pvalue[tb3$Pvalue == "0.000"] = "<0.001"
      
      tb3Order <- tb3[match(descOrder[1:10], desc),]
      colnames(tb3Order) <- c(
        "Variable Description",
        "Odds Ratio",
        "95% CI for Odds Ratio",
        "P-Value"
      )
      
      if(num_summary$resp.pos.n[match(iOutcome, outcomeList)] >= 100 | BpartOutput==F){
        write.table(tb3Order
                    , paste0(resultCohDir
                             , 'Table4'
                             , lookupTabel$flag[match(iOutcome, lookupTabel$oldNm)]
                             , '_'
                             , lookupTabel$newNm[match(iOutcome, lookupTabel$oldNm)]
                             , '_old.csv')
                    , sep=','
                    # , col.names=NA
                    , row.names = F
        )
        
      }
        
    }
    
  }  
  
  
#   4.	Table type 4: Actual outcomes by quintile of predicted risk score based on EN LR model for top ~10 variables 
#   a.	Columns:
#     i.	Quintile group [goes from 1 to 5]
#     ii.	Outcome 1
#     iii.	Outcome 2
#     iv.	Outcome 3
#     v.	Outcome 4
#     vi.	Outcome 5
#     vii.	Outcome 6
#     viii.	Each cell shows the column % to the nearest 1decimal place. Also put a row at the bottom showing 100% for each column so it is clear that this refers to column percents.
#     

  # read in the prediction and label of the full data
  
  
  subVarsEnetQuintile <- lapply(outcomeList
                                      ,getQuintile
                                      ,dir=subVarsCohDir
                                      ,coh=coh
                                      ,inFile="_probsAndLabels.csv"
                                      ,n.bkt=5)  
  
  tb4 <- as.data.frame(t(ldply(lapply(subVarsEnetQuintile, function(X)X[, 2]), rbind)))
#   rownames(tb4) <- paste0("Group", 1:nrow(tb4))
  # colnames(tb4) <- outcomeList
  tb4 <- cbind(Group=1:nrow(tb4),tb4)
  # reverse the order of group number (decreasing)
  tb4 <- tb4[order(tb4$Group, decreasing = T),]
  # name back the group number
  tb4$Group <- 1:nrow(tb4)
  
  colnames(tb4) <- c('Group', outcomeList)
  
  tb4[, -1][, num_summary$resp.pos.n < 100] <- NA
  # rename the outcome names
  colnames(tb4) <- c('Quintile Group'
                     , as.character(lookupTabel$newNm[match(outcomeList, lookupTabel$oldNm)])
                     )
    
  write.table(tb4
              , paste0(resultCohDir, 'Table2_old.csv')
              , sep=','
              , row.names=F
              # , col.names = NA
  )
  return(varsWithoutDescLst)
  
}

tempLst <- list()
i=1
for(iRepeat in 1:n.repeat){
  cat(iRepeat, '\n\n')
  for(coh in cohList){
    cat(coh, '\n\n')
    tempLst[[i]] <- generateTables(coh, iRepeat)  
    i <- i+1
  }
  
}

unique(unlist(lapply(tempLst, function(elm){
  x <- unique(unlist(elm))
  return(x)
})))



