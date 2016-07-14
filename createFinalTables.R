library(dplyr)
library(plyr)
ENetModelDir <- "F:\\Jie\\MS\\02_Code\\MS_InitModel\\Results\\2016-07-12 14.54.21\\"
subVarsENetModelDir <- "F:\\Jie\\MS\\02_Code\\MS_InitModel\\Results\\2016-07-14 07.17.09\\"

GlmModelDir <- "F:\\Jie\\MS\\02_Code\\MS_NonRegularisedGLM\\Results\\2016-07-14 06.10.58\\"

varDescDir <- "F:\\Jie\\MS\\01_Data\\ModelData\\data4Model\\"
varDescFile <- "lookup_20160714.csv"

outcomeList <- c("relapse_fu_any_01", "edssprog", "edssconf3",
                           "relapse_or_prog", "relapse_and_prog", "relapse_or_conf")

cohList <- c("Cmp")

n.repeat <- 1

timeStamp <- as.character(Sys.time())
timeStamp <- gsub(":", ".", timeStamp)  # replace ":" by "."
resultDir <- paste("./Results/", timeStamp, "/", sep = '')
dir.create(resultDir, showWarnings = TRUE, recursive = TRUE, mode = "0777")



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
  
  testAuc_glm <- ldply(lapply(outcomeList, function(iOutcome){
    outcomeDir <- paste0(GlmModelDir, coh, '\\', iOutcome, '\\')
    aucCi <- read.table(paste0(outcomeDir, 'ci_auc_glm.csv')
                        , sep=','
                        , header = F)
    return(aucCi)
  }), rbind)
  
  
  # calculate the acu confidence interval using formula
  # SE (AUC)= √((AUC ( 1-AUC)+(N_1-1)  (Q_1- 〖AUC〗^2 )+(N_2-1)(Q_2-〖 AUC〗^2) )/(N_1 N_2 ))
  getAucCi <- function(iOutcome, dir){
    resp <- read.table(paste0(dir, '\\', iOutcome, '\\Cmp_data_for_model.csv')
                       , sep=','
                       , header = T)[, "y"]
    n.pos <- sum(resp)
    n.neg <- sum(1-resp)
    if(dir==cohDir){
      auc <- EnetTestAuc[EnetTestAuc$Outcome==iOutcome, 'AUC on test']
      
    }else if(dir == subVarsCohDir){
      auc <- subVarsEnetTestAuc[subVarsEnetTestAuc$Outcome==iOutcome, 'AUC on test']
    }
    q1 <- auc/(2-auc)
    q2 <- (2*auc^2)/(1+auc)
    
    SE <- sqrt((auc*(1-auc)+(n.pos-1)*(q1-auc^2)+(n.neg-1)*(q2-auc^2))/(n.pos*n.neg))
    
    aucCiDiff <- 1.96*SE
    return(aucCiDiff)
  }
  EnetAucCi <- unlist(lapply(outcomeList, getAucCi, cohDir))
  subVarsEnetAucCi <- unlist(lapply(outcomeList, getAucCi, subVarsCohDir))
  GlmAucCi <- (testAuc_glm[3]-testAuc_glm[1])/2
  
  tb1 <- data.frame(Cohort=EnetTestAuc$Cohort
                    , Outcome=EnetTestAuc$Outcome
                    , Elastic_net=paste0(round(EnetTestAuc$`AUC on test`, 2), '+/-', round(EnetAucCi, 2))
                    , Elastic_net_10_variables=paste0(round(subVarsEnetTestAuc$`AUC on test`, 2), '+/-', round(subVarsEnetAucCi, 2))
                    , Unconstrained_LR=paste0(round(testAuc_glm[, 2], 2), '+/-', round(GlmAucCi[, 1], 2))
                    )
  write.table(tb1
              , paste0(resultCohDir, 'Table1.csv')
              , sep=','
              , row.names=F
               )
  
  for(iOutcome in outcomeList){
    cat(iOutcome, '\n\n')
    outcomeDir <- paste0(cohDir, iOutcome, '\\')
    
    # 2.	Table type 2: Variable importance for full EN model based on ~80 variables (sorted by descending variable importance)
    #   a.	Note there are separate tables for each outcome i.e. Tables 2A-2F.
    #   b.	Columns:
    #     i.	Variable rank [start at 1 and go to ~80]
    #     ii.	Variable description 
    #     iii.	Score 
    #     iv.	Number of times retained 
    #     v.	Average coefficient 
    #     vi.	Average odds ratio [footnote describing meaning!] 
    
    avgCoef <- read.table(paste0(outcomeDir, "av_coefs_Cmp.csv")
                          , sep=','
                          , header = T
                          , stringsAsFactors = F)
    
    avgRank <- read.table(paste0(outcomeDir, "av_ranking_Cmp.csv")
                          , sep=','
                          , header = T
                          , stringsAsFactors = F)
    
    varDesc <- read.table(paste0(varDescDir, varDescFile)
                          , sep=','
                          , header = T
                          , stringsAsFactors = F)
    
    coefAllFoldAlpha <- read.table(paste0(outcomeDir, "coefs_Cmp.csv")
                                   , sep=','
                                   , header = T
                                   , stringsAsFactors = F)
    
    # avgRank <- cbind(vars=rownames(avgRank), avgRank)
    avgCoefRank <- cbind(score=avgRank, Coef=avgCoef[match(rownames(avgRank), rownames(avgCoef)),])
    avgCoefRank <- avgCoefRank[order(avgCoefRank$x),]
    rank <- 1:nrow(avgCoefRank)
    OR <- exp(avgCoefRank$Coef)
    desc <- varDesc[match(rownames(avgRank), varDesc[, 1]), 2]
    n.retained <-apply(apply(coefAllFoldAlpha, 2, function(x)x!=0), 1, sum)
    n.retained <- paste0(n.retained, '/', ncol(coefAllFoldAlpha))
    
    tb2 <- cbind(avgCoefRank, Rank=rank, Desc=desc, OR=OR, n_retained=n.retained)
    names(tb2)[1] <- c("Score")
    
    
    
    tb2 <- tb2[1:10, c('Rank', 'Desc', 'Score', 'n_retained', 'Coef', 'OR')]
    
    write.table(tb2
                , paste0(resultCohDir, 'Table2_', iOutcome, '.csv')
                , sep=','
                , row.names = T
                )
    
    #   3.	Table type 3: Odds ratio for unconstrained LR based on most important ~10 variables
    #     a.	One table for each outcome e.g. Table 3A-3F.
    #     b.	Columns:
    #       i.	Variable description 
    #     ii.	Odds ratio 
    #     iii.	95% CI for odds ratio 
    #     iv.	P-value 
    cohDir_glm <- paste0(GlmModelDir, coh, '\\', iOutcome, '\\')
    coefInf_GLM <- read.table(paste0(cohDir_glm, 'coef_info.csv')
                              , sep=','
                              , header = T
                              , stringsAsFactors = F
                              , check.names =F)[-1, ]
    
    desc <- varDesc[match(coefInf_GLM[, 1], varDesc[, 1]), 2]
    
    # orDiff <- coefInf_GLM[, grepl("odds_97.5", names(coefInf_GLM))]-coefInf_GLM[, grepl("odds_2.5", names(coefInf_GLM))]
    
    # orInf <- paste0(coefInf_GLM[, 'odds'], "+/-", orDiff/2)
    tb3 <- data.frame(Variable=coefInf_GLM[, 1]
                      , Desc=desc
                      , OR=coefInf_GLM$odds
                      , OR_2.5=coefInf_GLM$`odds_2.5%`
                      , OR_97.5=coefInf_GLM$`odds_97.5%`
                      , Pvalue=coefInf_GLM$`Pr(>|z|)`)
    
    write.table(tb3
                , paste0(resultCohDir, 'Table3_', iOutcome, '.csv')
                , sep=','
                , row.names=F)
    
  }  
  
  
  
  
}


for(iRepeat in 1:n.repeat){
  cat(iRepeat, '\n\n')
  for(coh in cohList){
    cat(coh, '\n\n')
    generateTables(coh, iRepeat)  
  }
  
}




