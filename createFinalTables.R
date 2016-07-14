library(dplyr)
ENetModelDir <- "F:\\Jie\\MS\\02_Code\\MS_InitModel\\Results\\2016-07-12 14.54.21\\"
GlmModelDir <- "F:\\Jie\\MS\\02_Code\\MS_NonRegularisedGLM\\Results\\2016-07-14 03.08.07\\"

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


# 1.	Table type 1: AUC by model type 
# a.	Columns:
#   i.	Elastic-net (EN) LR model with ~80 variables
#   ii.	EN LR model with ~10 variables
#   iii.	Unconstrained LR model with ~10 variables
#   iv.	[AUC and 95% CI for each model]. 
# b.	Rows: One row per outcome 
# c.	Total: 3 columns * 6 outcomes. Each cell (i.e. row / column combination has AUC and 95%ci).
# d.	Variable selection: 
#   i.	Initial ~80 variables to be based on domain knowledge
#   ii.	Top ~10 variables to be based on variable importance 

generateTables <- function(coh, iRepeat){
  iRepeatDir <- paste0(ENetModelDir, iRepeat, '\\')
  cohDir <- paste0(iRepeatDir, coh, '\\')
  
  EnetTestAuc <- read.table(paste0(iRepeatDir, 'AvAUCs_test.csv')
                        , sep=','
                        , header = F
                        , stringsAsFactors = F)
  names(EnetTestAuc) <- c('Cohort', 'Outcome', 'AUC on test')
  
  
  
# 2.	Table type 2: Variable importance for full EN model based on ~80 variables (sorted by descending variable importance)
#   a.	Note there are separate tables for each outcome i.e. Tables 2A-2F.
#   b.	Columns:
#     i.	Variable rank [start at 1 and go to ~80]
#     ii.	Variable description 
#     iii.	Score 
#     iv.	Number of times retained 
#     v.	Average coefficient 
#     vi.	Average odds ratio [footnote describing meaning!] 
  
  for(iOutcome in outcomeList){
    cat(iOutcome, '\n\n')
    outcomeDir <- paste0(cohDir, iOutcome, '\\')
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
    rank <- 1:nrow(avgCoefRank)
    OR <- exp(avgCoefRank$Coef)
    desc <- varDesc[match(rownames(avgRank), varDesc[, 1]), 2]
    n.retained <-apply(apply(coefAllFoldAlpha, 2, function(x)x!=0), 1, sum)
    n.retained <- paste0(n.retained, '/', ncol(coefAllFoldAlpha))
    
    tb2 <- cbind(avgCoefRank, Rank=rank, Desc=desc, OR=OR, n_retained=n.retained)
    names(tb2)[1] <- c("Score")
    
    
    
    tb2 <- tb2[1:10, c('Rank', 'Desc', 'Score', 'n_retained', 'Coef', 'OR')]
    
    resultCohDir <- paste0(resultDir, iRepeat, '\\', coh, '\\')
    if(!dir.exists(resultCohDir))
      dir.create(resultCohDir, recursive = T)
    write.table(tb2
                , paste0(resultCohDir, 'Table2_', iOutcome, '.csv')
                , sep=','
                , row.names = T
                )
    
  }  
  
  
  
  
}







