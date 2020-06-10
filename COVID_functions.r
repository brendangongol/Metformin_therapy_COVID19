

#### Frequency counting ############################################################################
####################################################################################################
freqsdt <- function(DTstr, xstr, percent=TRUE)
{    
    if (percent)
        return(eval(parse(text=sprintf('%s[,.(frequency=.N),.(%s)]', DTstr, xstr)))[
            order(-frequency)][,percent:=100*frequency/sum(frequency)]) 
    else
        return(eval(parse(text=sprintf('%s[,.(frequency=.N),.(%s)]', DTstr, xstr)))[
            order(-frequency)]) 
}  

#### Column Combination generator ##################################################################
####################################################################################################
col_combs <- function(ncols){
    comb <- do.call(CJ, replicate(ncols, 0:1, FALSE))
    names <- as.numeric(gsub("V", "",colnames(comb)))
	comb2 <- data.table()
    cols <- colnames(comb)
        for(i in 1:length(cols)){
        comb1 <- comb[,i, with = FALSE]
        comb1[[cols[i]]][comb1 == 1] <- names[i]
        comb2 <-cbind(comb2, comb1)}
    comb2[, key := do.call(paste, c(.SD, sep = "_")), .SDcols = names(comb2)]
    k <- comb2$key
    return(k)}

#### Function that generates all combinations or models for a partucilar set of columns #########################
#### The data table must contain only columns one wants to run and the dependent variable needs to ##############
#### be in a separate column called "outcome" ###################################################################
#### If significant = T, only components of a model that are significant will be returned #######################
glmcompileR <- function(DT, key, significant){
  compiled_dt <- data.table()
  pb <- txtProgressBar(min = 0, max = length(key), style = 3)
  for(m in 2:length(key)){
    #### subset data table to contain only required columns
    spl <- strsplit(key[m], split = "_")[[1]]
    spl <- spl[!spl == "0"]
    newdt <- data.table()
    for(i in 1:length(spl)){
      temp <- DT[,c(as.numeric(spl[i])), with = FALSE]
      newdt <- cbind(newdt, temp)
      }
    #### add outcome variable to last column 
    newdt$outcome <- DT$outcome
    
    #### generate formula 
    nam <- colnames(newdt[, -(ncol(newdt)), with = FALSE])
    if(length(nam) > 1){
      plus <- "+"
      p1 <- paste(nam[-length(nam)], plus, sep = " ")
      p2 <- paste(p1, nam[length(nam)])
      fmla <- as.formula(paste0(colnames(newdt)[ncol(newdt)], " ~ ", p2))
    }else{
      fmla <- as.formula(paste0(colnames(newdt)[ncol(newdt)], " ~ ", nam))
      }
    
    #### fit glm model
    glm.fit <- glm(fmla, data=newdt, family = binomial)
    
    #### create data table 
    dtf <- data.table(coef(summary(glm.fit)))#[,"Pr(>|z|)"]
    dtf$names <- rownames(summary(glm.fit)$coefficients)
    chr <- as.character(glm.fit$formula)
    dtf$formula <- paste(chr[2], chr[1], chr[3:length(chr)])
    if(significant == "T"){
        #### subset data table to retain only significant indepenent variables 
        dtf <- dtf[dtf$`Pr(>|z|)` < 0.05 & (!dtf$names == "(Intercept)"),]
    }
    if(nrow(dtf) > 0){
      compiled_dt <- rbind(compiled_dt, dtf)
      }
    setTxtProgressBar(pb, m)
  }
  close(pb)
  compiled_dt <- compiled_dt[!duplicated(compiled_dt[,c(1:3)]),]
  compiled_dt <- compiled_dt[(!compiled_dt$names == "(Intercept)"),]
 
  compiled_dt$Significance <- "NA"
  compiled_dt[`Pr(>|z|)` < 0.05,]$Significance <- "Significant"
  compiled_dt[`Pr(>|z|)` > 0.05,]$Significance <- "Not-Significant"
  return(compiled_dt)
}

#### function that runs a logistic regression model on each independent variable in a data table ####
#### The key is a numeric string indicating the columns to use ######################################
#### e.g. key <- 2:length(colnames(dt_sub))-1 #######################################################
univariantglmR <- function(DT, key, significant){
  compiled_dt <- data.table()
  pb <- txtProgressBar(min = 0, max = length(key), style = 3)
  for(m in 2:length(key)){
    #### subset data table to contain only required columns
    newdt <- DT[,c(as.numeric(key[m])), with = FALSE]
    #### add outcome variable to last column 
    newdt$outcome <- DT$outcome
    
    #### generate formula 
    nam <- colnames(newdt[, -(ncol(newdt)), with = FALSE])
    if(length(nam) > 1){
      plus <- "+"
      p1 <- paste(nam[-length(nam)], plus, sep = " ")
      p2 <- paste(p1, nam[length(nam)])
      fmla <- as.formula(paste0(colnames(newdt)[ncol(newdt)], " ~ ", p2))
    }else{
      fmla <- as.formula(paste0(colnames(newdt)[ncol(newdt)], " ~ ", nam))
      }
    
    #### fit glm model
    glm.fit <- glm(fmla, data=newdt, family = binomial)
    
    #### create data table 
    dtf <- data.table(coef(summary(glm.fit)))#[,"Pr(>|z|)"]
    dtf$names <- rownames(summary(glm.fit)$coefficients)
    chr <- as.character(glm.fit$formula)
    dtf$formula <- paste(chr[2], chr[1], chr[3:length(chr)])
    if(significant == "T"){
        #### subset data table to retain only significant indepenent variables 
        dtf <- dtf[dtf$`Pr(>|z|)` < 0.05 & (!dtf$names == "(Intercept)"),]
    }
    if(nrow(dtf) > 0){
      compiled_dt <- rbind(compiled_dt, dtf)
      }
    setTxtProgressBar(pb, m)
  }
  close(pb)
  compiled_dt <- compiled_dt[!duplicated(compiled_dt[,c(1:3)]),]
  compiled_dt <- compiled_dt[(!compiled_dt$names == "(Intercept)"),]
  
  compiled_dt$Significance <- "NA"
  compiled_dt[`Pr(>|z|)` < 0.05,]$Significance <- "Significant"
  compiled_dt[`Pr(>|z|)` > 0.05,]$Significance <- "Not-Significant"
  return(compiled_dt)
}




#### function that runs a logistic regression model on an independent variable to determine what it predicts in a data table ####
#### The key is a numeric string indicating the columns to use ##################################################################
#### all dependent variables that are going to be tested must be a numeric column consisting of 0 or 1. #########################
#### the independent variable must be in a column called "test". ################################################################
#### e.g. key <- 2:length(colnames(dt_sub))-1 ##########################################
dependentglmR <- function(DT, key, independent_column, significant){
  compiled_dt <- data.table()
  pb <- txtProgressBar(min = 0, max = length(key), style = 3)
  for(m in 2:length(key)){
    #### subset data table to contain only required columns
    newdt <- DT[,c(as.numeric(key[m])), with = FALSE]
    #### add outcome variable to last column 
    newdt <- cbind(newdt, DT[,independent_column, with = FALSE])
    
    #### generate formula 
    nam <- colnames(newdt[, -(ncol(newdt)), with = FALSE])
    if(length(nam) > 1){
      plus <- "+"
      p1 <- paste(nam[-length(nam)], plus, sep = " ")
      p2 <- paste(p1, nam[length(nam)])
      fmla <- as.formula(paste0(p2, " ~ ", colnames(newdt)[ncol(newdt)]))
    }else{
      fmla <- as.formula(paste0(nam, " ~ ", colnames(newdt)[ncol(newdt)]))
    }
    
    #### fit glm model
    glm.fit <- glm(fmla, data=newdt, family = binomial)
    
    #### create data table 
    dtf <- data.table(coef(summary(glm.fit)))#[,"Pr(>|z|)"]
    dtf$names <- rownames(summary(glm.fit)$coefficients)
    chr <- as.character(glm.fit$formula)
    dtf$formula <- paste(chr[2], chr[1], chr[3:length(chr)])
    if(significant == "T"){
      #### subset data table to retain only significant indepenent variables 
      dtf <- dtf[dtf$`Pr(>|z|)` < 0.05 & (!dtf$names == "(Intercept)"),]
    }
    if(nrow(dtf) > 0){
      compiled_dt <- rbind(compiled_dt, dtf)
    }
    setTxtProgressBar(pb, m)
  }
  close(pb)
  compiled_dt <- compiled_dt[!duplicated(compiled_dt[,c(1:3)]),]
  compiled_dt <- compiled_dt[(!compiled_dt$names == "(Intercept)"),]
  
  compiled_dt$Significance <- "NA"
  compiled_dt[`Pr(>|z|)` < 0.05,]$Significance <- "Significant"
  compiled_dt[`Pr(>|z|)` > 0.05,]$Significance <- "Not-Significant"
  return(compiled_dt)
}











#### This function runs a normality test on selected columns ####
#### DT is the data table to analyze ############################
#### cols are the columns (in numeric fashion) to analyze #######
#################################################################
#cols <- c(6, 9, 10:12, 18, 20:25, 27)
normalityfindR <- function(DT, cols){
    finDT <- data.table()
    for(i in 1:length(cols)){

        test <- shapiro.test(as.numeric(DT[[cols[i]]]))

        if(test$p.value < 0.05){
          distribution <- "Not normally distributed"
        }else{
          distribution <- "normally distributed"
        }

    DT2 <- data.table(method = test$method,
                     pvalue = test$p.value,
                     variable = colnames(DT[,cols[i], with = FALSE]),
                     normality = distribution
                     )
    finDT <- rbind(finDT,DT2)
    }
    return(finDT)
}









