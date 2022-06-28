library(ggplot2)
library(fdrame)
library(samr)
library(Rcpp)
library(BH)
theme_set(theme_classic())
Sys.setenv("PKG_CXXFLAGS"="-std=c++17")
sourceCpp("PermFDRAdjust.cpp")

#### FUNCTIONS USED IN THE SIMPLYCORRECT APP.

# This function simulates a quantitative omics experiment in which no genes/or proteins are changing in truth. It returns the resulting data, colors for volcano plotting, and the threshold used.
simulateData = function(seed = 1,
                        setSeed = T,
                        m = 400,
                        avg = 22,
                        sd = 3,
                        nc = 3,
                        nt = 3,
                        expVar = 0.01,
                        threshold = 0.05,
                        adjMeth = "None",
                        nperms = 100,
                        model = "NoChange",
                        effectMagnitude = 0.5,
                        effectSD = 1,
                        bioVar = 0.05, 
                        bioVGamma = c(3, 0.5),
                        s0 = 0,
                        s0.perc = -1,
                        permMeth = "SRS") {
  if (setSeed) {
    set.seed(seed)
  }
  
  #Generate proteins
  data = data.frame(rnorm(m, avg, sd)) # Make control group True Quantites
  names(data)[1] = "TrueQuantC" # Rename control group True Quantity column
  data$Protein = "Protein" # Make Protein column with "protein names"
  data$Protein = paste(data$Protein, row.names(data)) # Index protein names
  
  ## Generate intensities across samples
  
  #### MODEL 1: NO CHANGE ####
  if (model == "NoChange") {
    for (each in 3:(2+nc+nt)) {
      data[,each]= rnorm(nrow(data), data$TrueQuantC, data$TrueQuantC * expVar)
      colnames(data)[each] = paste("Int",each-2,sep="")
    }
    
    #### MODEL 2: ALL CHANGE ####
  } else if (model == "AllChange") {
    #We will randomly change some proteins down and some proteins up
    sign = sample(c(-1,1),nrow(data), replace = T)
    #make replicates within the control group
    for (each in 3:(2+nc)) {
      data[,each]= rnorm(nrow(data), data$TrueQuantC, data$TrueQuantC * expVar)
      colnames(data)[each] = paste("Int",each-2,sep="")
    }
    #make replicates in the test group
    for (each in (nc+3):(2+nc+nt)) {
      #Add or subtract effectMagnitude to this and generate random values.
      data[,each]= rnorm(nrow(data), data$TrueQuantC + sign * effectMagnitude, (data$TrueQuantC + sign * effectMagnitude) * expVar)
      colnames(data)[each] = paste("Int",each-2,sep="")
    }
    
    #### MODEL 3: REALISTIC ####
  } else if (model == "Realistic") {
    #For each protein, generate a true effect from the normal distribution.
    realEffects = rnorm(nrow(data), 0, effectSD)
    
    #Each protein will have its own biological variation as a percentage of the protein intensity. Sample this from a gamma dist. 
    #The gamma dist with shape = k and scale = theta has a mean of k*theta. 
    #The user specifies the mean = (scaling factor) * k * theta. 
    #We produce a scaled gamma dist with (scaling factor) = mean/(k * theta).
    bioVars = rgamma(nrow(data), 
                     shape = bioVGamma[1], 
                     scale = bioVGamma[2]) * bioVar / (bioVGamma[1] * bioVGamma[2])
    
    #make replicates within the control group
    for (each in 3:(2+nc)) {
      #The variation will be the sum of the biol. and exptl. variation
      data[,each]= rnorm(nrow(data), data$TrueQuantC, 
                         (data$TrueQuantC * bioVars + data$TrueQuantC * expVar))
      colnames(data)[each] = paste("Int",each-2,sep="")
    }
    
    #In biology, some may consider a protein "unchanged" if the biological variability is greater than the effect (underlying shift in mean). We'll make these log2fc's actually 0, so that the null is true.
    
    # if (real Effect < biological variability)
    nullIsTrue = (abs(realEffects) < (bioVars * data$TrueQuantC))
    # then set the real effect to 0 for those proteins
    realEffects[nullIsTrue] = 0
    
    #make replicates in the treatment group
    for (each in (2+nc+1):(2+nc+nt)) {
      data[,each]= rnorm(nrow(data), data$TrueQuantC + realEffects, 
                         (data$TrueQuantC * bioVars + data$TrueQuantC * expVar))
      colnames(data)[each] = paste("Int",each-2,sep="")
    }
    
    #save expVar and bioVars
    data$bioVar = bioVars
    data$expVar = expVar
    
    #Label changing or not
    data$nullIsTrue = nullIsTrue
    
  } else {
    stop("Invalid model name passed to simulateData function. Must be 'NoChange', 'AllChange', or 'Realistic' ")
  }
  
  #Compare group to group
  beginC = 3
  endC = 2 + nc
  beginT = nc + 3
  endT = 2 + nc + nt
  data$logFC =  rowMeans(data[,beginT:endT]) - rowMeans(data[,beginC:endC])
  data$P.Val = NA
  for (each in 1:nrow(data)) {
    test = t.test(data[each,beginC:endC], data[each,beginT:endT])
    data$P.Val[each] = test$p.value
  }
  
#### Adjust p-values / rejection threshold ####
  if (adjMeth == "BH") {threshold = BH.adjust(data, threshold)}
  if (adjMeth == "Bonf") {threshold = threshold / nrow(data)}
  #For permutation...
  if (adjMeth == "Perm") {
    #Generate design and expression matrix
    myDesign = c(rep(1, nc), rep(2, nt))
    intOnly = as.matrix(data[,grep("Int", names(data))])
    #Perform perm FDR
    set.seed(seed)
    nPerms = nperms
    if (permMeth == "samr") {
      listForSamr = list(x = intOnly, y = myDesign, genenames = data$Protein, geneid = data$nullIsTrue, logged2 = T)
      permList = samr(data = listForSamr, resp.type = "Two class unpaired",
                      assay.type = "array", s0 = s0, s0.perc = s0.perc,
                      nperms = nperms, testStatistic = "standard", 
                      random.seed = seed)
      FDRTable = samr.compute.delta.table(permList)
      adjList = samr.compute.siggenes.table(permList, del = NULL, data = listForSamr, delta.table = FDRTable, all.genes = T)
      adjData = as.data.frame(rbind(adjList$genes.up, adjList$genes.lo) )
      
      # Rename samr output columns to match what we had before
      data = adjData[,c("Gene ID", "Gene Name", "Fold Change", "q-value(%)")]
      names(data) = c("Protein", "nullIsTrue", "logFC", "P.Val")
      # Convert columns to proper types
      data$Protein = as.character(data$Protein)
      data$nullIsTrue = as.logical(data$nullIsTrue)
      data$logFC = log2(as.numeric(as.character(data$logFC)))
      data$P.Val = as.numeric(as.character(data$P.Val)) / 100
      
      #If the adj p-values are NA, give these p = 1
      data$P.Val[is.na(data$P.Val)] = 1
      #If the adj p-values are 0, give them the next lowest or threshold / 2, whichever is lower
      if (min(data$P.Val[data$P.Val != 0]) > threshold / 2) {
        data$P.Val[data$P.Val == 0] = threshold / 2
      } else {
        data$P.Val[data$P.Val == 0] = min(data$P.Val[data$P.Val != 0])
      }
    }
    else {
      #threshold = perm.FDR.adjust(data, threshold, myDesign, intOnly, nperms, nc, nt)
      threshold = permFDRAdjustRcpp(data, threshold, myDesign, intOnly, nperms, nc, nt)
    }
  }
  
  #Label false positives and true negatives
  if (model == "NoChange") {
    data$ResultType = factor(ifelse(data$P.Val < threshold, "FalsePos", "TrueNeg"),
                             levels = c("FalsePos", "TruePos", "FalseNeg", "TrueNeg"))
  } else if (model == "AllChange") {
    data$ResultType = factor(ifelse(data$P.Val < threshold, "TruePos", "FalseNeg"),
                             levels = c("FalsePos", "TruePos", "FalseNeg", "TrueNeg"))
  } else if (model == "Realistic") {
    data$ResultType = factor(
      ifelse(
        data$nullIsTrue,
        ifelse(data$P.Val < threshold, "FalsePos", "TrueNeg"),
        ifelse(data$P.Val < threshold, "TruePos", "FalseNeg")
      ),
      levels = c("FalsePos", "TruePos", "FalseNeg", "TrueNeg")
    )
  }
  
  ####Color coding ####
  #Choose colors
  colorList = list(
    FalsePoscolor = col2rgb("red")[,1],
    TruePoscolor = col2rgb("blue")[,1],
    FalseNegcolor = col2rgb("cyan")[,1],
    TrueNegcolor = col2rgb("pink")[,1]
  )
  #Convert colors from name to hex
  myColors = c()
  for (each in colorList  ) {
    each = unlist(each)
    each = each / 255
    each = rgb(each[1], each[2], each[3])
    myColors = c(myColors, each)
  }
  #Tweak the colors
  myColors[1] = "#D00030"
  myColors[2] = "#0000BB"
  myColors[3] = "#0090FF"
  myColors[4] = "#FFB08B"
  
  #### Generate results ####
  #How many of each?
  results = data.frame( Metric = c(
    "Total Hits",
    "Total False Positives",
    "Total True Positives",
    "Total False Negatives",
    "Total True Negatives"
  ),
  ValueNum = c(
    sum(data$ResultType == "FalsePos" | data$ResultType == "TruePos"),
    sum(data$ResultType == "FalsePos"),
    sum(data$ResultType == "TruePos"),
    sum(data$ResultType == "FalseNeg"),
    sum(data$ResultType == "TrueNeg")
  )
  )
  
  #What percent of the dataset is each? 
  results = rbind(results, data.frame(Metric = c(
    "Hits (% of Data)",
    "False Positives (% of Data)",
    "True Positives (% of Data)",
    "False Negatives (% of Data)",
    "True Negatives (% of Data)"
  ),
  ValueNum = c(
    results[1,"ValueNum"] / nrow(data) * 100,
    results[2,"ValueNum"] / nrow(data) * 100,
    results[3,"ValueNum"] / nrow(data) * 100,
    results[4,"ValueNum"] / nrow(data) * 100,
    results[5,"ValueNum"] / nrow(data) * 100
  )
  ))
  
  ### FDP, FNP, Sensitivity, Specificity
  ### First avoid dividing by 0
  ### Then calculate the fraction
  
  # FDP
  if (is.na(results[2,"ValueNum"]) | results[2,"ValueNum"] == 0) {
    fdp = 0
  } else {
    fdp = results[2,"ValueNum"] / (results[2,"ValueNum"] + results[3,"ValueNum"]) * 100
  }
  
  # FNP
  if (is.na(results[4,"ValueNum"]) | results[4,"ValueNum"] == 0) {
    fnp = 0
  } else {
    fnp = results[4,"ValueNum"] / (results[4,"ValueNum"] + results[5,"ValueNum"]) * 100
  }
  
  # Sensitivity
  if (results[3,"ValueNum"] + results[4,"ValueNum"] == 0) {
    sens = NA
  } else {
    sens = results[3,"ValueNum"] / (results[3,"ValueNum"] + results[4,"ValueNum"]) * 100
  }
  
  # Specificity
  if (results[2,"ValueNum"] + results[5,"ValueNum"] == 0) {
    spec = NA
  } else {
    spec = results[5,"ValueNum"] / (results[5,"ValueNum"] + results[2,"ValueNum"]) * 100
  }
  
  #Save results
  results = rbind(results, data.frame(Metric = c(
    "False Discovery Proportion (FDP)",
    "False Negative Proportion (FNP)",
    "Observed Sensitivity",
    "Observed Specificity"
  ),
  ValueNum = c(fdp, fnp, sens, spec)
  ))
  
  #Convert results to character strings
  results$Value = as.character(round(results$ValueNum, 1))
  #Add "%"s where appropriate
  results$Value[grep("%|F[DN]R|Sens|Spec", results$Metric)] =
    paste0(results$Value[grep("%|F[DN]R|Sens|Spec", results$Metric)], "%")
  
  #Reorder results and cut out ValueNum
  finalResults = rbind(results[11:14,c(1,3)],
                       results[1:10,c(1,3)])
  
  return(list(data, myColors, finalResults, threshold))
  
}

# This function draws a volcano plot showing the false & true positives and negatives. 
volcanoPlot = function(data, colors, threshold, title = "No Correction" , adj = F) {
  p = ggplot(data,aes(x=logFC, y=-log10(P.Val))) + 
    geom_hline(yintercept=-log10(threshold), colour='black', linetype='solid') +
    geom_point(aes(color = ResultType), 
               size = 3,
               alpha = 0.7) +
    discrete_scale("color", "myScale", 
                   palette = colorRampPalette(colors), drop = F) +
    ggtitle(title) +
    xlab("Observed Effect") +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 13),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 13)) +
    scale_y_continuous(limits=c(0, max(c(-log10(data$P.Val), -log10(threshold)))))
  
  if (adj)
    p = p + ylab("-log10(Q.Val)")
  
  p
}

# This function controls FDR using the manual method described in Benjamini & Hochberg, 1995. It doesn't "adjust" p-values per se, but raises the threshold for rejection. This results in nicer-looking plots than the stats::p.adjust(), which uses a cumulative minimum function resulting in ugly plots, but that method doesn't require a pre-specified rejection threshold. 
# The two methods are equivalent in terms of what hypotheses end up accepted and rejection, and therefore equivalent in terms of FDR control.
# It returns the new threshold for P-value rejection.
BH.adjust = function(data, threshold) {
  
  #Rank p-values from lowest to highest, rank them
  data = data[order(data$P.Val),]
  data$rank = 1:nrow(data)
  
  #Multiply the p value by the number of analytes m, then divide by the rank i
  data$adjP = data$P.Val * nrow(data) / data$rank
  
  #Find k which is the largest rank i for which (adjP <= threshold) is true
  data$ineq = data$adjP <= threshold
  #If there are none...
  if (length(which(data$ineq)) == 0) {
    #... then k = 0
    k = 0
  } else { k = max(which(data$ineq)) }
  
  #Now go back to your P values and reject all H_i where i <=k
  #If you don't do this, you will accept some H_i you don't want to accept. Notice how adjP is non-monotonic.
  data$reject = (data$rank <= k)
  
  #The new threshold is after the last rejected H_i.
  #If none are rejected, we consider that if the highest-ranking p-value had satisfied the inequality (adjP = threshold), that would have been the last and only rejected H_i.
  #Solve for that P and set it as the new rejection threshold.
  if (sum(data$reject) == 0) {
    #Notice that this is equivalent to Bonferroni!
    newThreshold = threshold / nrow(data)
    #If all are rejected, then our previous threshold was already past the last rejected H_i, so keep it
  } else if (sum(!data$reject) == 0) {
    newThreshold = threshold
  } else {
    #Otherwise, use the midpoint of the log10 to draw the line on the volcano plot nicely.
    lastRejectedP = data$P.Val[max(which(data$reject))]
    firstAcceptedP = data$P.Val[max(which(data$reject)) + 1]
    newThreshold = 10 ** (mean(c(log10(lastRejectedP),
                                 log10(firstAcceptedP))))
  }
  
  return(newThreshold)
  
}

# This function controls FDR using the permutation method described in our manuscript. Like the BH method above, it corrects the rejection threshold rather than the p-values themselves.  
# It returns the new threshold for P-value rejection.
perm.FDR.adjust = function(data, threshold, myDesign, intOnly, nPerms, nc, nt) {
  
  # Generate [nPerms] random permutations of the samples.
  # For each permutation, perform t-tests and 
  # save the p-values.
  randPValList = list()
  for (i_perm in 1:nPerms) {
    # Generate random design
    randomDesign = randBalDesign(nc, nt)
    
    # Do t test on each protein and log the p value
    randPVals = c()
    for (i_prot in 1:nrow(intOnly)) {
      ttest = t.test(intOnly[i_prot, randomDesign == 1],
                     intOnly[i_prot, randomDesign == 2])
      randPVals = c(randPVals, ttest$p.value)
    }
    
    # Add p-values to list
    randPValList[[i_perm]] = randPVals
  }
  
  # At each p-value measured in our real experiment, calculate the estimated FDP = V / R where R is the number of real proteins with this p value or less (the rank).
  data = data[order(data$P.Val),]
  data$rank = 1:nrow(data)
  data$estFDP = NA
  
  for (rank in data$rank) {
    eachThresh = data[rank, "P.Val"]
    hitCounts = c()
    
    # Calculate the average # hits at this threshold
    for (pVec in randPValList) {
      hitCount = length(which(pVec <= eachThresh))
      hitCounts = c(hitCounts, hitCount)
    }
    v = mean(hitCounts)
    
    # Save the estimated FDP
    fdp = v / rank
    data$estFDP[rank] = fdp
  }
  
  # Return the highest p-value threshold for which est FDP <= threshold
  validThreshInds = which(data$estFDP <= threshold)
  if (length(validThreshInds) < 1) {
    return(data[1,"P.Val"] / 2)
  } else if (length(validThreshInds) == length(data$estFDP)) {
    if (data$P.Val[length(data$P.Val)] + 0.05 > 1)
      return((data$P.Val[length(data$P.Val)] + 1) / 2)
    else
      return(data$P.Val[length(data$P.Val)] + 0.05)
  } else {
    return(mean(c(data$P.Val[max(validThreshInds)],
                  data$P.Val[max(validThreshInds) + 1])))
  }
}

# This function controls FDR using the permutation method described in our manuscript. Like the BH method above, it corrects the rejection threshold rather than the p-values themselves.  
# It returns the new threshold for P-value rejection.
# It uses the Rcpp and BH packages to leverage fast C++ code.
permFDRAdjustRcpp = function(data, threshold, myDesign, intOnly, nPerms, nc, nt) {
  pVals = data$P.Val[order(data$P.Val)]
  intMatrix = as.matrix(intOnly)
  return(permFDRAdjustCpp(pVals, threshold, myDesign, intMatrix, nPerms, nc, nt))
}


# This function takes N's and returns a design vector of 1s and 2s that is randomized within the constraint of being maximally balanced.
randBalDesign = function(nc, nt) {
  # First find the maximally balanced allocation.
  possibleCinC = nc:(nc-nt)
  possibleCinT = nc - possibleCinC
  balance = abs(possibleCinC / nc - possibleCinT / nt)
  i_alloc = which(balance==min(balance))[1]
  
  # Next randomly allocate samples into each group at the maximally balanced proportions
  CinC = possibleCinC[i_alloc]
  CinT = possibleCinT[i_alloc]
  
  controlGroup = c(rep(1,CinC), rep(2,nc-CinC))
  controlGroup = controlGroup[sample(1:nc)]
  
  testGroup = c(rep(1,CinT), rep(2,nt-CinT))
  testGroup = testGroup[sample(1:nt)]
  
  return(c(controlGroup, testGroup))
}
