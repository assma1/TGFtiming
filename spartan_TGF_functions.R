#####################################################################
###   Targeted gene flow under a changing environment.            ###
###   Modified from generic TGF model by Kelly & Phillips 2018    ###  
###   Mods by Adam Smart                                          ###
###   For publication in Smart & Phillips 2020. "X"               ###
###                                                               ###
###   Updated 25.02.2020                                          ###
#####################################################################

### Other scripts required for full analysis: 'spartan_TGF_parameters.R; spartan_TGF_dataClean.R'

### State variables ###
# Vt = total variance of trait  
# Vg = genetic variace 
# Ve = env. variance
# h2 = heritability of trait
# eSize = effect size 
# nP = number of phenotypic loci (equal effects at each locus)
# nI = number of loci involved in D-M incompatibility 
# nN = number of neutral loci used to track prop. of recipeint genome intact post TGF
# nH = number of neutral loci used to track the heterozygosity of the loci post TGF
# f = frequency of '1' alleles at each locus
# lLoci =  random location of loci on each simulated 'chromosome'
# initMean = trait mean valaue at t = 0 

# Calculate input values for model

gSetup<- function(f, h2, VT, nP, nN, nI, nH){
  Ve<- (1-h2)*VT
  Vg<- h2*Ve/(1-h2)
  eSize<- sqrt(Vg/(2*nP*f*(1-f)))
  initMean<- f*2*nP*eSize
  lLoci<- runif(sum(nP, nN, nI, nH), 0, 1)
  initSD<- sqrt(Vg+Ve)
  list(eSize=eSize, Vg=Vg, Ve=Ve, initMean=initMean, lLoci=lLoci, initSD=initSD)
}

# Initialises the simulation with a array of individuals of dimensions [n, (nP+nI+nN), 2]
# recipient = setting nN to either 1, or 0 to track source vs. recipient genome
# Age is used to set the lifespan of agents in the population. Not used in base model.

initInds<-function(n , nP, nN, nI, nH, f, recipient=TRUE){
  Age<- rep(1, n)
  allelesP1 <- matrix(rbinom(n*nP, 1, f),
                      nrow=n, dimnames=list(NULL, rep("P", nP)))
  allelesP2 <- matrix(rbinom(n*nP, 1, f),
                      nrow=n, dimnames=list(NULL, rep("P", nP)))
  
  if (recipient) gCode <- 1 else gCode<- 20
  
  allelesH1 <- matrix(rbinom(n*nH, gCode, 0.5), 
                      nrow=n, dimnames=list(NULL, rep("H", nH)))
  allelesH2 <- matrix(rbinom(n*nH, gCode, 0.5), 
                      nrow=n, dimnames=list(NULL, rep("H", nH)))
  
  if (recipient) gCode<- 0 else gCode<- 1
  
  allelesI1 <- matrix(rep(gCode, n*nI), 
                      nrow=n, dimnames=list(NULL, rep("I", nI)))
  allelesI2 <- matrix(rep(gCode, n*nI), 
                      nrow=n, dimnames=list(NULL, rep("I", nI)))
  allelesN1 <- matrix(rep(gCode, n*nN), 
                      nrow=n, dimnames=list(NULL, rep("N", nN)))
  allelesN2 <- matrix(rep(gCode, n*nN), 
                      nrow=n, dimnames=list(NULL, rep("N", nN)))
  
  chrome1<- cbind(allelesP1, allelesI1, allelesN1, allelesH1, Age)
  chrome2<- cbind(allelesP2, allelesI2, allelesN2, allelesH2, Age)
  abind(chrome1, chrome2, along = 3)
}

# Calculates the phenotype for each row in the array, based on sum of nP alleles and effect size
# scaler is the initMean, used to center the initial train distribution around '0'

aPhen<- function(pop, eSize, Ve, scaler){
  expPhen<- apply(pop[,1:20,], 1, sum)*eSize - scaler 
  rnorm(length(expPhen), mean=expPhen, sd=sqrt(Ve))
}

# Create gametes for each individual in the population. This function simulates switching of alleles 
# based on a number of cuts to a simulated 'chromosome'.
# recombRate dictates how many cuts [nCut] (and resultant segment switches) of the simulated chromosome.

# C++ integration for defining nCuts
# Complite function for use in recombRate

cppFunction(
  "List cuts_cpp(NumericVector nCut) {
        int n = nCut.size();
        List ans(n);
        for (int i=0; i<n; i++){
          NumericVector vals = runif(nCut[i]);  
          std::sort(vals.begin(), vals.end());
          ans[i] = vals;
        }  
        return ans;
    }")

recom<- function(pop, lLoci, recombRate){
  nCut<- rpois(nrow(pop), recombRate)
  lCut<-cuts_cpp(nCut)
  if (length(lCut) == 0) return(pop) 
  damSet<- mapply(findInterval, lCut, MoreArgs = list(x = lLoci)) %% 2 == 0
  startSide<- rbinom(n = dim(pop)[1], size = 1, prob = 0.5) == 1
  damSet<- sweep(x = damSet, MARGIN = 2, STATS = startSide, FUN = "==")
  ifelse(t(damSet), pop[,,1], pop[,,2])
}

# Beverton Holt model of population growth 
# Rmax is maximum reproductive output
# Nstar is the carrying capacity of the system
# alpha [a] dictates the strength of the density dependance in the system

bevHolt<- function(N, Rmax, Nstar){
  a<-(Rmax-1)/(Nstar*1) 
  Rmax/(1+a*N)
}

# Defines a poisson curve to model the environmental shift in the model
# m dictates the width [slope] of the response aka. flattening term

optShiftLog<- function(initSD, m, gens){
  x<-seq(from = -gens, to = gens-1, by = 2)
  (initSD*2)/(1+exp(-m*(x-(x/2))))
}

# Defines stabilising selection, the magnitude of reduction in fitness for each increment away from the environmental optimum
# phen is drawn from generations mean phenotype
# opt is the environmental optimum at each timestep (calculated below)
# decay is the strength of the pentaly for each increment away from the environmental optimum

fitFun<- function(phen, opt, decay){
  exp(-decay*(phen-opt)^2)
}

# Reproduce individuals for each subsequent generation
# EW is the expected reproductive output for each individual
# W is reproductive output drawn from a poisson dist. with lambda = EW

reproduce<- function(pop, Rmax, Nstar, eSize, Ve, scaler, opt, decay, recombRate, lLoci, nP, nI, nN, nH){
  if(nrow(pop)<=2) return(pop)
  
  N<- dim(pop)[1]
  phen<- aPhen(pop, eSize, Ve, scaler)
  fitness<- fitFun(phen, opt, decay)
  EW<- bevHolt(N, Rmax, Nstar)*fitness
  W<- rpois(N, lambda = EW)
  
  offDam<- recom(pop[rep(1:N, W),,, drop = FALSE], lLoci, recombRate)
  offDam<- cbind(offDam, rep(0, nrow(offDam)))
  if (sum(offDam)<=1) return(pop)
  colnames(offDam) <- rep(c("P", "I", "N", "H", "Age"), times = c(nP, nI, nN, nH, 1))
  
  offSire<- recom(pop[sample(1:N, sum(W), replace = TRUE), , , drop=FALSE], lLoci, recombRate)
  offSire<- cbind(offSire, rep(0, nrow(offSire)))
  colnames(offSire) <- rep(c("P", "I", "N", "H", "Age"), times = c(nP, nI, nN, nH, 1))
  
  if (length(offDam) != length(offSire)) return (pop)
  
  abind(offDam, offSire, along = 3)
}


# Survival based on DM incompatibilty.
# outbreeding = is.logical = include outbreeding depression

survival<- function(pop, eSize, Ve, outbreeding, nLoci, h2, alpha, H2){
  pSurv<- (pop[,"Age",]==0) * 1 + (pop[,"Age",]==1) * 0
  if (outbreeding == TRUE) pSurv<- pSurv*Vmod(pop, H2, alpha, nLoci) else pSurv<- pSurv
  surv<- as.logical(rbinom(nrow(pop), size = 1, prob = pSurv)) #survival based on pSurv
  pop<- pop[surv, , ,drop=FALSE]
  pop[,"Age",]<- pop[,"Age",] + 1
  pop
}

# Function to score proportion of genome that is from recipient population
# r.current = proportion of genome from recipient across the population
# r.new = proportion of loci with recipient alleles

neutral.gTracker<-function(pop, nLoci){
  if (nrow(pop)<=2) return (pop)
  gCols<- grepl("N", colnames(pop)) #gentoype columns "N"
  gCols<- pop[, gCols, , drop = FALSE]
  r.current<- 1-sum(gCols)/length(gCols) # proportion of genome from recipeint across population #list
  r.current
  }

# Tracks various ratios and gene frequencies throughout the simulation. 
# Including recording outputs of homozygosity scores for Turelli and Orr model of incompatibility 
# P1 = homozygous at damAllele '1'
# P2 = homozygous at sireAllee '2'
# PH = hetrozygous
# freqs = frequency of '1' allele at each locus across the population

gTracker<- function(pop, gVars, nLoci){
  freqs<- apply(pop, 2, sum)/(dim(pop)[1]*2) #frequencies at each locus
  
  #Simpsons Diversity Index
  gCols_h<- grepl("H", colnames(pop))
  gCols_h<- pop[, gCols_h, , drop = FALSE] #
  SDI_freqs<- apply(gCols_h, 2, FUN=function(x){sum(x[x==1])/2})
  SDI_sum<- sum(SDI_freqs)
  temp<- SDI_freqs-1
  SDI_freqs<- SDI_freqs*temp
  SDI_freqs<- sum(SDI_freqs)
  SDI<- 1-(SDI_freqs)/(SDI_sum*(SDI_sum-1)) # 1 - SDI to make 1 equal to high diveristy
  
  #scores for 'I' incompatibility allels
  gCols<- grepl("I", colnames(pop))
  nCols<- sum(gCols) #number of alleles
  nLoci<- nCols
  gCols<- pop[, gCols, , drop = FALSE] #genotypes
  damAllele<- gCols[, , 1, drop = FALSE] #relevant subset
  sireAllele<- gCols[, , 2, drop = FALSE]
  homozygous<- damAllele==sireAllele #conditions
  homP1<- homozygous & damAllele == 0
  homP2<- homozygous & damAllele == 1
  P1<- rowSums(homP1)/nLoci
  P2<- rowSums(homP2)/nLoci
  PH<- 1-(P1+P2)
  list(P1 = P1, P2 = P2, PH = PH, freqs=freqs,  SDI=SDI)
}

# Provides a breakdown score for incompatibility bewteen hetrozygous loci [H0], heterozygous 
# and a homozygous loci [H1] and homozygous loci [H2] respectivly.

hWeight<-function(H2){
  H1<- 0.5*H2
  H0<- 0.25*H2
  list(H2, H1, H0)
}

# Turelli and Orr model of  allele incompatibility. Computes a 'hybrid breakdown score' for each individual
# based on thier composition of nI alleles.
#ExpS is the probilibty of survival from outbreeding depression
#nLoci is number of alleles involved in DM incomptibility (nLoci = nI)

TOmod<- function(pop, nLoci, H2){
  pList<- gTracker(pop, gVars, nLoci)[1:3]
  hList<- hWeight(H2)
  ExpS<- nLoci*(pList[[1]]*pList[[2]]*hList[[1]]+(pList[[1]]+pList[[2]])*pList[[3]]*hList[[1]]+pList[[3]]^2*hList[[3]])
  ExpS
}

# Simple exponential function to link hybrid breakdown scores to fitness.
# alpha is some positive constant, can vary to manipulate strength of outbreeding depression

Vmod<- function(pop, H2, alpha, nLoci){
  exp(-alpha*TOmod(pop, nLoci, H2))
}

# Top level function to set up simulation with targeted gene flow and a variable shift in the environment
# TGFtime = time to introduce individuals
# burnIn = simulation burn in time 

evolve<- function(initN, nP, nN, nI, nH, initFreq, h2, VT, gens, Rmax, Nstar, decay, recombRate, m, nIntro, TGFtime, nLoci, H2, outbreeding, alpha, burnIn = 10){
  
  #initilise the population
  pop<- initInds(initN, nP, nN, nI, nH, f=initFreq)
  
  #global parameters
  gVars<- gSetup(initFreq, h2, VT, nP, nN, nI, nH)
  opt<<- optShiftLog(initSD=gVars$initSD, m, gens)
  
  #find 'f' that corresponds to opt @ gg=50 #0.1=Vg
  fint<<- (opt[1:50]^2)/((2*nP*0.1)+opt[1:50]^2)
  
  # bins for collecting stuff
  introductees<-c()
  genome<-c()
  IntProp<-c()
  pheno<-c()
  simpsonsDI<-c()
  popSize<-c()
  
  # function for implementing TGF - placing new individuals in recipeint population
  inject<- function(nIntro, nP, nN, nI, nH, f, pop, recipient = FALSE){
    if (nIntro == 0) return(pop)
    nIntro<- ceiling(nrow(pop)*nIntro)
    injectSample<- initInds(n=nIntro, nP, nN, nI, nH, f, recipient = recipient)
    pop<- abind(pop, injectSample, along = 1)
    return(pop)
  }
  
  #function for iterating head function over generations, given a start and finish time
  iterate<-function(start = 0, finish = gens, collect = TRUE, outbreeding){
    for (gg in start:finish){
     
      if (nrow(pop)<=2) {break}
      if (gg == TGFtime & nrow(pop) >1){pop<<- inject(nIntro, nP, nN, nI, nH, f=fint[gens], pop)} # translocate TGF inds. into recip. pop.
      if (gg == TGFtime & nrow(pop) >1){IntProp<<- c(IntProp, nIntro/nrow(pop))} # record propotions of TGF inds. to recip. inds.
      
      pop<<- reproduce(pop, Rmax, Nstar, eSize=gVars$eSize, Ve=gVars$Ve, scaler=gVars$initMean, opt=opt[gg], decay, recombRate, lLoci=gVars$lLoci, nP, nI, nN, nH)
      
      if (nrow(pop)<=2) {break}
      
      pop<<- survival(pop, gVars$eSize, Ve, outbreeding, nLoci, h2, alpha, H2)
      
      if (collect){

        popSize[gg]<<-nrow(pop)
        temp<<- gTracker(pop, gVars, nLoci) # track gene frequencies 
        genome<<- neutral.gTracker(pop, nLoci)
        pheno[gg]<<- mean(aPhen(pop, eSize=gVars$eSize, Ve=gVars$Ve, scaler=gVars$initMean)) # mean phenotype
        simpsonsDI[gg]<<- gTracker(pop, gVars, nLoci)$SDI # average heterozygosity across H alleles each generation
        introductees[gg]<<- ceiling(nrow(pop)*nIntro)
      }
    }
  }
  
  # iterate over generations
  iterate(1, burnIn, collect = TRUE, outbreeding) # burn in stage
  iterate(burnIn+1, gens, collect = TRUE, outbreeding) # focal stage
  
  # collect and store the following in [output] >> passed to 'spartan_TGF_paramters.R'
  list(popSize=popSize, genome=genome, introductees=introductees, pheno=pheno, simpsonsDI=simpsonsDI)
}
