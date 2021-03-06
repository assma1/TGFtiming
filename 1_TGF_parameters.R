#####################################################################
###   Targeted gene flow under a changing environment.            ###
###   Set up parameter space and run simulations                  ###
###   For publication in Smart & Phillips 2020. "X"               ###
###                                                               ###
###   Updated 25.02.2020                                          ###
#####################################################################

### Other scripts required for full analysis: '2_TGF_functions.R; 3_TGF_figures.R'

###########################################################################
### STEP 1:  Command line arguments for SPARTAN runs                    ###
###########################################################################

# Set seed
set.seed(1000)

# Create comnand line for job_submission.slurm script
command_args<- commandArgs(trailingOnly = TRUE)
parameter_index<- as.numeric(command_args[1])

# Check for package dependices and install/load in packages
.libPaths("/home/asmart1/R/lib/3.6")
lib<- .libPaths()[1]
repo<- "https://cran.ms.unimelb.edu.au/"
pkgs = c("abind", "Rcpp")
new.pkgs <- pkgs[!(pkgs %in% installed.packages(lib=lib)[,"Package"])]
if(length(new.pkgs)) install.packages(new.pkgs, lib=lib, repos=repo)
inst = lapply(pkgs, library, character.only = TRUE)


###########################################################################
### STEP 2:  Load in functions from sibling script                      ###
###########################################################################

source("src/2_TGF_functions.R")

##########################################################################
### STEP 3:  Set input paramters for simulation runs                   ###
##########################################################################

## Setup for outbreeding scenarios
H2<-1
h1<-H2*0.5
h0<-H2*0.25

# Set up loci values (all hetrozygotes)
p1<-0 # e.g. (1,0)
p2<-0 # e.g. (0,1)
pH<-1 # e.g. (1,0) or (0,1)

# Calclate ExpS
ExpS<- 20*(p1*p2*H2+(p1+p2)*pH*h1+pH*h0) # Numeric for # loci
alpha50<- (-log(0.5, base=exp(1))/ExpS) 
alpha10<- (-log(0.9, base=exp(1))/ExpS)

## Possible ID options for each parameter in managment space.
initFreq_options = c(0.05, 0.3) #initial frequency of alleles
recomb_options = c(10, 50) 
h2_options = c(0.1, 0.3)
Rmax_options = c(2, 5)
Nstar_options = c(500, 2000)
m_options = c(0.15, 1)
decay_options = 1
alpha_options = c(0, alpha10, alpha50)
change_type_options = c(0, 1)

# Create grid of all possible parameter combinations
testGrid<-expand.grid(initFreq_options, # intitail frequency of alleles
                      recomb_options,   # chromosome recombination rate
                      h2_options,       # trait heritability 
                      Rmax_options,     # maximum reproductive output
                      Nstar_options,    # carrying capcity 
                      m_options,        # paramter controlling env. shape
                      decay_options,    # penetaly of trait mismatch
                      alpha_options,
                      change_type_options)    # strength of outbreeding depression

initFreq<-testGrid[parameter_index, 1] 
recomb<- testGrid[parameter_index, 2]
h2<-testGrid[parameter_index, 3]  
Rmax <- testGrid[parameter_index, 4]    
Nstar <- testGrid[parameter_index, 5]  
m <- testGrid[parameter_index, 6]       
decay<- testGrid[parameter_index, 7]
alpha<- testGrid[parameter_index, 8]
change<- testGrid[parameter_index, 9]

##########################################################################
### STEP 5:  Run TGF simulations                                       ###
##########################################################################

# Simulation head function
manageSpace<- function(initFreq, recomb, h2, Rmax, Nstar, m, decay, alpha, change, props=seq(from=0, to=0.5, by=0.025), yr=seq(from=0, to=50, by=2), simulations = 100){

  # Output matrix
  TGFstack<- matrix(ncol=6, nrow= 0)

  for (k in yr){
    
    # Collection bins
    extinctRisk<- vector()
    r.pop<- vector() 
    simpsonsDI<- vector()
    meanPheno<- vector()
    
    # Diagnostic outputs
    for (j in props){
      cat("\n Calculating extinction probability for generic model at:",
          "\n initFreq = ", initFreq,
          "\n recombRate =", recomb,
          "\n h2 = ", h2,
          "\n Rmax = ", Rmax,
          "\n Nstar = ", Nstar,
          "\n m = ", m,
          "\n decay = ", decay,
          "\n alpha = ", alpha, 
          "\n type = ", change,
          "\n Percentage introduced = ", j,
          "\n Year of Introduction =", k)
      
      # Collection bins
      extinct<- vector() 
      genomePop<- vector() 
      heterozygosity<- vector()
      pheno<- vector()

      # Run simulation
      for (i in 1:simulations){

        output<- evolve(initN = Nstar, 
                        nP = 20, # pheno alleles
                        nN = 20, # neutral alleles
                        nI = 20, # incompatibility alleles
                        nH = 20, # heterozygosity alleles
                        initFreq = initFreq, #frequency of "1" alleles initially
                        h2 = h2, #heritability 
                        Ve = 1, # total variation in system
                        Rmax = Rmax, # max reproductive output
                        Nstar = Nstar, # carrying capacity 
                        decay = decay, # rate of environmental decay
                        recombRate = recomb, # rate of recombination
                        gens = 50, # number of generations
                        m = m, # width of environemtnal optimum
                        nIntro = j, # number of TGF individuals
                        TGFtime = k, # time to introduce TGF individuals
                        nLoci = 20, # number of loci 
                        H2 = 1, # outbreeding value
                        outbreeding = TRUE, # outbreeding logical
                        alpha = alpha, # reduction in popsize due to outbreeding depression
                        type = change
        )
        extinct[i]<- ifelse(length(output$popSize) < 50, TRUE, FALSE) # if extinct TRUE
        genomePop[i]<- ifelse(is.null(output$genome), NA, tail(output$genome, 1))  # record % recipient genome per ind.
        heterozygosity[i]<- ifelse(length(output$simpsonsDI) < 50, 0, tail(output$simpsonsDI, 1))# record the proportion of introduced inds. to recipient inds.
        pheno[i]<- ifelse(is.null(output$pheno), NA, tail(output$pheno, 1))
      }
      extinctRisk<- c(extinctRisk, length(extinct[extinct==TRUE])/simulations) # proportion of extinctions over sim. runs
      r.pop<- c(r.pop, mean(as.numeric(genomePop), na.rm=T)) # record proportion of recipient genome
      simpsonsDI<- c(simpsonsDI, mean(as.numeric(heterozygosity), na.rm=T)) # recond Gini-Simpson indec
      meanPheno<- c(meanPheno, mean(as.numeric(pheno), na.rm=T)) # mean phenotype 
    }
    TGFstack<- rbind(TGFstack, cbind(k, extinctRisk, r.pop, props, simpsonsDI, meanPheno)) # create output array >> used in 'spartan_TGF_dataClean.R'
  }
  saveRDS(TGFstack, sprintf("out/output_%04d.rds", parameter_index))
  }


# Run sensitivity analysis across entire testGrid space [dim: 972  8]
manage<-manageSpace(initFreq, recomb, h2, Rmax, Nstar, m, decay, alpha, change)

