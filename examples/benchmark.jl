# Benchmarking

using MonotoneDecomposition

# Candidate functions
fs = ["x^2" "x^3" "exp(x)" "sigmoid" "SE_1" "SE_0.1" "Mat12_1" "Mat12_0.1" "Mat32_1" "Mat32_0.1" "RQ_0.1_0.5" "Periodic_0.1_4"]

benchmarking.(
    fs[1:2];
    σs = [0.1, 0.2, 0.4, 0.5, 1.0, 1.5, 2.0], # noise level to be surveyed
    jplot = false, # μerr vs σs
    nrep = 1, # NB: for fast auto-generation procedure, only use nrep = 1; in the paper, use nrep = 100
    competitor = "ss_single_lambda", 
    nfold = 10, # number of folds
    one_se_rule=false, 
    nλ = 2, # the number of λ to be searched
    rλ = 0.1, # the search region of λ, (1-rλ, 1+rλ)*λ
)



