# Benchmarking

using MonotoneDecomposition

# Candidate functions
fs = ["x^2" "x^3" "exp(x)" "sigmoid" "SE_1" "SE_0.1" "Mat12_1" "Mat12_0.1" "Mat32_1" "Mat32_0.1" "RQ_0.1_0.5" "Periodic_0.1_4"]

nrep = 1
nλ = 2
nfold = 2
competitor = "ss_single_lambda"
nλ = ifelse(occursin("single_lambda", competitor), 1, nλ)
one_se_rule = false
resfolder0 = "/tmp"
if length(ARGS) > 0
    @info "Use args passed from CLI"
    competitor = ARGS[1]
    resfolder0 = ARGS[2]
    if !isdir(resfolder0)
        mkdir(resfolder0)
    end    
    nλ = parse(Int, ARGS[3])
    nrep = parse(Int, ARGS[4])
    nfold = parse(Int, ARGS[5])
    one_se_rule = parse(Bool, ARGS[6])
end
timestamp = replace(strip(read(`date -Iseconds`, String)), ":" => "_")
resfolder = joinpath(resfolder0, "nrep$nrep-nfold$nfold-nlam$nλ-1se$(one_se_rule)-$competitor-$timestamp")
mkdir(resfolder)
@info "Results are saved into $resfolder"
benchmarking.(
    fs;
    σs = [0.1, 0.2, 0.4, 0.5, 1.0, 1.5, 2.0], # noise level to be surveyed
    jplot = false, # μerr vs σs
    nrep = nrep, # NB: for fast auto-generation procedure, only use nrep = 1; in the paper, use nrep = 100
    competitor = competitor, 
    nfold = nfold, # number of folds
    one_se_rule = one_se_rule, 
    nλ = nλ, # the number of λ to be searched
    rλ = 0.5, # the search region of λ, (1-rλ, 1+rλ)*λ
    resfolder = resfolder,
    verbose = false
)

# !!! tip "run from command line"
#     ```
#     julia examples/benchmark.jl ss_single_lambda /tmp 2 1 2 false
#     ```
#     
#     You can also enable the debug mode to print more internal steps as follows
#     ```
#     JULIA_DEBUG=MonotoneDecomposition julia examples/benchmark.jl ss_single_lambda /tmp 2 1 2 false
#     ```

# summary the results
MonotoneDecomposition.summary(resfolder = resfolder, format = "tex")
