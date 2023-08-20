using Distributed # not necessary when calling with `-p n` option
using Serialization
using MonotoneDecomposition

gurobi(1)

function postfix()
    timestamp = replace(strip(read(`date -Iseconds`, String)), ":" => "_")
    ## get current (if modified are not saved, parent) commit
    pcommit = strip(read(`git log --format="%h" -1`, String))
    return "$(pcommit)_$timestamp"    
end

function experiment_control_typeI_error(; resfolder = "/tmp", # "../res/mono_test/"
                                          nrep = 100, kw...)
    res = pmap(x->single_test_compare_mono(; kw...), 1:nrep)
    serialize(joinpath(resfolder, "res_mono_test_mono_sup_nrep$(nrep)_$(postfix()).sil"), res)
end

function experiment_bowman_example(; resfolder = "/tmp", nrep = 100, kw...)
    res = pmap(x->single_test_compare_bowman(; kw...), 1:nrep)
    serialize(joinpath(resfolder, "res_mono_test_bowman_sup_nrep$(nrep)_$(postfix()).sil"), res)
end

function experiment_ghosal_example(; resfolder = "/tmp", nrep = 100, kw...)
    res = pmap(x->single_test_compare_ghosal(; kw...), 1:nrep) 
    serialize(joinpath(resfolder, "res_mono_test_ghosal_sup_nrep$(nrep)_$(postfix()).sil"), res)   
end

