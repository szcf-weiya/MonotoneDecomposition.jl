module MonotoneDecomposition

include("utils.jl")
include("genfunc.jl")
include("mono_decomp.jl")
include("mono_test.jl")

export gen_data,
       gen_kern,
       gp,

       mono_decomp,

       recover,
       div_into_folds,
       
       # CIs
       coverage_prob,
       conf_band_width,

       predict, # extend StatsBase.predict (The extended method can be also exported, see https://github.com/JuliaStats/StatsBase.jl/blob/master/src/StatsBase.jl)

       mono_decomp_cs,
       mono_decomp_ss,
       cv_mono_decomp_ss,
       cv_mono_decomp_cs,
       cvfit,
       cvfit_gss,
       cubic_spline,
       smooth_spline,
       cv_cubic_spline,
       cvplot,

       ecos,
       gurobi,

       build_model,
       
       benchmarking,
       benchmarking_cs,
       benchmarking_ss,

       single_test_compare_bowman,
       single_test_compare_ghosal,
       single_test_compare_mono,
       single_test_compare_desc
end 