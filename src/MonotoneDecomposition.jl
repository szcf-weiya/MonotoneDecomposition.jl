module MonotoneDecomposition

include("utils.jl")
include("genfunc.jl")
include("mono_decomp.jl")

export gen_data,
       mono_decomp,
       gp,

       recover,
       pick_knots,
       div_into_folds,
       
       # CIs
       coverage_prob,
       conf_band_width,

       sol_as_optim,

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
       benchmarking_ss

end 