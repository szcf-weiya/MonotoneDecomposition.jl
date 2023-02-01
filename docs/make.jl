ENV["PLOTS_TEST"] = "true"
ENV["GKSwstype"] = "100"
using Documenter, MonotoneDecomposition
using Literate
indir = joinpath(@__DIR__, "..", "examples")
outdir = joinpath(@__DIR__, "src", "examples")
for file in ["gp.jl", "md_SE.jl"]
    Literate.markdown(joinpath(indir, file), outdir; credit = false)
end
# using Pkg
# Pkg.activate("..")

makedocs(sitename="MonotoneDecomposition.jl",
        pages = [
            "Home" => "index.md",
            "Examples" => [
                "Gaussian Process" => "examples/gp.md",
                "Monotone Decomposition" => "examples/md_SE.md",
            ],
            "API" => "api.md"
        ]
)