ENV["PLOTS_TEST"] = "true"
ENV["GKSwstype"] = "100"
using Documenter, MonotoneDecomposition
using Literate
indir = joinpath(@__DIR__, "..", "examples")
outdir = joinpath(@__DIR__, "src", "examples")
for file in ["gp.jl", "md_SE.jl", "sample_size.jl", "anyJ,jl", "anylam.jl", "benchmark.jl", "benchmark_parallel.jl"]
    Literate.markdown(joinpath(indir, file), outdir; credit = false)
end
# using Pkg
# Pkg.activate("..")

makedocs(sitename="MonotoneDecomposition.jl",
        pages = [
            "Home" => "index.md",
            "Examples" => [
                "Gaussian Process (GP)" => "examples/gp.md",
                "Monotone Decomposition (MD)" => [
                    "On A GP Random Function" => "examples/md_SE.md",
                    "Effects of Sample Size" => "examples/sample_size.md",
                    "CS vs MDCS under any J" => "examples/anyJ.md",
                    "SS vs MDSS under any λ" => "examples/anylam.md"
                ],
                "Benchmarking" => "examples/benchmark.md",
                "Benchmarking in Parallel" => "examples/benchmark_parallel.md"
            ],
            "API" => "api.md"
        ]
)

deploydocs(
    repo = "github.com/szcf-weiya/MonotoneDecomposition.jl"
)