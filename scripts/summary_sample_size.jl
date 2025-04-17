using Serialization
using LaTeXStrings

## effect of sample sizes (fixJ)
resfolders = [
    "xps/nrep100-nfold10-nlam20-1se_false-bspl-2025-03-15T02_11_08-04_00-n20-snr_true/",
    "t460p/nrep100-nfold10-nlam20-1se_false-bspl-2025-03-21T10_31_21-04_00-n50-snr_true/",
    "t460p/nrep100-nfold10-nlam20-1se_false-bspl-2025-03-21T12_48_32-04_00-n100-snr_true/",
    "t460p/nrep100-nfold10-nlam20-1se_false-bspl-2025-03-21T12_48_35-04_00-n200-snr_true/",
    "t460p/nrep100-nfold10-nlam20-1se_false-bspl-2025-03-21T12_48_40-04_00-n500-snr_true/"
]

## single_lambda
resfolders = [
    "rocky/nrep100-nfold10-nlam1-1se_false-ss_single_lambda-2025-03-21T06_45_10+08_00-n20-snr_true",
    "rocky/nrep100-nfold10-nlam1-1se_false-ss_single_lambda-2025-03-21T06_45_26+08_00-n50-snr_true",
    "rocky/nrep100-nfold10-nlam1-1se_false-ss_single_lambda-2025-03-15T10_38_50+08_00-n100-snr_true",
    "rocky/nrep100-nfold10-nlam1-1se_false-ss_single_lambda-2025-03-21T06_46_54+08_00-n200-snr_true/",
    "rocky/nrep100-nfold10-nlam1-1se_false-ss_single_lambda-2025-03-21T06_47_57+08_00-n500-snr_true"
]

outputfolder = "~/Overleaf/paperMonoDecomp2/res/"

function summary_sample_size(resfolders, method = "bspl")
    @assert occursin(method, resfolders[1])
    # Candidate functions
    fs = ["x^2" "x^3" "exp(x)" "sigmoid" "SE_1" "SE_0.1" "Mat12_1" "Mat12_0.1" "Mat32_1" "Mat32_0.1" "RQ_0.1_0.5" "Periodic_0.1_4"]

    ns = [20, 50, 100, 200, 500]
    nn = length(ns)
    nf = length(fs)
    RES = Array{Array}(undef, nn, nf)
    for i in 1:nn
        for j in 1:nf
            RES[i, j] = deserialize(joinpath("../res/", resfolders[i], fs[j] * ".sil"))
        end
    end
    label = ["MDCS" "CS"]
    if occursin("ss", method)
        label .= ["MDSS" "SS"]
    end
    # mean error along sample size (the 3rd SNR)
    for j in 1:nf 
        μerr = vcat([mean(RES[i, j][:, [2, 4], 3], dims = 1) for i in 1:nn]...)
        σerr = vcat([std(RES[i, j][:, [2, 4], 3], dims = 1) / 10 for i in 1:nn]...)
        title = fs[j]
        if j <= 2
            title = latexstring(title)
        end
        fig = plot(ns, μerr, label = label, yerror = σerr, legend = :topright, xlab = L"n", ylab = "error", xticks = (ns, string.(ns)), title = title)
        savefig(fig, "$outputfolder/sample_size/$(method)_$(fs[j]).pdf")
    end    
end
