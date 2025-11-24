using Serialization
using LaTeXStrings
using StatsBase

## effect of sample sizes (fixJ)
resfolders = [
    "xps/nrep100-nfold10-nlam20-1se_false-bspl-2025-03-15T02_11_08-04_00-n20-snr_true/",
    "t460p/nrep100-nfold10-nlam20-1se_false-bspl-2025-03-21T10_31_21-04_00-n50-snr_true/",
    "t460p/nrep100-nfold10-nlam20-1se_false-bspl-2025-03-21T12_48_32-04_00-n100-snr_true/",
    "t460p/nrep100-nfold10-nlam20-1se_false-bspl-2025-03-21T12_48_35-04_00-n200-snr_true/",
    "t460p/nrep100-nfold10-nlam20-1se_false-bspl-2025-03-21T12_48_40-04_00-n500-snr_true/"
]

## with random seed
resfolders = [
    "zju/nrep100-nfold10-nlam20-1se_false-bspl-2025-10-04T12_15_00+08_00-n20-snr_true",
    "zju/nrep100-nfold10-nlam20-1se_false-bspl-2025-10-04T12_15_08+08_00-n50-snr_true",
    "zju/nrep100-nfold10-nlam20-1se_false-bspl-2025-10-04T12_15_12+08_00-n100-snr_true",
    "zju/nrep100-nfold10-nlam20-1se_false-bspl-2025-10-04T12_15_16+08_00-n200-snr_true",
    "zju/nrep100-nfold10-nlam20-1se_false-bspl-2025-10-04T12_15_20+08_00-n500-snr_true"
]

## single_lambda
resfolders = [
    "rocky/nrep100-nfold10-nlam1-1se_false-ss_single_lambda-2025-03-21T06_45_10+08_00-n20-snr_true",
    "rocky/nrep100-nfold10-nlam1-1se_false-ss_single_lambda-2025-03-21T06_45_26+08_00-n50-snr_true",
    "rocky/nrep100-nfold10-nlam1-1se_false-ss_single_lambda-2025-03-15T10_38_50+08_00-n100-snr_true",
    "rocky/nrep100-nfold10-nlam1-1se_false-ss_single_lambda-2025-03-21T06_46_54+08_00-n200-snr_true/",
    "rocky/nrep100-nfold10-nlam1-1se_false-ss_single_lambda-2025-03-21T06_47_57+08_00-n500-snr_true"
]

# fix_ratio
resfolders = [
    "rocky/nrep100-nfold10-nlam20-1se_false-ss_fix_ratio-2025-03-25T00_52_07+08_00-n20-snr_true",
    "rocky/nrep100-nfold10-nlam20-1se_false-ss_fix_ratio-2025-03-25T00_53_43+08_00-n50-snr_true",
    "rocky/nrep100-nfold10-nlam20-1se_false-ss_fix_ratio-2025-03-25T00_54_49+08_00-n100-snr_true",
    "rocky/nrep100-nfold10-nlam20-1se_false-ss_fix_ratio-2025-03-25T00_55_29+08_00-n200-snr_true",
    "rocky/nrep100-nfold10-nlam20-1se_false-ss_fix_ratio-2025-03-25T00_55_21+08_00-n500-snr_true"
]

resfolders = [
    "zju/nrep100-nfold10-nlam20-1se_false-ss_fix_ratio-2025-10-03T11_32_34+08_00-n20-snr_true",
    "zju/nrep100-nfold10-nlam20-1se_false-ss_fix_ratio-2025-10-03T11_32_38+08_00-n50-snr_true",
    "zju/nrep100-nfold10-nlam20-1se_false-ss_fix_ratio-2025-10-03T11_30_34+08_00-n100-snr_true",
    "zju/nrep100-nfold10-nlam20-1se_false-ss_fix_ratio-2025-10-03T11_32_26+08_00-n200-snr_true",
    "zju/nrep100-nfold10-nlam20-1se_false-ss_fix_ratio-2025-10-03T11_32_30+08_00-n500-snr_true"
]

# cvbspl2
resfolders = [
    "rocky/nrep100-nfold10-nlam20-1se_false-cvbspl2-2025-04-18T07_11_59+08_00-n20-snr_true",
    "rocky/nrep100-nfold10-nlam20-1se_false-cvbspl2-2025-04-18T07_19_05+08_00-n50-snr_true",
    "rocky/nrep100-nfold10-nlam20-1se_false-cvbspl2-2025-04-18T07_19_13+08_00-n100-snr_true",
    "rocky/nrep100-nfold10-nlam20-1se_false-cvbspl2-2025-04-18T07_20_02+08_00-n200-snr_true",
    "rocky/nrep100-nfold10-nlam20-1se_false-cvbspl2-2025-04-18T07_20_18+08_00-n500-snr_true"
]

# !! cvbspl2-2025-10-10

resfolders = [
    "zju/nrep100-nfold10-nlam20-1se_false-cvbspl2-2025-10-10T21_04_20+08_00-n20-snr_true",
    "zju/nrep100-nfold10-nlam20-1se_false-cvbspl2-2025-10-10T21_04_24+08_00-n50-snr_true",
    "zju/nrep100-nfold10-nlam20-1se_false-cvbspl2-2025-10-10T21_04_29+08_00-n100-snr_true",
    "zju/nrep100-nfold10-nlam20-1se_false-cvbspl2-2025-10-10T21_04_32+08_00-n200-snr_true",
    "zju/nrep100-nfold10-nlam20-1se_false-cvbspl2-2025-10-10T21_04_36+08_00-n500-snr_true"
]

# double_grid as single_lambda
resfolders = [
    "zju/nrep100-nfold10-nlam20-1se_false-ss_double_grid-2025-08-10T13_35_03+08_00-n20-snr_true",
    "zju/nrep100-nfold10-nlam20-1se_false-ss_double_grid-2025-08-10T13_35_07+08_00-n50-snr_true",
    "zju/nrep100-nfold10-nlam20-1se_false-ss_double_grid-2025-08-09T13_49_34+08_00-n100-snr_true",
    "zju/nrep100-nfold10-nlam20-1se_false-ss_double_grid-2025-08-10T13_35_11+08_00-n200-snr_true",
    "zju/nrep100-nfold10-nlam20-1se_false-ss_double_grid-2025-08-10T09_12_01+08_00-n500-snr_true"
]

# double_grid
resfolders = [
    "zju/nrep100-nfold10-nlam20-1se_false-ss_double_grid-2025-10-02T22_44_13+08_00-n20-snr_true",
    "zju/nrep100-nfold10-nlam20-1se_false-ss_double_grid-2025-10-02T22_44_17+08_00-n50-snr_true",
    "zju/nrep100-nfold10-nlam20-1se_false-ss_double_grid-2025-10-02T22_43_57+08_00-n100-snr_true",
    "zju/nrep100-nfold10-nlam20-1se_false-ss_double_grid-2025-10-02T22_44_06+08_00-n200-snr_true",
    "zju/nrep100-nfold10-nlam20-1se_false-ss_double_grid-2025-10-02T22_44_10+08_00-n500-snr_true"
]

# double_grid as single_lambda
resfolders = [
    ""
]

outputfolder = "~/Overleaf/paperMonoDecomp2/res/"
snrs = [0.1, 0.5, 1, 2, 10]

# by default, show the results when idx_snr = 3, i.e., snr = 1.0
function summary_sample_size(resfolders, method0 = "bspl", idx_snr = 3, method1 = "")
    @assert occursin(method0, resfolders[1])
    method = method0 * method1
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
        μerr = vcat([mean(RES[i, j][:, [2, 4], idx_snr], dims = 1) for i in 1:nn]...)
        σerr = vcat([std(RES[i, j][:, [2, 4], idx_snr], dims = 1) / 10 for i in 1:nn]...)
        title = fs[j]
        if j <= 2
            title = latexstring(title)
        end
        fig = plot(ns, μerr, label = label, yerror = σerr, legend = :topright, xlab = L"n", ylab = "MSPE", xticks = (ns, string.(ns)), title = title,
                    labelfontsize = 14, legendfontsize = 14, tickfontsize = 12, titlefontsize = 16)
        savefig(fig, "$outputfolder/sample_size/$(method)_$(fs[j])_snr$(snrs[idx_snr]).pdf")
    end
end


##  summary_sample_size(resfolders, "ss_double_grid", 2, "_as_single_lambda")
##  summary_sample_size(resfolders, "cvbspl2", 2, "")
## summary_sample_size(resfolders, "bspl-2025-10-04", 2)
## summary_sample_size(resfolders, "cvbspl2-2025-10-10", 2)
## summary_sample_size(resfolders, "cvbspl2-2025-10-10", 1) # to use
## summary_sample_size(resfolders, "ss_fix_ratio-2025-10-03", 1)
## summary_sample_size(resfolders, "ss_double_grid-2025-10-02", 1) # to use

