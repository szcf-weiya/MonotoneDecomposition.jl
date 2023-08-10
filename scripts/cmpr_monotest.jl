using Distributed # not necessary when calling with `-p n` option
using LaTeXStrings
using Serialization
using MonotoneDecomposition

gurobi()

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

function summary_res()
    res = deserialize("../res/mono_test/res_mono_test_bowman.sil")
    # too many lines
    plot(a[:,:,1]', ls = [:dash :dot :dashdot :solid], markershape = :star5)
    plot!(a[:,:,2]', ls = [:dash :dot :dashdot :solid], markershape = :dtriangle)
    plot!(a[:,:,3]', ls = [:dash :dot :dashdot :solid], markershape = :circle)
    μ = mean(res) # size (4, 3, 5, 5)
    A = Array{Matrix, 1}(undef, 5)
    for i = 1:5
        # cat along sigma noise level, then the format is (50, 100, 200, 50 100 200, 50 100 200)
        A[i] = hcat([μ[:, :, j, i] for j = [1,2,5]]...)
    end
    print2tex(A, ["MCS", "MSS", "Meyer", "Ghosal", "Bowman"], ["σ = 0.001", "σ = 0.01", "σ = 0.1"], subcolnames = ["n = 50", "100", "200"], subrownames = string.([0,0.15,0.25,0.45]), colnames_of_rownames = ["Methods", "a"], file = "res_mono_test_bowman.tex")

    plot(μ[:, 3, 5, :], markershape = [:star5 :diamond :dtriangle :rect :circle], legend = :topleft, xticks = (1:5, string.([0, 0.15, 0.25, 0.45])), alpha = 0.75, label = ["MD(CS)" "MD(SS)" "Meyer" "Ghosal" "Bowman"], title = "n = 200, σ = 0.1")
    hline!([0.05], ls = :dash, label = L"\alpha = 0.05", xlab = "a (Bowman's curves)")
    savefig("../res/mono_test/fig_n200_sigma0.1_bowman.pdf")

    # 
    res2 = deserialize("../res/mono_test/res_mono_test_ghosal.sil")
    μ2 = mean(res2) # 3x4x5
end

function summary_res4x5_3()
    res = deserialize("res_mono_test_bowman_4x5+3_no1se.sil")
    μ = mean(res)
    A = Array{Matrix, 1}(undef, 7)
    for (i, k) in enumerate([5, 10, 15, 20, 21, 22, 23])
        # cat along sigma noise level, then the format is (50, 100, 200, 50 100 200, 50 100 200)
        A[i] = hcat([μ[:, :, j, k] for j = [1,2,3]]...)
    end
    print2tex(A, ["v1: MD (CS)", "v1: MD (SS)", "v2: MD (CS)", "v2: MD (SS)", "Meyer", "Ghosal", "Bowman"], ["σ = 0.001", "σ = 0.01", "σ = 0.1"], subcolnames = ["n = 50", "100", "200"], subrownames = string.([0,0.15,0.25,0.45]), colnames_of_rownames = ["Methods", "a"], file = "res_mono_test_bowman_0511v2.tex", format = "raw")


    # res2 = deserialize("res_mono_test_ghosal_4x5+3_no1serule.sil")
    res2 = deserialize("res_mono_test_ghosal_4x5+3.sil")
    μ2 = mean(res2)
    A2 = Array{Matrix, 1}(undef, 7)
    for (i, k) in enumerate([5, 10, 15, 20, 21, 22, 23])
        # cat along sigma noise level, then the format is (50, 100, 200, 50 100 200, 50 100 200)
        A2[i] = hcat([μ2[1, j, :, k] for j = [1,2,3]]...)'
    end
    print2tex(A2, ["v1: MD (CS)", "v1: MD (SS)", "v2: MD (CS)", "v2: MD (SS)", "Meyer", "Ghosal", "Bowman"], ["n = 200"], subcolnames = "m" .* string.(1:4), subrownames = string.([0.001, 0.01,0.1]), colnames_of_rownames = ["Methods", "noise"], file = "res_mono_test_ghosal_0511v2_1se.tex", format = "raw")
end

## before the commit https://github.com/szcf-weiya/Clouds/commit/81ab460a6932c3e294b36b650093987a9776bff3#diff-d3b2f4590eb36d701d869fb64569fce96d45c5a4ca6f2750eff89ee6066798f4
## where the results for _sup and _supss become 3-vector instead of 2-vector
function summary_res_sup()
    resfolder = "../res/mono_test"
    res = deserialize(joinpath(resfolder, "res_mono_test_bowman_sup7.sil"))
    res2 = deserialize(joinpath(resfolder, "res_mono_test_ghosal_sup7.sil"))
    res3 = deserialize(joinpath(resfolder, "res_mono_test_mono_sup.sil"))

    # res = deserialize("../res/mono_test/res_mono_test_bowman_sup7_500.sil")
    # res2 = deserialize("../res/mono_test/res_mono_test_ghosal_sup7_500.sil")

    μ = mean(res) # 4 x 3 x 3 x 7
    μ2 = mean(res2) # 3 x 3 x 4 x 7
    μ3 = mean(res3) # 3 x 3 x 5 x 7
    A = Array{Matrix, 1}(undef, 7)
    A2 = Array{Matrix, 1}(undef, 7)
    A3 = Array{Matrix, 1}(undef, 7)
    for (i, k) in enumerate(vcat(4:7, 1:3))
        # cat along sigma noise level, then the format is (50, 100, 200, 50 100 200, 50 100 200)
        A[i] = hcat([μ[:, :, j, k] for j = [1,2,3]]...)
    end
    for (i, k) in enumerate(vcat(4:7, 1:3))
        # cat along sigma noise level, then the format is (50, 100, 200, 50 100 200, 50 100 200)
        A2[i] = hcat([μ2[:, j, :, k]' for j = [1,2,3]]...)
    end

    for (i, k) in enumerate(vcat(4:7, 1:3))
        # cat along sigma noise level, then the format is (50, 100, 200, 50 100 200, 50 100 200)
        A3[i] = hcat([μ3[:, j, :, k]' for j = [1,2,3]]...)
    end
    print2tex(A, ["MD (CS) sup", "MD (CS)", "MD (SS) sup", "MD (SS)", "Meyer", "Ghosal", "Bowman"], ["\$\\sigma = 0.001\$", "\$\\sigma = 0.01\$", "\$\\sigma = 0.1\$"], subcolnames = ["n = 50", "100", "200"], subrownames = string.([0,0.15,0.25,0.45]), colnames_of_rownames = ["Methods", "a"], file = "../res/mono_test/res_mono_test_bowman_sup_500.tex", format="raw")
    print2tex(A2, ["MD (CS) sup", "MD (CS)", "MD (SS) sup", "MD (SS)", "Meyer", "Ghosal", "Bowman"], ["\$\\sigma = 0.001\$", "\$\\sigma = 0.01\$", "\$\\sigma = 0.1\$"], subcolnames = ["n = 50", "100", "200"], subrownames = "m" .* string.(1:4), colnames_of_rownames = ["Methods", "curve"], file = "../res/mono_test/res_mono_test_ghosal_sup_500.tex", format="raw")

    print2tex(A3, ["MD (CS) sup", "MD (CS)", "MD (SS) sup", "MD (SS)", "Meyer", "Ghosal", "Bowman"],  ["\$\\sigma = 0.001\$", "\$\\sigma = 0.01\$", "\$\\sigma = 0.1\$"], subcolnames=["n = 50", "100", "200"],subrownames = ["\$x\$", "\$x^3\$", "\$x^{1/3}\$", "\$e^x\$", "\$1/(1+e^{-x})\$"], colnames_of_rownames = ["Methods", "curve"], file = "../res/mono_test/res_mono_test_mono_sup.tex", format="raw")

    plt = plot(μ[:, 3, 1, vcat(4:7, 1:3)], markershape = [:star5 :star5 :diamond :diamond :dtriangle :rect :circle], legend = :topleft, xticks = (1:4, string.([0, 0.15, 0.25, 0.45])), alpha = [0.9 0.2 0.9 0.2 0.9 0.9 0.9], label = ["MD (CS) sup" "MD (CS)" "MD (SS) sup" "MD (SS)" "Meyer" "Ghosal" "Bowman"], title = "n = 200, σ = 0.001")
    hline!(plt, [0.05], ls = :dash, xlab = "Bowman's curves", label = "")
    old_ticks = yticks(plt[1])
    merged_ticks = (vcat(old_ticks[1], 0.05), vcat(old_ticks[2], "0.05"))
    yticks!(merged_ticks)
    savefig("../res/mono_test/fig_n200_sigma0.001_bowman_sup.pdf")

    plt = plot(μ2[3, 3, :, vcat(4:7, 1:3)], markershape = [:star5 :star5 :diamond :diamond :dtriangle :rect :circle], legend = :topleft, xticks = (1:4, "m" .* string.(1:4)), alpha = [0.9 0.2 0.9 0.2 0.9 0.9 0.9], label = ["MD (CS) sup" "MD (CS)" "MD (SS) sup" "MD (SS)" "Meyer" "Ghosal" "Bowman"], title = "n = 200, σ = 0.1")
    hline!(plt, [0.05], ls = :dash, xlab = "Ghosal's curves", label = "")
    old_ticks = yticks(plt[1])
    merged_ticks = (vcat(old_ticks[1], 0.05), vcat(old_ticks[2], "0.05"))
    yticks!(merged_ticks)

    # https://stackoverflow.com/questions/51651364/how-to-add-specific-x-ticks-labels-and-vertical-lines-to-a-plot-using-plots-jl
    savefig("../res/mono_test/fig_n200_sigma0.1_ghosal_sup.pdf")


end

function summary_res4x3_3()
    res = deserialize("res_mono_test_bowman_4x3+3_all.sil")
    res2 = deserialize("res_mono_test_ghosal_4x3+3_all.sil")
    # too many lines
    plot(a[:,:,1]', ls = [:dash :dot :dashdot :solid], markershape = :star5)
    plot!(a[:,:,2]', ls = [:dash :dot :dashdot :solid], markershape = :dtriangle)
    plot!(a[:,:,3]', ls = [:dash :dot :dashdot :solid], markershape = :circle)
    μ = mean(res) # size (4, 3, 5, 15)
    μ2 = mean(res2) # size (3, 3, 4, 15)
    A = Array{Matrix, 1}(undef, 5)
    for (i, k) in enumerate([1, 4, 13, 14, 15])
        # cat along sigma noise level, then the format is (50, 100, 200, 50 100 200, 50 100 200)
        A[i] = hcat([μ[:, :, j, k] for j = [1,2,5]]...)
    end

    A2 = Array{Matrix, 1}(undef, 5)
    for (i, k) in enumerate([1, 4, 13, 14, 15])
        # cat along sigma noise level, then the format is (50, 100, 200, 50 100 200, 50 100 200)
        A2[i] = hcat([μ2[:, :, j, k] for j = [1,2,5]]...)
    end

    print2tex(A, ["MD (CS)", "MD (SS)", "Meyer", "Ghosal", "Bowman"], ["σ = 0.001", "σ = 0.01", "σ = 0.1"], subcolnames = ["n = 50", "100", "200"], subrownames = string.([0,0.15,0.25,0.45]), colnames_of_rownames = ["Methods", "a"], file = "res_mono_test_bowman_0511.tex")

    plot(μ[:, 3, 5, [1, 4, 13, 14, 15]], markershape = [:star5 :diamond :dtriangle :rect :circle], legend = :topleft, xticks = (1:5, string.([0, 0.15, 0.25, 0.45])), alpha = 0.75, label = ["MD(CS)" "MD(SS)" "Meyer" "Ghosal" "Bowman"], title = "n = 200, σ = 0.1")
    hline!([0.05], ls = :dash, label = L"\alpha = 0.05", xlab = "a (Bowman's curves)")
    savefig("../res/mono_test/fig_n200_sigma0.1_bowman_0511.pdf")

    # sigma = 0.001
    plot(μ[:, 3, 1, [1, 4, 13, 14, 15]], markershape = [:star5 :diamond :dtriangle :rect :circle], legend = :topleft, xticks = (1:5, string.([0, 0.15, 0.25, 0.45])), alpha = 0.75, label = ["MD(CS)" "MD(SS)" "Meyer" "Ghosal" "Bowman"], title = "n = 200, σ = 0.001")
    hline!([0.05], ls = :dash, label = L"\alpha = 0.05", xlab = "a (Bowman's curves)")
    savefig("../res/mono_test/fig_n200_sigma0.001_bowman_0511.pdf")

    # sigma = 0.01
    plot(μ[:, 3, 2, [1, 4, 13, 14, 15]], markershape = [:star5 :diamond :dtriangle :rect :circle], legend = :topleft, xticks = (1:5, string.([0, 0.15, 0.25, 0.45])), alpha = 0.75, label = ["MD(CS)" "MD(SS)" "Meyer" "Ghosal" "Bowman"], title = "n = 200, σ = 0.01")
    hline!([0.05], ls = :dash, label = L"\alpha = 0.05", xlab = "a (Bowman's curves)")
    savefig("../res/mono_test/fig_n200_sigma0.01_bowman_0511.pdf")

    # sigma = 0.001, ghosal
    plot(μ2[3, 1, :, [1, 4, 13, 14, 15]], markershape = [:star5 :diamond :dtriangle :rect :circle], legend = :topleft, xticks = (1:4, "m" .* string.(1:4)), alpha = 0.75, label = ["MD(CS)" "MD(SS)" "Meyer" "Ghosal" "Bowman"], title = "n = 200, σ = 0.001")
    hline!([0.05], ls = :dash, label = L"\alpha = 0.05", xlab = "(Ghosal's curves)")
    savefig("../res/mono_test/fig_n200_sigma0.001_ghosal_0511.pdf")

    # sigma = 0.01, ghosal
    plot(μ2[3, 2, :, [1, 4, 13, 14, 15]], markershape = [:star5 :diamond :dtriangle :rect :circle], legend = :topleft, xticks = (1:4, "m" .* string.(1:4)), alpha = 0.75, label = ["MD(CS)" "MD(SS)" "Meyer" "Ghosal" "Bowman"], title = "n = 200, σ = 0.01")
    hline!([0.05], ls = :dash, label = L"\alpha = 0.05", xlab = "(Ghosal's curves)")
    savefig("../res/mono_test/fig_n200_sigma0.01_ghosal_0511.pdf")

    # sigma = 0.1
    plot(μ2[3, 3, :, [1, 4, 13, 14, 15]], markershape = [:star5 :diamond :dtriangle :rect :circle], legend = :topleft, xticks = (1:4, "m" .* string.(1:4)), alpha = 0.75, label = ["MD(CS)" "MD(SS)" "Meyer" "Ghosal" "Bowman"], title = "n = 200, σ = 0.1")
    hline!([0.05], ls = :dash, label = L"\alpha = 0.05", xlab = "(Ghosal's curves)")
    savefig("../res/mono_test/fig_n200_sigma0.1_ghosal_0511.pdf")

    # 
    res2 = deserialize("../res/mono_test/res_mono_test_ghosal.sil")
    μ2 = mean(res2) # 3x4x5
end