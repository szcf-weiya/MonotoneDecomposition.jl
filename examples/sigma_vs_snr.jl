snr = 10
σs = Float64[]
for f in [x->x^2 x->x^3 x->exp(x) x -> 1 / (1 + exp(-5x)) "SE_1" "SE_0.1" "Mat12_1" "Mat12_0.1" "Mat32_1" "Mat32_0.1" "RQ_0.1_0.5" "Periodic_0.1_4"]
    append!(σs, mean([gen_data(100, nothing, f, snr = snr)[5] for _ in 1:100]))
end

serialize("../res/sigma_vs_snr/snr$snr.sil", σs)


σs01 = deserialize("../res/sigma_vs_snr/snr0.1.sil")
σs05 = deserialize("../res/sigma_vs_snr/snr0.5.sil")
σs1 = deserialize("../res/sigma_vs_snr/snr1.sil")
σs2 = deserialize("../res/sigma_vs_snr/snr2.sil")
σs10 = deserialize("../res/sigma_vs_snr/snr10.sil")

ntitle = 12
σs_from_snr = Array{AbstractMatrix{Float64}}(undef, ntitle)

for i = 1:ntitle
    σs_from_snr[i] = reshape([σs01[i], σs05[i], σs1[i], σs2[i], σs10[i]], (5, 1))
end

serialize("../res/sigma_vs_snr/all.sil", σs_from_snr)

σs_from_snr = deserialize("../res/sigma_vs_snr/all.sil")


######### 
arr_σs = zeros(100, 5)
for i in 1:100
    x, m1, m2, m3, m4, m5, σs = MonotoneDecomposition.gen_mono_data(n = 1000, σ = nothing, snr = 10000, equidistant = true)
    arr_σs[i, :] .= σs
end
μ_σs = mean(arr_σs, dims = 1)[:]
serialize("../res/sigma_vs_snr/test_snr10000.sil", μ_σs)

σs1w = deserialize("../res/sigma_vs_snr/test_snr10000.sil")
σs100 = deserialize("../res/sigma_vs_snr/test_snr100.sil")
σs1 = deserialize("../res/sigma_vs_snr/test_snr1.sil")