using MonotoneDecomposition
using Plots
# B splines with different J

f = x->x^5
x, y, xnew = gen_data(200, 0.1, f)
p = scatter(x, y)
plot!(p, x, f.(x))
for J = 4:6
    yhat, yhatnew = cubic_spline(J)(x, y, xnew)
    plot!(p, x, yhat, label = "Bspl (J = $J)")
end
p