### This is a sketch, needs updating


T(x, y, n, p) = (q = (p[:μ]-1/p[:K]^(2-p[:k])); 
            sum( 
                n.^2 ./ ((x .- q*n .- (p[:k] - 1).*n.^(p[:k] - 1)).^2 .+ y^2)
            )
            )

function boundary(p)
    a = interaction_matrix(p)
    n = equilibrium(p, a)
    println(n)
    if all(n .> p[:cutoff])
        s = spectrum(p, a)
        
        s = s[real.(s) .> .99*minimum(real.(s))]
        # println(T(0, 0, n, p) - 1/(p[:σ])^2)
        # println( p[:S]/((1-p[:k])p[:μ]p[:S])^2 - 1/(p[:σ])^2 )
        X = range(minimum(real.(s)), maximum(real.(s)); length = 200)
        Y = range(1.1*minimum(imag.(s)), 1.1*maximum(imag.(s)); length =  200)

        scatter(s, alpha = .3, legend = false, dpi = 500,
        aspect_ratio = 1)

        Plots.contour!(X, Y, (x,y) -> T(x, y, n, p), levels = [1/(p[:σ])^2],
        linewidth = 2)
    else
        print("The system is not feasible")
    end
    # println((1-p[:k])maximum(a*n))
    # vline!([- (1-p[:k])maximum(a*n)], color = "green")
    # vline!([- (1-p[:k])minimum(a*n)], color = "green")
end