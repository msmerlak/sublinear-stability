### This is a sketch, needs updating


function spectrum(p)
    @assert p[:converged]
    return eigen(J(p[:equilibrium], p)).values
end


T(x, y, n, p) = (q = (p[:μ]-1/p[:K]^(2-p[:k])); 
            sum( 
                n.^2 ./ ((x .- q*n .- (p[:k] - 1).*n.^(p[:k] - 1)).^2 .+ y^2)
            )
            )



function boundary(p)

    if !haskey(p, :a) random_interactions!(p) end
    if !haskey(p, :equilibrium) evolve!(p) end
    s = spectrum(p)


    n = p[:equilibrium]
    if p[:converged] #&& all(n .> p[:n0]) 
        s = spectrum(p)
        
        s = s[real.(s) .> .8*minimum(real.(s))]
        # println(T(0, 0, n, p) - 1/(p[:σ])^2)
        # println( p[:S]/((1-p[:k])p[:μ]p[:S])^2 - 1/(p[:σ])^2 )
        X = range(minimum(real.(s)), maximum(real.(s)); length = 100)
        Y = range(1.1*minimum(imag.(s)), 1.1*maximum(imag.(s)); length =  100)

        scatter!(s, alpha = .3, legend = false,
        aspect_ratio = 1)

        Plots.contour!(X, Y, (x,y) -> T(x, y, n, p), levels = [1/(p[:σ])^2],
        linewidth = 2, color = :auto, colorbar = false)
    else
        print("The system is not feasible")
    end
    # println((1-p[:k])maximum(a*n))
    # vline!([- (1-p[:k])maximum(a*n)], color = "green")
    # vline!([- (1-p[:k])minimum(a*n)], color = "green")
end