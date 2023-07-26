### Plots the spectrum of the community matrix 
### and countour lines

function spectrum(p)
    @assert p[:converged]
    return eigen(J(p[:equilibrium], p)).values
end

T(x, y, n, p) = (q = (p[:μ]-1/p[:K]^(2-p[:k])); 
            sum( 
                n.^2 ./ ((x .- q*n .- (p[:k] - 1).*n.^(p[:k] - 1)).^2 .+ y^2)
            )
            )

function boundary(p; overprint = false)

    if !haskey(p, :a) random_interactions!(p) end
    if !haskey(p, :equilibrium) evolve!(p) end
    if p[:scaled] p[:σ] = p[:σ]/sqrt(p[:S]) end
    if p[:scaled] p[:μ] = p[:μ]/p[:S] end

    s = spectrum(p)

    n = p[:equilibrium]
    if p[:converged]
        s = spectrum(p)
        
        s = s[real.(s) .> .8*minimum(real.(s))]

        X = range(1.2*minimum(real.(s)), .8*maximum(real.(s)); length = 100)
        Y = range(1.2*minimum(imag.(s)), 1.2*maximum(imag.(s)); length =  100)

        if overprint
            scatter!(s, alpha = .3, legend = false,
            aspect_ratio = 1, grid = false, color = COLOR_SUB35)
        else
            scatter(s, alpha = .3, legend = false,
            aspect_ratio = 1, grid = false, color = COLOR_SUB35)
        end

        Plots.contour!(X, Y, (x,y) -> T(x, y, n, p), levels = [1/(p[:σ])^2],
        linewidth = 2, color = :auto, colorbar = false)
    else
        print("The system is not feasible")
    end
end