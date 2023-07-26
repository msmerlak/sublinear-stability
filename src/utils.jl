#= Defines of useful functions and constants =#

ppart(x) = max(x, 0)

function offdiag(A::AbstractMatrix)
    [A[ι] for ι in CartesianIndices(A) if ι[1] ≠ ι[2]]
end 

function expand(c::AbstractDict)

    dict = typeof(c)

    iterable_fields = filter(k -> typeof(c[k]) <: AbstractVector, keys(c))
    non_iterables = setdiff(keys(c), iterable_fields)

    iterable_dict = dict(iterable_fields .=> getindex.(Ref(c), iterable_fields))
    non_iterable_dict = dict(non_iterables .=> getindex.(Ref(c), non_iterables))

    vec(
        map(Iterators.product(values(iterable_dict)...)) do vals
            dd = [k=>convert(eltype(c[k]),v) for (k,v) in zip(keys(iterable_dict),vals)]
            if isempty(non_iterable_dict)
                Dict(dd)
            elseif isempty(iterable_dict)
                non_iterable_dict
            else
                # We can't use merge here because it promotes types.
                # The uniqueness of the dictionary keys is guaranteed.
                dict(dd..., collect(non_iterable_dict)...)
            end
        end
    )
end

function xlogx(x)
    if x <= 0.
        return 0.
    else
        return x*log(x)
    end
end

function Ω(c::Vector)
    return exp(-sum(xlogx.(c./norm(c, 1))))/length(c)
end

function uniquetol(A; kws...)
    S = []
    for a in A
         if !any(s -> isapprox(s, a; kws...), S)
             push!(S, a)
         end
    end
    return S
end

COLOR_LOG = :turquoise3
COLOR_SUB = :darkgoldenrod1

COLOR_LOG49 = "#009af9"
COLOR_SUB49 = "#fb8d00"

COLOR_LOG35 = "#0070b2"
COLOR_SUB35 = "#b26300" 

K = 100.
CUTOFF = 1.