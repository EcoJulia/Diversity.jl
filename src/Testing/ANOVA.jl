
# Compute F-statistic, given vector of vectors
#ANOVA
function f(R)
    X = vcat(R...)
    N = length(X)
    nj = length.(R)
    X̅ = mean(X)
    X̅j = mean.(R)
    k = length(R)
    top = sum(nj .*(X̅j .- X̅) .^2) /(k-1)
    bottom = sum([sum((R[i] .- X̅j[i]) .^2) for i in 1:k])/(N-k)
    return top/bottom
end
function f_slow(R,N,nj,X̅,k)
    X̅j = mean.(R)
    top = sum(nj .*(X̅j .- X̅) .^2) /(k-1)
    bottom = sum([sum((R[i] .- X̅j[i]) .^2) for i in 1:k])/(N-k)
    return top/bottom
end

function f(R,N,nj,X̅,k)
    X̅j = mean.(R)
    top = 0.0
    bottom = 0.0
    @inbounds for i in 1:k
        top += nj[i] * (X̅j[i] - X̅)^2
        for j in 1:nj[i]
            bottom += (R[i][j] - X̅j[i])^2
        end
    end
    top/= (k-1)
    bottom /= (N-k)
    return top/bottom
end

#Pairwise ANOVA
function f_pairs(R) 
    n = length(R)
    F = zeros(n,n)
    for j in 1:n-1
        for i in j+1:n
            F[i,j] = f([R[i],R[j]])
        end
    end
    return F
end
function f_pairs(R ::Vector,N ::Array{Int},nj ::Matrix{Tuple},X̅ ::Matrix{Float64}) 
    n = length(R)
    F = zeros(n,n)
    for j in 1:n-1
        for i in j+1:n
            F[i,j] = f([R[i],R[j]],N[i,j],nj[i,j],X̅[i,j],2)
        end
    end
    return F
end
# get parameter values corresponding to each pairwise combination of groups
function get_pars(R)
    n = length(R)
    N_p = Array{Int}(undef,n,n)
    nj_p= Array{Tuple}(undef,n,n)
    X̅_p = Array{Float64}(undef,n,n)

    for j in 1:n-1
        for i in j+1:n
            N_p[i,j] = length(R[i]) + length(R[j])
            nj_p[i,j] = (length(R[i]) , length(R[j]))
            X̅_p[i,j] = mean(vcat(R[i]...,R[j]...))
        end
    end
    return N_p,nj_p,X̅_p
end
