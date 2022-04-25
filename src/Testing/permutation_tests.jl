
function permutest(disp ::Disp,n_perm = 10000)
    """
    This function permutes the residuals returned by 'dispersion' to provide global and pairwise
    P-values. Currently, no corrections are made for multiple testing. MultipleTesting.jl provides
    this functionality, if desired

    The function accepts:
        disp: The named tuple returned by 'dispersion'
        n_perm: the number of permutations

    This function returns a named tuple containing:
        P = Global P-value
        pairwise_P = pairwise P-value
        F = Global F-Statistic 
        pairwise_F = pairwise F-statistics

    https://github.com/juliangehring/MultipleTesting.jl
    """
    # read in values from `disp`
    F = disp.F
    R = disp.residuals
    F_pairs = disp.pairwise_F
    group = disp.group
    levels = disp.levels
    level_names = string.(levels)
    inds = [group .== level for level in levels]
    n = length(levels)

    # pre-calculate values for ANOVA function
    r = vcat(R...)
    X = vcat(R...)
    N = length(X)
    nj = Tuple(length.(R))
    X̅ = mean(X)
    k = length(R)
    # generate empty containers for F values
    perm = Vector{Float64}(undef,n_perm)  
    #rs = copy(R)#Vector{Vector}(undef,n)  
    # run permutation
    @inbounds for p in 1:n_perm
        shuffle!(r)
        rs = [view(r,inds[i]) for i in 1:n]
   
        perm[p] =f(rs,N,nj,X̅,k)
    end
    ppairs =zeros(k,k)

    @inbounds for i in 1:k-1
        for j in i+1:k
            r_1 = r[group .== levels[i]]
            r_2 = r[group .== levels[j]]
            pstat =permutest(r_1,r_2, F_pairs[j,i],n_perm)
            ppairs[j,i] = pstat
        end
    end
    # calculate P-values and return P and F values
    P = sum(perm .> F)/n_perm
    P_pairs = AxisArray(ppairs,level_names,level_names)
    F_pairs = AxisArray(LowerTriangular(F_pairs), level_names,level_names)
    
    return (P =P,pairwise_P = P_pairs,F = F, pairwise_F =F_pairs)
end

function permutest(r_1,r_2 ,F,n_perm = 1000)
    r = vcat(r_1,r_2)  
    nj = (length(r_1),length(r_2))
    n = 2
    N = sum(nj)
    # pre-calculate values for ANOVA function
    X̅ = mean(r)
    k = 2
    inds = (1:nj[1],1:nj[2])
    # generate empty containers for F values
    perm = Vector{Float64}(undef,n_perm) 
    # run permutation
    @inbounds for p in 1:n_perm
        shuffle!(r)
        rs = (view(r,inds[1]) ,view(r,inds[2]) )
        perm[p] =f(rs,N,nj,X̅,k)
    end  
    P = sum(perm .> F)/n_perm
    return P
end


get_tri(x) = [x[i,j] for i in 2:size(x)[1] for j in 1:i]
function mantel(x_1,x_2, n_perm = 1000)
    """
    A permutation based test of the correlation between two matrices. Given two dissimilarity/distance matrices of equal size, calculates the P-value as
    the probability of obtaining the same or greater correlation coefficient after permuting one of the matrices.
    """
    x = get_tri(x_1)
    y = get_tri(x_2)
    perm_c = Vector{Float64}(undef,n_perm)

    c = cor(x,y)
    @inbounds for i in 1:n_perm
        shuffle!(y)
        perm_c[i] = cor(x,y)
    end
    
    sum(perm_c .> c)/n_perm
end



