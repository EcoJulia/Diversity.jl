# structs

struct Disp
    F ::Float64
    pairwise_F ::Array{Float64} 
    medians ::Union{Tuple,Vector}
    residuals ::Vector{Matrix{Float64}} 
    means ::Vector{Float64}
    group  ::Vector
    levels  ::Vector
    group_residuals ::Vector
end
#Helper functions

# calculate distances from medians
function get_residuals(x,c)
    d = x .-c
    sum(d .^2,dims =2) 
end

# doube center matrix
function double_center(D)
    A = D .- mean(D,dims = 1)
    return A .- mean(A,dims = 2) 
end

# calculate medians 
# Get spatial median as a row vector
geo_median(X) = transpose(DirectionalStatistics.geometric_median([row for row in eachrow(X)]))
#apply `geo_median` function to get positive_inds and negative_inds median for each group
function geo_median(pco_space, group, positive_inds)

    spMedPos = [geo_median(pco_space[group .== g,positive_inds]) for g in unique(group)]
    spMedNeg = [geo_median(pco_space[group .== g, .!positive_inds]) for g in unique(group)]

    return (spMedPos,spMedNeg )
end
function geo_median(pco_space, group)

    return [geo_median(pco_space[group .== g,:]) for g in unique(group)]
end

function dispersion(D,group;metric = false)
    """
    Quantify the dispersion around the group spatial median for each group in D, according to the given distance metric.
    Once X is converted to a distance matrix, it is then transformed by principal coordinate analysis before medians and 
    distances are calculated. The algorithm is based on Anderson (2006). 
    
    The function accepts:
        D: A dissimiarity or distance matrix
        group: A vector containing group labels
        metric: (key word) boolean value specifying whether or not the distance used used metric
    
    This function returns a named tuple containing:
        F = Global F-Statistic 
        pairwise_F = pairwise F-statistics
        medians = spatial (aka geometric) median of each group in the transformed coordinates
        residuals = vector containing each observations distance from its group median
        means = vector containing the means distance to median for each group
        group = the original grouping vector
        levels = unique items from 'group'
    
    Anderson, M.J. (2006) Distance-based tests for homogeneity of multivariate dispersions. Biometrics 62(1), 245--253.
    """
    levels = sort(unique(group))
    
    if metric == false
        # Double center the matrix.
        A = double_center(D .^2 ) 
        vals, vecs = eigen(Hermitian(-A/2))
        pco_space = vecs * diagm(sqrt.(abs.(vals))) 
        # record indices corresponding to postitive eigenvalues
        positive_inds = real.(vals) .> 0.0 
        med = geo_median(pco_space,ones(size(pco_space)[1]),positive_inds)
        medians = geo_median(pco_space, group, positive_inds)
        # calculate distances (to median) for "positive_inds" and "negative_inds" vectors separately. in orer to 
        # subtract "negative_inds" from "positive_inds". See Anderson (2006) for details.
        dis_pos = [get_residuals(pco_space[group .== levels[i],positive_inds], medians[1][i]) for i in 1:length(levels)]
        dis_neg = [get_residuals(pco_space[group .== levels[i], .!positive_inds], medians[2][i]) for i in 1:length(levels)]
        # Where dis_neg > dis_pos, set residual = 0 by taking only the real part 
        # of the square root of a negative. 
        residuals = [(real.(sqrt.(Complex.(dis_pos[i] .- dis_neg[i]))) )   for i in 1:length(dis_pos)] 
        group_dis_pos = get_residuals(vcat(medians[1]...), med[1][1]) 
        group_dis_neg = get_residuals(vcat(medians[2]...), med[2][1])
        group_residuals = [(real.(sqrt.(Complex.(group_dis_pos[i] .- group_dis_neg[i]))) )   for i in 1:length(group_dis_pos)] 

    else
        medians = geo_median(D, group)
        med = geo_median(D,ones(size(D)[1]))
        residuals = [get_residuals(D[group .== levels[i],:], medians[i]) for i in 1:length(levels)]
        group_residuals = vec(get_residuals(vcat(medians...), med[1])) 
    end
    F = f(residuals)
    F_pairs = f_pairs(residuals)
    means = AxisArray(mean.(residuals),string.(levels))
    return Disp(F, F_pairs, medians, residuals, means,group, levels,group_residuals)
end

