using Pkg
Pkg.activate(".")
Pkg.instantiate()
using DataFrames
using Diversity
using Plots
pop = DataFrame(SpeciesA = [1, 1, 0], SpeciesB = [2, 0, 0]; SpeciesC =[3, 1, 4])
pop = [1 1 0; 2 0 0; 3 1 4];
names(pop)
pop = pop / sum(pop);
meta = Metacommunity(pop);
Z = [1.0 1 1; 1 1 1; 1 1 1]
meta_z = Metacommunity(pop, Z)

rho_z = RawAlpha(meta_z)
rho = RawAlpha(meta)

redundancies = subdiv(rho, 0)
redundancies_z = subdiv(rho_z, 0)

f(x) = x^0.25

using Diversity.Ecology
using Diversity.API
communitydata = [10 20 30 20 0; #5 species (columns) and 6 sites (rows)
                10 0 50 80 10;
                60 10 90 0 0; 
                10 10 10 10 10;
                70 70 70 70 70;
                10 0 0 90 0]
communitydata /= sum(communitydata)
bar = (Gamma(Metacommunity(communitydata)))

a = [1 2 0; 2 0 0;3 1 1]
sum(x->x>0, communitydata, dims=2)

function bubba(x)
    print(x+1)
end

function generalisedpielou(level::DiversityLevel,
    proportions::AbstractArray,
    sim::AbstractTypes)
if (level == subcommunityDiversity)
dm = ᾱ
elseif (level == metacommunityDiversity)
dm = Gamma
else
error("Can't calculate richness for $level")
end
ns = sum(x->x>0, proportions, dims=2) 
gp = level(dm(Metacommunity(proportions, sim)), 1)
gp[!,:diversity] .= log.(gp[!,:diversity])./log.(ns)
gp[!,:measure] .= "Pielou"
select!(gp, Not(:q))
return gp
end

community = [10, 20, 20]

community /= sum(community)

ecosystem = [2 2 0.; 0 2 2]'

ecosystem /= sum(ecosystem)

pielou(proportions::AbstractVecOrMat) =
    generalisedpielou(subcommunityDiversity, proportions,
                       UniqueTypes(size(proportions, 1)))
function pielou(asm::EcoBase.AbstractAssemblage)
    hassimilarity(asm) && error("function cannot run with $(typeof(gettypes(asm))) types as contains similarity")
    return pielou(occurrences(asm))
end

pielou(ecosystem)