using Revise, Oscar

includet( "Multiindices.jl" )
includet( "Misc.jl" )
includet( "AlgebraicGroups.jl" )

using Main.Misc
using Main.Representations


r = 4
p = 2

G = Representations.AlgGroupRep( r )

S = Representations.subgroup( G, ideal(G.R, det(G.A)-1))

H = Representations.representation_on_ext_power( S, p )

R, J = null_cone(H)
