using Revise, Oscar

includet( "Multiindices.jl" )
includet( "Misc.jl" )
includet( "AlgebraicGroups.jl" )

using Main.Misc
using Main.Representations


r = 2
p = 2
d = 4

G = Representations.AlgGroupRep( r )

S = Representations.subgroup( G, ideal(G.R, det(G.A)-1))

H = Representations.representation_on_sym_power( S, d )

R, J = null_cone(H)
