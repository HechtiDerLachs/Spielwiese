
using Oscar

include( "./Multiindices.jl" )

using Main.Multiindices

alpha = OrderedMultiindex(5, 3)

for i in OrderedMultiindex( 6, 4) 
  print( i )
end
