
#using Oscar

include( "./Multiindices.jl" )

using Main.Multiindices

alpha = HomogMultiindex(4, 3)

for i in alpha
  print( i )
  println( power( [1,2,3,5], i ))
end
