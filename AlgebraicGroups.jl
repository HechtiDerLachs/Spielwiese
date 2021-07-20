module Representations

using Oscar
using Main.Multiindices
using Main.Misc

import AbstractAlgebra.Field
import AbstractAlgebra.Generic.MatSpaceElem

export AlgGroupRep, subgroup, representation_on_ext_power, null_cone

export apply_diff_op


#abstract type AlgebraicGroupElem end

#abstract type Field = Union{ Nemo.FlintRationalField } end

mutable struct AlgGroupRep
  R::MPolyRing 		# The underlying polynomial ring
  I::MPolyIdeal		# The ideal describing the closed subgroup
  A::MatSpaceElem	# The matrix for the representation

  # default constructor for GL(n;QQ)
  function AlgGroupRep( n::Int )
    #R_t,t = PolynomialRing(QQ,"t")
    R, A,t = PolynomialRing(QQ, "a" => (1:n, 1:n), "t" => (1:1))
    A = matrix(A)
    t = t[1]
    I = ideal(R, [1-det(A)*t] )
    return new( R, I, A )
  end

  function AlgGroupRep( R::MPolyRing, I::MPolyIdeal, A::MatSpaceElem )
    return new( R,I,A )
  end

end

function subgroup( G::AlgGroupRep, I::MPolyIdeal )
  return AlgGroupRep( G.R, G.I + I, G.A )
end

mutable struct AlgGroupElem
  parent::AlgGroupRep
  coord::Matrix{fmpq}
  function AlgGroupElem( G::AlgGroupRep, A::Matrix{fmpq} ) 
    # sanity check 
    if !(nrows(A) == nrows( G.A ) && ncols(A) == ncols( G.A ))
      println( "matrix sizes are incompatible" )
      return
    end
    return new(G,A)
  end
end


function mult( X::AlgGroupElem, Y::AlgGroupElem ) 
  if ncols(X) != nrows(Y)
    error( "Incompatible sizes of matrices" )
    return
  end
  if X.parent != Y.parent 
    error( "Elements do not live in the same group" )
    return
  end
  return Z = AlgGroupElem( X.parent, X.A*Y.A )
end

function representation_on_sym_power( G::AlgGroupRep, d::Int )
  # return the representation induced on the d-th symmetric power 
  # of the vector space V on which G is defined.
  @show "call to representation on symmetric power" 
  n = ncols( G.A )
  N = binomial( n+d-1, d )
  @show N
  S, y = PolynomialRing( G.R, [ "y$k" for k in (1:n) ])
  M = MatrixSpace( G.R, N, N )
  B = zero(M)
  A = matrix( S.(G.A) )
  @show A
  z = A*matrix(y)
  @show z
  for i in HomogMultiindex( n, d )
    @show [z[k,1] for k in (1:n)]
    f = power( [z[k,1] for k in (1:n)], i )
    for j in HomogMultiindex( n, d )
      B[ linear_index(j), linear_index(i) ] = coeff( f, power( y, j ))
    end
  end
  return AlgGroupRep( G.R, G.I, B )
end
      

function representation_on_ext_power( G::AlgGroupRep, p::Int )
  # return the representation on the p-th exterior power of the vector 
  # space V on which G is defined. 

  n = ncols( G.A )
  N = binomial( n, p )
  M = MatrixSpace( G.R, N, N )
  B = M(0)
  for i in OrderedMultiindex( n, p )
    for j in OrderedMultiindex( n, p )
      B[ linear_index( i ), linear_index(j) ] = index_signature(i)*index_signature(j)*det( G.A[i.index,j.index] )
    end
  end

  return AlgGroupRep( G.R, G.I, B )
end

function null_cone( H::AlgGroupRep )
  # Assemble the ring for the homomorphism needed for the computation
  n = ncols( H.A )
  M, phi, v = add_variables( H.R, vcat( [ "x$i" for i in (1:n) ], [ "y$i" for i in (1:n) ]))
  A = matrix(phi.(H.A))
  h = matrix(v[n+1:2*n]) - A*matrix(v[1:n])
  I = ideal(M,vec(collect(h))) + ideal( M, [ phi(f) for f in gens(H.I) ] )
  b = eliminate( I, gens(M)[1:length(gens(H.R))] )
  R, x = PolynomialRing( base_ring(M), [ "x$i" for i in (1:n) ] ) 
  psi = AlgebraHomomorphism( M, R, vcat( [ zero(R) for i in (1:length(gens(M))-n)], x ))
  J = [ psi(f) for f in gens(b) ]
  return R, J
end
  
# Interpret f as a differential operator by substitution of 
# the variables by their respective partial derivatives and 
# apply it to the polynomial g.
function apply_diff_op( f::MPolyElem, g::MPolyElem )
  result = 0
  c = vec( collect( coefficients( f )))
  m = vec( collect( monomials( f )))
  e = vec( collect( exponent_vectors( f )))
  N = length( c )
  x = gens( parent( f ))
  n = length( x )
  for i in (1:N) 
    tmp = g
    for k in (1:n)
      a = e[i][k]
      for l in (1:a)
        tmp = derivative( tmp, x[k])
      end
    end
    result = result + c[i]*tmp
  end
  return result
end
    
end
