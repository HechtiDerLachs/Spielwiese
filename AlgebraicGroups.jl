module Representations

using Oscar
using Main.Multiindices
using Main.Misc

import AbstractAlgebra.Field
import AbstractAlgebra.Generic.MatSpaceElem

export AlgGroup, subgroup, representation_on_ext_power, null_cone
export LinearRep

export apply_diff_op


#abstract type AlgebraicGroupElem end

#abstract type Field = Union{ Nemo.FlintRationalField } end

mutable struct AlgGroup
  R::MPolyRing 		# The underlying polynomial ring
  I::MPolyIdeal		# The ideal describing the closed subgroup
  Z::MatSpaceElem	# The matrix of coordinates in GL_n(k). 
  			# Note that the polynomial ring might have additional 
			# variables.

  # default constructor for GL(n;QQ)
  function AlgGroup( n::Int )
    #R_t,t = PolynomialRing(QQ,"t")
    R, Z,t = PolynomialRing(QQ, "a" => (1:n, 1:n), "t" => (1:1))
    Z = matrix(Z)
    t = t[1]
    I = ideal(R, [1-det(Z)*t] )
    return new( R, I, Z )
  end

  function AlgGroup( R::MPolyRing, I::MPolyIdeal, Z::MatSpaceElem )
    return new( R,I,Z )
  end

end

function subgroup( G::AlgGroup, I::MPolyIdeal )
  return AlgGroup( G.R, G.I + I, G.Z )
end

mutable struct AlgGroupElem
  parent::AlgGroup
  coord::Matrix{fmpq}
  function AlgGroupElem( G::AlgGroup, A::Matrix{fmpq} ) 
    # sanity check 
    if !(nrows(A) == nrows( G.Z ) && ncols(A) == ncols( G.Z ))
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

mutable struct LinearRep
  G::AlgGroup
  A::MatSpaceElem

  function LinearRep( G::AlgGroup, A::MatSpaceElem )
    # TODO: Add plausibility checks.
    return new( G, A )
  end

  # The identity representation
  function LinearRep( G )
    return new( G, G.Z )
  end
end

function representation_on_sym_power( rep::LinearRep, d::Int )
  # return the representation induced on the d-th symmetric power 
  # of the vector space V on which the representation rep is 
  # defined.
  n = ncols( rep.A )	# the rank of the original representation
  N = binomial( n+d-1, d )	# the rank of the induced representation
  # A polynomial ring containing variables for the coordinates of 
  # the original representation:
  S, y = PolynomialRing( rep.G.R, [ "y$k" for k in (1:n) ])
  # An empty matrix for the induced representation
  M = MatrixSpace( rep.G.R, N, N )
  B = zero(M)
  # The matrix of the original representation, but promoted to the 
  # new ring S
  A = matrix( S.(rep.A) )
  # images of the variables under the action of the group
  z = A*matrix(y)
  # Populate the matrix B for the induced representation...
  for i in HomogMultiindex( n, d )
    # ...by going through the multiindices of the monomials 
    # in the given degree and taking the respective powers...
    f = power( [z[k,1] for k in (1:n)], i )
    for j in HomogMultiindex( n, d )
      # ...and extracting the coefficients again.
      B[ linear_index(j), linear_index(i) ] = coeff( f, power( y, j ))
    end
  end
  return LinearRep( rep.G, B )
end
      

function representation_on_ext_power( G::AlgGroup, p::Int )
  # return the representation on the p-th exterior power of the vector 
  # space V on which G is defined. 

  n = ncols( G.Z )
  N = binomial( n, p )
  M = MatrixSpace( G.R, N, N )
  B = M(0)
  for i in OrderedMultiindex( n, p )
    for j in OrderedMultiindex( n, p )
      B[ linear_index( i ), linear_index(j) ] = index_signature(i)*index_signature(j)*det( G.Z[i.index,j.index] )
    end
  end

  return AlgGroup( G.R, G.I, B )
end

function null_cone( rep::LinearRep )
  # Assemble the ring for the homomorphism needed for the computation
  n = ncols( rep.A )
  M, phi, v = add_variables( rep.G.R, vcat( [ "x$i" for i in (1:n) ], [ "y$i" for i in (1:n) ]))
  B = matrix(phi.(rep.A))
  h = matrix(v[n+1:2*n]) - B*matrix(v[1:n])
  I = ideal(M,vec(collect(h))) + ideal( M, [ phi(f) for f in gens(rep.G.I) ] )
  b = eliminate( I, gens(M)[1:length(gens(rep.G.R))] )
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
    
# Apply the Omega process for the representation rep to 
# the polynomial f. Here f is a polynomial in the coordinate 
# functions of the vector space V on which the representation 
# is defined. 
function omega_process( rep::LinearRep, f::MPolyElem ) 

end
  
end
