module Multiindices

export OrderedMultiindex
export linear_index
export index_signature

abstract type Multiindex end

mutable struct OrderedMultiindex <: Multiindex 
  n::Int
  p::Int
  index::Vector{Int}
  function OrderedMultiindex(n::Int, p::Int)
    v = Vector{Int}(undef,p) 
    for i in 1:p
      v[i] = i
    end
    return new( n, p, v )
  end
end

function Base.show( io::Base.TTY, i::OrderedMultiindex ) 
  outstr= "0"
  for k in 1:length(i.index) 
    outstr = outstr * " < " * string(i.index[k]) 
  end
  outstr = outstr * " <= " * string( i.n ) * "\n"
  #Base.show( io, "n = " * string(i.n) * "; p = " * string(i.p) * "; index = " * string( i.index ) )
  Base.print( io, outstr )
end

function Base.iterate( i::OrderedMultiindex ) 
  for k in 1:i.p
    i.index[k] = k
  end
  return i, i.index
end

function Base.iterate( i::OrderedMultiindex, state ) 
  # determine the first entry that has not yet reached its 
  # limit 
  k = i.p
  while ( k>0 ) && ( i.index[k] == i.n - (i.p-k) )
    k = k-1
  end

  # whenever this statement is true, the iteration is exhausted
  if k == 0 
    return nothing
  end

  # otherwise iterate by one
  i.index[k] = i.index[k]+1
  # and adjust all the follow-up indices
  for j in k+1:i.p
    i.index[j] = i.index[j-1]+1
  end
  return i, i.index
end

# return the linear index in the enumeration pattern 
# implemented in the iterator 
function linear_index( i::OrderedMultiindex )
  N = binomial( i.n, i.p )
  for k in i.p:-1:1 
    N = N - binomial( i.n-i.index[k], i.p-k+1 )
  end
  return N
end

function ordered_multiindex( l::Int, n::Int, p::Int )
  i = OrderedMultiindex( n, p )
  if l > binomial( n, p ) 
    return nothing
  end
  for k in 1:p
    while l > binomial( n-i.index[k], p-k )
      l = l - binomial( n-i.index[k], p-k )
      i.index[k] = i.index[k]+1
    end
    for r = k+1:p
      i.index[r] = i.index[r-1]+1
    end
  end
  return i
end

function index_signature( alpha::OrderedMultiindex )
  sign = 1
  for k in alpha.p:-1:1
    if mod( k-alpha.index[k], 2 ) != 0
      sign = sign*(-1)
    end
  end
  return sign
end

end



