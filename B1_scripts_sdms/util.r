


clipper <- function(vec, min_value=-Inf, max_value=Inf){
  
  # This function clips the values of received vector vec to lie in [min_value, 
  # max_value].
  #
  # arguments:
  #   - vec:                    vector of values to clip.
  #   - min_value:              lower bound for clipping vec.
  #   - max_value:              upper bound for clipping vec
  #
  # returns:
  #   - the clipped version of vec
  
  return(pmax(rep(min_value, length(vec)), pmin(vec, rep(max_value, length(vec)))))
}

