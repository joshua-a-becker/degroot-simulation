simDegroot = function(A, rounds, TRUTH=0, return.distribution=F) {
  weighted_deg_cent = apply(A,2,sum) %>% DescTools::Gini

  output = data.frame()
  new_guess = (A %^% rounds) %*% belief
  for(i in 0:GAME_LENGTH) {
    output_row = data.frame(
      N=vcount(g),
      deg.cent=g$deg.cent,
      truth=TRUTH,
      mean = mean(V(g)$guess),
      median = median((V(g)$guess)),
      mean_ind_err = mean(abs(V(g)$guess-TRUTH)),
      core = V(g)[1]$guess,
      sd.pool=sd(V(g)$guess),
      round = i,
      stringsAsFactors=F
    )
    output = rbind(output, output_row)
    this_guess = V(g)$guess
    V(g)$guess = V(g)$guess*(1-V(g)$alpha) + V(g)$alpha*(V(g)$guess %*% as_adj(g, sparse=F) / rep(1, vcount(g)) %*% as_adj(g, sparse=F))
  }
  output
  if(return.distribution) {
    return(list(
      distribution=this_guess,
      output=output
    ))
  } else { return(output) }
}