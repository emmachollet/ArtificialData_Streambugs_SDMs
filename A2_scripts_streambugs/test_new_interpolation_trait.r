

# need the following
par.env.global.orig <- par.env.global
par.env.global.update <- par.env.global
par.invtraits.orig <- par.invtraits
par.invtraits.update <- par.invtraits
# Invertebrates
no.class.new <- 20

# retrieve name environmental factors
# names.par.env.global <- names(par.env.global)
# names.env.fact <- unique(sub("\\_.*", "", names.par.env.global))

acronym.env.fact <- c("tempmax", "current", "sapro")
                      # , "orgmicropoll")

for (acro in acronym.env.fact) {
  # acro <- acronym.env.fact[1]
  # retrieve original classes of environmental factor
  ind.par.env <- which(grepl(acro, names(par.env.global.update)))
  class.env.orig <- par.env.global.update[ind.par.env]
  name.env.par <- unique(sub("\\_.*", "",names(class.env.orig)))
  no.class.orig <- length(class.env.orig)
  
  # create new (more) classes
  val.class.env.new <- round(seq(min(class.env.orig), max(class.env.orig), length.out = no.class.new), digits = 2)
  class.env.new <- val.class.env.new
  names.class.new <- paste("class", 1:no.class.new, sep = "")
  names(class.env.new) <- paste(name.env.par, names.class.new, sep = "_")
  
  for(taxon in Invertebrates){
    # taxon <- Invertebrates[5]  
    # retrieve original preference trait extracted from database
    ind.tax.env.score <- which(grepl(taxon, names(par.invtraits.update)) & grepl(acro, names(par.invtraits.update)))
    scores.tax.env <- par.invtraits.update[ind.tax.env.score]
    name.score <- str_split(names(class.env.orig), "_")[[1]][1]

    # create "new" preference trait from a polynomial interpolation of the original classes and traits
    linear.interp <- approx(class.env.orig, scores.tax.env,xout=class.env.orig,rule=2)$y
    polynomial.interp <- pracma::pchip(class.env.orig, scores.tax.env, class.env.new)
    
    # plot
    plot(class.env.orig, scores.tax.env, col='red', pch=13, main = paste(taxon, name.score, sep = "_"))
    points(class.env.new, polynomial.interp, col='blue', pch=5)
    lines(class.env.orig, linear.interp, col='green', lwd=1)
    lines(class.env.new, polynomial.interp, col='purple', lwd=1)

    scores.tax.new <- round(polynomial.interp, digits = 2)
    names(scores.tax.new) <- paste(paste0(taxon, "_", name.score, "_"), names.class.new, sep = "_")
    
    # remove old scores and append new one
    par.invtraits.update <- par.invtraits.update[-ind.tax.env.score]
    par.invtraits.update <- append(par.invtraits.update, scores.tax.new)
    
  }
  # remove old classes and append new one
  par.env.global.update <- par.env.global.update[-ind.par.env]
  par.env.global.update <- append(par.env.global.update, class.env.new)
  
}

