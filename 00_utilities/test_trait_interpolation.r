

# test approx function in get_par_trait_cond.r

s.x <- list.df.preferences$tempmaxtolval$Values
s.y <- list.df.preferences$tempmaxtolval$Baetisalpinus

xout <- seq(1, 23, by = 1)

# install.packages("pracma")

linear.interp <- approx(s.x,s.y,xout=xout,rule=2)$y
polynomial.interp <- pracma::pchip(s.x, s.y, xout)
 
# f.sap.i.untrans <- pracma::pchip(xi = s.x,yi = s.y, x = xout)

plot(s.x, s.y, col='red', pch=13)
lines(xout, linear.interp, col='green', lwd=1)
lines(xout, polynomial.interp, col='blue', lwd=1)


# test interpolation before get_par_trait_cond.r

y.names.par.catch <- list.variables.par.catch$Aare

y.names <- y.names.par.catch$y.names
par.fixc <- y.names.par.catch$par.fixc
reaches.orig <- y.names.par.catch$y.names$reaches
inp <- NA

if ( !is.list(y.names) ) y.names <- decode.statevarnames(y.names)

sys.def <- streambugs.get.sys.def(y.names=y.names,par=par.fixc,inp=NA)

for (i in 1:length(sys.def)) {
  print(head(sys.def[[i]]))
}



