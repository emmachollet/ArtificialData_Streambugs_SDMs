#test rhs:
  
# file to run the rhs in R manually, these preprocessing steps are needed then.

  
times <- tout
par <- par.fix
#inp <- NA

file.def=NA
file.res=NA
file.add=NA

return.res.add = FALSE
tout.add=NA


if ( !is.list(y.names) ) y.names <- decode.statevarnames(y.names)

# extract structured system definition from parameter vector and input list:
# --------------------------------------------------------------------------

sys.def <- streambugs.get.sys.def(y.names=y.names,par=par,inp=inp)

#file.def <- paste("output/Sys.def.file_",name.run,"3.dat",sep="")

if ( !is.na(file.def) ) streambugs.write.sys.def(sys.def=sys.def,file=file.def)

# update (time-dependent) parameters to initial time:
# ---------------------------------------------------

sys.def <- streambugs.updatepar(sys.def,times[1])

# compile initial state:
# ----------------------

w     <- sys.def$par.envcond.reach$parvals[,"w"]
fA    <- sys.def$par.envcond.habitat$parvals[,"fA"]
D.ini <- sys.def$par.initcond$parvals[,"Dini"]
y.ini <- D.ini * w * fA
names(y.ini) <- y.names$y.names

y     = y.ini
parms = sys.def
t=1

add.out=TRUE
