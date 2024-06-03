#######

plot.foodweb <- function(y.names,pars,file=NA,cex=1,col=1,font=1,title="",
                         lwd=1,bg=colors()[1],lcol=colors()[555],ncrit=8,lcrit=20,
                         survivals=NA,observed=NA,texts=TRUE,pointcol=FALSE,...)
  
  # plot foodweb
  #
  # last update: 16.04.2013
  #
  # call:                           
  # plot.foodweb(y.names,pars,file=NA,...)
  #
  # arguments:
  # y.names:       list of decoded states or vector with names of states
  # pars:          definition of parameters as used for rhs.benthos
  # file:          name of optional file to plot to
  # texts:         if texts=FALSE plot points else taxa names
  # pointcol:      if TRUE color of taxanames or points are varicolored (eawagfarben)
  # ncrit:         number of inverts in one line at which they are shifted up and down alternately
  # lcrit:         number of letters/characters of state names which are plottet in level 2
  #
  #
  # ...:           attributes passed to the pdf command to open the (optional)
  #                output file
# output:
# plot of the foodweb

{
  
  #analyse structure
  # ------------------
  
  # decode states:
  
  if ( !is.list(y.names) ) y.names <- decode.statevarnames(y.names)
  
  state_names <- c("SusPOM",y.names$taxa)
  
  # get levels and food of each state variable:
  
  stoich.Cons  <- get.par.stoich.web  ("Cons" ,pars)
  
  level         <-  rep(1,length(state_names))
  names(level)  <-  state_names
  foods         <-  list()
  
  # identify foods
  # loop over state variables:
  
  for ( i in 1:length(state_names) )
  {
    state_name <- state_names[i]
    
    foods[[state_name]] <- names(stoich.Cons[[state_name]])
    
  }  #loop over state variables
  
  # set levels of consumers    
  
  count = 1
  repeat 
  {
    #loop over state variables
    for (i in 1:length(foods) )
    {
      
      state_name <- names(foods[i])
      
      #update own level
      
      level[state_name] <- max( (max( level[paste(foods[[state_name]])] ) +1) ,level[state_name] )
      
      #update level of my predators
      
      predators<-NA
      for (  j in 1:length(foods))
      {
        pind <- which(foods[[j]]==state_name)
        if (sum(pind) >0) {predators <- c( predators, names(foods[j]) )}
      }
      predators    <- unique(predators[!is.na(predators)] )
      
      
      if (length(predators)>0)
      {
        for (p in 1:length(predators))
        {
          predator <- predators[p]
          level[predator] <- max( (level[state_name]+1),level[predator] )
        }
      }
    }
    
    count <- count+1
    
    if(count>10){break} 
    
  } #repeat loop
  
  
  
  # ( level )
  # ( foods  )
  
  structure           <- matrix(NA,nrow=length(state_names),ncol=4)
  colnames(structure) <- c("level","width","xpos","ypos")
  rownames(structure) <- state_names
  structure           <- as.data.frame(structure)
  structure[,"level"] <- level
  structure[,"width"] <- 1
  
  ind.topdown         <- order(structure$level,decreasing=F)
  names(ind.topdown)  <- state_names[ind.topdown]
  ind.bottomup        <- order(structure$level,decreasing=T)
  names(ind.bottomup) <- state_names[ind.bottomup]
  
  nr_levels <- NA           #nr of states at each level
  for (l in 1:max(level) )
  {
    nr_levels[l] <-  length(structure[structure[,"level"]==l,"level"])
  }
  
  oldlevel <- 0
  offset   <- 0
  for ( stat in ind.bottomup )
  {
    #stat <- ind.bottomup[2]
    
    newlevel <- structure$level[stat]
    n        <- nr_levels[newlevel]
    
    if ( newlevel != oldlevel )
    {
      offset   <- 0
      structure$xpos[stat] <- offset + 1/(n+1)
      offset   <- structure$xpos[stat]
      oldlevel <- newlevel
    } else
    {
      structure$xpos[stat] <- offset + 1/(n+1)
      offset <- structure$xpos[stat]
    }
  }
  
  structure$ypos <- structure$level /(max(level)+1)
  
  ypos <- unique(structure$ypos)
  for (i in 1:length(ypos))
  {
    yposi <- ypos[i]
    indy <- which(structure$ypos==yposi)
    if( length(indy)> ncrit)
    {
      indx <- grep(paste(sort(structure$xpos[indy]),collapse="|"),structure$xpos[indy])
      
      indh <- which(indx%%2==0)
      indl <- which(indx%%2==1)
      structure$ypos[indy][indh] <- yposi*1.1
      structure$ypos[indy][indl] <- yposi*0.9
    }
  }
  
  
  
  
  # open pdf file:
  
  if ( !is.na(file) ) pdf(file,...)
  
  # initialize plot:
  
  # par.def <- par(no.readonly=TRUE)
  par(mar=c(0,0,0,0),bg=bg)
  plot(numeric(0),numeric(0),xlim=c(1,0),ylim=c(0,1),xlab="",ylab="",axes=F)# ,...
  
  
  
  
  #check survivals
  if ( length(survivals)==1 ) 
  {
    no.survivals <- TRUE
    survivals <- rep(NA,length(state_names))
    names(survivals) <- state_names 
  } else no.survivals <- FALSE
  
  #check observed
  if ( length(observed)==1 ) 
  {
    no.observed <- TRUE
    observed <- rep(NA,length(state_names))        
    names(observed) <- state_names 
  } else no.observed <- FALSE
  
  
  # draw food connections
  
  for ( s in 1:length(foods) )
  {
    food <- foods[[s]]
    cons <- names(foods[s])
    if (!is.na(food)[1])
    {
      for (f in 1:length(food))
      {
        from <- cons
        to   <- food[f]
        x <- numeric(2)
        x[1]<- structure[paste(from),"xpos"]
        x[2]<- structure[paste(to),"xpos"]
        y   <- numeric(2)
        y[1]<- structure[paste(from),"ypos"]
        y[2]<- structure[paste(to),"ypos"]
        
        indcons <-  which(from==names(survivals))
        indconsm <- which(from==names(observed))
        
        indfood <- which(to==names(survivals))
        indfoodm <- which(to==names(observed))
        
        if( sum(indfoodm)==0 ) warning (to, " not in observed")
        if( sum(indfood)==0 )  warning (to, " not in survivals")
        
        cond.1 <- (survivals[indcons]=="extinct"|
                     survivals[indfood]=="extinct"|
                     observed[indconsm]=="never"  |
                     observed[indfoodm]=="never")  
        
        cond.2 <- ( observed[indconsm]=="sometimes" |
                      observed[indfoodm]=="sometimes" |
                      observed[indconsm]=="NA"        | 
                      observed[indfoodm]=="NA" )
        
        cond.3 <- !no.observed # length(observed)>1 & !is.na(observed)
        
        
        clf <- lcol
        
        if( !is.na(cond.1) & cond.1 ) { clf=NA } else {
          if( !is.na(cond.2) & cond.2 & cond.3) {clf=gray(0.5) } }
        
        lines(x,y,col=clf,lwd=lwd)
        
      }
    }
  }
  
  
  
  # plot text or points
  
  for ( l in 1 : max(level) )
  {
    
    
    label = rownames(structure)[structure$level == l]
    
    
    
    cl <- rep(1,length(label))
    
    for (m in 1:length(label))
    {
      
      if (!pointcol)
      {
        #cl[m] <- gray(0.5)
        
        if(!no.survivals)
        {
          inds <- which(label[m]==names(survivals))
          if( !is.na(survivals[inds]) )
          {
            if(survivals[inds]=="extinct")  { cl[m] <- gray(0.6) }
          }
        }
        
        if(!no.observed)
        {
          indm <- which(label[m]==names(observed))
          if( !is.na( observed[indm]) )
          {
            if(observed[indm]=="notobserved" | observed[indm]=="never") 
            {cl[m] <- gray(0.7)} else {
              if(observed[indm]== "sometimes" ) 
              {cl[m] <- gray(0.45)} else {
                if (observed[indm]== "always" | observed[indm]=="observed") 
                {cl[m] <- 1} else {
                  if (observed[indm]== "NA")
                  { cl[m] <- gray(0.5) }
                  
                  
                }
              } 
            }
          }
        }
        
      } else 
      {
        if (pointcol)
        {  
          
          if(!no.survivals)
          {      
            inds <- which(label[m]==names(survivals))
            if( !is.na(survivals[inds]) )
            {
              ifelse(survivals[inds]=="extinct", 
                     cl[m] <- rgb(240,123,0,maxColorValue = 255), 
                     cl[m] <- rgb(0,173,221,maxColorValue = 255))
            } 
          }
          
          if(!no.observed)
          {
            indm <- which(label[m]==names(observed))
            if( !is.na( observed[indm] ) )
            {
              if(observed[indm]=="notobserved" | observed[indm]=="never") 
              {cl[m] <- rgb(240,123,0,maxColorValue = 255)} else {
                if(observed[indm]== "sometimes" ) 
                {cl[m] <- rgb(188,208,45,maxColorValue = 255)} else {
                  if (observed[indm]== "always" | observed[indm]=="observed") 
                  {cl[m] <- rgb(0,173,221,maxColorValue = 255)} else {
                    if (observed[indm]=="NA") 
                    { cl[m] <-gray(0.5)}
                  }
                }
              } 
            }
          } 
          
        } 
      }
    } #end for m
    
    if(texts)
    {
      if(l==2) {label = substr(rownames(structure)[structure$level == l],1,lcrit)}
      
      if( length(label)!=length(unique(label))  ) {warning("lcrit: ",lcrit," to small,labels not unique")}
      
      text(x  = structure$xpos[structure$level == l],
           y  = structure$ypos[structure$level == l],
           labels = label,
           col    = cl,
           font   = font,
           cex    = cex
      )
      
    }  else
    {
      points(x   = structure$xpos[structure$level == l],
             y   = structure$ypos[structure$level == l],
             pch = 19,
             col = cl,
             cex = cex*3
      )
    }
    
  } #end for l
  
  mtext(title,3,line=-2,cex=cex,font=2)
  
  # reset plot parameters:
  #par(par.def)
  
  # close pdf file:
  
  if ( !is.na(file) ) dev.off()
  
  # return(list(structure=structure[,c(1,4)])) #foods=foods,products=products
  
} # end function