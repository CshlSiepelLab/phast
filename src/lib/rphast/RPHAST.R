#Alexandra Denby
#RPHAST.R
#Alpha version
#5/10/08

phast.init <- function(){
#   dyn.load("/home/ajd45/phast/lib/librphast.so")
   .C("init",PACKAGE="librphast")
   return()
}

#creates a new alignment object
msa.makeobj <- function(){
  a=list()
  class(a)="msa"
  return(a)
}

#spec is a character vector containing species names
msa.new <- function(nseqs, spec=rep("",nseqs),seqs=rep("",nseqs)){

   return=.C("rph_msa_new",nseqs=as.integer(nseqs),spec=as.character(spec), seqs=as.character(seqs), address=as.numeric(0), numberSpecies=as.integer(0), length=as.integer(0), alphabet=as.character(""), error=as.integer(0), errstr=as.character(""), PACKAGE="librphast")

   if (return$error!=0){
      print(return$errstr)
      return()
   }

   msa=msa.makeobj()
   attr(msa,"address")=return$address
   attr(msa,"format")="FASTA"
   attr(msa,"numberSpecies")=return$numberSpecies
   attr(msa,"species")=readspecies(msa)
   attr(msa,"length")=return$length
   attr(msa,"alphabet")=return$alphabet
   attr(msa,"loaded")=1
   return(msa)

}

#read an alignment file, return the alignment object associated with it
msa.read <- function(fname, format="FASTA", ss=FALSE, gff=0, do4d=FALSE){

   if(class(gff)=="gff"){
	gffad=attr(gff,"address")
   }else{
	gffad=0
   }

   return=.C("rph_msa_read",fname=as.character(fname), format=as.character(format), address=as.numeric(0), gffAddress=as.numeric(gffad), fourD=as.integer(do4d), SS=as.integer(ss), numberSpecies=as.integer(0),length=as.integer(0),alphabet=as.character(""),error=as.integer(0), errstr=as.character(""),PACKAGE="librphast")

   if (return$error!=0){
      print(return$errstr)
      return()
   }

   msa=msa.makeobj()
   attr(msa,"address")=return$address
   attr(msa,"format")=format
   attr(msa,"numberSpecies")=return$numberSpecies
   attr(msa,"length")=return$length
   attr(msa,"alphabet")=return$alphabet
   attr(msa,"species")=readspecies(msa)
   attr(msa,"loaded")=1
   return(msa)
}

readspecies <- function(msa){
   return=.C("readSpecies", address=as.numeric(attr(msa, "address")), numberSpecies=as.integer(attr(msa, "numberSpecies")), species=as.character(rep("",attr(msa, "numberSpecies"))), error=as.integer(0), errstr=as.character(""),PACKAGE="librphast")
   if (return$error!=0){
      print(return$errstr)
      return()
   }
   return(return$species)
}


#writes an alignment to a file in specified format
#default is format which was read in as
msa.write <- function(fname, msa, format=attr(msa,"format")){
   if(class(msa)!="msa"){
	cat("Error: align must be of class \"msa\"\n")
        return()
   }
   if(attr(msa,"loaded")==0){
	cat("Error: alignment was freed from memory\nPlease reload\n")
	return()
   }
   return=.C("rph_msa_print",fname=as.character(fname), format=as.character(format), address=as.numeric(attr(msa,"address")), error=as.integer(0), errstr=as.character(""),PACKAGE="librphast")
   if (return$error!=0){
      print(return$errstr)
      return()
   }

}

#free an alignment from memory
#TODO: get this to also remove from the R side
msa.free <- function(msa){
   
   if(class(msa)!="msa"){
	cat("Error: align must be of class \"msa\"\n")
        return()
   }

   if(attr(msa,"loaded")==0){
	cat("Error: alignment was freed from memory\nPlease reload\n")
	return()
   }

   return=.C("rph_msa_free",address=as.numeric(attr(msa,"address")), error=as.integer(0), errstr=as.character(""),PACKAGE="librphast")

   if (return$error!=0){
      print(return$errstr)
      return()
   }

   attr(msa,"loaded")=0
   return(msa)
}

#print function, alignments
print.msa<- function(align){
   if(attr(align,"loaded")==0){
	cat("Alignment unloaded\n")
	return()
   }
   cat("\nMultiple Sequence Alignment\n")
   cat(attr(align,"numberSpecies"),"species,",attr(align,"length"),"bases\n\n")
}

#summary function, alignments
summary.msa<- function(align){
   if(attr(align,"loaded")==0){
	cat("Error: alignment was freed from memory\nPlease reload\n")
	return()
   }
   print(align)
   cat("Species:\n")
   print(attr(align,"species"))
   cat("\n")
}

msa.append <- function(base, new){

   if(class(base)!="msa"){
	cat("Error: base must be of class \"msa\"\n")
        return()
   }
   if(class(new)!="msa"){
	cat("Error: base must be of class \"msa\"\n")
        return()
   }
   if(attr(base,"address")==attr(new,"address")){
	cat("Error: msa.append does not support duplication\n")
	return()
   }

   return=.C("rph_msa_concatenate",baseAddress=as.numeric(attr(base,"address")), newAddress=as.numeric(attr(new,"address")),error=as.integer(0), errstr=as.character(""),PACKAGE="librphast")

   if (return$error!=0){
      print(return$errstr)
      return()
   }

   attr(base,"length")=attr(base,"length")+attr(new,"length")
   return(base)
}


####################TREE MODELS######################


tm.new <- function(){
  a=list()
  class(a)="tm"
  return(a)
}


#reads a tree model from a file
tm.read <- function(fname){
   return=.C("rph_tm_read",fname=as.character(fname), treeAddress=as.numeric(0), modAddress=as.numeric(0), error=as.integer(0), errstr=as.character(""), PACKAGE="librphast")

   if (return$error!=0){
      print(return$errstr)
      return()
   }

   tm=tm.new()
   attr(tm, "address")=return$modAddress
   attr(tm, "tree")=return$treeAddress
   return(tm)
}

#scales a tree model's tree based on a specified constant
tm.scale <- function(treemod, scale){
   if (class(treemod)!="tm"){
      cat("Error: treemod must be of class \"tm\"\n")
      return()
   }
   return=.C("rph_tr_scale",modAddress=as.numeric(attr(treemod, "address")), treeAddress=as.numeric(0), scale=as.numeric(scale), error=as.integer(0), errstr=as.character(""),PACKAGE="librphast")
   
   if (return$error!=0){
      print(return$errstr)
      return()
   }

   tm=tm.new()
   attr(tm, "address")=return$modAddress
   attr(tm, "tree")=return$treeAddress
   return(tm)
}

#Writes a tree model to a file
tm.write <- function(fname, treemod){
   if (class(treemod)!="tm"){
      cat("Error: treemod must be of class \"msa\"\n")
      return()
   }
   return=.C("rph_tm_print",fname=as.character(fname), address=as.numeric(attr(treemod, "address")), error=as.integer(0), errstr=as.character(""), PACKAGE="librphast")

   if (return$error!=0){
      print(return$errstr)
      return()
   }

}

#Frees a tree model
#TODO: free from R memory as well
tm.free <- function(tm){
   if(class(tm)!="tm"){
	cat("Error: object must be of class \"tm\"\n")
	return()
   }
   return=.C("rph_tm_free",address=as.numeric(attr(tm,"address")), error=as.integer(0), errstr=as.character(""),PACKAGE="librphast")

   if (return$error!=0){
      print(return$errstr)
      return()
   }

}


#print function, tree model
print.tm <- function(tm){
   cat("Tree model object\n")
}

summary.tm <- function(tm){
   print(tm)
}


#Does not yet support parameters scale/subtree
tm.fit <- function(msa, substmod="REV", tree="", treemod=0, precision="HIGH", use_em=FALSE){
   
   if(treemod!=0 && class(treemod)!="tm"){
	cat("Error: treemod must be of class \"tm\"\n")
	return()
   }

   if(class(msa)!="msa"){
	cat("Error: msa must be of class \"msa\"\n")
	return()
   }

   if(treemod==0 && tree==""){
	cat("Error: must enter either an initial model or a tree to build from\n");
	cat("Please use parameter tree or treemod\n");
	return();
   }
   
   if(precision=="HIGH"){
	prec=3
   }
   else if (precision=="MED"){
	prec=2
   }
   else if (precision=="LOW"){
	prec=1
   }

   return=.C("preprocessMSA",msaAddress=as.numeric(attr(msa,"address")), error=as.integer(0), errstr=as.character(""), PACKAGE="librphast")
   if (return$error!=0){
      print(return$errstr)
      return()
   }

   return=.C("rph_tm_fit",modAddress=as.numeric(attr(treemod, "address")), msaAddress=as.numeric(attr(msa,"address")), substModel=as.character(substmod), treestr=as.character(tree), prec=as.integer(prec), em=as.integer(use_em), fittedAddress=as.numeric(0),error=as.integer(0), errstr=as.character(""), PACKAGE="librphast")   

   if (return$error!=0){
      print(return$errstr)
      return()
   }

   tm=tm.new()
   attr(tm, "address")=return$fittedAddress
   return(tm)
}

tm.computeLogLikelihood <- function(treemod, msa, colScores=FALSE){

   if(class(treemod)!="tm"){
	cat("Error: treemod must be of class \"tm\"\n")
	return()
   }

   if(class(msa)!="msa"){
	cat("Error: msa must be of class \"msa\"\n")
	return()
   }

   if(colScores){
	cs=rep(0,attr(msa,"length"))
   }else{
	cs=-1
   }	

   return=.C("rph_tl_compute_log_likelihood", modAddress=as.numeric(attr(treemod, "address")), msaAddress=as.numeric(attr(msa, "address")), columnScores=as.numeric(cs), likelihood=as.numeric(0), error=as.integer(0), errstr=as.character(""),PACKAGE="librphast")

   if (return$error!=0){
      print(return$errstr)
      return()
   }

   if (colScores){
        l=list()
	l$colScores=return$columnScores
	l$likelihood=return$likelihood
	return(l)
   }else{
	return(return$likelihood)
   }

}



###################GFF####################
gff.new <- function(){
  a=list()
  class(a)="gff"
  return(a)
}

gff.read <- function(fname){

   return=.C("rph_gff_read", fname=as.character(fname), gffAddress=as.numeric(0), error=as.integer(0), errstr=as.character(""), PACKAGE="librphast")

   if (return$error!=0){
      print(return$errstr)
      return()
   }

   gff=gff.new()
   attr(gff, "address")=return$gffAddress
   
   return(gff)
}

gff.writeSet <- function(gff, fname){

   if(class(gff)!="gff"){
	print("Error: gff must be of class \"gff\"\n")
	return()
   }

   return=.C("rph_gff_print_set",address=as.numeric(attr(gff, "address")), fname=as.character(fname), error=as.integer(0), errstr=as.character(""), PACKAGE="librphast")

   if (return$error!=0){
      print(return$errstr)
      return()
   }
   
}

gff.writeFeatures <- function(gff, fname){
   if(class(gff)!="gff"){
	print("Error: gff must be of class \"gff\"\n")
	return()
   }

   return=.C("rph_gff_print_features",address=as.numeric(attr(gff, "address")), fname=as.character(fname), error=as.integer(0), errstr=as.character(""), PACKAGE="librphast")

   if (return$error!=0){
      print(return$errstr)
      return()
   }

}

gff.free <- function(gff){

  if(class(gff)!="gff"){
	print("Error: gff must be of class \"gff\"\n")
	return()
   }

  return=.C("rph_gff_free", address=as.numeric(attr(gff,"address")), error=as.integer(0), errstr=as.character(""))

   if (return$error!=0){
      print(return$errstr)
      return()
   }

}

print.gff <- function(gff){
   cat("GFF Object\n")
}

summary.gff <- function(gff){
   print(gff)
}


###########FROM PHYLOBOOT#################

tm.generateMSA<-function(treemod, nsites){
   if (class(treemod)!="tm"){
      cat("Error: treemod must be of class \"tm\"\n")
      return()
   }
   return=.C("rph_tm_generate_msa",modAddress=as.numeric(attr(treemod, "address")), numSites=as.integer(nsites), msaAddress=as.numeric(0), numberSpecies=as.integer(0), length=as.integer(0), alphabet="", error=as.integer(0), errstr=as.character(""), PACKAGE="librphast")

   if (return$error!=0){
      print(return$errstr)
      return()
   }
   
   msa=msa.makeobj()
   attr(msa,"address")=return$msaAddress
   attr(msa,"format")="FASTA"
   attr(msa,"numberSpecies")=return$numberSpecies
   attr(msa,"length")=return$length
   attr(msa,"alphabet")=return$alphabet
   attr(msa,"species")=readspecies(msa)
   attr(msa,"loaded")=1
   return(msa)
}

################################PHYLOP###############################
phylop.baseScores <- function(msa, treemod, method="SPH", mode="CON", fitModel=FALSE, epsilon=0, stats=FALSE){
 
   if(class(msa)!="msa"){
	cat("Error: align must be of class \"msa\"\n")
        return()
   }

   if (class(treemod)!="tm"){
      cat("Error: treemod must be of class \"tm\"\n")
      return()
   }

   if (mode=="CON"){
	modeInt=1
   }else if (mode=="NNEUT"){
	modeInt=2
   }else if (mode=="ACC"){
	modeInt=3
   }else{
	cat("ERROR: unrecognized mode\n")
	cat("Please use one of: CON, NNEUT, ACC\n")
	return()
   }

   return=.C("preprocessMSA",msaAddress=as.numeric(attr(msa,"address")), error=as.integer(0), errstr=as.character(""), PACKAGE="librphast")
   if (return$error!=0){
      print(return$errstr)
      return()
   }

   if (method != "SPH" && (epsilon !=0 || fitModel)){
	cat("WARNING: fitModel/epsilon only supported in method SPH\n")
	cat("         discarding these values\n")
   }

   if (method=="SPH"){
	if (epsilon==0){
		epsilon=1e-6
		cat("Using default epsilon value 1e-6\n")	  
	}
	return=.C("phyloP_SPH",msaAddress=as.numeric(attr(msa,"address")), tmAddress=as.numeric(attr(treemod, "address")), fm=as.integer(fitModel), eps=as.numeric(epsilon), stats=as.integer(stats), modeNum=as.integer(modeInt), pVals=as.numeric(rep(0,attr(msa,"length"))), postMeans=as.numeric(rep(0,attr(msa,"length"))), postVars=as.numeric(rep(0,attr(msa,"length"))), error=as.integer(0), errstr=as.character(""),PACKAGE="librphast")
	
   	if (return$error!=0){
      		print(return$errstr)
      		return()
   	}
	
	l=list()
	l$msa=msa
	l$treemod=treemod
	l$method="SPH"
	l$stats=stats
	l$pvals=return$pVals
	if (stats){
		l$postMeans=return$postMeans
		l$postVars=return$postVars
	}
	class(l)="bScore"
	return(l)
   }	

   if (method=="LRT"){
	return=.C("phyloP_LRT",msaAddress=as.numeric(attr(msa,"address")), tmAddress=as.numeric(attr(treemod, "address")), stats=as.integer(stats), modeNum=as.integer(modeInt), pVals=as.numeric(rep(0,attr(msa,"length"))), llrs=as.numeric(rep(0,attr(msa,"length"))), scales=as.numeric(rep(0,attr(msa,"length"))), error=as.integer(0), errstr=as.character(""),PACKAGE="librphast")

   	if (return$error!=0){
      		print(return$errstr)
      		return()
   	}

	l=list()
	class(l)="bScore"
	l$msa=msa
	l$tm=treemod
	l$method="LRT"
	l$stats=stats
	l$pvals=return$pVals
	if (stats){
	   l$llrs=return$llrs
	   l$sscales=l$scales
	}
	return(l)
   }
   if (method=="SCORE"){
	return=.C("phyloP_SCORE",msaAddress=as.numeric(attr(msa,"address")), tmAddress=as.numeric(attr(treemod, "address")), stats=as.integer(stats), modeNum=as.integer(modeInt), pVals=as.numeric(rep(0,attr(msa,"length"))), teststats=as.numeric(rep(0,attr(msa,"length"))), derivs=as.numeric(rep(0,attr(msa,"length"))), error=as.integer(0), errstr=as.character(""),PACKAGE="librphast")

   	if (return$error!=0){
      		print(return$errstr)
      		return()
   	}
	
	l=list()
	l$msa=msa
	l$tm=treemod
	lmethod="SCORE"
	l$stats=stats
	l$pvals=return$pVals
	if (stats){
	   l$teststats=return$teststats
	   l$derivs=return$derivs
	}
	return(l)
  }

  if (method == "GERP"){
	return=.C("phyloP_GERP",msaAddress=as.numeric(attr(msa,"address")), tmAddress=as.numeric(attr(treemod, "address")), stats=as.integer(stats), modeNum=as.integer(modeInt), nrejected=as.numeric(rep(0,attr(msa,"length"))), nneut=as.numeric(rep(0,attr(msa,"length"))), nobs=as.numeric(rep(0,attr(msa,"length"))), nspec=as.numeric(rep(0,attr(msa,"length"))), error=as.integer(0), errstr=as.character(""),PACKAGE="librphast")

   	if (return$error!=0){
      		print(return$errstr)
      		return()
   	}
	
	l=list()
	l$msa=msa
	l$tm=treemod
	l$method="GERT"
	l$stats=stats
	l$nrejected=return$nrejected
	if (stats){
	    l$nneut=return$nneut
	    l$nobs=return$nobs
	    l$nspec=return$nspec
      	}
	return(l)
  }

  cat("ERROR: unsupported method\n")
  cat("Use one of {LRT, SPH, SCORE, GERP}\n")
  return()

}

print.bScore<- function(bScore){
   cat("Base score object\n")
   cat("Method used: ", bScore$method, "\n")
}

summary.bScore <- function(bScore){
   print(bScore)
   if (bScore$method=="SPH"){
	if (bScore$stats){
	   cat("Contains: pvals, postMeans, postVars\n")
	}
	cat("Contains: pvals\n")
   }
   if (bScore$method=="LRT"){
	if (bScore$stats){
	   cat("Contains: pvals, llrs, scales\n")
	}
	cat("Contains: pvals\n")
   }
   if (bScore$method=="SCORE"){
	if (bScore$stats){
	   cat("Contains: pvals, teststats, derivs\n")
	}
	cat("Contains: pvals\n")
   }
   if (bScore$method=="GERT"){
	if (bScore$stats){
	   cat("Contains: nrejected, nobs, nneut, nspec\n")
	}
	cat("Contains: nrejected\n")
   }

}

