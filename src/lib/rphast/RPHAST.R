#Alexandra Denby
#RPHAST.R
#First test case for RPHAST
# 4/12/08

phast.init <- function(){
   dyn.load("/home/alex/ResearchProject/phast/lib/librphast.so")
}

#creates a new alignment object
msa.new <- function(){
  a=list()
  class(a)="msa"
  return(a)
}


#read an alignment file, return the alignment object associated with it
#set maxseq>100 if >100 species in file
msa.read <- function(fname, format="FASTA", maxseq=100){
   return=.C("rph_msa_read",fname=as.character(fname), format=as.character(format), address=as.integer(0), numberSpecies=as.integer(0),length=as.integer(0),alphabet=as.character(""),error=as.integer(0), errstr=as.character("")),PACKAGE="librphast")

   if (return$error!=0){
      print(return$errstr)
      return()
   }

   msa=msa.new()
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
   return=.C("readSpecies", address=as.integer(attr(msa, "address")), numberSpecies=as.integer(attr(msa, "numberSpecies")), species=as.character(rep("",attr(msa, "numberSpecies"))), error=as.integer(0), errstr=as.character(""),PACKAGE="librphast")
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
   return=.C("rph_msa_print",fname=as.character(fname), format=as.character(format), address=as.integer(attr(msa,"address")), error=as.integer(0), errstr=as.character(""),PACKAGE="librphast")
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

   .C("rph_msa_free",address=as.integer(attr(msa,"address")), error=as.integer(0), errstr=as.character(""),PACKAGE="librphast")

   if (return$error!=0){
      print(return$errstr)
      return()
   }

   attr(msa,"loaded")=0
   return(msa)
}

#print function, alignments
print.alignment<- function(align){
   if(attr(align,"loaded")==0){
	cat("Alignment unloaded\n")
	return()
   }
   cat("\nMultiple Sequence Alignment\n")
   cat(attr(align,"numberSpecies"),"species,",attr(align,"length"),"bases\n\n")
}

#summary function, alignments
summary.alignment<- function(align){
   if(attr(align,"loaded")==0){
	cat("Error: alignment was freed from memory\nPlease reload\n")
	return()
   }
   print(align)
   cat("Species:\n")
   print(attr(align,"species"))
   cat("\n")
}




####################TREE MODELS######################


tm.new <- function(){
  a=list()
  class(a)="tm"
  return(a)
}


#reads a tree model from a file
tm.read <- function(fname){
   return=.C("rph_tm_read",fname=as.character(fname), treeAddress=as.integer(0), modAddress=as.integer(0), error=as.integer(0), errstr=as.character(""), PACKAGE="librphast")

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
      cat("Error: align must be of class \"tm\"\n")
      return()
   }
   return=.C("rph_tr_scale",modAddress=as.integer(attr(treemod, "address")), treeAddress=as.integer(0), scale=as.numeric(scale), error=as.integer(0), errstr=as.character(""),PACKAGE="librphast")
   
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
      cat("Error: align must be of class \"msa\"\n")
      return()
   }
   return=.C("rph_tm_print",fname=as.character(fname), address=as.integer(attr(treemod, "address")), error=as.integer(0), errstr=as.character(""), PACKAGE="librphast")

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
   .C("rph_tm_free",address=as.integer(attr(tm,"address")), error=as.integer(0), errstr=as.character(""),PACKAGE="librphast")

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


###########MORE FUNCTIONS#################

#Again, used max species
tm.generateMSA<-function(treemod, nsites, maxsp=100){
   if (class(treemod)!="tm"){
      cat("Error: align must be of class \"tm\"\n")
      return()
   }
   return=.C("rph_tm_generate_msa",modAddress=as.integer(attr(treemod, "address")), numSites=as.integer(nsites), msaAddress=as.integer(0), numberSpecies=as.integer(0), length=as.integer(0), alphabet="", species=rep("",maxsp), error=as.integer(0), PACKAGE="librphast")
   
   msa=msa.new()
   attr(msa,"address")=return$msaAddress
   attr(msa,"format")="FASTA"
   attr(msa,"numberSpecies")=return$numberSpecies
   attr(msa,"length")=return$length
   attr(msa,"alphabet")=return$alphabet
   attr(msa,"species")=readspecies(msa)
   attr(msa,"loaded")=1
   return(msa)
}