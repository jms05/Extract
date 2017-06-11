pocessaMarixReges <- function(matrix)
{
  
  names = c()
  ontology = c()
  index = 1
  regEX = "GO:[0-9]+"
  regexOnto = ".Ontology"
  nomes = row.names(matrix)
  for(i in 1:length(nomes)){
    val_regex = str_extract(nomes[i],regEX)
    if(!is.na(val_regex)){
      val_ontology = str_extract(nomes[i],regexOnto)
      if(!is.na(val_ontology)){
        names[index] = (val_regex)
        ontology[index] = (matrix[nomes[i],1])
        index = index+1
      }
    }
  }
  new.df= data.frame(GO = names,Ontology = ontology)
  return(new.df)
}

install.packages("stringr")
library(stringr)
resproc = pocessaMarixReges(as.matrix(go1))
