
library(bibtex)
library(RefManageR)

Bib <- ReadBib("Ref.bib", check = FALSE)

ref<-function(refkey){
  fname = Bib[refkey]$note
  system(paste0("open ", fname))
}

