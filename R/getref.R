ref<-function(refkey){
  library(bibtex)
  library(RefManageR)
  Bib <- ReadBib("Ref.bib", check = FALSE)
  fname = Bib[refkey]$note
  system(paste0("open ", fname))
}

