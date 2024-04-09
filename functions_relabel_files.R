
grabLaCieBAR <- function(date, times, destn){
  
  path = "/Volumes/LaCie5/"
  times.for.weird <- c(times[1]-1, times, times[length(times)]+1)
  
  fls <- list.files(path, recursive = T, pattern=date)
  searchcode <- paste0(date, "T", times.for.weird, collapse="|")
  flsIntime <- fls[grep(searchcode, fls)]
  shortnames <- gsub(".*/", "", flsIntime)
  ids <- sub(".*/","",sub("_.*", "", flsIntime))
  idNames <- paste0(ids, "_", shortnames)
  print(idNames)
  cont <- readline("...Do thse names look right?")
  res <- c()
  if (cont == "T"){
     res1 <- file.copy(from = paste0(path,flsIntime), to = paste0(destn, "/", idNames), copy.date = TRUE)
     res <- c(res,res1)
  }
  return(res)
}



relabelBARfiles <- function(type="file", PAM=NULL){
  
  if (type=="file"){
    (fl <- list.files(".", pattern="S20*")[!list.files(".", pattern="S20*") %in% list.files(".", pattern="_S20*")])
    if (length(fl)!=0){
      for (fli in fl){
        file.rename(fli, paste0(PAM,"_",fli))
      }
    }
    cat("Renamed", length(fl), "file(s) to PAM", PAM)
  }
  
  if (type=="folder"){
    dir = list.dirs(recursive = F)
    if (length(list.files(".", pattern="S20.*.wav"))==0){
      for (i in dir){
        getwd()
        setwd(i);getwd()
        dirs <- list.dirs(recursive = T)
        for (j in 2:length(dirs)){
          setwd(dirs[j]);print(getwd())
          fl <- list.files(pattern="S202*.")
          pamID=substr(dirs[j], start = 3, stop = 3)
          for (fli in fl){
            file.rename(fli, paste0("../../",pamID,"_",fli)) #change to file.rename when you're confident!
          }
          setwd("../")
        }
        setwd("../")
      }
    }
  }
  
  if ((!type %in% c("file", "folder"))){
    warning("type must be either file or folder")
  }
  
}
