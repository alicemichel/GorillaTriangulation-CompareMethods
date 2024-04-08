

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
