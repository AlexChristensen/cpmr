.onload <- function(libname, pkgname)
{library.dynam("cpmr",package=pkgname,lib.loc=libname)}

.onAttach <- function(libname, pkgname)
{
    temp <- packageDescription("cpmr")
    msg <- paste("Package: ",temp$Package,": ",temp$Title,"\n",
               "Version: ",temp$Version,"\n",
               "Updated on: ",
               temp$Date,"\n", sep="")
    msg <- paste(msg,"Maintainer: Alexander P. Christensen, University of Pennsylvania\n",sep="")
    msg <- paste(msg,"Maintainer: Christopher N. Wahlheim, University of North Carolina at Greensboro\n",sep="")
    msg <- paste(msg,'For citation information, type citation("cpmr")\n')
    msg <- paste(msg,"For bugs and errors, submit an issue to <https://github.com/AlexChristensen/cpmr/issues>")
    packageStartupMessage(msg)
}