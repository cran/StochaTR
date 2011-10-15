solveSLP <-
function (fname = "SLPparams.dat", fdir = getwd()) 
{
    .C("SolveSLP_WMTR", as.character(fname), as.character(paste(fdir, 
        "/", sep = "")), PACKAGE = "StochaTR")
    "The solution has been successfully calculated!"
}

