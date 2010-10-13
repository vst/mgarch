.First.lib <-
function(lib, pkg)
{
    library.dynam("mgarch", pkg, lib)
    mylib <- dirname(.path.package("mgarch"))
    ver <- packageDescription("mgarch", lib = mylib)["Version"]
    vertxt <- paste("\n\t`mgarch' version:", ver, "\n")
    introtxt <-
        paste("\n\t`mgarch' is a package for simulating, estimating and diagnosing MGARCH processes\n",
              "\t See `library (help=mgarch)' for details.\n\n",
              sep = "")
    if(interactive() || getOption("verbose"))
        cat(paste(vertxt, introtxt))
}
