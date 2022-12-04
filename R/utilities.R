.pretty <-
    function(x, label = "")
{
    txt <- paste0(label, ": ", paste(x, collapse = " "))
    paste0(strwrap(txt, indent = 0, exdent = 2), collapse = "\n")
}
