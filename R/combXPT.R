combXPT = function(folders, domain)
{
# Author: Kyun-Seop Bae k@acr.kr
# Read and row-combine one SDTM domain (e.g. EX, PC) across several folders.
  nFolder = length(folders)
  if (nFolder == 0) stop("You did not specify any folder!")
  if (length(domain) != 1) stop("One domain name needs.")
  folders = Trim(folders)
  domain = UT(domain)

  for (i in seq_len(nFolder)) {
    cFolder = folders[i]
    if (!(substr(cFolder, nchar(cFolder), nchar(cFolder)) %in% c("/", "\\"))) folders[i] = paste0(cFolder, "/")
  }

  XPT = NULL
  colNames = NULL
  for (i in seq_len(nFolder)) {
    cFileName = paste0(folders[i], domain, ".XPT")
    if (file.exists(cFileName)) {
      cXPT = foreign::read.xport(cFileName)
      for (j in seq_len(ncol(cXPT))) {
        cXPT[,j] = as.character(cXPT[,j])
        cXPT[is.na(cXPT[,j]),j] = ""
      }
      if (is.null(XPT)) {  # first file actually found (not necessarily folders[1])
        colNames = colnames(cXPT)
        XPT = cXPT
      } else {
        if (!setequal(colNames, colnames(cXPT))) warning(paste0(folders[i], domain, " does not have the same column!"))
        colNames = intersect(colNames, colnames(cXPT))
        XPT = rbind(XPT[colNames], cXPT[colNames])
      }
    }
  }
  if (is.null(XPT)) stop(paste0("No folder contained ", domain, ".XPT !"))
  return(XPT)
}
