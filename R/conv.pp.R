conv.pp = function(nca)
{
# Author: Kyun-Seop Bae k@acr.kr
# Convert an NCA result matrix into a data.frame, coercing the known PP
# parameter columns to numeric and leaving the identifier columns as character.
# INPUT
#   nca: matrix or data.frame of NCA results, with PPTESTCD column names
# RETURNS
#   data.frame with an "NCA" attribute of "pp"
  numericCols = c("Dose", "Dur", "b0", "CMAX", "CMAXD", "TMAX", "TLAG", "CLST", "CLSTP", "TLST", "LAMZHL", "LAMZ",
                  "LAMZLL", "LAMZUL" , "LAMZNPT","CORRXY", "R2", "R2ADJ", "AUCLST", "AUCALL", "AUCIFO",
                  "AUCIFOD", "AUCPEO", "AUCIFP", "AUCIFPD", "AUCPEP", "AUCPBEO", "AUCPBEP", "AUMCLST", "AUMCIFO",
                  "AUMCPEO", "AUMCIFP", "AUMCPEP",
                  "C0", "MRTIVLST", "MRTIVIFO", "MRTIVIFP", "MRTEVLST", "MRTEVIFO", "MRTEVIFP",
                  "VZFO", "VZFP", "CLFO", "CLFP", "VZO", "VZP", "CLO", "CLP", "VSSO", "VSSP")

  if (is.null(dim(nca))) stop("nca should be a matrix or a data.frame!")
  colNames = colnames(nca)
  if (is.null(colNames)) stop("nca should have column names!")
  nca = matrix(as.character(unlist(nca)), nrow=nrow(nca), ncol=ncol(nca), dimnames=list(NULL, colNames))

# Build the columns as a list first. The old loop started at column 2 and so
# failed outright on a single-column result, and grew the frame with cbind(),
# which re-coerced the already converted columns back to character.
  Cols = lapply(colNames, function(cCol) {
    if (cCol %in% numericCols) suppressWarnings(as.numeric(nca[,cCol])) else nca[,cCol]
  })
  names(Cols) = colNames
# row.names is given explicitly because a one-row input yields length-1 columns
# that carry the column name, which data.frame() would otherwise adopt as the
# row name (a single-row result was labelled "STUDYID").
  Res = data.frame(Cols, stringsAsFactors=FALSE, check.names=FALSE,
                   row.names=seq_len(nrow(nca)))
  attr(Res, "NCA") = "pp"
  return(Res)
}
