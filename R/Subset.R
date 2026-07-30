Subset = function(Tbl, key, Cond)
{
# Author: Kyun-Seop Bae k@acr.kr
# Last modification: 2026.7.29
# Input
#  Tbl: table to be subsetted
#  key: column names for the condition
#  Cond: subsetting condition values
# Returns
#   subsetted table

  nCol = length(key)
  if (length(Cond) != nCol) stop("key and Cond should be same length!")
  missingKey = setdiff(key, colnames(Tbl))
  if (length(missingKey) > 0) stop(paste("No such column:", paste(missingKey, collapse=", ")))

# One vectorised comparison per key column instead of a scalar comparison per
# cell: the old double loop ran nrow(Tbl) * length(key) R-level iterations.
  vRow = rep(TRUE, nrow(Tbl))
  for (j in seq_len(nCol)) {
    vRow = vRow & (as.character(Tbl[,key[j]]) == as.character(Cond[j]))
  }
  vRow[is.na(vRow)] = FALSE
  return(Tbl[vRow, , drop=FALSE])
}
