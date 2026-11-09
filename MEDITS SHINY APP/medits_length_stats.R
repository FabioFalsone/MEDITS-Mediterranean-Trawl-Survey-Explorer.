# Standalone robust yearly LFD summaries (zero/one supported bin safe).
medits_length_stats_yearly <- function(LF_matrix) {
  LF_matrix <- as.matrix(LF_matrix)
  len <- suppressWarnings(as.numeric(rownames(LF_matrix)))
  if (any(is.na(len))) len <- seq_len(nrow(LF_matrix))
  yrs <- colnames(LF_matrix)
  
  res <- lapply(seq_along(yrs), function(j) {
    freq <- suppressWarnings(as.numeric(LF_matrix[,j]))
    # Compute only from finite positive frequencies at finite lengths.
    # A one-bin year has an observable mean and L95 at that exact length.
    # A zero-bin year remains NA; it must NOT break the other years' LFDs.
    ok <- is.finite(len) & is.finite(freq) & freq > 0
    if (!any(ok))
      return(data.frame(year=yrs[j],mean=NA_real_,l95=NA_real_))
    ll <- len[ok]; ff <- freq[ok]
    ord <- order(ll); ll <- ll[ord]; ff <- ff[ord]
    mean_len <- sum(ll*ff)/sum(ff)
    if (length(ll)==1L) {
      l95 <- ll[[1]]
    } else {
      cum <- cumsum(ff)/sum(ff)
      l95 <- stats::approx(x=cum,y=ll,xout=0.95,ties='ordered',rule=2)$y
    }
    data.frame(year=yrs[j],mean=mean_len,l95=l95)
  })
  out <- do.call(rbind, res)
  out$year <- factor(out$year, levels = yrs)
  out
}
