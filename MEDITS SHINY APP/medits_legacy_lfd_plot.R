# MEDITS Ocean v9 – true v4-style LFD charts.
# Reuses the exact ggplotFL::plotLengths() + ggpubr::ggarrange() rendering
# approach from the previous working application, with data already produced
# by MEDITS Ocean v8 (including optional pooled-LOESS sex redistribution).
# No re-estimation of length frequencies or survey indices takes place here.
# An emergency fallback to the v8 renderer is provided only for exceptional
# plotLengths failures, avoiding a crash of the whole Shiny session.

medits_lfd_legacy_plot <- function(a, view='total', sex=c('C','F','M','U'),
                                   corrected=TRUE, fallback=TRUE) {
  sex <- match.arg(sex)
  alloc <- medits_lfd_allocate(a,view)
  d <- switch(sex,
    C=alloc$total,
    F=if(isTRUE(corrected)) alloc$F_corrected else alloc$F_observed,
    M=if(isTRUE(corrected)) alloc$M_corrected else alloc$M_observed,
    U=alloc$U_observed)
  caption <- switch(sex,C='total',F='female',M='male',U='undetermined')
  shiny::validate(shiny::need(is.data.frame(d) && nrow(d)>0 &&
     all(c('Year','Length','Frequency') %in% names(d)) &&
     any(is.finite(d$Frequency) & d$Frequency>0),
     paste0('No positive ',caption,' length frequencies in the selected area.')))
  d <- d[is.finite(d$Year) & is.finite(d$Length) & d$Length>0 &
           is.finite(d$Frequency) & d$Frequency>=0,,drop=FALSE]
  yrs <- sort(unique(d$Year[d$Frequency>0]))
  d <- d[d$Year %in% yrs,,drop=FALSE]

  # The original v4 figure for undetermined sex was a faceted bar chart,
  # not the ggplotFL two-block display used for total/female/male.
  if(identical(sex,'U')) {
    return(ggplot2::ggplot(d,ggplot2::aes(x=Length,y=Frequency)) +
      ggplot2::geom_col(fill='#766E9C') +
      ggplot2::facet_wrap(~Year,scales='free_y') +
      ggplot2::theme_bw() +
      ggplot2::labs(title=paste('Undetermined LFD -',view),
                    x='Length (mm)',y='n/km2'))
  }

  render_original <- function() {
    # The v4 plotting code first split the available years into two groups.
    n <- length(yrs)
    i <- seq_len(ceiling(n/2))
    y1 <- yrs[i]
    y2 <- if(length(i)<n) yrs[-i] else yrs[0]
    all_lengths <- sort(unique(d$Length))

    make_block <- function(years) {
      z <- d[d$Year %in% years,,drop=FALSE]
      z <- stats::aggregate(Frequency ~ Year + Length,data=z,FUN=sum)
      # A complete length x year grid avoids fragile as.FLQuant reshaping and
      # preserves valid one-bin years and genuine zeros in the raw estimates.
      mat <- matrix(0,nrow=length(all_lengths),ncol=length(years),
        dimnames=list(as.character(all_lengths),as.character(years)))
      ix <- cbind(match(as.character(z$Length),rownames(mat)),
                  match(as.character(z$Year),colnames(mat)))
      mat[ix] <- z$Frequency
      mat[!is.finite(mat) | mat<0] <- 0
      len_stats <- medits_length_stats_yearly(mat)
      flq <- FLCore::FLQuant(array(as.numeric(mat),
        dim=c(nrow(mat),ncol(mat),1L,1L,1L,1L),
        dimnames=list(len=rownames(mat),year=colnames(mat),
                      unit='all',season='all',area='all',iter='1')),
        units='mm')
      # EXACT chart family from v4: plotLengths, blue mean, red L95,
      # white theme and Length (mm) x axis. Not a geom_col replacement.
      plot <- ggplotFL::plotLengths(flq) +
        ggplot2::geom_vline(data=len_stats,
          ggplot2::aes(xintercept=l95),linewidth=1,colour='red') +
        ggplot2::geom_vline(data=len_stats,
          ggplot2::aes(xintercept=mean),linewidth=1,colour='blue') +
        ggplot2::xlab('Length (mm)')
      # ggplot objects can fail during rendering rather than construction:
      # catch such failures now, while the v8 fallback remains available.
      ggplot2::ggplot_build(plot)
      plot
    }
    p1 <- make_block(y1)
    if(length(y2)) {
      p2 <- make_block(y2)
      ggpubr::ggarrange(p1,p2,ncol=1)
    } else p1
  }
  if(!isTRUE(fallback)) return(render_original())
  tryCatch(render_original(),error=function(e) {
    warning('The original ggplotFL plotLengths chart could not be rendered for ',
      caption,' / ',view,': ',conditionMessage(e),
      '. Showing the safe v8 display for this selection only.',call.=FALSE)
    medits_lfd_display_plot(a,view,sex,corrected=corrected)
  })
}