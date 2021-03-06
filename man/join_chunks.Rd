% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NewPhaseParent.R
\name{join_chunks}
\alias{join_chunks}
\alias{link_dad_haps}
\title{\code{Joining phased haplotype chunks. }}
\usage{
join_chunks(GBS.array, chunklist, join_length, verbose, OR)

link_dad_haps(GBS.array, dad_haps_lofl, hapidx, join_length, OR)
}
\arguments{
\item{GBS.array}{A GBS.array object generated from \code{impute_parent} or\code{sim.array} functions.}

\item{chunklist}{A list of haplotype chunks generated from \code{phase_chunk}.}

\item{join_length}{The length of each neighboring chunks used to connect them into a longer one.}

\item{verbose}{Writing verbose messages. Default=TRUE.}
}
\value{
return A list of haplotype chunks.
}
\description{
Extending/joining the phased haplotype chunks generated from \code{phase_chunk}.
}
\examples{
GBS.array <- sim.array(size.array=50, numloci=1000)
win_length = 10
haps <- setup_haps(win_length)
probs <- error_mx(hom.error=0.02, het.error=0.8, imiss=0.2) 
chunklist <- phase_chunk(GBS.array, win_length, haps, verbose=TRUE)
res <- join_chunks(GBS.array, chunklist, join_length, verbose=TRUE)

}

