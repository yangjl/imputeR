% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim.array.R
\name{sim.array}
\alias{add_error}
\alias{copy.mom}
\alias{kid}
\alias{missing.idx}
\alias{ran.hap}
\alias{sim.array}
\title{\code{A progeny array simulator.}}
\usage{
sim.array(size.array, numloci, hom.error = 0.02, het.error = 0.8,
  rec = 0.25, selfing = 0.5, imiss = 0.5, misscode = 3)

kid(p1, p2, het.error, hom.error, rec, imiss, misscode)

ran.hap(numloci, p)

add_error(diploid, hom.error, het.error)

copy.mom(mom, co_mean)

missing.idx(nloci, imiss)
}
\arguments{
\item{size.array}{Size of the progeny array.}

\item{numloci}{Number of loci to simulate.}

\item{hom.error}{Homozygous error rate, default=0.02.}

\item{het.error}{Heterozygous error rate, default=0.8.}

\item{rec}{Recombination rate, default=0.25.}

\item{selfing}{Proportion of selfed progeny, default=0.5.}

\item{imiss}{Individual missing rate, default=0.5.}

\item{misscode}{Missing code, default=3.}

\item{p1}{first parent (a list)}

\item{p2}{second parent (a list)}

\item{p}{a vector of allele freq with length numloci.}

\item{diploid}{a vector of diploid genotype, for example, c(0, 1, 1, 0, 2, 2).}

\item{mom}{a list of mom's haplotype.}

\item{co_mean}{mean number of crossovers per chromosome.}

\item{numloci}{number of loci.}

\item{hom.error}{homozygous error rate.}

\item{het.error}{heterozygous error rate.}
}
\value{
Return GBS.array object. 

Slot1: true_parents, a list of data.frame(hap1, hap2).
Slot2: gbs_parents, a list of genotypes. For example, c(1, 2, 2, 0, 3).
Slot3: true_kids, a list of data.frame(hap1, hap2).
Slot4: gbs_kids, a list of kid genotypes. For example, c(1, 1, 3, 1, 2).
Slot5: pedigree, a data.frame (kid, p1, p2). Note, p1 is the focal parent.
Slot6: freq, a vector of reference allele freq for all SNPs.

  See \url{https://github.com/yangjl/imputeR/blob/master/vignettes/imputeR-vignette.pdf} for more details.

a list of true [[1]] and observed [[2]] kid

a vector of genotype,

a vector of diploid genotype.

a vector of diploid genotype.
}
\description{
Simulate a S4 GBS.array object with true_parents, gbs_parents, true_kids, gbs_kids and pedigree.
}
\examples{
test <- sim.array(size.array=50, numloci=10)
class(test)

ran.hap(3, c(0.5, 0.8, 0.1))

add_error(c(0, 1, 1, 0, 2, 2), 0.02, 0.8)

copy.mom(list(c(0, 1, 1, 0, 1, 1)), 0.25)
}

