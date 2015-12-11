# AllClasses.R

#' A S4 class stores genomic data.
#' 
#' Stores GBS data of progeny array and parental data, as well as pedigree information.
#' In the GBS data, 0 ref, 1 het, 2 alt and 3 missing.
#'
#' @slot true_parents List of parents' read data. List index must consistent with pedigree.
#' @slot gbs_parents List of parents' GBS data. 
#' @slot true_kids List of progeny's readl data.
#' @slot gbs_kids List of progeny's GBS data.
#' @slot pedigree Pedigree information.
#' 
#' @exportClass GBS.array
#' 
setClass("GBS.array",
         representation=representation(
             true_parents = "list", # list of data.frame(hap1, hap2)
             gbs_parents = "list", # for phased parents, list of data.frame(chunk, hap1, hap2)
             true_kids = "list",
             gbs_kids = "list",
             pedigree = "data.frame",
             snpinfo = "data.frame"),             
         prototype=prototype(
             true_parents = list(), # list of data.frame(hap1, hap2)
             gbs_parents = list(),
             true_kids = list(),
             gbs_kids = list(),
             pedigree = data.frame(),
             snpinfo = data.frame())
         )
