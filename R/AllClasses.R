# AllClasses.R

#' An S4 class that stores links to Geno4imputeR objects, genotypes, loci missing rate, and MAF
#'
#' @slot genomx a matrix of bialleic genotypes, 0 ref, 1 het, 2 alt and 3 missing
#' @slot info SNP information, including SNP ID, chr, position, ref and alt alleles, loci missing rate.
#' @slot imiss individual plant missing rate
#' 
#' @exportClass Geno4imputeR
setClass("Geno4imputeR",
         slots=list(#filename="character",
                    genomx = "matrix",
                    info = "data.frame",
                    imiss= "data.frame"
                    #ref="integer",
                    #alt="IntegerList",
                    #genotypes="matrix",
                    #samples="character",
                    ))

# AllClasses.R
# Copyright (C) 2014 Vince Buffalo <vsbuffaloAAAAA@gmail.com>
# Distributed under terms of the BSD license.

#' An S4 class that stores links to Tassel HDF5 objects, loci positions, and genotypes
#'
#' @slot filename path to HDF5 file
#' @slot ranges a \code{GenomicRanges} object of loci
#' @slot ref a Tassel-encoded integer vector of reference alleles
#' @slot alt a Tassel-encoded \code{IntegerList} of alternate alleles
#' @slot genotypes a matrix of bialleic genotypes
#' @slot samples sample names
#' @slot version string Tassel version
#'
#' @exportClass TasselHDF5
setClass("TasselHDF5",
         slots=list(filename="character",
                    #seqnames="Rle",
                    #positions="Rle",
                    ranges="GRanges",
                    #ref="character",
                    #alt="CharacterList",
                    ref="integer",
                    alt="IntegerList",
                    genotypes="matrix",
                    samples="character",
                    version="character"))