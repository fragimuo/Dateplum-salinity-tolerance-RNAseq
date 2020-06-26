#!/usr/bin/env Rscript
#' title: "Extract the sequences for B2G"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/storage/reference/Picea-abies/v1.1")
#' ```

#' Load libraries
suppressPackageStartupMessages(library(Biostrings))
#suppressPackageStartupMessages(library(pander))
#suppressPackageStartupMessages(library(VennDiagram))
#suppressPackageStartupMessages(library(org.At.tair.db))
#suppressPackageStartupMessages(library(genomeIntervals))

#' Source utilities
source("/mnt/picea/home/fmunoz/Git/UPSCb/src/R/blastUtilities.R")

#' Read the AA sequence
prot <- readAAStringSet ("/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/trinity/Trinity.format.fasta.transdecoder.pep")

#' # Export
#' Prepare the blast input
chk <- breakInChunks(length(prot),400)


dev.null <- mclapply(1:length(chk),function(i){
  writeXStringSet(prot[start(chk)[i]:end(chk)[i]],
                  file=paste0("/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/annotation/hmmer/tmp/Trinity.fasta.transdecoder.pep.",i))  
})

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```

