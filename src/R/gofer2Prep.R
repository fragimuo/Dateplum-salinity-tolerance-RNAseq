#' ---
#' title: "Picea abies gopher input file generation"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Working directory
setwd("/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance")
#' ```
#' Libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

#' Out dir
dir.create("gofer2",showWarnings = FALSE)

#' Callback function
cb <- function(x,pos){
  x %>% select(Row.names, 
               matches("GO.ID"),
               Enzyme.Code,
               InterPro.Signatures,
               Division) %>% 
    filter(!(is.na(Annotation.GO.ID) &
               is.na(Enzyme.Code) &
               is.na(InterPro.Signatures)) &
               !is.na(Division))
}

#' GO Black List (Top parents synonyms)
BL <- c("GO:0000004","GO:0007582","GO:0044699","GO:0008150",
        "GO:0005554","GO:0003674","GO:0008372","GO:0005575",
        # and obsolete terms
        "GO:0006803","GO:0006118")

#' # Data
b2g <- read_csv_chunked("annotation/B2GO_annot.txt",
                        callback = DataFrameCallback$new(cb),
                        chunk_size = 1e5)

#' # Extraction
#' ## GO
#' ### Viridiplantae
go.plant <- b2g %>% filter(Division=="Viridiplantae",!is.na(Annotation.GO.ID)) %>% 
  select("Row.names","Annotation.GO.ID")

#' Remove the black list
go.plant[,2] <- sapply(strsplit(go.plant[,2,drop=TRUE],"\\|"),function(g){
  paste(g[!g %in% BL],collapse="|")
})

#' Sanity check
stopifnot(nrow(go.plant %>% filter(Annotation.GO.ID == "")) == 0)

#' Export
write_delim(go.plant,"gofer2/persimmon_gene_to_go.tsv",col_names = FALSE,delim="\t")
            
#' ### Fungi
go.fungi <- b2g %>% filter(Division=="Fungi",!is.na(Annotation.GO.ID)) %>% 
  select("Row.names","Annotation.GO.ID")

#' Remove the black list
go.fungi[,2] <- sapply(strsplit(go.fungi[,2,drop=TRUE],"\\|"),function(g){
  paste(g[!g %in% BL],collapse="|")
})

#' Sanity check
stopifnot(nrow(go.fungi %>% filter(Annotation.GO.ID == "")) == 0)

#' Export
write_delim(go.fungi,"gofer2/fungi_gene_to_go.tsv",col_names = FALSE,delim="\t")

#' ## KEGG EC codes
#' ### Plant
kegg.plant <- b2g %>% filter(Division=="Viridiplantae",!is.na(Enzyme.Code)) %>% 
  select(Row.names,Enzyme.Code) %>% mutate(EC=gsub("EC:","",Enzyme.Code))
write_delim(kegg.plant[,c(1,3)],"gofer2/persimmon_gene_to_kegg.tsv",col_names = FALSE,delim="\t")

#' ### Fungi
kegg.fungi <- b2g %>% filter(Division=="Fungi",!is.na(Enzyme.Code)) %>% 
  select(Row.names,Enzyme.Code) %>% mutate(EC=gsub("EC:","",Enzyme.Code))
write_delim(kegg.fungi[,c(1,3)],"gofer2/fungi_gene_to_kegg.tsv",col_names = FALSE,delim="\t")

#' ## PFAM
#' ### Plant
pfam.plant <- b2g %>% filter(Division=="Viridiplantae",!is.na(InterPro.Signatures),grepl("PFAM",InterPro.Signatures)) %>% 
  select(Row.names,InterPro.Signatures) %>% mutate(PFAM=unlist(mclapply(strsplit(InterPro.Signatures,"\\|"),
                                                                        function(sig){
                                                                          paste(sub(" .*","",sig[grepl("PFAM",sig)]),collapse="|")
                                                                        },mc.cores=16L)))
write_delim(pfam.plant[,c(1,3)],"gofer2/persimmon_gene_to_pfam.tsv",col_names = FALSE,delim="\t")

#' ### Fungi
pfam.fungi <- b2g %>% filter(Division=="Fungi",!is.na(InterPro.Signatures),grepl("PFAM",InterPro.Signatures)) %>% 
  select(Row.names,InterPro.Signatures) %>% mutate(PFAM=unlist(mclapply(strsplit(InterPro.Signatures,"\\|"),
                                                                        function(sig){
                                                                          paste(sub(" .*","",sig[grepl("PFAM",sig)]),collapse="|")
                                                                        },mc.cores=16L)))
write_delim(pfam.fungi[,c(1,3)],"gofer2/fungi_gene_to_pfam.tsv",col_names = FALSE,delim="\t")

#' # Session Info
#' ```{r sessionInfo, echo=FALSE}
#' sessionInfo()
#' ```
