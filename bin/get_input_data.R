#!/usr/bin/env Rscript 

# extract reference and query data to be used in the control workflow
# required data include:
# - reference expression matrix (10X dir)
# - query expression matrix 
# - reference SDRF file (need cell ids, inferred cell types, and CL terms)
# - marker genes file for reference dataset 

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(R.utils))

option_list = list(
    make_option(
        c("-a", "--reference-mat-url"),
        action = "store",
        default = NA,
        type = 'character',
        help = ""
    ),
    make_option(
        c("-b", "--reference-barcodes-url"),
        action = "store",
        default = NA,
        type = 'character',
        help = ""
    ),
    make_option(
        c("-c", "--reference-genes-url"),
        action = "store",
        default = NA,
        type = 'character',
        help = ""
    ),
    make_option(
        c("-d", "--query-matrix-url"),
        action = "store",
        default = NA,
        type = 'character',
        help = ""
    ),
    make_option(
        c("-e", "--query-barcodes-url"),
        action = "store",
        default = NA,
        type = 'character',
        help = ""
    ),
    make_option(
        c("-f", "--query-genes-url"),
        action = "store",
        default = NA,
        type = 'character',
        help = ""
    ),
    make_option(
        c("-g", "--reference-metadata"),
        action = "store",
        default = NA,
        type = 'character',
        help = ""
    ),
    make_option(
        c("-i", "--marker-genes-file"),
        action = "store",
        default = NA,
        type = 'character',
        help = ""
    )
)

opt = wsc_parse_args(option_list, mandatory = c("reference_mat_url", "reference_barcodes_url", "reference_genes_url",
                                                 "query_matrix_url", "query_barcodes_url", "query_genes_url",
                                                 "reference_metadata", "marker_genes_file"))

urls = unlist(opt[1:8])
# standard file names
names = c(rep(c("matrix.mtx", "barcodes.tsv",  "genes.tsv"), 2), "ref_metadata.txt", "ref_marker_genes.txt")
# set up directories
ref_dir = "reference_10X_dir/"
query_dir = "query_10X_dir/"
dir.create(ref_dir)
dir.create(query_dir)

for(idx in 1:length(names)){
    if(idx <=3) d = paste0(ref_dir, names[idx], sep = "")
    else if(idx <=6) d = paste0(query_dir, names[idx], sep = "")
    else d = names[idx]

    zipped = FALSE
    if(endsWith(urls[idx], ".gz")){ 
        d = paste0(d, ".gz", sep="")
        zipped = TRUE
    }

    download.file(url = urls[idx], destfile=d)
    if(!file.exists(d)) stop(paste0("Failed to download file: ", d, sep=""))
    # gunzip, if necessary
    if(zipped) gunzip(d, overwrite = TRUE, remove = TRUE)
}
