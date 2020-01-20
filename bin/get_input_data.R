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
        c("-t", "--data-type"),
        action = "store",
        default = NA,
        type = 'character',
        help = "Type of dataset being downloaded (must be 'query' or 'reference')"
    ),
    make_option(
        c("-a", "--mat-url"),
        action = "store",
        default = NA,
        type = 'character',
        help = "Link to download reference expression matrix in .mtx format"
    ),
    make_option(
        c("-b", "--barcodes-url"),
        action = "store",
        default = NA,
        type = 'character',
        help = "Link to download reference barcodes file in .tsv format"
    ),
    make_option(
        c("-c", "--genes-url"),
        action = "store",
        default = NA,
        type = 'character',
        help = "Link to download reference genes file in .tsv format"
    ),
    make_option(
        c("-g", "--ref-metadata"),
        action = "store",
        default = NA,
        type = 'character',
        help = "Link to reference metadata file in text format. Must be provided if data-type is 'reference'"
    ),
    make_option(
        c("-i", "--marker-genes-file"),
        action = "store",
        default = NA,
        type = 'character',
        help = "Link to reference marker genes file in .tsv format. Must be provided if data-type is 'reference'"
    ),
    make_option(
        c("-j", "--output-10x-dir"),
        action = "store",
        default = NA,
        type = 'character',
        help = "Output path for reference 10X directory"
    ), 
    make_option(
        c("-k", "--ref-metadata-path"),
        action = "store",
        default = "ref_metadata.txt",
        type = 'character',
        help = "Output path for reference metadata file in text format"
    ),
    make_option(
        c("-l", "--markers-path"),
        action = "store",
        default = "ref_marker_genes.txt",
        type = 'character',
        help = "Output path for marker gene table in text format"
    )
)

opt = wsc_parse_args(option_list, mandatory = c("data_type", "output_10x_dir", 
                                                "mat_url", "barcodes_url", "genes_url"))

out_dir = opt$output_10x_dir
dir.create(out_dir)
data_type = opt$data_type
# standard 10x direcoty structure 
std_10x_files = c("matrix.mtx", "barcodes.tsv",  "genes.tsv")

# stage urls and file names for downloading
if(data_type == "query"){
    urls = c(opt$mat_url, opt$barcodes_url, opt$genes_url)
    dest = paste(out_dir, std_10x_files, sep="/")
} else if(data_type == "reference") {
    if(is.na(opt$marker_genes_file) | is.na(opt$ref_metadata)){ 
        stop("Necessary arguments not provided. Run get_input_data.R --help to see documentation")
    }
    urls = c(opt$mat_url, opt$barcodes_url, opt$genes_url, opt$ref_metadata, opt$marker_genes_file)
    dest = c(paste(out_dir, std_10x_files, sep="/"), opt$ref_metadata_path, opt$markers_path)
} else {
    stop(paste("Incorrect data type provided", data_type))
}

# download data
for(idx in 1:length(urls)){
    url = urls[idx]
    d = dest[idx]
    # manage archived files
    zipped = FALSE
    if(endsWith(url, ".gz")){ 
        d = paste0(d, ".gz", sep="")
        zipped = TRUE
    }

    download.file(url = url, destfile=d)
    if(!file.exists(d)) stop(paste0("Failed to download file:", d))
    # gunzip, if necessary
    if(zipped) gunzip(d, overwrite = TRUE, remove = TRUE)
}
