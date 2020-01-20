#!/usr/bin/env nextflow 

// Meta-workflow that controls execution of steps in cell type 
// prediction tools evaluation framework


// TODO: should we add cross-validation indices generation??
// TODO: add control statements on which tools need to be run ?


// after input data are obtained, run individual workflows 


// extract reference and query data
if(params.data_download.run == "True"){
    process fetch_query_data{
        publishDir "${baseDir}/data", mode: 'copy'
        conda 'envs/load_data.yaml'

        output:
            file("query_10X_dir") into QUERY_10X_DIR 

        """
        get_input_data.R\
                --data-type "query"\
                --mat-url ${params.data_download.query_matrix_url}\
                --barcodes-url ${params.data_download.query_barcodes_url}\
                --genes-url ${params.data_download.query_genes_url}\
                --output-10x-dir query_10X_dir
        """
    }
    process fetch_ref_data{
        publishDir "${baseDir}/data", mode: 'copy'
        conda 'envs/load_data.yaml'

        output:
            file("reference_10X_dir") into REFERENCE_10X_DIR 
            file("ref_metadata.txt") into REFERENCE_METADATA
            file("ref_marker_genes.txt") into REF_MARKER_GENES

        """
        get_input_data.R\
                --data-type "reference"\
                --mat-url ${params.data_download.reference_mat_url}\
                --barcodes-url ${params.data_download.reference_barcodes_url}\
                --genes-url ${params.data_download.reference_genes_url}\
                --ref-metadata ${params.data_download.reference_metadata}\
                --marker-genes-file ${params.data_download.marker_genes_file}\
                --output-10x-dir reference_10X_dir
        """
    }
}

// make channels re-usable
REFERENCE_10X_DIR = REFERENCE_10X_DIR.first()
QUERY_10X_DIR = QUERY_10X_DIR.first()
REFERENCE_METADATA = REFERENCE_METADATA.first()
REF_MARKER_GENES = REF_MARKER_GENES.first()

// run garnett 
if(params.garnett.run == "True"){
    process run_garnett_workflow {
        publishDir "${params.tool_outputs_dir}", mode: 'copy'
        conda 'envs/nextflow.yaml'

        input:
            file(reference_10X_dir) from REFERENCE_10X_DIR
            file(query_10X_dir) from QUERY_10X_DIR
            //file(ref_metadata) from REFERENCE_METADATA
            file(ref_marker_genes) from REF_MARKER_GENES

        output:
             file("garnett_output.txt") into GARNETT_OUTPUT

        """
        RESULTS_DIR=\$PWD
        SUBDIR="garnett"
        mkdir -p $WORK_DIR/\$SUBDIR

        nextflow run $GARNETT_BASE_DIR/main.nf\
                            -work-dir $WORK_DIR/\$SUBDIR\
                            --results_dir \$RESULTS_DIR\
                            --ref_10x_dir ${reference_10X_dir}\
                            --query_10x_dir ${query_10X_dir}\
                            --marker_genes ${ref_marker_genes}
        """
    } 
}

// run scmap-cell
if(params.scmap_cell.run == "True"){ 
    process run_scmap_cell_workflow {
        publishDir "${params.tool_outputs_dir}", mode: 'copy'
        conda 'envs/nextflow.yaml'

        input:
            file(reference_10X_dir) from REFERENCE_10X_DIR
            file(query_10X_dir) from QUERY_10X_DIR
            file(ref_metadata) from REFERENCE_METADATA

        output: 
            file("scmap-cell_output.txt") into SCMAP_CELL_OUTPUT

        """
        RESULTS_DIR=\$PWD
        SUBDIR="scmap_cell"
        mkdir -p $WORK_DIR/\$SUBDIR     

        nextflow run $SCMAP_BASE_DIR/main.nf\
                            -work-dir $WORK_DIR/\$SUBDIR\
                            --results_dir \$RESULTS_DIR\
                            --projection_method ${params.scmap_cell.projection_method}\
                            --query_10x_dir ${query_10X_dir}\
                            --reference_10x_dir ${reference_10X_dir}\
                            --reference_metadata ${ref_metadata}
        """
    }
}


// run scmap-cluster 
if(params.scmap_cluster.run == "True"){
    process run_scmap_cluster_workflow {
        publishDir "${params.tool_outputs_dir}", mode: 'copy'
        conda 'envs/nextflow.yaml'

        input:
            file(reference_10X_dir) from REFERENCE_10X_DIR
            file(query_10X_dir) from QUERY_10X_DIR
            file(ref_metadata) from REFERENCE_METADATA

        output:
            file("scmap-cluster_output.txt") into SCMAP_CLUST_OUTPUT

        """
        RESULTS_DIR=\$PWD
        SUBDIR="scmap_clust"
        mkdir -p $WORK_DIR/\$SUBDIR 

        nextflow run $SCMAP_BASE_DIR/main.nf\
                            -work-dir $WORK_DIR/\$SUBDIR\
                            --results_dir \$RESULTS_DIR\
                            --projection_method ${params.scmap_cluster.projection_method}
                            --query_10x_dir ${query_10X_dir}\
                            --reference_10x_dir ${reference_10X_dir}\
                            --reference_metadata ${ref_metadata}
        """
    }
}

// run scpred 
if(params.scpred.run == "True"){
    process run_scpred_workflow {
        publishDir "${params.tool_outputs_dir}", mode: 'copy'
        conda 'envs/nextflow.yaml'

        input:
            file(reference_10X_dir) from REFERENCE_10X_DIR
            file(query_10X_dir) from QUERY_10X_DIR
            file(ref_metadata) from REFERENCE_METADATA

        output:
            file("scpred_output.txt") into SCPRED_OUTPUT

        """
        RESULTS_DIR=\$PWD
        SUBDIR="scpred"
        mkdir -p $WORK_DIR/\$SUBDIR 

        nextflow run $SCPRED_BASE_DIR/main.nf\
                            -work-dir $WORK_DIR/\$SUBDIR\
                            --results_dir \$RESULTS_DIR\
                            --method ${params.scpred.method}\
                            --training_10x_dir ${reference_10X_dir}\
                            --prediction_10x_dir ${query_10X_dir}\
                            --metadata_file ${ref_metadata}
        """
    }
}
// run the analysis of predicted labels 
TOOL_OUTPUTS_DIR = Channel.fromPath(params.tool_outputs_dir)
//REF_LAB_FILE = Channel.fromPath(params.data_download.reference_metadata)
if(params.label_analysis.run == "True"){
    process run_label_analysis {
        conda 'envs/nextflow.yaml' 
        publishDir "${params.label_analysis.output_dir}", mode: 'copy'

        input:
            file(tool_outputs_dir) from TOOL_OUTPUTS_DIR
            file(ref_lab_file) from REFERENCE_METADATA

        output:
            file("${params.label_analysis.cell_anno_table}") into CELL_ANNO_TABLE
            file("${params.label_analysis.tool_perf_table}") into TOOL_PERF_TABLE
            file("${params.label_analysis.tool_table_pvals}") into TOOL_TABLE_PVALS

        """
        RESULTS_DIR=\$PWD
        SUBDIR="label_analysis"
        mkdir -p $WORK_DIR/\$SUBDIR 

        nextflow run $LABEL_ANALYSIS_BASE_DIR/main.nf\
                            -work-dir $WORK_DIR/\$SUBDIR\
                            --results_dir \$RESULTS_DIR\
                            --input_dir ${tool_outputs_dir}\
                            --ref_labels_file ${ref_lab_file}\
                            --tool_perf_table ${params.label_analysis.tool_perf_table}\
                            --cell_anno_table ${params.label_analysis.cell_anno_table}\
                            --tool_table_pvals ${params.label_analysis.tool_table_pvals}
        """

    }
}
