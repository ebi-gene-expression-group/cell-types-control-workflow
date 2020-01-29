#!/usr/bin/env nextflow 

// Meta-workflow that controls execution of steps in cell type 
// prediction tools evaluation framework


//TODO: add metadata files for both datasets

// extract reference and query data
if(params.data_download.run == "True"){
    process fetch_query_data{
        publishDir "${baseDir}/data", mode: 'copy'
        conda 'envs/load_data.yaml'

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }

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

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }

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

// run garnett 
// -work-dir $WORK_DIR/\$SUBDIR\
if(params.garnett.run == "True"){
    process run_garnett_workflow {
        publishDir "${params.tool_outputs_dir}", mode: 'copy'
        conda 'envs/nextflow.yaml'

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }
        
        input:
            file(reference_10X_dir) from REFERENCE_10X_DIR
            file(query_10X_dir) from QUERY_10X_DIR
            file(ref_marker_genes) from REF_MARKER_GENES

        output:
             file("garnett_output.txt") into GARNETT_OUTPUT

        """
        RESULTS_DIR=\$PWD

        nextflow run $GARNETT_GIT\
                            -r $GARNETT_GIT_BRANCH\
                            --results_dir \$RESULTS_DIR\
                            --ref_10x_dir ${reference_10X_dir}\
                            --query_10x_dir ${query_10X_dir}\
                            --marker_genes ${ref_marker_genes}\
                            --ref_cds_gene_id_type ${params.garnett.ref_cds_gene_id_type}\
                            --query_cds_gene_id_type ${params.garnett.query_cds_gene_id_type}\
                            --database ${params.garnett.database}\
                            --marker_gene_id_type ${params.garnett.marker_gene_id_type}\
                            --classifier_gene_type ${params.garnett.classifier_gene_type}\
                            --n_outgroups ${params.garnett.n_outgroups}\
                            --cell_id_field ${params.garnett.cell_id_field}\
                            --predicted_cell_type_field ${params.garnett.predicted_cell_type_field}
        """
    } 
} else{
    GARNETT_OUTPUT = Channel.empty()
}

// run scmap-cell
if(params.scmap_cell.run == "True"){ 
    process run_scmap_cell_workflow {
        publishDir "${params.tool_outputs_dir}", mode: 'copy'
        conda 'envs/nextflow.yaml'

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }

        input:
            file(reference_10X_dir) from REFERENCE_10X_DIR
            file(query_10X_dir) from QUERY_10X_DIR
            file(ref_metadata) from REFERENCE_METADATA

        output: 
            file("scmap-cell_output.txt") into SCMAP_CELL_OUTPUT

        """
        RESULTS_DIR=\$PWD    

        nextflow run $SCMAP_GIT\
                            -r $SCMAP_GIT_BRANCH\
                            --results_dir \$RESULTS_DIR\
                            --projection_method ${params.scmap_cell.projection_method}\
                            --query_10x_dir ${query_10X_dir}\
                            --reference_10x_dir ${reference_10X_dir}\
                            --reference_metadata ${ref_metadata}\
                            --output_dir_cell ${params.scmap_cell.output_dir_cell}\
                            --col_names ${params.scmap_cell.col_names}\
                            --cell_id_col ${params.metadata.barcode_col_name}\
                            --cluster_col ${params.metadata.ref_label_col_name}\
                            --plot_file ${params.scmap_cell.plot_file}\
                            --threshold ${params.scmap_cell.threshold}
        """
    }
} else {
    SCMAP_CELL_OUTPUT = Channel.empty()
}


// run scmap-cluster 
if(params.scmap_cluster.run == "True"){
    process run_scmap_cluster_workflow {
        publishDir "${params.tool_outputs_dir}", mode: 'copy'
        conda 'envs/nextflow.yaml'

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }

        input:
            file(reference_10X_dir) from REFERENCE_10X_DIR
            file(query_10X_dir) from QUERY_10X_DIR
            file(ref_metadata) from REFERENCE_METADATA

        output:
            file("scmap-cluster_output.txt") into SCMAP_CLUST_OUTPUT

        """
        RESULTS_DIR=\$PWD

        nextflow run $SCMAP_GIT\
                            -r $SCMAP_GIT_BRANCH\
                            --results_dir \$RESULTS_DIR\
                            --projection_method ${params.scmap_cluster.projection_method}\
                            --query_10x_dir ${query_10X_dir}\
                            --reference_10x_dir ${reference_10X_dir}\
                            --reference_metadata ${ref_metadata}\
                            --output_dir_cluster ${params.scmap_cluster.output_dir_cluster}\
                            --col_names ${params.scmap_cluster.col_names}\
                            --cell_id_col ${params.metadata.barcode_col_name}\
                            --cluster_col ${params.metadata.ref_label_col_name}\
                            --plot_file ${params.scmap_cluster.plot_file}\
                            --threshold ${params.scmap_cluster.threshold}
        """
    }
} else{
    SCMAP_CLUST_OUTPUT = Channel.empty()
}

// run scpred 
if(params.scpred.run == "True"){
    process run_scpred_workflow {
        publishDir "${params.tool_outputs_dir}", mode: 'copy'
        conda 'envs/nextflow.yaml'

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }
        
        input:
            file(reference_10X_dir) from REFERENCE_10X_DIR
            file(query_10X_dir) from QUERY_10X_DIR
            file(ref_metadata) from REFERENCE_METADATA

        output:
            file("scpred_output.txt") into SCPRED_OUTPUT

        """
        RESULTS_DIR=\$PWD

        nextflow run $SCPRED_GIT\
                            -r $SCPRED_GIT_BRANCH\
                            --results_dir \$RESULTS_DIR\
                            --method ${params.scpred.method}\
                            --training_10x_dir ${reference_10X_dir}\
                            --prediction_10x_dir ${query_10X_dir}\
                            --metadata_file ${ref_metadata}\
                            --eigenvalue_plot_path ${params.scpred.eigenvalue_plot_path}\
                            --train_probs_plot_path ${params.scpred.train_probs_plot_path}\
                            --prediction_probs_path ${params.scpred.prediction_probs_path}\
                            --model_predictions_path ${params.scpred.model_predictions_path}\
                            --confusion_table_path ${params.scpred.confusion_table_path}\
                            --normalised_counts_slot ${params.scpred.normalised_counts_slot}\
                            --cell_id_col_name ${params.metadata.barcode_col_name}\
                            --cell_types_col_name ${params.metadata.ref_label_col_name}\
                            --col_names ${params.scpred.col_names}\
                            --log_transform ${params.scpred.log_transform}\
                            --model ${params.scpred.model}
        """
    }
} else {
    SCPRED_OUTPUT = Channel.empty()
}

// Combine method outputs into single channel
 ALL_RESULTS = 
    GARNETT_OUTPUT
    .concat(SCMAP_CLUST_OUTPUT)
    .concat(SCMAP_CELL_OUTPUT)
    .concat(SCPRED_OUTPUT)
// place tool outputs into single dir
process combine_results{
    input:
        file(method_outputs) from ALL_RESULTS.collect()

    output:
        file('results_dir') into COMBINED_RESULTS_DIR

    """
    mkdir -p results_dir/
    for file in ${method_outputs}
    do
        mv \$file results_dir
    done
    """
}

// run analysis of predicted labels 
if(params.label_analysis.run == "True"){
    process run_label_analysis {
        conda 'envs/nextflow.yaml' 
        publishDir "${params.label_analysis.output_dir}", mode: 'copy'

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }

        input:
            file(tool_outputs_dir) from COMBINED_RESULTS_DIR
            file(ref_lab_file) from REFERENCE_METADATA

        output:
            file("${params.label_analysis.cell_anno_table}") into CELL_ANNO_TABLE
            file("${params.label_analysis.tool_perf_table}") into TOOL_PERF_TABLE
            file("${params.label_analysis.tool_table_pvals}") into TOOL_TABLE_PVALS

        """
        RESULTS_DIR=\$PWD 

        nextflow run $LABEL_ANALYSIS_GIT\
                            -r $LABEL_ANALYSIS_GIT_BRANCH\
                            --results_dir \$RESULTS_DIR\
                            --input_dir ${tool_outputs_dir}\
                            --ref_labels_file ${ref_lab_file}\
                            --tool_perf_table ${params.label_analysis.tool_perf_table}\
                            --cell_anno_table ${params.label_analysis.cell_anno_table}\
                            --tool_table_pvals ${params.label_analysis.tool_table_pvals}\
                            --num_iter ${params.label_analysis.num_iter}\
                            --num_cores ${params.label_analysis.num_cores}\
                            --cell_ontology_col ${params.metadata.ref_CL_col_name}\
                            --barcode_col_ref ${params.metadata.barcode_col_name}\
                            --label_column_ref ${params.metadata.ref_label_col_name}\
                            --semantic_sim_metric ${params.label_analysis.semantic_sim_metric}\
                            --ontology_graph ${params.label_analysis.ontology_graph}\
                            --empirical_dist ${params.label_analysis.empirical_dist}
        """
    }
}
