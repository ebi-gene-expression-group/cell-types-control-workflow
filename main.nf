#!/usr/bin/env nextflow 

// Meta-workflow that controls execution of steps in cell type 
// prediction tools evaluation framework

// extract reference and query data
if(params.data_download.run == "True"){
    query_data = params.data_download.query_output_dir
    query_n_clust = params.data_download.query_num_clust.toString()
    query_markers = "marker_genes_" + query_n_clust + ".tsv"

    process fetch_query_data{
        publishDir "${baseDir}/data", mode: 'copy'
        conda "${baseDir}/envs/load_data.yaml"

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }

        output:
            file("${query_data}/query_10x_data") into QUERY_10X_DIR
            file("${query_data}/query_sdrf.txt") into CONDENSED_SDRF_QUERY
            file("${query_data}/query_${query_markers}") into QUERY_MARKERS

        """
        get_experiment_data.R\
                --accesssion-code ${params.data_download.query_acc_code}\
                --expr-data-type ${params.data_download.expr_data_type}\
                --normalisation-method ${params.data_download.normalisation_method}\
                --output-dir-name ${params.data_download.query_output_dir}\
                --get-condensed-sdrf\
                --get-marker-genes\
                --number-of-clusters ${params.data_download.query_num_clust}

        # rename files to avoid name collisions in subsequent processes
        mv ${query_data}/10x_data ${query_data}/query_10x_data
        mv ${query_data}/condensed-sdrf.tsv ${query_data}/query_sdrf.txt
        mv ${query_data}/${query_markers} ${query_data}/query_${query_markers}
        """
    }
    ref_data = params.data_download.ref_output_dir
    ref_n_clust = params.data_download.ref_num_clust.toString()
    ref_markers = "marker_genes_" + ref_n_clust + ".tsv"

    // condensed sdrf files need 'un-melting' 
    process unmelt_sdrf_query {
        conda "${baseDir}/envs/exp_metadata.yaml"
        input:
            file(condensed_sdrf) from CONDENSED_SDRF_QUERY

        output:
            file("query_metadata.tsv") into UNMELTED_SDRF_QUERY

        """
        unmelt_condensed.R\
                -i ${condensed_sdrf}\
                -o query_metadata.tsv\
                --retain-types\
                --has-ontology                 
        """
    }

    process fetch_ref_data{
        publishDir "${baseDir}/data", mode: 'copy'
        conda "${baseDir}/envs/load_data.yaml"

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }

        output:
            file("${ref_data}/ref_10x_data") into REF_10X_DIR
            file("${ref_data}/ref_sdrf.txt") into CONDENSED_SDRF_REF
            file("${ref_data}/ref_${ref_markers}") into REF_MARKERS

        """
        get_experiment_data.R\
                --accesssion-code ${params.data_download.ref_acc_code}\
                --expr-data-type ${params.data_download.expr_data_type}\
                --normalisation-method ${params.data_download.normalisation_method}\
                --output-dir-name ${params.data_download.ref_output_dir}\
                --get-condensed-sdrf\
                --get-marker-genes\
                --number-of-clusters ${params.data_download.ref_num_clust}

        # rename files to avoid name collisions in subsequent processes
        mv ${ref_data}/10x_data ${ref_data}/ref_10x_data
        mv ${ref_data}/condensed-sdrf.tsv ${ref_data}/ref_sdrf.txt
        mv ${ref_data}/${ref_markers} ${ref_data}/ref_${ref_markers}

        """
    }
    
    process unmelt_sdrf_ref {
        conda "${baseDir}/envs/exp_metadata.yaml"
        input:
            file(condensed_sdrf) from CONDENSED_SDRF_REF

        output:
            file("ref_sdrf_proc.tsv") into UNMELT_SDRF_REF

        """
        unmelt_condensed.R\
                -i ${condensed_sdrf}\
                -o ref_sdrf_proc.tsv\
                --retain-types\
                --has-ontology                 
        """
    }
    
}else{
    REF_10X_DIR = Channel.fromPath(params.cv_ref_10x_dir).first()
    QUERY_10X_DIR = Channel.fromPath(params.cv_query_10x_dir).first()
    UNMELT_SDRF_REF = Channel.fromPath(params.cv_unmelt_sdrf_ref).first()
    UNMELTED_SDRF_QUERY = Channel.fromPath(params.cv_unmelt_sdrf_query).first()
    REF_MARKERS = Channel.fromPath(params.cv_ref_markers).first()
}

// run garnett 
if(params.garnett.run == "False"){
    process run_garnett_workflow {
        publishDir "${params.tool_outputs_dir}", mode: 'copy'
        conda "${baseDir}/envs/nextflow.yaml"

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }
        
        input:
            file(reference_10X_dir) from REF_10X_DIR
            file(query_10X_dir) from QUERY_10X_DIR
            file(ref_marker_genes) from REF_MARKERS

        output:
             file("garnett_output.txt") into GARNETT_OUTPUT

        """
        RESULTS_DIR=\$PWD

        nextflow run $EVAL_WORKFLOWS/garnett-eval-workflow/main.nf\
                            -profile cluster\
                            --results_dir \$RESULTS_DIR\
                            --ref_10x_dir ${reference_10X_dir}\
                            --query_10x_dir ${query_10X_dir}\
                            --marker_genes ${ref_marker_genes}\
                            --pval-col ${params.garnett.pval_col}\
                            --groups-col ${params.garnett.groups_col}\
                            --gene-names ${params.garnett.gene_names}\
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
        conda "${baseDir}/envs/nextflow.yaml"

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }

        input:
            file(reference_10X_dir) from REF_10X_DIR
            file(query_10X_dir) from QUERY_10X_DIR
            file(ref_metadata) from UNMELT_SDRF_REF

        output: 
            file("scmap-cell_output.txt") into SCMAP_CELL_OUTPUT

        """
        RESULTS_DIR=\$PWD    

        nextflow run $EVAL_WORKFLOWS/scmap-eval-workflow/main.nf\
                            -profile cluster\
                            --results_dir \$RESULTS_DIR\
                            --projection_method ${params.scmap_cell.projection_method}\
                            --query_10x_dir ${query_10X_dir}\
                            --reference_10x_dir ${reference_10X_dir}\
                            --reference_metadata ${ref_metadata}\
                            --output_dir_cell ${params.scmap_cell.output_dir_cell}\
                            --col_names ${params.scmap_cell.col_names}\
                            --cell_id_col ${params.metadata.ref_barcode_col_name}\
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
        conda "${baseDir}/envs/nextflow.yaml"

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }

        input:
            file(reference_10X_dir) from REF_10X_DIR
            file(query_10X_dir) from QUERY_10X_DIR
            file(ref_metadata) from UNMELT_SDRF_REF

        output:
            file("scmap-cluster_output.txt") into SCMAP_CLUST_OUTPUT

        """
        RESULTS_DIR=\$PWD

        nextflow run $EVAL_WORKFLOWS/scmap-eval-workflow/main.nf\
                            -profile cluster\
                            --results_dir \$RESULTS_DIR\
                            --projection_method ${params.scmap_cluster.projection_method}\
                            --query_10x_dir ${query_10X_dir}\
                            --reference_10x_dir ${reference_10X_dir}\
                            --reference_metadata ${ref_metadata}\
                            --output_dir_cluster ${params.scmap_cluster.output_dir_cluster}\
                            --col_names ${params.scmap_cluster.col_names}\
                            --cell_id_col ${params.metadata.ref_barcode_col_name}\
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
        conda "${baseDir}/envs/nextflow.yaml"

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }
        
        input:
            file(reference_10X_dir) from REF_10X_DIR
            file(query_10X_dir) from QUERY_10X_DIR
            file(ref_metadata) from UNMELT_SDRF_REF

        output:
            file("scpred_output.txt") into SCPRED_OUTPUT

        """
        RESULTS_DIR=\$PWD

        nextflow run $EVAL_WORKFLOWS/scpred-eval-workflow/main.nf\
                            -profile cluster\
                            --results_dir \$RESULTS_DIR\
                            --method ${params.scpred.method}\
                            -latest\
                            --training_10x_dir ${reference_10X_dir}\
                            --prediction_10x_dir ${query_10X_dir}\
                            --metadata_file ${ref_metadata}\
                            --eigenvalue_plot_path ${params.scpred.eigenvalue_plot_path}\
                            --train_probs_plot_path ${params.scpred.train_probs_plot_path}\
                            --prediction_probs_path ${params.scpred.prediction_probs_path}\
                            --model_predictions_path ${params.scpred.model_predictions_path}\
                            --confusion_table_path ${params.scpred.confusion_table_path}\
                            --normalised_counts_slot ${params.scpred.normalised_counts_slot}\
                            --cell_id_col_name ${params.metadata.ref_barcode_col_name}\
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
        conda "${baseDir}/envs/nextflow.yaml"
        publishDir "${params.label_analysis.output_dir}", mode: 'copy'

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }

        input:
            file(tool_outputs_dir) from COMBINED_RESULTS_DIR
            // NB: use query labels as 'true' labels; ref labels were used for model training
            file(query_lab_file) from UNMELTED_SDRF_QUERY

        output:
            file("${params.label_analysis.tool_perf_table}") into TOOL_PERF_TABLE
            file("${params.label_analysis.tool_table_pvals}") into TOOL_TABLE_PVALS

        """
        RESULTS_DIR=\$PWD 

        nextflow run $EVAL_WORKFLOWS/label-analysis-eval-workflow/main.nf\
                            -profile cluster\
                            --results_dir \$RESULTS_DIR\
                            --input_dir ${tool_outputs_dir}\
                            --ref_labels_file ${query_lab_file}\
                            --tool_perf_table ${params.label_analysis.tool_perf_table}\
                            --cell_anno_table ${params.label_analysis.cell_anno_table}\
                            --tool_table_pvals ${params.label_analysis.tool_table_pvals}\
                            --num_iter ${params.label_analysis.num_iter}\
                            --num_cores ${params.label_analysis.num_cores}\
                            --cell_ontology_col ${params.metadata.query_CL_col_name}\
                            --barcode_col_ref ${params.metadata.query_barcode_col_name}\
                            --label_column_ref ${params.metadata.query_label_col_name}\
                            --semantic_sim_metric ${params.label_analysis.semantic_sim_metric}\
                            --ontology_graph ${params.label_analysis.ontology_graph}\
                            --empirical_dist ${params.label_analysis.empirical_dist}
        """
    }
}
