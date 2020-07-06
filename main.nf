#!/usr/bin/env nextflow 

// Meta-workflow that controls execution of steps in cell type prediction tools evaluation framework

// download data
if(params.data_import.run == "True"){

	// extract matrix types required by the tools, conditioned on them being 'on' 
	tool_switch = ["True":0, "False":1]
	garnett_matrix_type = [params.garnett.matrix_type, null][tool_switch[params.garnett.run]]
	scmap_cluster_matrix_type = [params.scmap_cluster.matrix_type, null][tool_switch[params.scmap_cluster.run]]
	scmap_cell_matrix_type = [params.scmap_cell.matrix_type, null][tool_switch[params.scmap_cell.run]]
	scpred_matrix_type = [params.scpred.matrix_type, null][tool_switch[params.scpred.run]]
	
	UNIQUE_MATRIX_TYPES = Channel.of(
				garnett_matrix_type,
	                        scmap_cluster_matrix_type,
	                        scmap_cell_matrix_type,
	                        scpred_matrix_type)
	                    .filter{ it != null }
	                    .unique()

	// parse txt file with dataset accessions into a queue channel; build combinations with matrix types 
	IMPORT_PARAMS = Channel
	                .fromPath(params.data_import.training_datasets)
	                .splitCsv(header:false, sep:",")
	                .combine(UNIQUE_MATRIX_TYPES)

process fetch_training_datasets {
    conda "${baseDir}/envs/load_data.yaml"

    input:
        tuple val(dataset_id), val(seq_method), val(num_clust), val(barcode_col), val(cell_type_col), val(matrix_type) from IMPORT_PARAMS

    output:
        tuple file(dataset_id), val(dataset_id), val(barcode_col), val(cell_type_col), val(matrix_type) into TRAINING_DATA
        val(num_clust) into N_CLUST

    """
    if [ ${seq_method} ==  "droplet" ]; then 
        MATRIX_TYPE_UPD="CPM"
    else
        MATRIX_TYPE_UPD=${matrix_type}
    fi
    get_experiment_data.R\
                --accesssion-code ${dataset_id}\
                --config-file ${params.data_import.scxa_import_config_file}\
                --matrix-type \$MATRIX_TYPE_UPD\
                --output-dir-name ${dataset_id}\
                --get-sdrf ${params.data_import.get_sdrf}\
                --get-condensed-sdrf ${params.data_import.get_cond_sdrf}\
                --get-idf ${params.data_import.get_idf}\
                --get-marker-genes ${params.data_import.get_marker_genes}\
                --number-of-clusters ${num_clust}
    """
}

// to avoid problems with sdrf-barcode matching, unmelt condensed SDRF and use it in downstream processes
if(params.data_import.unmelt_sdrf.run == "True"){
    process unmelt_condensed_sdrf {
        conda "${baseDir}/envs/exp_metadata.yaml"

        input:
            tuple file(data), val(dataset_id), val(barcode_col), val(cell_type_col), val(matrix_type) from TRAINING_DATA

        output:
            tuple file("${data}.query"), file("${data}.ref"), val(dataset_id), val(barcode_col), val(cell_type_col), val(matrix_type) into TRAINING_DATA_UNMELT
	    file("${data}/unmelted_sdrf.tsv") into UNMELTED_SDRF_QUERY
        """
        unmelt_condensed.R\
                -i ${data}/condensed-sdrf.tsv\
                -o ${data}/unmelted_sdrf.tsv\
                --has-ontology\
                --retain-types ${params.data_import.unmelt_sdrf.retain_types}
	# rename data dir name to avoid downstream file name collision
	parallel cp -R ${data} ::: ${data}.query ${data}.ref 
        """
    } 
    TRAINING_DATA_UNMELT.set{ TRAINING_DATA_PROCESSED }
}else{
    TRAINING_DATA.set{ TRAINING_DATA_PROCESSED }
	}
	
	// fork queue channel contents into channels for corresponding tools
	TRAINING_DATA_PROCESSED.into{
	    GARNETT_DATA
	    SCMAP_CELL_DATA
	    SCMAP_CLUSTER_DATA
	    SCPRED_DATA
	}
	
	// add number of clusters to Garnett data 
	GARNETT_FULL_DATA = GARNETT_DATA.merge(N_CLUST)

}else{
	// parent cross-validation data input 
	
	MANUAL_INPUT_DATA = Channel.fromPath("${params.input_data}/*")
	BARCODE_COL = Channel.from(params.metadata.query_barcode_col_name).first()
	CELL_LABEL_COL = Channel.from(params.metadata.query_label_col_name).first()

	// add dataset id and matrix type to tuple	
	MANUAL_INPUT_DATA.map{it -> tuple(it, it.getBaseName().toString().split('\\.')[0], it.getBaseName().toString().split('\\.')[2] ) }.set{ZIP_FILES}
	//group by matrix type
	ZIP_FILES.groupTuple(by:[1, 2]).map{it -> tuple(it[0][0], it[0][1], it[1], it[2]) }.set{GROUPED_DATA}
	
	// unzip cross-validation data
	process unzip_data {
	
	input: 
	tuple file(test_zip), file(train_zip), val(dataset_id), val(matrix_type) from GROUPED_DATA
	val(barcode_col) from BARCODE_COL
	val(cell_label_col) from CELL_LABEL_COL	
	
	output: 
	tuple file("*.test.*"), file("*.train.*"), val(dataset_id), val(barcode_col), val(cell_label_col), val(matrix_type) into DATA
	file("*.test.*/marker_genes_*") into MARKERS
	file("*.test.*/unmelted_sdrf.tsv") into UNMELTED_SDRF_QUERY
	"""
	# unzip test 
	test_zipdir=\$(unzip -qql $test_zip | head -n1 | tr -s ' ' | cut -d' ' -f5- | sed 's|/.*||')
	unzip $test_zip
	# unzip train 
	unzip $train_zip
	"""
	}
	
	// extract marker's number of clusters value	
	MARKERS.map{it -> it.getBaseName().toString().split('\\_')[2] }.first().set{N_CLUST}

	// send data to different tool's channel 
	DATA.into{
	    GARNETT_DATA
	    SCMAP_CELL_DATA
	    SCMAP_CLUSTER_DATA
	    SCPRED_DATA
	}
		
	// add number of clusters to Garnett data
	GARNETT_FULL_DATA = GARNETT_DATA.merge(N_CLUST)
}

// keep only relevant version of the dataset 
GARNETT_FILTERED_DATA = GARNETT_FULL_DATA.filter{ it[5] == params.garnett.matrix_type }

//run garnett 
if(params.garnett.run == "True"){
    process run_garnett_workflow {
        publishDir "${params.tool_outputs_dir}", mode: 'copy'
        conda "${baseDir}/envs/nextflow.yaml"

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }
        
        input:
            tuple file(test_data), file(training_data), val(dataset_id), val(barcode_col), val(cell_label_col), val(matrix_type), val(num_clust) from GARNETT_FILTERED_DATA
        output:
             file("garnett_output.txt") into GARNETT_OUTPUT

        """
        RESULTS_DIR=\$PWD
        nextflow run $EVAL_WORKFLOWS/garnett-eval-workflow/main.nf\
                            -profile ${params.profile}\
                            --results_dir \$RESULTS_DIR\
                            --ref_10x_dir "${training_data}/10x_data"\
                            --query_10x_dir "${test_data}/10x_data"\
                            --marker_genes ${training_data}/marker_genes_${num_clust}.tsv\
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
}else{
    GARNETT_OUTPUT = Channel.empty()
}

// keep only relevant version of the dataset 
SCMAP_CELL_FILTERED_DATA = SCMAP_CELL_DATA.filter{ it[5] == params.scmap_cell.matrix_type }

// run scmap-cell
if(params.scmap_cell.run == "True"){
    process run_scmap_cell_workflow {
        publishDir "${params.tool_outputs_dir}", mode: 'copy'
        conda "${baseDir}/envs/nextflow.yaml"

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }

        input:
		tuple file(test_data), file(training_data), val(dataset_id), val(barcode_col), val(cell_label_col), val(matrix_type) from SCMAP_CELL_FILTERED_DATA
        output: 
            file("scmap-cell_output.txt") into SCMAP_CELL_OUTPUT

        """
        RESULTS_DIR=\$PWD    
        nextflow run $EVAL_WORKFLOWS/scmap-eval-workflow/main.nf\
                            -profile ${params.profile}\
                            --results_dir \$RESULTS_DIR\
                            --projection_method ${params.scmap_cell.projection_method}\
                            --query_10x_dir "${test_data}/10x_data"\
                            --reference_10x_dir "${training_data}/10x_data"\
                            --reference_metadata ${training_data}/unmelted_sdrf.tsv\
                            --output_dir_cell ${params.scmap_cell.output_dir_cell}\
                            --col_names ${params.scmap_cell.col_names}\
                            --cell_id_col ${params.metadata.ref_barcode_col_name}\
                            --cluster_col ${params.metadata.ref_label_col_name}\
                            --plot_file ${params.scmap_cell.plot_file}\
                            --threshold ${params.scmap_cell.threshold}
        """
    }
}else{
    SCMAP_CELL_OUTPUT = Channel.empty()
}

// keep only relevant version of the dataset 
SCMAP_CLUSTER_FILTERED_DATA = SCMAP_CLUSTER_DATA.filter{ it[5] == params.scmap_cluster.matrix_type }

// run scmap-cluster 
if(params.scmap_cluster.run == "True"){
    process run_scmap_cluster_workflow {
        publishDir "${params.tool_outputs_dir}", mode: 'copy'
        conda "${baseDir}/envs/nextflow.yaml"

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }

        input:
		tuple file(test_data), file(training_data), val(dataset_id), val(barcode_col), val(cell_label_col), val(matrix_type) from SCMAP_CLUSTER_FILTERED_DATA

        output:
            file("scmap-cluster_output.txt") into SCMAP_CLUST_OUTPUT

        """
        RESULTS_DIR=\$PWD
        nextflow run $EVAL_WORKFLOWS/scmap-eval-workflow/main.nf\
                            -profile ${params.profile}\
                            --results_dir \$RESULTS_DIR\
                            --projection_method ${params.scmap_cluster.projection_method}\
                            --query_10x_dir "${test_data}/10x_data"\
                            --reference_10x_dir "${training_data}/10x_data"\
                            --reference_metadata  ${training_data}/unmelted_sdrf.tsv\
                            --output_dir_cluster ${params.scmap_cluster.output_dir_cluster}\
                            --col_names ${params.scmap_cluster.col_names}\
                            --cell_id_col ${params.metadata.ref_barcode_col_name}\
                            --cluster_col ${params.metadata.ref_label_col_name}\
                            --plot_file ${params.scmap_cluster.plot_file}\
                            --threshold ${params.scmap_cluster.threshold}
        """
    }
}else{
    SCMAP_CLUST_OUTPUT = Channel.empty()
}

// keep only relevant version of the dataset 
SCPRED_FILTERED_DATA = SCPRED_DATA.filter{ it[5] == params.scpred.matrix_type }

// run scpred 
if(params.scpred.run == "True"){
    process run_scpred_workflow {
        publishDir "${params.tool_outputs_dir}", mode: 'copy'
        conda "${baseDir}/envs/nextflow.yaml"

        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }
        
        input:
		tuple file(test_data), file(training_data), val(dataset_id), val(barcode_col), val(cell_label_col), val(matrix_type) from SCPRED_FILTERED_DATA

        output:
            file("scpred_output.txt") into SCPRED_OUTPUT

        """
        RESULTS_DIR=\$PWD
        nextflow run $EVAL_WORKFLOWS/scpred-eval-workflow/main.nf\
                            -profile ${params.profile}\
                            --results_dir \$RESULTS_DIR\
                            --method ${params.scpred.method}\
                            -latest\
                            --training_10x_dir "${training_data}/10x_data"\
                            --prediction_10x_dir "${test_data}/10x_data"\
                            --metadata_file ${training_data}/unmelted_sdrf.tsv\
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
}else{
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
        publishDir "${params.label_analysis_outdir}", mode: 'copy'
        errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
        maxRetries 5
        memory { 16.GB * task.attempt }

        input:
            file(tool_outputs_dir) from COMBINED_RESULTS_DIR
            // NB: use query labels as 'true' labels; ref labels were used for model training
            file(query_lab_file) from UNMELTED_SDRF_QUERY.first()

        output:
            file("${params.label_analysis.tool_perf_table}") into TOOL_PERF_TABLE
            file("${params.label_analysis.tool_table_pvals}") into TOOL_TABLE_PVALS

        """
        RESULTS_DIR=\$PWD 
        nextflow run $EVAL_WORKFLOWS/label-analysis-eval-workflow/main.nf\
                            -profile cluster\
			    --results_dir \$RESULTS_DIR\
                            --input_dir ${tool_outputs_dir}\
                            --condensed_sdrf ${params.label_analysis.condensed_sdrf}\
		   	    --parallel ${params.label_analysis.parallel}\
                            --ontology_dict ${params.label_analysis.ontology_dict}\
                            --ontology_graph ${params.label_analysis.ontology_graph}\
			    --tool_perf_table ${params.label_analysis.tool_perf_table}\
                            --cell_anno_table ${params.label_analysis.cell_anno_table}\
                            --tool_table_pvals ${params.label_analysis.tool_table_pvals}\
                            --ref_labels_file ${query_lab_file}\
                            --empirical_dist ${params.label_analysis.empirical_dist}\
                            --num_iter ${params.label_analysis.num_iter}\
                            --num_cores ${params.label_analysis.num_cores}\
                            --cell_ontology_col ${params.metadata.ref_CL_col_name}\
                            --barcode_col_ref ${params.label_analysis.barcode_col_ref}\
                            --barcode_col_pred ${params.label_analysis.barcode_col_pred}\
                            --label_column_ref ${params.metadata.ref_label_col_name}\
                            --label_column_pred ${params.label_analysis.label_column_pred}\
                            --semantic_sim_metric ${params.label_analysis.semantic_sim_metric}
        """
    }
}
