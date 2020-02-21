# cell-types-eval-control-workflow
A meta-workflow for evaluation of multiple scRNAseq cell type classification tools. It controls the Nextflow pipelines for multiple tools in a scalable manner. The following tools are included: 
* [garnett-eval-workflow](https://github.com/ebi-gene-expression-group/garnett-eval-workflow)
* [scmap-eval-workflow](https://github.com/ebi-gene-expression-group/scmap-eval-workflow) - 2 versions: 'cluster' and 'cell'
* [scpred-eval-workflow](https://github.com/ebi-gene-expression-group/scpred-eval-workflow)
* [label-analysis-eval-workflow](https://github.com/ebi-gene-expression-group/label-analysis-eval-workflow)

Schematic representation of the workflow is shown below:

![](https://github.com/ebi-gene-expression-group/cell-types-control-workflow/blob/develop/control_pipeline.png)

Rounded blocks represent processes or stand-alone pipelines, whereas rectangles represent output tables. Note that for each tool, independent pipeline is run, and additional tools can be easily incorporated if necessary. Each workflow/process are run in isolated environment and in parallel (when run on compute cluster).  

### Configuring the pipeline
This workflow is highly configurable, allowing to control both general and tool-scpecific parameters. To configure the pipeline, open `nextflow.config` file and make necessary changes. Config file is organised in a nested structure, with blocks corresponding to individual processes/tools. Most of the blocks can be disabled by setting `run = 'False'`.
Key parameters to configure are:
* `data_download`: this block controls the query and reference expression data sets that will be privided as inputs to the analysed tools.
    * `ref_acc_code` and `query_acc_code` correspond to Expression Atlas accession codes of reference and query data sets, respectively. 
    * `expr_data_type` controls whether the raw, filtered or normalised versions of the data sets are imported. 
    * `normalisation_method` indicates whether TPM- or CPM-normalised data should be imported. 
    * `query(ref)_num_clust` indicates the number of clusters in the marker gene files. See documentation of the `atlas-data-import` repository for more details. 
* `metadata`: this block sets the necessary metadata fields for downstream analysis. For both query and control metadata files, the following fields need to be provided:
    * Barcode: ID of the individual cells in the data set
    * Label: inferred cell type
    * Cell Ontology (CL): cell ontology term corresponding to each label
* See documentation for individual sub-workflows to explore which tool parameters can be modified. 

### Triggering the pipeline 
Obviously, you will need Nextflow installed to run analyses. The best way to install it is via Conda. It is recommended to create a clean environment to avoid dependency conflicts. Run the following commands:
```
conda create --name nextflow
conda activate nextflow
conda install nextflow
```

Before running pipeline, make sure you are within the `cell-types-eval-control-workflow` directory.

Then issue the following command:
```
./bin/fetch-tool-eval-workflows.sh
```
This will pull the latest version of repo containing individual pipelines from GitHub. To trigger pipeline exectution, run 
```
nextflow run main.nf --profile cluster
```
To run locally, remove `--profile cluster` parameter from the command. 

### Results 
When the pipeline execution finishes, tool analysis table will be located in `data/label_analysis_output/` directory. 








