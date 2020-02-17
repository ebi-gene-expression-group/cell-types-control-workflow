# cell-types-control-workflow
A meta-workflow that controls the Nextflow pipelines for multiple tools in a scalable manner. It consists of: 
* [garnett-workflow](https://github.com/ebi-gene-expression-group/garnett-workflow)
* [scmap-workflow](https://github.com/ebi-gene-expression-group/scmap-workflow) - 2 versions: 'cluster' and 'cell'
* [scpred-workflow](https://github.com/ebi-gene-expression-group/scpred-workflow)
* [cell-types-analysis-workflow](https://github.com/ebi-gene-expression-group/cell-types-analysis-workflow)

Schematic representation of the workflow is shown below:

![](https://github.com/ebi-gene-expression-group/cell-types-control-workflow/blob/add_zooma_mapping/control_pipeline.png)

Rounded blocks represent processes or stand-alone pipelines, whereas rectangles represent output tables. Note that for each tool, independent pipeline is run, and additional tools can be easily incorporated if necessary. Each workflow/process are run in isolated environment and in parallel (when run on compute cluster).  

### Configuring the pipeline
This workflow is highly configurable, allowing to control both general and tool-scpecific parameters. To configure the pipeline, open `nextflow.config` file and make necessary changes. Config file is organised in a nested structure, with blocks corresponding to individual processes/tools. 
Key parameters to configure are:
- `data_download`: this block controls the query and reference expression data sets that will be privided as inputs to the analysed tools. Thus, `ref_acc_code` and `query_acc_code` correspond to Expression Atlas accession codes of reference and query data sets, respectively. `expr_data_type` controls whether the raw, filtered or normalised versions of the data sets are imported. `normalisation_method` indicates whether TPM- or CPM-normalised data should be imported. `query(ref)_num_clust` indicates the number of clusters in the marker gene files. See documentation of the `atlas-data-import` repository for more details. 



### Triggering the pipeline 
Obviously, you will need Nextflow installed to run analyses. The best way to install it is via Conda. It is recommended to create a clean environment to avoid dependency conflicts. Run the following commands:
```
conda create --name nextflow
conda activate nextflow
conda install nextflow
```
To run the pipeline, make sure you are within the `cell-types-control-workflow`directory. Then issue the following command:
```
nextflow run main.nf
```
To enable parallel execution, add `--profile cluster` parameter to the command. 








