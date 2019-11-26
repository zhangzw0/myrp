###SAVERg

#### Overview

SAVERg is designed to carry out two experiments (e.g. cell clustering and pseudotime trajectory analysis) to assess the imputation accuracy of SAVER.

_Yiqiu Tang 2019-09-19_



#### Instruction

The package depends the R package TSCAN which is archived in the Bioconductor repository. This package can be be install with BiocManager using the code below.

```R
BiocManager::install("TSCAN")
```



#### Contributions



#### Examples

##### Cell clustering

We take the iPSC dataset as an example. The rows correspond to genes and the columns correspond to cells. The cell labels are known and can be treated as the gold standard.
You can use the raw count matrix as the input matrix:

```bash
data('ipsc_count_data')
cell_clustering(ipsc_count_data, need.imputation=TRUE, imputed.data = FALSE)
```
To save time, you can just input the imputed result of SAVER and carry out cell clustering analysis:

```bash
data('ipsc_saver')
cell_clustering(ipsc_saver, need.imputation=FALSE, imputed.data = TRUE)
```
You can also input the raw count matrix and not perform SAVER imputation to compare the results with imputation and without imputation:

```bash
data('ipsc_count_data')
cell_clustering(ipsc_count_data, need.imputation=FALSE, imputed.data = FALSE)
```

#####Pseudotime trajectory analysis

We take the deng dataset as an example. The rows correspond to genes and the columns correspond to cells. The cell labels are known and can be treated as the gold standard.
You can use the raw count matrix as the input matrix:

```bash
data("deng_count_data") 
deng_cellLabels <- factor(colnames(deng_count_data),
                          levels=c('zygote', 'early 2-cell', 'mid 2-cell', 'late 2-cell',
                                   '4-cell', '8-cell', '16-cell', 'early blastocyst',
                                   'mid blastocyst', 'late blastocyst'))
trajectory_analysis(deng_count_data, cellLabels=deng_cellLabels,
                    need.imputation=FALSE, imputed.data = TRUE)
```
To save time, you can just input the imputed result of SAVER and carry out pseudotime trajectory analysis:

```bash
data("deng_saver") 
deng_cellLabels <- factor(colnames(deng_saver),
                          levels=c('zygote', 'early 2-cell', 'mid 2-cell', 'late 2-cell',
                                   '4-cell', '8-cell', '16-cell', 'early blastocyst',
                                   'mid blastocyst', 'late blastocyst'))
trajectory_analysis(deng_saver, cellLabels=deng_cellLabels,
                    need.imputation=TRUE, imputed.data=FALSE)
```
You can also input the raw count matrix and not perform SAVER imputation to compare the results with imputation and without imputation:

```bash
data("deng_count_data") 
deng_cellLabels <- factor(colnames(deng_count_data),
                          levels=c('zygote', 'early 2-cell', 'mid 2-cell', 'late 2-cell',
                                   '4-cell', '8-cell', '16-cell', 'early blastocyst',
                                   'mid blastocyst', 'late blastocyst'))
trajectory_analysis(deng_count_data, cellLabels=deng_cellLabels,
                    need.imputation=FALSE, imputed.data = FALSE)
```


For detialed usages, please refer to "SAVERg-manual.pdf".







































































