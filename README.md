<h3 align="center">CAFCOMPARER</h3> 


<!-- ABOUT THE PROJECT -->
## About The Project
This repository includes the private CRCT06 Tool to allow the Easy exploration 
of our CAF collection.


## Functionalities
Within this tool 3 main functionalities are fulfilled

### Analyze the enrichment of each CAF line in the different CAF subtypes

"data/CAF_signatures.xlsx" contains a sheet per article, including the subtypes 
defined innit + the genes given to characterize these subtypes.
Number of genes, and the way they were obtention varies, and should be considered, 
as it affect the enrichment obtained.

To add new signatures, just create a new sheet in "data/CAF_signatures.xlsx", and 
fill it using the same structure as the rest.

This tool works with Humans as well as mice signatures, although the culture CAFs are human.
For that, using ontology the closest human equivalent signature is built.

### Perform a comparative analysis of the different CAF signatures.

When new signatures are added it will be of iterest to explore how they relate to
the old ones. For that one could look at the enrichments produced in the first function. 
Yet, we also offer the option of performing this comparison Both Quantitatively and Qualitatively
via a Correlogram of the signature enrichments, and an Upsetplot respectively


### Explore the enrichment and expression of other genes

In the same structure as with the CAF subtypes, this tool offers the possibility 
of analyzing the enrichment of custom list of genes, as for example it is Tyrosyne 
kynase receptors ("data/Receptors_HGNC.csv"). For analyzing new list of genes just
imitate the structure of this file in a new CSV, and save it in the "data/" folder.

In more simple cases, This tool also allows to rank the CAF cell lines based in a 
single gene expression

## How does it work?

Runing SignatureTesting.R (click "source" in the upright corner) will produce 
several plots, which if specified, are generated into a file in "output/". We can classify according to the funcionalities described 
above. In this section I will list this resulting plots as well as specific details.

### Analyze the enrichment of each CAF line in the different CAF subtypes
- output/GSEA_heatmap.pdf -> general heatmap
- output/author_journal_year.pdf -> heatmap comparing the enrichment in the top with the genes used to determined this enrichment. Usefull to see how relevant these differences are.

### Perform a comparative analysis of the different CAF signatures.
- output/corrplot_CAF_subtypes.pdf -> Correlogram of the pairwise comparisons between the CAF signatures


### Explore the enrichment and expression of other genes
- output/gene_name.pdf: barplot ranking the CAFs based in the expression of the given name. 
  This is specified and can be easily changed in the Parameters section

- output/filename_heatmap.pdf -> Heatmap comparing the enrichment of the list (or lists) of other genes, together with their Zscore gen expression. Usefull to see how relevant these differences are.
- output/filename_boxplot.png -> grouped boxplot summarizing the expression of the list (or lists) of genes in each individual CAF line. Useful to choose a CAF expressing for instance a lot of Tyrosine Kynase receptors.



