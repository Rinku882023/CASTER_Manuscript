### CASTER_Manuscript


**Summary**

**"Systems Genomics Reveals microRNA Regulation of ICS Response in Childhood Asthma" (submitted: under review)**

The variability and difficulty in quantifying responses to inhaled corticosteroids (ICS) among asthmatic patients have been previously addressed by the development of a measure known as the Cross-sectional Asthma STEroid Response (CASTER) (Kho _et al._, 2020). MicroRNAs (miRNAs) have been found to have significant impacts on asthma and inflammatory processes. The objective of this study was to identify key association between circulating miRNAs and ICS response in childhood asthma.

To identify miRNAs associated with ICS response, small RNA sequencing was conducted on peripheral blood serum obtained from 580 asthmatic children under ICS treatment as part of The Genetics of Asthma in Costa Rica Study (GACRS). Generalized linear models were employed for this purpose. To validate the findings, replication was performed on children receiving ICS from the Childhood Asthma Management Program (CAMP) cohort. The relationship between the replicated miRNAs and the transcriptome of Lymphoblastoid Cell Lines in response to a glucocorticoid was evaluated.

**Files in Repo**
```
Code_for_Table1.R                                 # Code for Table1 generation
LCL_geneExpressionDataAnalysis.R                  # Code for LCL differential gene Expression data analysis
PearsonCorrelationMatrix_pairedDesignExpressionData.R # Code for generating Pearson Correlation matrix for paired design expression data
Rep_3miRNAs_LCLgeneAssociationAnalysis.R              # Code for miRNA and gene expression association analysis
WGCNA_Analysis.R                                      # Code for weighted gene coexpression network analysis: gene module identification & miRNA- gene module 
                                                         eigenvalue association test
miRNA_CASTER_AssociationAnalysis.R                    # Code for miRNA-CASTER phenotype association analysis

```
