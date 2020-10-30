# BrugiaPopulationGenomics
Code Used in Analysis of "Extreme loss of heterozygosity on chromosome X in natural and laboratory populations of Brugia nematodes"


All shell scripts with the Download suffix are used for downloading SRA sequencing data. SRA numbers can be found in the paper. 

The shell scripts Brugia_malayi_population_genomics_shell.sh and Brugia_pahangi_population_genomics_shell.sh contain the bulk of the commands for
mapping, sorting and de-duplicating reads for individuals from those respective species. They also contain the commands for variant calling and filtering
and the VCFTools commands used to merge VCF files and call Tajima's D, relatedness, runs of homozygosity, and all other pop genome results presented in the 
manuscript. It also generated Tables 1-4 in .tsv form, which were formatted in excel in their final forms. 

The R script Final_PopGenome_Figures.R contains all of the R code used to plot every figure presented in the manuscript. In addition, it also contains the 
code for the calculation of SNV density in all species (B. pahangi, B. malayi, W. bancrofti, B. timori) as well as the assignment of contigs for W. bancrofti
and B. timori to B. malayi chromosomes. 

Any questions about the code can be directed to John Mattick at jmattick@som.umaryland.edu, or to jeanmattique@gmail.com. 
