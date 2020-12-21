#' Dataset preparation
#'
#' This step performs all preparation necessary to perform feamiR analysis, taking a set of mRNAs, a set of miRNAs and an interaction dataset and creating corresponding positive and negative datasets for ML modelling. 
#' PLEASE NOTE:
#' This analysis is run in Python so python must be installed and location specified if not on PATH. 
#' Both sreformat and PaTMaN must also be installed and path specified if not on PATH. 
#' Python >= 3.6 is required to use the neccesary packages.
#' The Python component required the following libraries: os, Bio, gtfparse, pandas, numpy, math, scipy.stats, matplotlib.pyplot, seaborn as sns, statistics, logging. Please ensure these are installed for the verison of Python you supply.
#' 
#' 
#' The function saves various files (using specified output_prefix) and if you wish to start preparation using one of these pre-output files then these can be specified and preparation will skip to that point (this should only be done with files output by the function).
#' @param mRNA_3pUTR Fasta file of only 3'UTRs, with gene name as name attribute (e.g. Serpinb8)
#' @param miRNA_full Fasta file of full mature miRNA hairpins, with miRNA ID as name attribute (e.g. hsa-miR-576-3p)
#' @param interactions CSV file containing only validated interactions between miRNA and mRNA (e.g. from miRTarBase). Must have columns miRNA (e.g. hsa-miR-576-3p), Target Gene (e.g. Serpinb8) and optionally Experiments (e.g. qRT-PCR) and/or Support Type (with values Functional MTI, Functional MTI (Weak), Non-Functional MTI, Non-Functional MTI (Weak))
#' @param annotations GTF file (e.g. from Ensembl) with attributes seqname (chromosome), feature (with 3'UTRs labelled exactly 'three_prime_utr'), transcript_id, gene_id and gene_name matching fullchromosomes and interactions
#' @param fullchromosomes Fasta file (e.g. top level file from Ensembl) containing full sequence for each chromosome with name as chromosome (e.g. 1, matching seqname from annotations)
#' @param nonseed_miRNA Binary, 1 if full miRNA features should be included in statistical analysis. Seed features are always included. Default: 0.
#' @param flanking_mRNA Binary, 1 if flanking region mRNA features should be included in statistical analysis. Seed features are always included. Default: 0.
#' @param UTR_output String. File name 3'UTR fasta file should be saved as (when annotations and full chromosomes files are supplied)
#' @param chr Number of chromosomes for species in question.
#' @param o Output prefix for any files created and saved.
#' @param positiveset CSV file containing validated pairs of miRNAs and mRNAs as output by initial stage of analysis. If positiveset and negative set are input, analysis begins at final statistical analysis stage.
#' @param negativeset CSV file containing non-validated pairs of miRNAs and mRNAs as output by initial stage of analysis. If positiveset and negative set are input, analysis begins at final statistical analysis stage.
#' @param pythonversion File path for installed Python version (default: python)
#' @param sreformatpath File path for installed sreformat (default: sreformat)
#' @param patmanpath File path for installed patman (default: patman)
#' @param minvalidationentries Minimum number of entries for a validation category to be considered separately in statistical analysis (default: 40)
#' @param patmanout TXT file containing patman output (saved as output_prefix + patman_seed.txt). If supplied, analysis begins at patman output processing stage.
#' @param num_runs Number of subsamples to create (default: 100)
#' @return CSV containing full positive and negative sets. Folder statistical_analysis of heatmaps showing significance of various features under Fisher exact and Chi-squared tests. Seed analysis will always be run, full miRNA and flanking analysis if the respective parameters are set to 1. Folder subsamples containing CSVs for 100 subsamples with positive and negative samples equal for use in classifiers and feature selection.
#' @export
#' @examples
#' preparedataset(pythonversion='python3',miRNA_full='miRNA_full_mmu.fasta',annotations='Mus_musculus.GRCm38.100.chr.gtf',fullchromosomes='Mus_musculus.GRCm38.dna.toplevel.fa',interactions='mmu_MTI.csv',chr=20,o='feamiR_ignore')
#' preparedataset(pythonversion='python3',interactions='mmu_MTI.csv',positiveset='feamiR_seed_positive.csv',negativeset='feamiR_seed_negative.csv',o='feamiR_')
preparedataset <- function(pythonversion='python',mRNA_3pUTR='',miRNA_full='',interactions='',annotations='',fullchromosomes='',seed=1,nonseed_miRNA=0,flankingmRNA=0,UTR_output='',chr='',o='feamiR_',positiveset='',negativeset='',sreformatpath='sreformat',patmanpath='patman',patmanoutput='',minvalidationentries=40,num_runs=100){
  if (!missing(mRNA_3pUTR)){mRNA_3pUTR = paste('-mRNA_3pUTR ',mRNA_3pUTR,sep='')}
  if (!missing(miRNA_full)){miRNA_full = paste('-miRNA_full ',miRNA_full,sep='')}
  if (!missing(interactions)){interactions = paste('-interactions ',interactions,sep='')}
  if (!missing(annotations)){annotations = paste('-annotations ',annotations,sep='')}
  if (!missing(fullchromosomes)){fullchromosomes = paste('-fullchromosomes ',fullchromosomes,sep='')}
  seed = paste('-seed ',seed,sep='')
  nonseed_miRNA = paste('-nonseed_miRNA ',nonseed_miRNA,sep='')
  flankingmRNA = paste('-flankingmRNA ',flankingmRNA,sep='')
  if (!missing(UTR_output)){UTR_output = paste('-UTR_output ',UTR_output,sep='')}
  if (!missing(chr)){chr = paste('-chr ',chr,sep='')}
  o = paste('-o ',o,sep='')
  if (!missing(positiveset)){positiveset = paste('-positiveset ',positiveset,sep='')}
  if (!missing(negativeset)){negativeset = paste('-negativeset ',negativeset,sep='')}
  patmanpath = paste('-patmanpath ',patmanpath,sep='')
  sreformatpath = paste('-sreformatpath ',sreformatpath,sep='')
  minvalidationentries = paste('-minvalidationentries ',minvalidationentries,sep='')
  num_runs = paste('-num_runs ',num_runs,sep='')
  if (!missing(patmanoutput)){patmanoutput = paste('-patmanoutput ',patmanoutput,sep='')}
  print(paste(pythonversion,'-W ignore ',system.file('python', "feamiR_prepare_dataset.py", package = "feamiR"),' ',mRNA_3pUTR,miRNA_full,interactions,annotations,fullchromosomes,seed,nonseed_miRNA,flankingmRNA,UTR_output,chr,o,positiveset,negativeset,sreformatpath,patmanpath,patmanoutput,minvalidationentries,num_runs))
  path=paste(pythonversion,'-W ignore ',system.file('python', "feamiR_prepare_dataset.py", package = "feamiR"),' ',mRNA_3pUTR,miRNA_full,interactions,annotations,fullchromosomes,seed,nonseed_miRNA,flankingmRNA,UTR_output,chr,o,positiveset,negativeset,sreformatpath,patmanpath,patmanoutput,minvalidationentries,num_runs)
  system(stringr::str_squish(path))
}
