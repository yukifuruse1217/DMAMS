# DMAMS

DMAMS is a program written in Ruby ver. 2.4.

The program detects genetic signature associated with selective sweep within the population.

Referene: Identifying potentially beneficial genetic mutations associated with monophyletic selective sweep and a proof-of-concept study with viral genetic data. Furuse Y. mSystems; 2021 Feb 23;6(1):e01151-20.

DOI: 10.1128/mSystems.01151-20



To run the script, the followings are required:

'fileutils'

   "gem install mechanize" for installation

'spreadsheet'

   "gem install spreadsheet" for installation

'machanize'

   "gem install mechanize"

'daru'

   "gem install daru -v0.1.6" for installation

'statsample'

   "gem install statsample" for installation


To run the DMAMS program, five arguments are required:

1) fasta_sequnce_file_name,

2) nwk_tree_fine_name,

3) "nuc" OR "pro" for clustering,

4) minimum_size_of_sequences_for_clustering (proportion), and

5) output_file_name.


Example command:

**ruby DMAMS_ver1.21.rb example.fas example.nwk pro 0.05 test.csv**


Sequence file should be in single-line fasta format; the first sequence should be outgroup.

Tree file should be in newick format; the tree must be rooted using the single outgroup.

In the two input files, there must not be any symbols including "space" for sequence name except "underscore".

Genetic signature can be either nucleotide or deduced amino acid. Please select the option either "nuc" or "pro". To detected genetic signature of deduced amino acid, sequence should be in-frame. (To detected genetic signature at single nucleotide, sequence can be untranslated region or out-frame coding region.)

Minimum_size_of_strains_for_clustering is a proportion of strains of all sequences to determine subpopulation. The value should be smaller than 0.5

Output file will be CSV.

In the output, genetic signatures associated with significant selective sweep are marked '1' in the "detection" column. Using example files, there should be '1' for the row of position 10 amino acid L. The mutation was designated as a beneficial mutation to generate the sample data.


*Please contact Yuki FURUSE (furusey.kyoto@gmail.com) for any inquiry.*
