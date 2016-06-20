#####################################
ALL CODES AND DATA FOR THE PAPER 'Inference of cellular level signaling networks using single-cell gene expression data in C. elegans reveals mechanisms of cell fate specification' by Xiao-Tai Huang, Yuan Zhu, Leanne Lai Hang Chan, Zhongying Zhao and Hong Yan.
#####################################


DATA:

  1. All single cell raw data are in '/data_and_codes/data/single_cell_data/' folder. The sub-folder 'w' includes all wild-type data, while the sub-folder 'm' includes all mutant data. The 'experiment_info.csv' file records single cell experiment information.

  2. The background network consists of PPI and PDI. All PPI and PDI data are contained in folder '/data_and_codes/data/background_network/'.

  3. The cause-effect pairs information are in '/data_and_codes/data/pathway_infer_include/Cause_effect_paris.txt' file. It is integrated by genetic interaction ('/data_and_codes/data/pathway_infer_include/GI.csv') and gene regulatory effect data ('/data_and_codes/data/pathway_infer_include/GR.csv').


CODES:
  All programs are coded by the integration of Perl5, SQL and R. If run the program, user should firstly install perl, SQLite, and R.

  1. Import2db_sqlite.pl
     This program can import all single cell data into tables of database (SQLite). Run the program as default or set options like this:

     $ perl Import2db_sqlite.pl -i=../data/single_cell_data/ -o=../data/single_cell_data/cdfiles_nhr25.db3 -x=../data/single_cell_data/experiment_info.csv

  2. Infer_effect_in_cells.pl
     This program can infer knockdown genes' effects in every cell of C. elegans embryo. R should be installed firstly before running this program. Run the program as default or set options like this:

     $ perl Infer_effect_in_cells.pl -d=../data/single_cell_data/cdfiles_nhr25.db3 -o=../results/single_cell_effect/ -c=1466.8

  3. Infer_effect_under_founders.pl
     This program can infer knockdown genes' effects under founder cells. R should be installed firstly before running this program. Run the program as default or set options like this:

     $ perl Infer_effect_under_founders.pl -i=../results/single_cell_effect/effect_matrix.csv -o=../results/founder_cell_effect/ -f='ABarp ABpla ABpra Caa Cpa'

  4. Find_connecting_paths.pl
     This program can find the connecting paths between the cause gene and the effect gene in the background network. To run the program, the perl module 'Paths::Graph' should be installed firstly. Then copy the file 'Graph2.pm' into the 'Paths::Graph' module folder, for example '/usr/local/share/perl5/Paths/' in Linux or C:\Perl\site\lib\Paths\ in Windows. The file 'Graph2.pm' is in the '/data_and_codes/codes/' folder, which is modified by Xiao-Tai Huang from the module 'Paths::Graph.pm'. The 'Graph2.pm' can limit the length of the finding paths. Run the program as default or set options like this:

     $ perl Find_connecting_paths.pl -l=5 -s=Y39E4B.1 -e=F11C1.6 -PPI=../data/background_network/PPIs.csv -PDI=../data/background_network/PDIs.csv
     $ perl Find_connecting_paths.pl -l=5 -s='Y39E4B.1 R10E11.1 C06G3.10 W09C2.1 K02B9.4' -e='F11C1.6 C38D4.6 T28H11.4'

  5. Confirm_candidate_paths_with_weight.pl
     This program can weight all the connecting paths and select some of them as the candidate paths. Run the program as default or set options like this:

     $ perl Confirm_candidate_paths_with_weight.pl -i=../results/connecting_paths/ -o=../results/candidate_paths/ -n=5

  6. Infer_pathway.pl
     This program can infer pathway in different founders based on candidate paths. Run the program as default or set options like this:

     $ perl Infer_pathway.pl -f='ABarp ABpla ABpra Caa Cpa'

  7. Draw_cell_lineage_tree.pl
     This program can visualize single cell data by drawing cell lineage tree into .svg file. To run this program, firstly all single cell data should be imported into SQLite database by 'Import2db_sqlite.pl'. To run the program, user should set the first parameter representing to draw which data set, for example 'w_t3' means the wild-type table 3 in the SQLite database. The other options can be set as default or like this:

     $ perl Draw_cell_lineage_tree.pl w_t3 --root P0 --cellstage=350 --autozoom
     $ perl Draw_cell_lineage_tree.pl w_t3 --root P0 --cellstage=350 -b --axis --brand --label --bottomlabel
     $ perl Draw_cell_lineage_tree.pl w_t3 --root P0 --endtime=140 -m --axis --brand --label --bottomlabel


RUN THE PROGRAMS:
  To only draw the cell lineage tree, user should run the programs in this order:
     1. run Import2db_sqlite.pl
     2. run Draw_cell_lineage_tree.pl

  To get the results of inferred pathways, user should run all the programs in this order:
     1. run Import2db_sqlite.pl
     2. run Infer_effect_in_cells.pl
     3. run Infer_effect_under_founders.pl
     4. run Find_connecting_paths.pl
     5. run Confirm_candidate_paths_with_weight.pl
     6. run Infer_pathway.pl


BINARY DISTRIBUTION:
  We also support the binary distribution of 'Import2db_sqlite.pl' and 'Draw_cell_lineage_tree.pl' both in Linux and Windows. The binary files of these two programs are also included in '/data_and_codes/codes/' folder.


RESULTS:
  All the results are included in '/data_and_codes/results/' folder. User can test all programs to get the results.
