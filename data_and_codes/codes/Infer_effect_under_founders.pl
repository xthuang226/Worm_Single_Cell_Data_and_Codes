### Infer knockdown gene effect in the selected founder cells by considering their progeny.
### Return 1(activating), -1(inhibiting) or 0(no effect) for the knockdown gene effect in each founder cell.

use strict;
use warnings;
use feature 'say';
use DBI;
use Statistics::R;
use Getopt::Long;


#
# Infer_effect_under_founders.pl
# Reqirements: Perl5
# Author: Xiao-Tai Huang
# Email: xthuang226@gmail.com
# Version: 1.0
# Date: 2016/06/19
#


my $file        = "Infer_effect_under_founders.pl";
my $program     = uc($&) if $file =~ /^\w+/;
my $version     = "1.0";
my $date        = "2016/06/19";
my $author      = "Xiao-Tai Huang";


### program identification
my $title = "$program $version, $date $author\n";


### option variables
$::opt_help         = 0;
$::opt_founders     = 'ABarp ABpla ABpra Caa Cpa';
$::opt_input        = '../results/single_cell_effect/effect_matrix.csv';
$::opt_output       = '../results/founder_cell_effect/';


my $usage = <<"END_OF_USAGE";
${title}Syntax:   \L$program\E [options]
Function: Infer knockdown gene effect under founder cells.
Options:                                                        (defaults:)
  -h, --help         print usage
  -f, --founders     set set founder cells to be infered        ($::opt_founders)
  -i, --input        set input single cell effect .csv file     ($::opt_input)
  -o, --output       set output folder for the results          ($::opt_output)
Examples:
  \L$program\E
  \L$program\E -i=../results/single_cell_effect/effect_matrix.csv -o=../results/founder_cell_effect/ -f='ABarp ABpla ABpra Caa Cpa'
END_OF_USAGE


### process options
my @OrgArgv = @ARGV;
GetOptions(
  "help|h!",
  "founders|f=s",
  "input|i=s",
  "output|o=s"
) or die $usage;
!$::opt_help or die $usage;


### Set parameters
my ($founders, $input, $output) = ($::opt_founders, $::opt_input, $::opt_output);



my $output_tmp = $output;
my $output_whole;

while ($output_tmp =~ s/^([\.\/]*[^\.\/]+\/?)//) {
	$output_whole .= $&;
	mkdir $output_whole;
}




my @founders = split (/[\s,]+/, $founders);

### Load single cell effect data into a matrix
my $SCE = &file2matrix($input);

my $all_cells_under_founders = &get_cells($SCE, \@founders);

my $all_knock_genes = &get_genes($SCE);


### main func
&main($all_knock_genes, $output);












### Load file into a matrix ###
sub file2matrix {
	my ($file) = @_;
	my @M;
	open (INPUT,$file)||die "Error!\n$!";
	my $i = 0;
	foreach (<INPUT>){
		chomp;
		my @array = split(/[\s,]+/,$_);
		$M[$i] = \@array;
		$i++;
	}
	close INPUT;
	return \@M;
}


### Get all founders and their progeny
sub get_cells {
	my ($SCE,$founders) = @_;
	
	my $cells_info;
	my $namelist = $$SCE[0];
	shift @$namelist;

	my %founders;

	foreach (@$founders) {
		my $founder = $_;
		my $cellpath_f = &get_cellpath($founder);
		my $cells_info_f;
		my $j = 1;
		foreach (@$namelist) {
			my $cell = $_;
			my $cellpath_l = &get_cellpath($cell);
			if ($cellpath_l =~ /$cellpath_f/){
				$$cells_info_f{$cell} = $j;
			}
			$j++;
		}
		$founders{$founder} = $cells_info_f;
	}

	return \%founders;
}


### Get cell path. Cell path represents a path from the zygote ('P0') to the specified cell. For example of cell 'ABala', its cell path is 'P0.AB.a.l.a';
sub get_cellpath {
	my ($cell) = @_;
	my @cellpath;

	while ($cell =~ s/([a-z])$//) {
		unshift @cellpath, $1;
	}
	unshift @cellpath, $cell;

	until($cellpath[0] eq 'P0'){
		if ($cellpath[0] eq 'AB') {
			unshift @cellpath, 'P0';
		}
		if ($cellpath[0] eq 'P1') {
			unshift @cellpath, 'P0';
		}
		if ($cellpath[0] eq 'EMS') {
			unshift @cellpath, 'P1';
		}
		if ($cellpath[0] eq 'P2') {
			unshift @cellpath, 'P1';
		}
		if ($cellpath[0] eq 'MS') {
			unshift @cellpath, 'EMS';
		}
		if ($cellpath[0] eq 'E') {
			unshift @cellpath, 'EMS';
		}
		if ($cellpath[0] eq 'C') {
			unshift @cellpath, 'P2';
		}
		if ($cellpath[0] eq 'P3') {
			unshift @cellpath, 'P2';
		}
		if ($cellpath[0] eq 'D') {
			unshift @cellpath, 'P3';
		}
		if ($cellpath[0] eq 'P4') {
			unshift @cellpath, 'P3';
		}
		if ($cellpath[0] eq 'Z2') {
			unshift @cellpath, 'P4';
		}
		if ($cellpath[0] eq 'Z3') {
			unshift @cellpath, 'P4';
		}
	}

	return join ('.', @cellpath);
}


### Get knockout gene list ###
sub get_genes {
	my ($DATA) = @_;
	my %genes_info;
	for (1..$#$DATA){
		$genes_info{$$DATA[$_][0]} = $_;
	}
	return \%genes_info;
}


### Extract data for one knockout gene in one founder cell
sub extract_data_gf {
	my ($ko_gene, $founder_cell, $SCE, $all_cells_under_founders, $all_knock_genes) = @_;

	my $cells_fp = $$all_cells_under_founders{$founder_cell};

	my @gfdata;
	
	my $row = $$all_knock_genes{$ko_gene};

	foreach (sort keys %$cells_fp){
		my $cell = $_;
		my $col = $$cells_fp{$cell};

		push @gfdata, $$SCE[$row][$col] unless $$SCE[$row][$col] =~ /_/;
	}
	return \@gfdata;
}


### Build contingency table
sub build_contingency_table {
	my ($gfdata) = @_;
	my $wt_sig = 0;
	my $wt_all = scalar @$gfdata;
	my $mu_sig = 0;
	my $mu_act = 0;
	my $mu_rep = 0;

	foreach (@$gfdata){
		$mu_sig += abs $_;
		if ($_ == 1){
			$mu_act++;
		}
		if ($_ == -1){
			$mu_rep++;
		}
	}

	my $mu_all = @$gfdata - $mu_sig;
	return ($wt_sig,$wt_all,$mu_sig,$mu_all,$mu_act,$mu_rep);
}


####### fisher exact test ########
sub fisher_exact_test {
	my ($wt_sig,$wt_all,$mu_sig,$mu_all) = @_;
	
	use Statistics::R;
	my $R = Statistics::R->new();
	$R->startR;
	$R->send('options(warn=-1)');
	
	$R->send("data <- matrix(c($wt_sig,$wt_all,$mu_sig,$mu_all),nrow=2,byrow=T)");
	$R->send('res <- fisher.test(data)');
	$R->send('res$p.value');
	my $pvalue = $R->read;
	$pvalue =~ s/^\[.+\]\s+//;
	$R->stopR;
	return $pvalue;
}


######## compute sig per cell per gene ########
sub compute_sig {
	my ($wt_sig,$wt_all,$mu_sig,$mu_all,$mu_act,$mu_rep) = @_;
	
	my $pvalue = &fisher_exact_test($wt_sig,$wt_all,$mu_sig,$mu_all);
	
	if ($pvalue <= 0.1){		### a setting for significance
		if ($mu_act >= 2*$mu_rep){
			return (1, $pvalue);
		}elsif ($mu_rep >= 2*$mu_act){
			return (-1, $pvalue);
		}else{
			return (0,1.1);
		}
	}else{
		return (0, $pvalue);
	}
}


### Infer effect for one knockout gene in one founder cell
sub infer_effect_gf {
	my ($ko_gene, $founder_cell) = @_;
	
	my $gfdata = &extract_data_gf($ko_gene, $founder_cell, $SCE, $all_cells_under_founders, $all_knock_genes);

	my ($wt_sig,$wt_all,$mu_sig,$mu_all,$mu_act,$mu_rep) = &build_contingency_table($gfdata);

	my ($sig,$pvalue) = &compute_sig($wt_sig,$wt_all,$mu_sig,$mu_all,$mu_act,$mu_rep);

	return ($sig,$pvalue,$wt_sig,$wt_all,$mu_sig,$mu_all,$mu_act,$mu_rep);
}


### main
sub main {
	my ($genes, $output) = @_;

	open (OUTPUT,">$output/founder_effect_matrix.csv")||die "Error!\n$!";
	open (OUTPUT2,">$output/founder_pvalue_matrix.csv")||die "Error!\n$!";
	open (OUTPUT3,">$output/log_file")||die "Error!\n$!";

	### print head line (cell name) #####
	print OUTPUT 'Gene\Cell,';
	say OUTPUT join (',', @founders);

	print OUTPUT2 'Gene\Cell,';
	say OUTPUT2 join (',', @founders);

	foreach (sort keys %$genes){
		my $ko_gene = $_;
		say "\n--->Working in $ko_gene";
		say OUTPUT3 "\n--->Working in $ko_gene";
		print OUTPUT "$ko_gene";
		print OUTPUT2 "$ko_gene";
		foreach (@founders){
			my $founder_cell = $_;
			print "=>$ko_gene $founder_cell\t";
			print OUTPUT3 "=>$ko_gene $founder_cell\t";
			my ($sig, $pvalue, $wt_sig, $wt_all, $mu_sig, $mu_all, $mu_act, $mu_rep) = &infer_effect_gf($ko_gene,$founder_cell);
			say "$wt_sig,$wt_all,$mu_sig,$mu_all,$mu_act,$mu_rep";
			say OUTPUT3 "$wt_sig,$wt_all,$mu_sig,$mu_all,$mu_act,$mu_rep";
			say "$sig , $pvalue";
			say OUTPUT3 "$sig , $pvalue";
			print OUTPUT ",$sig";
			print OUTPUT2 ",$pvalue";
		}
		print OUTPUT "\n";
		print OUTPUT2 "\n";
	}
	
	close OUTPUT;
	close OUTPUT2;
	close OUTPUT3;
}

