### Detect significant difference in every cell between wild-type and mutant datasets by using Wilcox statistic test.
### Return 1(activating), -1(inhibiting) or 0(no effect) for the gene expression difference in each cell.

use strict;
use warnings;
use feature 'say';
use DBI;
use Statistics::R;
use Getopt::Long;


#
# Infer_effect_in_cells.pl
# Reqirements: Perl5
# Author: Xiao-Tai Huang
# Email: xthuang226@gmail.com
# Version: 1.0
# Date: 2016/06/19
#


my $file        = "Infer_effect_in_cells.pl";
my $program     = uc($&) if $file =~ /^\w+/;
my $version     = "1.0";
my $date        = "2016/06/19";
my $author      = "Xiao-Tai Huang";


### program identification
my $title = "$program $version, $date $author\n";


### option variables
$::opt_help         = 0;
$::opt_cutoff       = 1466.8;
$::opt_indb         = '../data/single_cell_data/cdfiles_nhr25.db3';
$::opt_output       = '../results/single_cell_effect/';


my $usage = <<"END_OF_USAGE";
${title}Syntax:   \L$program\E [options]
Function: Infer knockdown gene effect in every cell.
Options:                                                  (defaults:)
  -h, --help        print usage
  -c, --cutoff      set cutoff value                      ($::opt_cutoff)
  -d, --indb        set input SQLite .db3 file path       ($::opt_indb)
  -o, --output      set output folder for the results     ($::opt_output)
Examples:
  \L$program\E
  \L$program\E -d=../data/single_cell_data/cdfiles_nhr25.db3 -o=../results/single_cell_effect/ -c=1466.8
END_OF_USAGE


### process options
my @OrgArgv = @ARGV;
GetOptions(
  "help|h!",
  "cutoff|c=f",
  "indb|d=s",
  "output|o=s"
) or die $usage;
!$::opt_help or die $usage;


### Set parameters
my ($indb, $output, $cutoff) = ($::opt_indb, $::opt_output, $::opt_cutoff);



my $output_tmp = $output;
my $output_whole;

while ($output_tmp =~ s/^([\.\/]*[^\.\/]+\/?)//) {
	$output_whole .= $&;
	mkdir $output_whole;
}


open (OUTPUT,">$output/effect_matrix.csv")||die "Error!\n$!";
open (OUTPUT2,">$output/pvalue_matrix.csv")||die "Error!\n$!";
open (OUTPUT3,">$output/log_file")||die "Error!\n$!";

### Connect to the database.
my $conn = DBI->connect("dbi:SQLite:dbname=$indb","","",{ RaiseError => 1 }) or die $DBI::errstr;
my $wt = &get_experiments_data('w');
my $ko_list = &get_ko_list;

say "Parameters: cutoff = $cutoff";
say OUTPUT3 "Parameters: cutoff = $cutoff";
say OUTPUT join(",", ('Gene\Cell', sort keys %$wt));
say OUTPUT2 join(",", ('Gene\Cell', sort keys %$wt));

my $i = 1;
foreach my $ko_gene (@$ko_list) {
	say "\n[$i/".(scalar @$ko_list)."]";
	say OUTPUT3 "\n[$i/".(scalar @$ko_list)."]";
	say "---> Working in $ko_gene ...";
	say OUTPUT3 "---> Working in $ko_gene ...";
	printf ("%-20s%-10s%s\n", ('Cell','Effect','p-value'));
	printf OUTPUT3 ("%-20s%-10s%s\n", ('Cell','Effect','p-value'));

	my ($line_sig, $line_pvalue) = &compute_sig_per_ko($wt, $ko_gene);

	say OUTPUT join(",", @$line_sig);
	say OUTPUT2 join(",", @$line_pvalue);
	$i++;
}


### Disconnect to the database.
$conn->disconnect;
$conn = undef;

close OUTPUT;
close OUTPUT2;
close OUTPUT3;












### Obtain gene expression data from either wild-type or mutant datasets.
sub get_experiments_data {
	my ($type, $ko_gene) = @_;
	my $sql;

	if ($type eq 'w') {			# Select wild-type tables and their endtime, exclued w_t1 because of low correlation coefficient.
		$sql = "
			SELECT * FROM table_info
			WHERE curr_name LIKE 'w%' AND curr_name != 'w_t1'
		";
	}
	if ($type eq 'm') {			# Select mutant tables and their endtime.
		$sql = "
			SELECT * FROM table_info
			WHERE ko_gene LIKE '$ko_gene'
		";
	}

	my $query = $conn->prepare($sql);
	$query->execute();

	my $merged_exp;
	while (my @row = $query->fetchrow_array) {
		my $table = $row[0];
		my $endtime = $row[7];

		my $cell_exp = &get_trial_data_by_timepoint($table, $endtime);

		$merged_exp = &merge_two_hash($cell_exp,$merged_exp);
	}
	undef $query;

	return $merged_exp;
}


### Obtain one trial data by experiment name and its endtime.
sub get_trial_data_by_timepoint {
	my ($table, $endtime) = @_;
	my $sql = 'SELECT * FROM '.$table.' WHERE time < '.$endtime;
	my $query = $conn->prepare($sql);
	$query->execute;
	
	my %cell_exp;
	while (my @row = $query->fetchrow_array) {
		my @exps = ($row[4]);
		$cell_exp{$row[0]}=\@exps;
	}
	undef $query;
	
	return \%cell_exp;
}


### Merge two hash together.
sub merge_two_hash {
	my ($hash1, $hash2) = @_;			# The hash is in form of $hash{\@array}. Merge two hash means extent @array for the same hash key. Or make new hash key with value \@array from either of two original hash value \@array.
	my %hash;

	foreach my $key (keys %$hash1) {
		unless (defined $hash{$key}) {
			if (defined $$hash2{$key}) {
				my $array1 = $$hash1{$key};
				my $array2 = $$hash2{$key};
				my @array = (@$array1,@$array2);
				$hash{$key} = \@array;
			}else{
				$hash{$key} = $$hash1{$key};
			}
		}
	}

	foreach my $key (keys %$hash2) {
		unless (defined $hash{$key}) {
			$hash{$key} = $$hash2{$key};
		}
	}
	
	return \%hash;
}


### Get knockdown gene list.
sub get_ko_list {
	my $sql = "
		SELECT ko_gene FROM table_info
		WHERE curr_name LIKE 'm%'
		GROUP BY ko_gene
		ORDER BY ko_gene;
	";
	my $query = $conn->prepare($sql);
	$query->execute();

	my @ko_list;
	while ( my @row = $query->fetchrow_array() ) {
		push @ko_list, $row[0];
	}

	undef $query;
	return \@ko_list;
}


### Divide array into subarrays.
sub dividearray {
	my ($array, $subnum) = @_;
	my $length = int ((scalar @$array)/$subnum);
	my @subarray;

	for (1..($subnum-1)){
		my $i = $_;
		my @sub = splice (@$array, 0, $length);
		$subarray[$i-1] = \@sub;
	}
	$subarray[$subnum-1] = $array;

	return @subarray;
}


### Compute significant difference per knockdown gene.
sub compute_sig_per_ko {
	my ($wt, $ko_gene) = @_;
	my $mu = &get_experiments_data('m', $ko_gene);
	my @line_sig;
	my @line_pvalue;

	push @line_sig, $ko_gene;
	push @line_pvalue, $ko_gene;

	foreach my $cn (sort keys %$wt) {
		my ($sig, $pvalue);

		if (defined $$wt{$cn} && defined $$mu{$cn}) {
			($sig, $pvalue) = &wilcoxon_test($$wt{$cn}, $$mu{$cn});

			if ($sig != 0) {
				my $mean_wt = &mean($$wt{$cn});
				my $mean_mu = &mean($$mu{$cn});
				if (($mean_wt > $cutoff && $mean_mu > $cutoff) || ($mean_wt < $cutoff && $mean_mu < $cutoff)) {			# The cutoff condition.
					($sig, $pvalue) = qw/0 1/;
				}
			}

		}elsif(defined $$wt{$cn}){
			($sig, $pvalue) = qw/0_ 1/;			# 0_ means no significant because of miss data in mutant. It can be replaced by 0.
		}

		push @line_sig, $sig;
		push @line_pvalue, $pvalue;
		printf ("%-20s%-10s%s\n", ("=>$cn",$sig,$pvalue));
		printf OUTPUT3 ("%-20s%-10s%s\n", ("=>$cn",$sig,$pvalue));
	}

	return (\@line_sig, \@line_pvalue);
}


### Wilcoxon rank-sum test
sub wilcoxon_test {
	my ($group1, $group2) = @_;
	
	my $wt;
	my $mu;

	my $R = Statistics::R->new() ;
	$R->startR;
	$R->run('options(warn=-1)');

	foreach (@$group1) {
		$wt .= "$_,";
	}
	foreach (@$group2) {
		$mu .= "$_,";
	}
	$wt =~ s/,$//;
	$mu =~ s/,$//;
	
	$R->run('wt = c('.$wt.')');
	$R->run('mu = c('.$mu.')');
	$R->run('reg <- wilcox.test(wt,mu)');
	$R->run('reg$p.value');
	my $pvalue = $R->read;
	$pvalue =~ s/^\[.+\]\s+//;
	if ($pvalue <= 0.05) {
		$R->run('great <- wilcox.test(wt,mu,alternative = \'g\')');
		$R->run('great$p.value');
		my $pvalue_g = $R->read;
		$pvalue_g =~ s/^\[.+\]\s+//;
		if ($pvalue_g <= 0.05){
			$R->stopR;
			return (1, $pvalue_g);			### wt > mu
		}else{
			$R->run('less <- wilcox.test(wt,mu,alternative = \'l\')');
			$R->run('less$p.value');
			my $pvalue_l = $R->read;
			$pvalue_l =~ s/^\[.+\]\s+//;
			$R->stopR;
			return (-1, $pvalue_l);			### wt < mu
		}
	}else{
		$R->stopR;
		return (0, $pvalue);				### wt == mu
	}
}


### Compute mean.
sub mean{
    my($vec)=@_;
    my $sum=0;
    foreach (@$vec) {
        $sum+=$_;
    }
    my $result=$sum/(scalar @$vec);
    return $result;
}

