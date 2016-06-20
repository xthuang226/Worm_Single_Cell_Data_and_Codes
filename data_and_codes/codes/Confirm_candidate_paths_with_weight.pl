use strict;
use warnings;
use feature 'say';
use Getopt::Long;



#
# Confirm_candidate_paths_with_weight.pl
# Reqirements: Perl5
# Author: Xiao-Tai Huang
# Email: xthuang226@gmail.com
# Version: 1.0
# Date: 2016/06/19
#


my $file        = "Confirm_candidate_paths_with_weight.pl";
my $program     = uc($&) if $file =~ /^\w+/;
my $version     = "1.0";
my $date        = "2016/06/19";
my $author      = "Xiao-Tai Huang";


### program identification
my $title = "$program $version, $date $author\n";


### option variables
$::opt_help         = 0;
$::opt_input        = '../results/connecting_paths/';
$::opt_output       = '../results/candidate_paths/';
$::opt_topn         = 5;


my $usage = <<"END_OF_USAGE";
${title}Syntax:   \L$program\E [options]
Function: Confirm candidate paths with weight.
Options:                                                     (defaults:)
  -h, --help      print usage
  -i, --input     set input folder of the connecting paths   ($::opt_input)
  -o, --output    set output folder for the results          ($::opt_output)
  -n, --topn      set top n-th of all candidate paths        ($::opt_topn)
Examples:
  \L$program\E
  \L$program\E -i=../results/connecting_paths/ -o=../results/candidate_paths/ -n=5
END_OF_USAGE



### process options
my @OrgArgv = @ARGV;
GetOptions(
  "help|h!",
  "input|i=s",
  "output|o=s",
  "topn|n=i"
) or die $usage;
!$::opt_help or die $usage;


### Set parameters
my ($input, $output, $topn) = ($::opt_input, $::opt_output, $::opt_topn);

my $output_tmp = $output;
my $output_whole;

while ($output_tmp =~ s/^([\.\/]*[^\.\/]+\/?)//) {
	$output_whole .= $&;
	mkdir $output_whole;
}


opendir (DIR, $input) ||die "Error!\n$!";

my @infiles = grep(/.csv$/, readdir DIR);

closedir DIR;



my $e1 = 1.192093e-007;


foreach (@infiles) {
	say;
	s/\.csv$//;
	&main($_);
}




sub main {
	my ($infile) = @_;

	my $max_length = 0;

	my %edge_with_weight;

	open (INPUT, "$input/$infile.csv") || die "Error!\n$!";
	foreach (<INPUT>) {
		chomp;
		my @elements = split(/,/,$_);
		my @nodes;
		my @score;
		foreach (@elements) {
			if (/./) {
				if (/[A-Z]/) {
					push @nodes, $_;
				}else{
					push @score, $_;
				}
			}
		}
		$max_length = scalar @nodes if scalar @nodes > $max_length;
		for (0..$#score) {
			$score[$_] = $e1 if $score[$_] == 0;
			my $node1 = $nodes[$_];
			my $node2 = $nodes[$_+1];

			if (defined $edge_with_weight{$node1}) {
				my $weight = $edge_with_weight{$node1};
				$$weight{$node2} = $score[$_];
			}else{
				my %weight;
				$weight{$node2} = $score[$_];
				$edge_with_weight{$node1} = \%weight;
			}
		}
	}
	close INPUT;

	say $max_length;

	foreach (keys %edge_with_weight) {
		my $node1 = $_;
		my $weight = $edge_with_weight{$node1};
		my $sum = 0;
		foreach (keys %$weight) {
			$sum += $$weight{$_};
		}
		foreach (keys %$weight) {
			my $node2 = $_;
			$$weight{$node2} *= 1/$sum;
		}
	}

	my %out;

	open (INPUT, "$input/$infile.csv") || die "Error!\n$!";
	foreach (<INPUT>) {
		chomp;
		my @elements = split(/,/,$_);
		my @nodes;
		my @score;
		foreach (@elements) {
			if (/./) {
				if (/[A-Z]/) {
					push @nodes, $_;
				}else{
					push @score, $_;
				}
			}
		}
		my @score_new;
		for (0..$#score) {
			my $node1 = $nodes[$_];
			my $node2 = $nodes[$_+1];
			my $weight = $edge_with_weight{$node1};
			push @score_new, $$weight{$node2};
		}
		my $multiply = 1;
		foreach (@score_new) {
			$multiply *= $_;
		}

		my $out_key = join(',', @nodes).","x(($max_length - scalar @nodes) + 1).join(',', @score_new).","x($max_length-1 - scalar @score_new);
		$out{$out_key} = $multiply;
	}
	close INPUT;

	my $i = 1;
	open (OUTPUT, ">$output/$infile\_weight.csv") || die "Error!\n$!";
	foreach (sort {$out{$b} <=> $out{$a}} keys %out) {
		if ($i>$topn) {
			last;
		}else{
			say OUTPUT "$_,$out{$_}";
		}
		$i++;
	}
	close OUTPUT;
}




