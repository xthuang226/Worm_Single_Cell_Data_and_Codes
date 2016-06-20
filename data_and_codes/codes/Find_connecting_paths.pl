use strict;
use warnings;
use feature 'say';
use Paths::Graph2;
use Getopt::Long;


#
# Find_connecting_paths.pl
# Reqirements: Perl5
# Author: Xiao-Tai Huang
# Email: xthuang226@gmail.com
# Version: 1.0
# Date: 2016/06/19
#

my $file        = "Find_connecting_paths.pl";
my $program     = uc($&) if $file =~ /^\w+/;
my $version     = "1.0";
my $date        = "2016/06/19";
my $author      = "Xiao-Tai Huang";


### program identification
my $title = "$program $version, $date $author\n";


### option variables
$::opt_help         = 0;
$::opt_startnodes   = 'Y39E4B.1 R10E11.1 C06G3.10 W09C2.1 K02B9.4 T05H4.14 M05B5.5 K08B4.1 C04F1.3 H12C20.3 C38D4.6 T28H11.4 T19E7.2 Y66A7A.8 Y47D3A.12 C24H11.3';
$::opt_endnodes     = 'F11C1.6';
$::opt_length       = '5';
$::opt_PPI          = '../data/background_network/PPIs.csv';
$::opt_PDI          = '../data/background_network/PDIs.csv';
$::opt_output       = '../results/connecting_paths/';


my $usage = <<"END_OF_USAGE";
${title}Syntax:   \L$program\E [options]
Function: Find paths connecting the cause gene and the effect gene in the background network.
Options:                                                     (defaults:)
  -h, --help           print usage
  -s, --startnodes     set startnodes of the path            ($::opt_startnodes)
  -e, --endnodes       set endnodes of the path              ($::opt_endnodes)
  -l, --length         set maximum path length               ($::opt_length)
  -o, --output         set output folder for the results     ($::opt_output)
  --PPI                set FILE path of PPI                  ($::opt_PPI)
  --PDI                set FILE path of PDI                  ($::opt_PDI)
Examples:
  \L$program\E
  \L$program\E -l=5 -s=Y39E4B.1 -e=F11C1.6 -PPI=../data/background_network/PPIs.csv -PDI=../data/background_network/PDIs.csv
  \L$program\E -l=5 -s='Y39E4B.1 R10E11.1 C06G3.10 W09C2.1 K02B9.4' -e='F11C1.6 C38D4.6 T28H11.4'
END_OF_USAGE



### process options
my @OrgArgv = @ARGV;
GetOptions(
  "help|h!",
  "startnodes|s=s",
  "endnodes|e=s",
  "length|l=i",
  "output|o=s",
  "PPI=s",
  "PDI=s"
) or die $usage;
!$::opt_help or die $usage;


### Set parameters
my ($max_path_length, $startnodes, $endnodes, $PPI_file, $PDI_file, $output) = ($::opt_length, $::opt_startnodes, $::opt_endnodes, $::opt_PPI, $::opt_PDI, $::opt_output);

### Program starting time
my $time_start = time;


my $output_tmp = $output;
my $output_whole;

while ($output_tmp =~ s/^([\.\/]*[^\.\/]+\/?)//) {
	$output_whole .= $&;
	mkdir $output_whole;
}


my @cause_list = split(/[\s,]+/, $startnodes);
my @effect_list = split(/[\s,]+/, $endnodes);


### Load Data
my ($orf2id,$id2orf,$ppgraph,$pdgraph) = &load_data($PPI_file, $PDI_file);

my $skgraph = &generate_skgraph($ppgraph, $pdgraph);

my $graph_run = $skgraph;

say "\n### Parameter: Max Path Length = $max_path_length\n";

my $total_path_number;		### Finded pathes number


foreach my $g_effect (@effect_list) {
	### Check effect genes in the background network
	if (defined $$orf2id{$g_effect}) {
		foreach my $g_cause (@cause_list) {
			### Check cause and effect genes in the background network
			if (defined $$orf2id{$g_cause}) {
				### Program starting time
				my $finding_start = time;

				$total_path_number = 0;

				open (OUTPUT,">$output/path_result_$max_path_length\_$g_cause\_to_$g_effect.csv")||die "Error!\n$!";

				my $start = $$orf2id{$g_cause};
				my $end = $$orf2id{$g_effect};

				say "$g_cause->...->$g_effect";
				&find_path($start, $end, $graph_run, $max_path_length-1);

				### Print totally finded pathes number
				say "Totally find $total_path_number pathes";

				### Ending time of each finding
				my $duration = time-$finding_start;
				my $print_time = &duration_time($duration);
				say "Finding duration time: $print_time";
			}else{
				say "The $g_cause is not defined in the background network";
			}
		}
	}else{
		say "The $g_effect is not defined in the background network";
	}
}



### Program duration time
my $duration = time-$time_start;
my $print_time = &duration_time($duration);
say "Program duration time: $print_time";


close OUTPUT;












sub load_data {
	my @infile = @_;

	my %name2id;
	my %id2name;
	my %ppgraph;
	my %pdgraph;

	my $i = 0;
	my $count = 0;
	foreach (@infile) {
		open (INPUT,$_)||die "Error!\n$!";
		foreach (<INPUT>){
			chomp;
			if (/^([^,\s]+)[,\s]+([^,\s]+)[,\s]+([^,\s]+)/){
				my @row = ($1, $2, $3);

				for (0..1) {
					unless (defined $name2id{$row[$_]}) {
						$name2id{$row[$_]} = $i;
						$id2name{$i} = $row[$_];
						$i++;
					}
				}

				if ($count == 0) {
					for (0..1) {
						my ($a, $b);

						if ($_ == 0) {
							$a = $row[0];
							$b = $row[1];
						}else{
							$a = $row[1];
							$b = $row[0];
						}

						if (defined $ppgraph{$name2id{$a}}) {
							my $b_nodes = $ppgraph{$name2id{$a}};
							$$b_nodes{$name2id{$b}} = $row[2];
							$ppgraph{$name2id{$a}} = $b_nodes;
						}else{
							my %b_nodes;
							$b_nodes{$name2id{$b}} = $row[2];
							$ppgraph{$name2id{$a}} = \%b_nodes;
						}
					}
				}
				if ($count == 1) {
					my ($a, $b) = ($row[0], $row[1]);

					if (defined $pdgraph{$name2id{$a}}) {
						my $b_nodes = $pdgraph{$name2id{$a}};
						$$b_nodes{$name2id{$b}} = $row[2];
						$pdgraph{$name2id{$a}} = $b_nodes;
					}else{
						my %b_nodes;
						$b_nodes{$name2id{$b}} = $row[2];
						$pdgraph{$name2id{$a}} = \%b_nodes;
					}
				}
			}
		}
		close INPUT;
		$count++;
	}
	return (\%name2id,\%id2name, \%ppgraph, \%pdgraph);
}


sub read_nodes_on_path {
	my ($file_path) = @_;
	my %nodes_on_path;
	open (INPUT,$file_path)||die "Error!\n$!";
	foreach (<INPUT>){
		chomp;
		$nodes_on_path{$_} = 1;
	}
	close INPUT;
	return %nodes_on_path;
}


sub generate_skgraph {
	my ($ppgraph, $pdgraph) = @_;

	my %skgraph;
	foreach (keys %$ppgraph) {
		$skgraph{$_} = &copy_interaction($$ppgraph{$_});
	}

	foreach (keys %$pdgraph) {
		my $a = $_;

		if (defined $$ppgraph{$a}) {
			my $b_nodes = $$pdgraph{$a};
			foreach (keys %$b_nodes) {
				my $b = $_;
				if (defined $$ppgraph{$b}) {
					%skgraph = &generate_direct_edge_in_exist_graph($a, $b, \%skgraph);
					%skgraph = &update_edge_value($a, $b, \%skgraph, $pdgraph);
				}else{
					my $b_nodes2 = $skgraph{$a};
					$$b_nodes2{$b} = $$b_nodes{$b};
					$skgraph{$a} = $b_nodes2;
				}
			}
		}else{
			$skgraph{$a} = &copy_interaction($$pdgraph{$a});
		}
	}

	return \%skgraph;
}


sub copy_interaction {
	my ($b_nodes) = @_;
	my %b_nodes;
	foreach (keys %$b_nodes) {
		my $b = $_;
		$b_nodes{$b} = $$b_nodes{$b};
	}
	return \%b_nodes;
}


sub generate_direct_edge_in_exist_graph {
	my ($a, $b, $graph) = @_;

	my $flag_1 = 0;				### edge existance from $a to $b
	my $flag_2 = 0;				### edge existance from $b to $a


	my $b_nodes_1 = $$graph{$a};
	$flag_1 = 1 if defined $$b_nodes_1{$b};

	my $b_nodes_2 = $$graph{$b};
	$flag_2 = 1 if defined $$b_nodes_2{$a};

	if ($flag_1 == 1 && $flag_2 == 1) {
		delete $$b_nodes_2{$a};
		$$graph{$b} = $b_nodes_2;
	}else{
		$$b_nodes_1{$b} = 1;
		$$graph{$a} = $b_nodes_1;
	}


	return %$graph;
}


sub update_edge_value {
	my ($a, $b, $graph1, $graph2) = @_;
	my $b_nodes = $$graph1{$a};
	my $b_nodes2 = $$graph2{$a};
	foreach (keys %$b_nodes) {
		if ($_ == $b) {
			$$b_nodes{$b} = $$b_nodes2{$b};
		}
	}
	$$graph1{$a} = $b_nodes;
	return %$graph1;
}


sub refine_graph {
	my ($cutoff, $graph) = @_;
	foreach (keys %$graph) {
		my $a = $_;
		my $b_nodes = $$graph{$a};
		foreach (keys %$b_nodes) {
			my $b = $_;
			if ($$b_nodes{$b} < $cutoff) {
				delete $$b_nodes{$b};
				if (defined $$graph{$b}) {
					my $b_nodes2 = $$graph{$b};
					delete $$b_nodes2{$a} if $$b_nodes2{$a} < $cutoff;
				}
			}
		}
	}
	return $graph;
}


sub find_path {
	my ($start, $end, $graph, $max_path_length)=@_;
	my $obj = Paths::Graph->new(-origin=>$start,-destiny=>$end,-graph=>$graph,-sub=>\&get_paths);
	$obj->free_path_event($start,$max_path_length);
}


sub get_paths {
    my ($self , @nodes) = @_;
	
	my @out_nodes;
	foreach (@nodes){
		my $node_orf = $$id2orf{$_};
		push @out_nodes, $node_orf;
	}

	my @out_nodes2 = @out_nodes;

	for (0..($#nodes-1)) {
		my $flag = &edge_in_graph($nodes[$_], $nodes[$_+1], $pdgraph);
		if ($flag == 1) {
			$out_nodes[$_] .= ' $a$';
			$out_nodes[$_+1] .= ' $b$';			
		}
	}

	my @edge_score = &get_edges_scores(@nodes);
	
	#### new print
	print OUTPUT join(",",@out_nodes2);
	for (1..(($max_path_length+1)-(scalar @nodes) + 1)){
		print OUTPUT ",";
	}
	say OUTPUT join(",",map{$_}@edge_score);

	$total_path_number++;
}


sub get_edges_scores {
	my @nodes = @_;
	my @nodes2 = @nodes;
	my @edge_score;
	for (1..(scalar @nodes - 1)){
		push @edge_score, &edge_score($nodes2[0], $nodes2[1], $graph_run);
		shift @nodes2;
	}
	return @edge_score;
}


sub edge_score {
	my ($a, $b, $graph) = @_;
	my $b_nodes = $$graph{$a};
	foreach (keys %$b_nodes) {
		if ($_ == $b) {
			return $$b_nodes{$b};
		}
	}
}


sub sum{
    my($vec)=@_;
    my $result=0;
    foreach(@$vec){
        $result+=$_;
        }
    return $result;
}
sub multiply{
    my($vec)=@_;
    my $result=1;
    foreach(@$vec){
        $result*=$_;
    }
    return $result;
}
sub mean{
    my($vec)=@_;
    my $sum=0;
    foreach(@$vec){
        $sum+=$_;
    }
    my $result=$sum/(scalar @$vec);
    return $result;
}
sub var{
    my($vec)=@_;
    my $m=&mean($vec);
    my $length=scalar @$vec;
    my @squre=map(($_-$m)*($_-$m),@$vec);
    my $result=1/($length-1)*&sum(\@squre);
    return $result;
}
sub sd {
    my($vec)=@_;
    my $result=sqrt(&var($vec));
    return $result;
}


sub duration_time {
	my ($duration) = @_;
	my $print_time;
	if ($duration<60){
		$print_time .= "$duration s\n\n";
	}elsif ($duration<3600){
		$print_time .= int($duration/60).'m '.$duration%(60).'s';
	}else{
		$print_time .= int($duration/3600).'h '.int($duration/60)%(60).'m '.$duration%(60).'s';
	}
	return $print_time;
}


sub graph_nodes {
	my ($graph) = @_;
	my %nodes;
	foreach (keys %$graph) {
		$nodes{$_} = 1;
		my $b_nodes = $$graph{$_};
		foreach (keys %$b_nodes) {
			$nodes{$_} = 1;
		}
	}
	return scalar keys %nodes;
}


sub graph_edges {
	my ($graph) = @_;
	my %exist;

	foreach (keys %$graph) {
		my $a = $_;
		my $b_nodes = $$graph{$a};
		foreach (keys %$b_nodes) {
			my $b = $_;

			$exist{$a."_".$b} = 1 unless defined $exist{$b."_".$a};
		}
	}

	return scalar keys %exist;
}


sub edge_in_graph {
	my ($a, $b, $graph) = @_;
	my $flag = 0;
	if (defined $$graph{$a}) {
		my $b_nodes = $$graph{$a};
		if (defined $$b_nodes{$b}) {
			$flag = 1;
			return $flag;
		}
		return $flag;
	}else{
		return $flag;
	}
}

