use strict;
use warnings;
use feature 'say';
use Storable qw(dclone);
use Cwd;
use Getopt::Long;



#
# Infer_pathway.pl
# Reqirements: Perl5
# Author: Xiao-Tai Huang
# Email: xthuang226@gmail.com
# Version: 1.0
# Date: 2016/06/19
#


my $file        = "Infer_pathway.pl";
my $program     = uc($&) if $file =~ /^\w+/;
my $version     = "1.0";
my $date        = "2016/06/19";
my $author      = "Xiao-Tai Huang";




### program identification
my $title = "$program $version, $date $author\n";


### option variables
$::opt_help           = 0;
$::opt_incandidates   = '../results/candidate_paths/';
$::opt_inothers       = '../data/pathway_infer_include/';
$::opt_output         = '../results/inferred_pathways/';
$::opt_founders       = 'ABarp ABpla ABpra Caa Cpa';


my $usage = <<"END_OF_USAGE";
${title}Syntax:   \L$program\E [options]
Function: Infer pathways in founder cells based on candidate paths.
Options:                                                       (defaults:)
  -h, --help        print usage
  --incandidates    set candidate paths file folder            ($::opt_incandidates)
  --inothers        set other input folder                     ($::opt_inothers)
  -o, --output      set output folder                          ($::opt_output)
  -f, --founders    set founders to be inferred pathways       ($::opt_founders)
Examples:
  \L$program\E
  \L$program\E -f='ABarp ABpla ABpra Caa Cpa'
END_OF_USAGE



### process options
my @OrgArgv = @ARGV;
GetOptions(
  "help|h!",
  "incandidates=s",
  "inothers=s",
  "output|o=s",
  "founders|f=s"
) or die $usage;
!$::opt_help or die $usage;



### Set parameters
my ($incandidates, $inothers, $output, $founders) = ($::opt_incandidates, $::opt_inothers, $::opt_output, $::opt_founders);





## Program starting time
my $time_start = time;
my $tp = $time_start;


opendir (DIR, $incandidates) || die "Cannot open dir : $!";
my @candidate_paths = readdir(DIR);
close DIR;


######## global setting
my $e1 = 1.192093e-007;
my $e2 = 1.192093e-007;


my %var_list_hash = (
	'k' => 2,
	'x' => 2,
	's' => 2,
	'd' => 2
);
my @var_list_ary = qw/k x s d/;


######## load K S D data
my $knockout = &load_K_S_D_data('k', $inothers);
my $regulation = &load_K_S_D_data('s', $inothers);
my $directions = &load_K_S_D_data('d', $inothers);


my @founders = split (/[,\s]+/, $founders);
my $top_line = 3;

foreach (@founders) {
	my $founder = $_;

	######## add ko data in each founder
	my $causes;
	($knockout,$causes) = &add_ko_per_founder($founder,$knockout);


	my $output_tmp = "$output/$founder";
	my $output_whole;

	while ($output_tmp =~ s/^([\.\/]*[^\.\/]+\/?)//) {
		$output_whole .= $&;
		mkdir $output_whole;
	}

	&main ($founder,$knockout,$regulation,$directions,$causes,\@candidate_paths,\@var_list_ary,\%var_list_hash,$e1,$e2,$top_line);
}



### Program ending time
my $print_time = &duration_time(time-$time_start);
say $print_time;














sub main {
	my ($founder,$knockout,$regulation,$directions,$causes,$candidate_paths,$var_list_ary,$var_list_hash,$e1,$e2,$top_line) = @_;

	foreach (@$causes) {
		my $cause = $_;
		foreach (@$candidate_paths) {
			my $cdpathfile = $_;
			$cdpathfile =~ s/\.[^\.]+$//;
			if ($cdpathfile =~ /$cause/) {
				open (INPUT,"$incandidates/$cdpathfile.csv")||die "Error!\n$!";
				my $i = 0;
				foreach (<INPUT>) {
					$i++;
					if ($i <= $top_line) {
						chomp;
						my @line = split(/,/,$_);

						open (OUTPUT1, ">>$output/$founder/log_$cdpathfile")||die "Error!\n$!";

						######## get path and its interaction info
						my ($interactions_on_path, $cdpath, $var_num) = &get_path_and_interaction_info(\@line);
						my ($path, $name_a_to_m, $name_m_to_a) = &path_name($cdpath);
						say "$founder";
						say "\n===> ".join(',',@$cdpath);
						say OUTPUT1 "\n===> ".join(',',@$cdpath);

						######## load KO info on the path
						my $ko_pairs = &load_KO_info_on_path($knockout,$name_a_to_m);
						print '->';
						print OUTPUT1 '->';
						foreach (sort keys %$ko_pairs) {
							print " $_ $$ko_pairs{$_};";
							print OUTPUT1 " $_ $$ko_pairs{$_};";
						}
						print "\n-> var: $var_num;  ko: ".(scalar keys %$ko_pairs).";";
						say "  table:".(2**(3*$var_num))*(2**(scalar keys %$ko_pairs));

						######## initialize all message
						my ($msg_var_to_func, $msg_func_to_var) = &init_all_msg($var_list_ary,$var_list_hash,$cdpath,$ko_pairs,$e1);

						######## associate var_to_func message
						foreach (keys %$msg_var_to_func) {
							my $msg = $$msg_var_to_func{$_};
							$msg = &associate_message($msg, $name_a_to_m, $name_m_to_a, $knockout, 'k', $e1) if $_ eq 'k';
							$msg = &associate_message($msg, $name_a_to_m, $name_m_to_a, $regulation, 's', $e1) if $_ eq 's';
							$msg = &associate_message($msg, $name_a_to_m, $name_m_to_a, $interactions_on_path, 'x', $e1) if $_ eq 'x';
							$msg = &associate_message($msg, $name_a_to_m, $name_m_to_a, $directions, 'd', $e1) if $_ eq 'd';
						}

						my $msg_var_to_func_new = dclone($msg_var_to_func);
						my $msg_func_to_var_new = dclone($msg_var_to_func);


						######## compute func_to_var from var_to_func
						say "-> generate conf tabel";
						my $conf_table = &configuration($var_list_ary,$var_list_hash,$ko_pairs,$var_num);

						say "-> update msg";
						&update_msg_func_to_var($msg_func_to_var_new,$conf_table,$var_num,$var_list_ary,$var_list_hash,$ko_pairs,$msg_var_to_func,$e2);

						print OUTPUT1 "\n";
						&show_message($msg_func_to_var_new);

						my ($path_value,$path_value2,$conf,$conf_value,$conf_value2) = &max_conf_and_value_msg($msg_func_to_var_new,$var_num,$line[$#line]);
						say "\t\t\tscore: $path_value\n";
						say "\t\t\tscore: $path_value2\n";
						say OUTPUT1 "score: $path_value\n";
						say OUTPUT1 "score: $path_value2\n";
						say OUTPUT1 join(',',@$conf);
						say OUTPUT1 join(',',@$conf_value);
						say OUTPUT1 join(',',@$conf_value2);

						open (OUTPUT2, ">>$output/$founder/$cdpathfile")||die "Error!\n$!";
						say OUTPUT2 join("\t", (join(',',@$cdpath), $path_value, join(',',@$conf), join(',',@$conf_value)));
						say OUTPUT2 join("\t", ("->",join(',',@$cdpath), $path_value2, join(',',@$conf), join(',',@$conf_value2)));
						close OUTPUT2;

						say &duration_time(time-$tp);
						$tp = time;
					}else{
						last;
					}
				}
				close INPUT;
			}
		}
	}
}



sub max_conf_and_value_msg {
	my ($msg_total,$var_num,$line_last) = @_;
	my @conf;
	my @conf_value;
	foreach (sort keys %$msg_total) {
		my $type = $_;
		next if $type eq 'k';
		my $msg = $$msg_total{$_};
		foreach (sort keys %$msg) {
			my $belief = $$msg{$_};
			map{s/^(.+\.\d{4})\d+$/$1/}@$belief;
			my ($var_max,$pos) = &max_array($$msg{$_});
			push @conf_value, $var_max;
			push @conf,join('\\', @$pos) if $type eq 's';
		}
	}
	my $path_value = &sum_array(\@conf_value);

	my @conf_value2;
	push @conf_value2, log($line_last);
	foreach (($var_num-1)..($var_num+2)) {
		push @conf_value2,$conf_value[$_];
	}
	my $path_value2 = &sum_array(\@conf_value2);

	return ($path_value,$path_value2,\@conf,\@conf_value,\@conf_value2);
}






##################
sub load_K_S_D_data {
	my ($type, $inothers) = @_;
	my %result;
	open (INPUT,"$inothers/Cause_effect_paris.txt")||die "Error!\n$!" if $type eq 'k';
	open (INPUT,"$inothers/Regulation.txt")||die "Error!\n$!" if $type eq 's';
	open (INPUT,"$inothers/Direction.txt")||die "Error!\n$!" if $type eq 'd';
	foreach (<INPUT>) {
		chomp;
		my @elements = split(/\s+/, $_);
		my @belief;
		if ($type eq 'd') {
			push @belief,$elements[2];
		}else{
			push @belief,($elements[2],$elements[3]);
		}
		$result{$elements[0].'_'.$elements[1]} = \@belief;
	}
	close INPUT;
	return \%result;
}
##################





##################
sub add_ko_per_founder {
	my ($founder,$knockout) = @_;
	my @cause;
	open (INPUT,"$inothers/KO_$founder.txt")||die "Error!\n$!";
	foreach (<INPUT>) {
		chomp;
		my @elements = split(/\s+/, $_);
		my @sign_belief = ($elements[2],$elements[3]);
		$$knockout{$elements[0].'_'.$elements[1]} = \@sign_belief;
		push @cause, $elements[0];
	}
	close INPUT;
	return ($knockout, \@cause);
}
##################





##################
sub get_path_and_interaction_info {
	my ($line) = @_;
	my %interactions_on_path;
	my @cdpath;

	my $i = 0;
	for (0..(scalar @$line/2 - 1)) {
		if ($$line[$_] =~ /./) {
			push @cdpath, $$line[$i];
		}
		$i++;
	}

	my $j = 1;
	for ($i..($i + (scalar @cdpath) - 2)) {
		my @belief = ($$line[$_]);
		$interactions_on_path{'v'.$j.'_v'.($j+1)} = \@belief;
		$j++;
	}

	my $var_num = (scalar @cdpath) - 1;
	return (\%interactions_on_path,\@cdpath,$var_num);
}
##################





##################
sub path_name {
	my ($old_path) = @_;
	my @new_path;
	my %name_a_to_m;
	my %name_m_to_a;
	my $i = 1;
	foreach my $old_name (@$old_path) {
		push @new_path, "v$i";
		$name_a_to_m{$old_name} = "v$i";
		$name_m_to_a{"v$i"} = $old_name;
		$i++;
	}
	return (\@new_path, \%name_a_to_m, \%name_m_to_a);
}
##################





##################
sub load_KO_info_on_path {
	my ($knockout,$name_a_to_m) = @_;
	my %ko_pairs;
	foreach (sort keys %$knockout) {
		my @nodes = split(/_+/, $_);
		my $flag = 0;
		foreach (@nodes) {
			$flag = 1 unless defined $$name_a_to_m{$_};
		}
		next if $flag == 1;
		my $s = $$name_a_to_m{$nodes[0]};
		my $e = $$name_a_to_m{$nodes[1]};
		foreach ($s,$e) {
			s/v//;
		}
		next if $s > $e;
		my $effct;
		my $info = $$knockout{$nodes[0].'_'.$nodes[1]};
		if ($$info[0] eq '-') {
			$effct = 0;
		}
		if ($$info[0] eq '+') {
			$effct = 1;
		}
		$ko_pairs{$$name_a_to_m{$nodes[0]}.'_'.$$name_a_to_m{$nodes[1]}} = $effct;
	}
	return \%ko_pairs;
}
##################





##################
sub init_all_msg {
	my ($var_list_ary,$var_list_hash,$cdpath,$ko_pairs,$e1) = @_;
	my $msg_var_to_func;
	foreach (@$var_list_ary) {
		my $var = $_;
		my $msg;
		if ($_ eq 'k') {
			foreach (sort keys %$ko_pairs) {
				my @belief;
				for (1..$var_list_hash{'k'}) {
					push @belief, $e1;
				}
				$$msg{$_} = \@belief;
			}
		}else{
			$msg = &init_msg(1,(scalar @$cdpath),$var_list_hash{$_},$e1);
		}
		$$msg_var_to_func{$var} = $msg;
	}
	my $msg_func_to_var = dclone($msg_var_to_func);
	return($msg_var_to_func, $msg_func_to_var);
}
##################





##################
sub init_msg {
	my ($s, $e, $col, $e1) = @_;
	my %message;
	my $i = $s;
	for ($s..$e-1) {
		my @belief;
		for (1..$col) {
			push @belief, $e1;
		}
		$message{"v$i\_v".($i+1)} = \@belief;
		$i++;
	}
	return \%message;
}
##################





##################
sub associate_message {
	my ($msg, $name_a_to_m, $name_m_to_a, $data, $flag, $e1) = @_;

	foreach (keys %$msg) {
		my $belief = $$msg{$_};

		### convert m name to a name
		my @nodes = split(/_+/, $_);
		my @new_nodes;
		my @new_nodes_rev;
		map{push @new_nodes, $$name_m_to_a{$_}}@nodes;
		map{push @new_nodes_rev, $$name_m_to_a{$_}}($nodes[1], $nodes[0]);
		my $key = join('_', @new_nodes);
		my $key_rev = join('_', @new_nodes_rev);

		if ($flag eq 'k' || $flag eq 's') {
			my $data_ary = $$data{$key};
			if (defined $data_ary) {
				if ($$data_ary[0] eq '-') {
					$$belief[0] = $$data_ary[1];
					$$belief[1] = 1-$$belief[0];
				}
				if ($$data_ary[0] eq '+') {
					$$belief[1] = $$data_ary[1];
					$$belief[0] = 1-$$belief[1];
				}
			}
		}elsif ($flag eq 'x') {
			my $data_ary = $$data{$key};
			if (defined $data_ary) {
				$$belief[1] = $$data_ary[0];
				$$belief[0] = 1-$$belief[1];
			}
		}elsif ($flag eq 'd') {
			if (defined $$data{$key_rev}) {
				my $ary = $$data{$key_rev};
				$$belief[0] = $$ary[0];
				$$belief[1] = 1-$$belief[0];
			}
			if (defined $$data{$key}) {
				my $ary = $$data{$key};
				$$belief[1] = $$ary[0];
				$$belief[0] = 1-$$belief[1];
			}
		}

		$$belief[0] = $e1 if $$belief[0] == 0;
		$$belief[1] = $e1 if $$belief[1] == 0;
		$$msg{$_} = $belief;
	}

	return ($msg);
}
##################





##################
sub configuration {
	my ($var_list_ary, $var_list_hash, $ko_pairs, $var_num) = @_;
	my $conf;
	my $flag = 0;
	foreach (@$var_list_ary) {
		my $var_type = $$var_list_hash{$_};
		my $var_num_new;
		if ($_ eq 'k') {
			$var_num_new = scalar keys %$ko_pairs;
		}else{
			$var_num_new = $var_num;
		}
		my $conf_sub = &gen_conf($var_num_new,$var_type);

		if ($flag == 0) {
			$conf = $conf_sub;
		}else{
			$conf = &multiply_array($conf,$conf_sub);
		}
		$flag = 1;
	}
	return $conf;
}


sub gen_conf {
	my ($var_num,$var_type) = @_;
	my @table;
	for (0..$var_type**$var_num-1) {
		my @result = &decimal_to_any($_,$var_type);
		my $value = sprintf "%0".$var_num."s",join('',@result);
		push @table, $value;
	}
	return \@table;
}


sub decimal_to_any {
	my ($number, $system) = @_;
	my $temp = $number;
	my @result;

	while () {
		my $head = int($temp/$system);

		if ($head < $system) {
			push @result, $temp % $system;
			push @result, $head if $head != 0;
			last;
		}else{
			push @result, $temp % $system;
		}
		$temp = $head;
	}

	my @result2;
	for (0..$#result) {
		push @result2, $result[$#result-$_];
	}
	return @result2;
}


sub multiply_array {
	my ($ary1, $ary2) = @_;
	my @result;
	foreach (@$ary1) {
		my $value1 = $_;
		foreach (@$ary2) {
			my $value2 = $_;
			push @result, $value1.$value2;
		}
	}
	return \@result;
}
##################





##################
sub update_msg_func_to_var {
	my ($msg_func_to_var_new,$conf_table,$var_num,$var_list_ary,$var_list_hash,$ko_pairs,$msg_var_to_func,$e2) = @_;
	say "  -> get tabel info";
	my ($table_with_info, $table_with_func) = &get_table($conf_table,$var_num,$var_list_ary,$ko_pairs,$msg_var_to_func,$e2);

	foreach (keys %$msg_func_to_var_new) {
		my $var = $_;
		if ($var ne 'k') {
			say "    -> update $var";
			my $var_type = $var_list_hash{$var};
			my $msg = $$msg_func_to_var_new{$var};
			foreach (sort keys %$msg) {
				my $edge_name = $_;
				for (0..$var_type-1) {
					my $sta = $_;
					&update_element_func_to_var($var, $edge_name, $sta, $msg_func_to_var_new, $table_with_info, $table_with_func, $var_num, $ko_pairs);
				}
			}
		}
	}
}


sub get_table {
	my ($conf_table,$var_num,$var_list_ary,$ko_pairs,$msg_total,$e2) = @_;
	my %table_with_info;
	my %table_with_func;
	foreach (@$conf_table) {
		my $conf = $_;
		my $func_value = &get_func_value($conf,$var_num,$var_list_ary,$ko_pairs,$e2);
		$table_with_func{$conf} = $func_value;
		my $belief_value = &get_belief_by_config($conf,$var_num,$msg_total,$var_list_ary);
		$table_with_info{$conf} = $belief_value;
	}
	return (\%table_with_info,\%table_with_func);
}


sub get_table_not_zero {
	my ($conf_table,$var_num,$var_list_ary,$ko_pairs,$msg_total,$e2) = @_;
	my %table_with_info;
	my %table_with_func;
	foreach (@$conf_table) {
		my $conf = $_;
		my $func_value = &get_func_value($conf,$var_num,$var_list_ary,$ko_pairs,$e2);
		unless ($func_value == 0) {
			$table_with_func{$conf} = $func_value;
			my $belief_value = &get_belief_by_config($conf,$var_num,$msg_total,$var_list_ary);
			$table_with_info{$conf} = $belief_value;
		}
	}
	return (\%table_with_info,\%table_with_func);
}


sub get_func_value {
	my ($conf, $var_num, $var_list_ary, $ko_pairs, $e2) = @_;
	my $indicator_X = 1;
	my $indicator_S = 1;
	my $indicator_D = 1;
	my @values = split (//, $conf);
	my %effct_on_path;
	foreach (@$var_list_ary) {
		if ($_ eq 'k') {
			foreach (sort keys %$ko_pairs) {
				$effct_on_path{$_} = 1 if $values[0] == 1;
				$effct_on_path{$_} = -1 if $values[0] == 0;
				shift @values;
			}
		}elsif ($_ eq 'x') {
			for (1..$var_num) {
				$indicator_X *= $values[0];
				shift @values;
			}
		}elsif ($_ eq 's') {
			my @s_conf;
			for (1..$var_num) {
				push @s_conf, shift @values;
			}
			foreach (sort keys %effct_on_path) {
				my $indicator_S_sub = 1;
				if (/^v(\d+)_v(\d+)$/) {
					my ($s,$e) = ($1,$2);
					for ($1..$2-1) {
						$indicator_S_sub *= 1 if $s_conf[$_-1] == 1;
						$indicator_S_sub *= -1 if $s_conf[$_-1] == 0;
					}
				}
				if ($indicator_S_sub == $effct_on_path{$_}) {
					$indicator_S *= 1;
				}else{
					$indicator_S *= 0;
				}
			}
		}elsif ($_ eq 'd') {
			for (1..$var_num) {
				$indicator_D *= $values[0];
				shift @values;
			}
		}
	}
	return (1 - $e2) * $indicator_X*$indicator_S*$indicator_D + $e2;
}


sub get_belief_by_config {
	my ($conf, $var_num, $msg_total, $var_list_ary) = @_;
	my @values = split (//, $conf);
	my @belief_value;
	foreach (@$var_list_ary) {
		my $msg = $$msg_total{$_};

		foreach (sort keys %$msg) {
			my $belief = $$msg{$_};
			my $value = shift @values;
			push @belief_value, $$belief[$value];
		}
	}
	return \@belief_value;
}


sub update_element_func_to_var {
	my ($var, $edge_name, $sta, $msg_total, $table_with_info, $table_with_func, $var_num, $ko_pairs) = @_;


	my $msg = $$msg_total{$var};
	my $loc = &loc_var_pos($var,$edge_name,$var_num,$ko_pairs);

	my ($table_with_info_part,$table_with_func_part) = &extract_table_part($loc,$sta,$table_with_info,$table_with_func);
	my @comp_func_belief;

	foreach (sort keys %$table_with_func_part) {
		my $belief_value = $$table_with_info_part{$_};
		my $value = $$table_with_func_part{$_} * &multiply_array_except_someone($belief_value,$loc);
		push @comp_func_belief, $value;
	}


	my $update_value = &sum_array(\@comp_func_belief);
	$update_value = log $update_value;


	my $belief_value = $$msg{$edge_name};
	$$belief_value[$sta] = $update_value;
}


sub loc_var_pos {
	my ($var,$edge_name,$var_num,$ko_pairs) = @_;
	my $ko_num = scalar keys %$ko_pairs;
	my $var_sub;
	if ($edge_name =~ /^v(\d+)_v/) {
		$var_sub = $1;
	}else{
		die "Error!";
	}
	my $loc = 0;
	if ($var eq 'x') {
		$loc += $ko_num - 1 + $var_sub;
	}elsif ($var eq 's') {
		$loc += $ko_num - 1 + $var_num + $var_sub;
	}elsif ($var eq 'd') {
		$loc += $ko_num - 1 + 2 * $var_num + $var_sub;
	}
	return $loc;
}


sub extract_table_part {
	my ($loc,$sta,$table_with_info,$table_with_func) = @_;
	my $table_with_info_part;
	my $table_with_func_part;

	my @conf_part;
	foreach (sort keys %$table_with_info) {
		my $conf = $_;
		if (/^\d{$loc}(\d)/) {
			if ($1 == $sta) {
				push @conf_part, $conf;
			}
		}
	}
	foreach (@conf_part) {
		$$table_with_info_part{$_} = $$table_with_info{$_};
		$$table_with_func_part{$_} = $$table_with_func{$_};
	}
	return ($table_with_info_part,$table_with_func_part);
}


sub multiply_array_except_someone {
	my ($array,$loc) = @_;
	my $result = 1;
	my $i = 0;
	foreach (@$array) {
		my $value = $_;
		$result *= $value unless $loc == $i;
		$i++;
	}
	return $result;
}


sub max_array {
	my ($array) = @_;
	my $max;
	my @pos;
	my $i = 0;
	foreach (@$array) {
		if ($i == 0) {
			$max = $_;
			push @pos,$i;
		}elsif ($_ > $max) {
			$max = $_;
			undef @pos;
			$pos[0] = $i;
		}elsif ($_ == $max) {
			push @pos, $i;
		}
		$i++;
	}
	return ($max,\@pos);
}


sub sum_array {
	my ($array) = @_;
	my $sum = 0;
	foreach (@$array) {
		$sum += $_;
	}
	return $sum;
}


sub times_array {
	my ($array) = @_;
	my $mul = 1;
	foreach (@$array) {
		$mul *= $_;
	}
	return $mul;
}
# ##################





##################
sub show_message {
	my ($msg_total, $var_list) = @_;

	foreach (sort keys %$msg_total) {
		say OUTPUT1;
		my $msg = $$msg_total{$_};
		foreach (sort keys %$msg) {
			my $ary = $$msg{$_};
			say OUTPUT1 join("\t",($_, @$ary));
		}
	}
}
##################



sub duration_time {
	my ($duration) = @_;
	my $print_time = "\nProgram duration time: ";
	if ($duration<60){
		$print_time .= "$duration s\n\n";
	}elsif ($duration<3600){
		$print_time .= int($duration/60).'m '.$duration%(60).'s';
	}else{
		$print_time .= int($duration/3600).'h '.int($duration/60)%(60).'m '.$duration%(60).'s';
	}
	return $print_time;
}

