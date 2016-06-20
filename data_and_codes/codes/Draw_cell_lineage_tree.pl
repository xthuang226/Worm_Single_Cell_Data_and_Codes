use strict;
use warnings;
use feature 'say';
use DBI;
use Tree::Binary2;
use SVG;
use Getopt::Long;


#
# Draw_cell_lineage_tree.pl
# Reqirements: Perl5
# Author: Xiao-Tai Huang
# Email: xthuang226@gmail.com
# Version: 1.0
# Date: 2016/06/19
#


my $file        = "Draw_cell_lineage_tree.pl";
my $program     = uc($&) if $file =~ /^\w+/;
my $version     = "1.0";
my $date        = "2016/06/19";
my $author      = "Xiao-Tai Huang";


### program identification
my $title = "$program $version, $date $author\n";


### option variables
my @bool = ("false", "true");
$::opt_help         = 0;
$::opt_root         = 'P0';
$::opt_endtime      = 130;
$::opt_cellstage    = 0;
$::opt_cutoff       = 1466.8;
$::opt_lefontsize   = 20;
$::opt_lafontsize   = 15;
$::opt_blafontsize  = 10;
$::opt_axfontsize   = 30;
$::opt_brfontsize   = 30;
$::opt_scale        = 0;
$::opt_lineinter    = 0;
$::opt_linewidth    = 0.1;
$::opt_width        = 6000;
$::opt_height       = 800;
$::opt_model        = 0;
$::opt_binary       = 0;
$::opt_label        = 0;
$::opt_bottomlabel  = 0;
$::opt_axis         = 0;
$::opt_brand        = 0;
$::opt_autozoom     = 0;
$::opt_indb         = '../data/single_cell_data/cdfiles_nhr25.db3';
$::opt_output       = '../results/cell_lineage_tree/';


my $usage = <<"END_OF_USAGE";
${title}Syntax:   \L$program\E <table> [options]
Function: Draw cell lineage tree into .svg format.
Options:                                                              (defaults:)
  -h, --help         print usage
  -r, --root         root of the tree                              ($::opt_root)
  --endtime          endtime of the tree                           ($::opt_endtime)
  --cellstage        cell stage of the tree                        ($::opt_cellstage)
  --cutoff           set cutoff for binarizing gene expression     ($::opt_cutoff)
  --lefontsize       set top leader cell font size                 ($::opt_lefontsize)
  --lafontsize       set cell label font size                      ($::opt_lafontsize)
  --blafontsize      set bottom cell label font size               ($::opt_blafontsize)
  --axfontsize       set axis font size                            ($::opt_axfontsize)
  --brfontsize       set brand font size                           ($::opt_brfontsize)
  --scale            control vertical line scale                   ($::opt_scale)
  --lineinter        set gap between two vertical lines            ($::opt_lineinter)
  --linewidth        set vertical line width                       ($::opt_linewidth)
  -w, --width        set graph width                               ($::opt_width)
  -h, --height       set graph height                              ($::opt_height)
  -i, --indb         set input SQLite database .db3 file           ($::opt_indb)
  -o, --output       set output tree graph folder path             ($::opt_output)
  -m, --model        draw model tree                               ($bool[$::opt_model])
  -b, --binary       draw tree with binary gene expression value   ($bool[$::opt_binary])
  --label            draw label for each cell                      ($bool[$::opt_label])
  --bottomlabel      draw label for bottom cells of the tree       ($bool[$::opt_bottomlabel])
  --axis             draw axis for the tree                        ($bool[$::opt_axis])
  --brand            draw brand for the tree                       ($bool[$::opt_brand])
  --autozoom         autozoom graph with window size               ($bool[$::opt_autozoom])
Examples:
  \L$program\E w_t3
  \L$program\E w_t3 --root P0 --cellstage=350 --autozoom
  \L$program\E w_t3 --root P0 --cellstage=350 -b --axis --brand --label --bottomlabel
  \L$program\E w_t3 --root P0 --endtime=140 -m --axis --brand --label --bottomlabel
END_OF_USAGE


### process options
my @OrgArgv = @ARGV;
GetOptions(
  "help|h!",
  "root|r=s",
  "endtime=f",
  "cellstage=i",
  "cutoff=f",
  "lefontsize=f",
  "lafontsize=f",
  "blafontsize=f",
  "axfontsize=f",
  "brfontsize=f",
  "scale=f",
  "lineinter=f",
  "linewidth=f",
  "width|w=f",
  "height|h=f",
  "model|m!",
  "binary|b!",
  "label!",
  "bottomlabel!",
  "axis!",
  "brand!",
  "autozoom!",
  "indb=s",
  "output=s"
) or die $usage;
!$::opt_help or die $usage;


### Set parameters
my $origin = 'P0';
my ($table, $root, $endtime, $cellstage, $indb, $output) = ($ARGV[0], $::opt_root, $::opt_endtime, $::opt_cellstage, $::opt_indb, $::opt_output);
die $usage unless defined $table;
my ($binary, $model, $cutoff) = ($::opt_binary, $::opt_model, $::opt_cutoff);
my ($label, $bottom_label, $axis, $brand) = ($::opt_label, $::opt_bottomlabel, $::opt_axis, $::opt_brand);
my ($width, $height, $linewidth, $lineinter, $scale, $autozoom) = ($::opt_width, $::opt_height, $::opt_linewidth, $::opt_lineinter, $::opt_scale, $::opt_autozoom);
my ($lefontsize, $lafontsize, $blafontsize, $axfontsize, $brfontsize) = ($::opt_lefontsize, $::opt_lafontsize, $::opt_blafontsize, $::opt_axfontsize, $::opt_brfontsize);


$binary = 1 if $model == 1;		##### If want to draw a model tree, automatically set as a binary tree.
if ($autozoom == 1) {
	$label = 0;
	$bottom_label = 0;
	$axis = 0;
	$brand = 0;
}


my ($marginleft, $marginright, $margintop, $marginbottom) = (1, 1, 1, 2);


my $axis_width = 0;
$axis_width =  2 if $axis == 1;
my $brand_width = 0;
$brand_width = 2 if $brand == 1;
my $rootheight = $lefontsize/10;
my $botlabheight = 0;
$botlabheight = $blafontsize*1 if $bottom_label == 1;


my $spaceleft = $marginleft + $axis_width;
my $spaceright = $marginright + $brand_width;
my $spacetop = $margintop + $rootheight;
my $spacebottom = $marginbottom + $botlabheight;


### Connect to DBMS
my $conn = DBI->connect("dbi:SQLite:dbname=$indb","","",{ RaiseError => 1 }) or die $DBI::errstr;


($endtime, my ($leafnumb, $treehash, $cell_exp_ary, $cell_exp_avg, $cell_length, $add_cells)) = &get_data_to_be_drew ($table, $endtime, $cellstage, $model, $origin);
map{say "Invalid data in cell: $_"}@$add_cells;

unless (defined $$cell_length{$root}) {
	die "\n[Warning: $root is invalid! Please input valid root name.]\n\n";
}


### Get sub tree info
my ($root_start_tp, $sub_tp, $sub_leafnumb) = &sub_tree_info($root, $treehash, $cell_length, $model, $add_cells, $endtime);


### Compute $scale and $lineinter if no input for them
if ($scale == 0 || $lineinter == 0) {
	$lineinter = (100 - $spaceleft - $spaceright) / ($sub_leafnumb - 1);
	$scale = (100 - $spacetop - $spacebottom) / $sub_tp;
}

my $treewidth = $lineinter * ($sub_leafnumb - 1);
my $treeheight = $sub_tp * $scale;

my $cellx = &fix_xaxis_position($treehash, $root, $spaceleft, $linewidth, $lineinter);


my $output_tmp = $output;
my $output_whole;

while ($output_tmp =~ s/^([\.\/]*[^\.\/]+\/?)//) {
	$output_whole .= $&;
	mkdir $output_whole;
}

my $outpath = "$output/$table\_$root\_tp$endtime\_c$leafnumb.svg";


&draw_tree ($root, $treehash, $cell_length, $cellx, $cell_exp_avg, $scale, $sub_tp, $model, $binary, $cutoff, $label, $bottom_label, $linewidth, $spacetop, $spacebottom, $spaceleft, $spaceright, $lefontsize, $lafontsize, $blafontsize, $axfontsize, $brfontsize, $autozoom, $width, $height, $root_start_tp, $treewidth, $treeheight, $marginleft, $endtime, $axis, $brand);


### Disconnect DBMS
$conn->disconnect();
$conn = undef;










sub get_data_to_be_drew {
	my ($table, $endtime, $cellstage, $model, $origin) = @_;

	my ($add_cells, $leafnumb, $treehash, $cell_exp_ary, $cell_exp_avg, $cell_length);

	if ($cellstage != 0) {
		$endtime = 130 unless defined $endtime;				##### If endtime is undefined, give a default endtime.
		($endtime, $add_cells, $leafnumb, $treehash, $cell_exp_ary, $cell_exp_avg, $cell_length) = &modify_endtime($table, $endtime, $cellstage, $model, $origin);	##### Get endtime with the given cellstage.
	}elsif (defined $endtime) {
		($cell_exp_ary, $cell_exp_avg, $cell_length) = &get_tree_data_by_endtime ($table, $endtime, $model);

		($add_cells, $cell_exp_ary, $cell_exp_avg, $cell_length) = &supplement_tree_data ($cell_exp_ary, $cell_exp_avg, $cell_length, $model);

		$treehash = &insert_treehash(keys %$cell_exp_avg);
		$leafnumb = &compute_cellstage($treehash, $origin);
		say "=> timepoint:$endtime\tcellstage:$leafnumb";
	}else{
		die 'Please give a specific $endtime or $cellstage !';
	}

	return ($endtime, $leafnumb, $treehash, $cell_exp_ary, $cell_exp_avg, $cell_length, $add_cells);
}


sub modify_endtime {
	my ($table, $endtime, $cellstage, $model, $origin) = @_;

	my ($cell_exp_ary, $cell_exp_avg, $cell_length) = &get_tree_data_by_endtime ($table, $endtime, $model);
	(my $add_cells, $cell_exp_ary, $cell_exp_avg, $cell_length) = &supplement_tree_data ($cell_exp_ary, $cell_exp_avg, $cell_length, $model);
	my $treehash = &insert_treehash(keys %$cell_exp_avg);
	my $leafnumb = &compute_cellstage($treehash, $origin);
	say "timepoint:$endtime\tcellstage:$leafnumb";

	my %record;

	while ($leafnumb != $cellstage){

		$record{$endtime} = abs($leafnumb-$cellstage);

		if ($leafnumb > $cellstage){
			$endtime--;
		}elsif($leafnumb < $cellstage){
			$endtime++;
		}

		if (defined $record{$endtime}) {
			my @record = sort {$record{$a} <=> $record{$b}} keys %record;
			$endtime = $record[0];
			($cell_exp_ary, $cell_exp_avg, $cell_length) = &get_tree_data_by_endtime ($table, $endtime, $model);
			($add_cells, $cell_exp_ary, $cell_exp_avg, $cell_length) = &supplement_tree_data ($cell_exp_ary, $cell_exp_avg, $cell_length, $model);
			$treehash = &insert_treehash(keys %$cell_exp_avg);
			$leafnumb = &compute_cellstage($treehash, $origin);
			say "timepoint:$endtime\tcellstage:$leafnumb";
			last;
		}else{
			($cell_exp_ary, $cell_exp_avg, $cell_length) = &get_tree_data_by_endtime ($table, $endtime, $model);
			($add_cells, $cell_exp_ary, $cell_exp_avg, $cell_length) = &supplement_tree_data ($cell_exp_ary, $cell_exp_avg, $cell_length, $model);
			$treehash = &insert_treehash(keys %$cell_exp_avg);
			$leafnumb = &compute_cellstage($treehash, $origin);
			say "timepoint:$endtime\tcellstage:$leafnumb";
		}

	}

	say "=> timepoint:$endtime\tcellstage:$leafnumb";
	return ($endtime, $add_cells, $leafnumb, $treehash, $cell_exp_ary, $cell_exp_avg, $cell_length);
}


sub get_tree_data_by_endtime {
	my ($table, $endtime, $model) = @_;
	my ($cell_exp_ary, $cell_exp_avg, $cell_length);

	my $sql = "SELECT * FROM $table WHERE time <= $endtime";
	my $query = $conn->prepare($sql);
	$query->execute();

	while ( my @row = $query->fetchrow_array() ) {
		my @exp_ary = split (/,/, $row[3]);
		$$cell_exp_ary{$row[0]} = \@exp_ary;
		$$cell_exp_avg{$row[0]} = $row[4];
		if ($model == 1) {
			$$cell_length{$row[0]} = ($row[5] =~ tr/././) + 1;
		}else{
			$$cell_length{$row[0]} = $row[2];
		}
	}

	undef $query;
	return ($cell_exp_ary, $cell_exp_avg, $cell_length);
}


sub supplement_tree_data {
	my ($cell_exp_ary, $cell_exp_avg, $cell_length, $model) = @_;
	my $add_cells;

	my $cell_level_sup_ref;
	$$cell_level_sup_ref{'P0'} = 1;
	$$cell_level_sup_ref{'AB'} = 2;
	$$cell_level_sup_ref{'P1'} = 2;
	$$cell_level_sup_ref{'ABa'} = 3;
	$$cell_level_sup_ref{'ABp'} = 3;
	$$cell_level_sup_ref{'EMS'} = 3;
	$$cell_level_sup_ref{'P2'} = 3;
	$$cell_level_sup_ref{'ABal'} = 4;
	$$cell_level_sup_ref{'ABar'} = 4;
	$$cell_level_sup_ref{'ABpl'} = 4;
	$$cell_level_sup_ref{'ABpr'} = 4;
	$$cell_level_sup_ref{'MS'} = 4;
	$$cell_level_sup_ref{'E'} = 4;
	$$cell_level_sup_ref{'C'} = 4;
	$$cell_level_sup_ref{'P3'} = 4;
	$$cell_level_sup_ref{'D'} = 5;
	$$cell_level_sup_ref{'P4'} = 5;
	$$cell_level_sup_ref{'Z2'} = 6;
	$$cell_level_sup_ref{'Z3'} = 6;

	# foreach (qw/P0 AB P1 ABa ABp EMS P2 ABal ABar ABpl ABpr MS E C P3 D P4 Z2 Z3/){
	foreach (qw/P0 AB P1 ABa ABp EMS P2/){
		my $add_cell = $_;
		unless (defined $$cell_exp_ary{$add_cell}){
			my @exp_ary = qw/0 0 0 0/;
			$$cell_exp_ary{$add_cell} = \@exp_ary;
			$$cell_exp_avg{$add_cell} = 0;

			if ($model == 1) {
				$$cell_length{$add_cell} = $$cell_level_sup_ref{$add_cell};
			}else{
				$$cell_length{$add_cell} = scalar @exp_ary;
			}

			push @$add_cells, $add_cell;
		}
	}
	return ($add_cells, $cell_exp_ary, $cell_exp_avg, $cell_length);
}


sub insert_treehash {
	my (@cells) = @_;
	my %treehash;

	### Put elements of @cells into a tree.
	foreach (@cells){
		$treehash{$_} = Tree::Binary2->new($_);
	}

	### Connect node to its parent.
	foreach (@cells){
		my $cell = $_;
		
		if ($cell =~ /^(A?B?M?S?E?C?D?.*)([aplrdv])$/){
			my $pres = $1;
			my $ends = $2;
			
			if ($ends =~ /[lad]/){
				$treehash{$pres} = Tree::Binary2->new( $pres ) unless defined $treehash{$pres};
				$treehash{$pres}->left( $treehash{$cell} );
			}
			if ($ends =~ /[rpv]/){
				$treehash{$pres} = Tree::Binary2->new( $pres ) unless defined $treehash{$pres};
				$treehash{$pres}->right( $treehash{$cell} );
			}
		}
		
		if ($cell =~ /^P1$/){
			$treehash{'P0'} = Tree::Binary2->new( 'P0' ) unless defined $treehash{'P0'};
			$treehash{'P0'}->right( $treehash{$cell} );
		}
		if ($cell =~ /^AB$/){
			$treehash{'P0'} = Tree::Binary2->new( 'P0' ) unless defined $treehash{'P0'};
			$treehash{'P0'}->left( $treehash{$cell} );
		}
		if ($cell =~ /^EMS$/){
			$treehash{'P1'} = Tree::Binary2->new( 'P1' ) unless defined $treehash{'P1'};
			$treehash{'P1'}->left( $treehash{$cell} );
		}
		if ($cell =~ /^P2$/){
			$treehash{'P1'} = Tree::Binary2->new( 'P1' ) unless defined $treehash{'P1'};
			$treehash{'P1'}->right( $treehash{$cell} );
		}
		if ($cell =~ /^P3$/){
			$treehash{'P2'} = Tree::Binary2->new( 'P2' ) unless defined $treehash{'P2'};
			$treehash{'P2'}->right( $treehash{$cell} );
		}
		if ($cell =~ /^C$/){
			$treehash{'P2'} = Tree::Binary2->new( 'P2' ) unless defined $treehash{'P2'};
			$treehash{'P2'}->left( $treehash{$cell} );
		}
		if ($cell =~ /^P4$/){
			$treehash{'P3'} = Tree::Binary2->new( 'P3' ) unless defined $treehash{'P3'};
			$treehash{'P3'}->right( $treehash{$cell} );
		}
		if ($cell =~ /^D$/){
			$treehash{'P3'} = Tree::Binary2->new( 'P3' ) unless defined $treehash{'P3'};
			$treehash{'P3'}->left( $treehash{$cell} );
		}
		if ($cell =~ /^Z3$/){
			$treehash{'P4'} = Tree::Binary2->new( 'P4' ) unless defined $treehash{'P4'};
			$treehash{'P4'}->left( $treehash{$cell} );
		}
		if ($cell =~ /^Z2$/){
			$treehash{'P4'} = Tree::Binary2->new( 'P4' ) unless defined $treehash{'P4'};
			$treehash{'P4'}->right( $treehash{$cell} );
		}
		if ($cell =~ /^E$/){
			$treehash{'EMS'} = Tree::Binary2->new( 'EMS' ) unless defined $treehash{'EMS'};
			$treehash{'EMS'}->right( $treehash{$cell} );
		}
		if ($cell =~ /^MS$/){
			$treehash{'EMS'} = Tree::Binary2->new( 'EMS' ) unless defined $treehash{'EMS'};
			$treehash{'EMS'}->left( $treehash{$cell} );
		}
	}
	
	return \%treehash;
}


sub compute_cellstage {
	my ($treehash, $origin) = @_;
	
	my $leafnumb;

	my $check_root = $$treehash{$origin};

	my $trav = $check_root->traverse($check_root->POST_ORDER);
	while ( my $node = $trav->() ) {
		unless ($node -> children){
			$leafnumb++;
		}
	}

	return $leafnumb;
}


sub compute_end_tp_for_cell {
	my ($cell, $treehash, $cell_length) = @_;
	my $tp = $$cell_length{$cell};

	my $parent = $$treehash{$cell} -> parent;
	while ($parent){
		$tp += $$cell_length{$parent->value};
		$parent = $parent -> parent;
	}
	return $tp;
}


sub compute_end_tp_and_leafnumb_for_tree {
	my ($root, $treehash, $cell_length) = @_;
	my $end_tp = 0;
	my $leafnumb = 0;

	my $check_root = $$treehash{$root};
	my $trav = $check_root->traverse($check_root->LEVEL_ORDER);
	while ( my $node = $trav->() ) {
		unless ($node -> children){
			my $tp = &compute_end_tp_for_cell($node->value, $treehash, $cell_length);
			$end_tp = $tp if $tp > $end_tp;
			$leafnumb++;
		}
	}
	return ($end_tp, $leafnumb);
}


sub sub_tree_info {
	my ($root, $treehash, $cell_length, $model, $add_cells, $endtime) = @_;

	my $root_start_tp = &compute_end_tp_for_cell($root, $treehash, $cell_length) - $$cell_length{$root};

	my ($sub_end_tp, $sub_leafnumb) = &compute_end_tp_and_leafnumb_for_tree($root, $treehash, $cell_length);

	my $sub_tp = $sub_end_tp - $root_start_tp;

	if ($model == 0) {
		my $sup_tp = 0;
		foreach (@$add_cells) {
			my $tp = &compute_end_tp_for_cell($_, $treehash, $cell_length);
			$sup_tp = $tp if $tp > $sup_tp;
		}
		if ($sup_tp < $endtime) {
			$sub_tp = $endtime + $sup_tp;
		}
	}

	return ($root_start_tp, $sub_tp, $sub_leafnumb);
}


sub fix_xaxis_position {
	my ($treehash, $root, $spaceleft, $linewidth, $lineinter) = @_;
	
	my %cellx;
	my $xleaf = $spaceleft;	

	my $check_root = $$treehash{$root};
	my $trav = $check_root->traverse($check_root->POST_ORDER);
	while ( my $node = $trav->() ) {
		
		my $cell = $node -> value;

		if (my @children = $node -> children){

			die "The data of cell $cell is not valid!\n" if scalar @children != 2;
			
			my $lx = $cellx{$children[0]->value};
			my $rx = $cellx{$children[1]->value};
			
			$cellx{$cell} = $lx + 0.5*($rx-$lx);
			
		}else{
			$cellx{$cell} = $xleaf;
			$xleaf += $lineinter;
		}
		
	}
	return \%cellx;
}


sub draw_tree {
	my ($root, $treehash, $cell_length, $cellx, $cell_exp_avg, $scale, $sub_tp, $model, $binary, $cutoff, $label, $bottom_label, $linewidth, $spacetop, $spacebottom, $spaceleft, $spaceright, $lefontsize, $lafontsize, $blafontsize, $axfontsize, $brfontsize, $autozoom, $width, $height, $root_start_tp, $treewidth, $treeheight, $marginleft, $endtime, $axis, $brand) = @_;

	### leader cells
	my %leader_cells;
	foreach (qw/P0 P1 P2 P3 P4 AB ABa ABp EMS E MS C D Z2 Z3/) {
		$leader_cells{$_} = 1;
	}

	undef %leader_cells if $autozoom == 1;
	
	my ($x1,$y1,$x2,$y2) = qw/0 0 0 0/;

	my $celly;
	$$celly{$root} = $spacetop;

	my $celly_start_tp;
	$$celly_start_tp{$root} = 0;
	
	my $gray = "rgb(160,160,160)";
	my @linecolors;
	for (0..9){
		push @linecolors, "rgb(".int($_*(255/9)).",0,0)";
	}

	my $svg;
	if ($autozoom == 1) {
		$svg = SVG->new();
	}else{
		$svg = SVG->new(width=>$width, height=>$height);
	}

	my $check_root = $$treehash{$root};
	my $trav = $check_root->traverse($check_root->LEVEL_ORDER);
	while ( my $node = $trav->() ) {
		my $cell = $node->value;

		$x1 = $$cellx{$cell};
		$y1 = $$celly{$cell};
		$x2 = $x1;
		$y2 = $y1 + ($$cell_length{$cell} * $scale);

		if (my @children = $node -> children){

			foreach (@children) {
				$$celly{$_->value} = $y2;
				$$celly_start_tp{$_->value} = $$celly_start_tp{$cell} + $$cell_length{$cell};
			}

			### Draw label
			my $x_label = $x1;
			my $y_label = $y1; 

			if (defined $leader_cells{$cell}) {
				if ($cell eq $root) {
					$x_label += -0.15;
					$y_label += -0.1;
				}elsif($cell eq 'E'){
					$x_label += -0.1;
					$y_label += -0.3;
				}elsif($cell eq 'C'){
					$x_label += -0.1;
					$y_label += -0.3;
				}elsif($cell eq 'D'){
					$x_label += -0.1;
					$y_label += -0.3;
				}else{
					$x_label += -0.2;
					$y_label += -0.3;
				}
				$svg->text(x => "$x_label%", y => "$y_label%", "fill"=>"balck", "font-family"=>"Times New Roman","font-size"=>$lefontsize, "-cdata" => $cell);
			}elsif ($label == 1) {
				$x_label += 0.1;
				$y_label += 0.5;
				my $cx = $width * $x_label / 100;
				my $cy = $height * $y_label / 100;
				$svg->text(x => "$x_label%", y => "$y_label%", fill => $gray, transform =>"rotate(90,$cx,$cy)", "font-family"=>"Times New Roman","font-size"=>$lafontsize, "-cdata" => $cell);
			}

		}else{
			$y2 = 100 - $spacebottom if $model == 1 || $y2 > (100 - $spacebottom);

			if ($bottom_label == 1) {
				my $x_bottom_label = $x2 - 0.05;
				my $y_bottom_label = 100 - $spacebottom + 1;
				my $cx = $width*$x_bottom_label/100;
				my $cy = $height*$y_bottom_label/100;
				$svg->text(x => "$x_bottom_label%", y => "$y_bottom_label%", transform =>"rotate(90, $cx, $cy)", "font-family"=>"Times New Roman","font-size"=>$blafontsize, "-cdata" => $cell);
			}
		}

		my $linecolor;
		if ($binary == 1) {
			if ($$cell_exp_avg{$cell} >= $cutoff){
				$linecolor = 'red';
			}else{
				$linecolor = 'black';
			}

			my $draw_method = 0;
			if ($draw_method == 0) {
				### Draw vertical line via line method
				$svg->line(x1 => "$x1%", y1 => "$y1%", x2 => "$x2%", y2 => "$y2%", stroke=>$linecolor,"stroke-width"=>"$linewidth%");
			}else{
				for (1..$$cell_length{$cell}) {
					my $cur_tp = $_;
					if (($$celly_start_tp{$cell}+$cur_tp) <= $sub_tp) {
						$svg->line(x1 => "$x1%", y1 => "$y1%", x2 => "$x2%", y2 => ($y1+$scale)."%", stroke=>$linecolor,"stroke-width"=>"$linewidth%");
						$y1 += $scale;
					}
				}
			}
		}else{

			### Draw vertical line via point method
			for (1..$$cell_length{$cell}) {
				my $cur_tp = $_;
				if (($$celly_start_tp{$cell}+$cur_tp) <= $sub_tp) {
					my $exps = $$cell_exp_ary{$cell};
					$linecolor = &set_linecolor($$exps[$cur_tp-1],\@linecolors);
					$svg->line(x1 => "$x1%", y1 => "$y1%", x2 => "$x2%", y2 => ($y1+$scale)."%", stroke=>$linecolor,"stroke-width"=>"$linewidth%");
					$y1 += $scale;
				}
			}

		}

		if (my @children = $node -> children) {
			### Draw branch
			$svg->line(x1 => $$cellx{$children[0]->value}."%", y1 => "$y2%", x2 => $$cellx{$children[1]->value}."%", y2 => "$y2%", stroke=>$linecolor,"stroke-width"=>"$linewidth%");
		}

	}


	### draw table name and cell stage
	if ($autozoom == 0) {
		$svg->text(x => "$spaceleft%", y => ($spacetop+1)."%", "font-family"=>"Times New Roman","font-size"=>$axfontsize, "-cdata" => "$table", "fill"=>"rgb(0,200,0)", "font-style"=>"italic");
		$svg->text(x => "$spaceleft%", y => ($spacetop+5)."%", "font-family"=>"Times New Roman","font-size"=>$axfontsize, "-cdata" => "$leafnumb-cell stage", "fill"=>"rgb(0,200,0)", "font-style"=>"italic");
	}


	### draw axis
	if ($axis == 1) {
		my $y1_axis = $spacetop;
		my $y4_axis = $spacetop + $treeheight;
		my $y2_axis = int ((1/3)*($y4_axis-$y1_axis)) + $y1_axis;
		my $y3_axis = int ((2/3)*($y4_axis-$y1_axis)) + $y1_axis;
		my $x1_axis = $marginleft;

		my $axislabel1 = $root_start_tp * 1.5;
		my $axislabel4 = $endtime * 1.5;
		my $axislabel2 = int ((1/3)*($axislabel4-$axislabel1));
		my $axislabel3 = int ((2/3)*($axislabel4-$axislabel1));

		$svg->line(x1 => "$x1_axis%", y1 => "$y1_axis%", x2 => "$x1_axis%", y2 => "$y4_axis%", stroke=>"black","stroke-width"=>"$linewidth%");
		my $x2_axis = $x1_axis + 0.5;
		my $x3_axis = $x2_axis + 0.2;

		$svg->line(x1 => "$x1_axis%", y1 => "$y1_axis%", x2 => "$x2_axis%", y2 => "$y1_axis%", stroke=>"black","stroke-width"=>"$linewidth%");
		$svg->text(x => "$x3_axis%", y => ($y1_axis+1)."%", "font-family"=>"Times New Roman","font-size"=>$axfontsize, "-cdata" => $axislabel1);

		$svg->line(x1 => "$x1_axis%", y1 => "$y2_axis%", x2 => "$x2_axis%", y2 => "$y2_axis%", stroke=>"black","stroke-width"=>"$linewidth%");
		$svg->text(x => "$x3_axis%", y => ($y2_axis+1)."%", "font-family"=>"Times New Roman","font-size"=>$axfontsize, "-cdata" => $axislabel2);

		$svg->line(x1 => "$x1_axis%", y1 => "$y3_axis%", x2 => "$x2_axis%", y2 => "$y3_axis%", stroke=>"black","stroke-width"=>"$linewidth%");
		$svg->text(x => "$x3_axis%", y => ($y3_axis+1)."%", "font-family"=>"Times New Roman","font-size"=>$axfontsize, "-cdata" => $axislabel3);

		$svg->line(x1 => "$x1_axis%", y1 => "$y4_axis%", x2 => "$x2_axis%", y2 => "$y4_axis%", stroke=>"black","stroke-width"=>"$linewidth%");
		$svg->text(x => "$x3_axis%", y => ($y4_axis+1)."%", "font-family"=>"Times New Roman","font-size"=>$axfontsize, "-cdata" => $axislabel4);

		my $cx = $width*($x3_axis-1)/100;
		my $cy = $height*(($y4_axis-$y1_axis)/2+7)/100;
		$svg->text(x => ($x3_axis-1)."%", y => (($y4_axis-$y1_axis)/2+7)."%", transform =>"rotate(-90, $cx, $cy)", "font-family"=>"Times New Roman","font-size"=>$axfontsize, "-cdata" => "Minute");
		# $svg->text(x => ($x3_axis-0.37)."%", y => ($y4_axis+4)."%", "font-family"=>"Times New Roman","font-size"=>$axfontsize, "-cdata" => "(min)");
	}


	## draw brand
	if ($brand == 1) {
		my $x_brand = $spaceright + $treewidth + 1;
		my $y1_brand = $spacetop;
		my $y2_brand = $spacetop + $treeheight;

		if ($binary == 1){
			$svg->line(x1 => "$x_brand%", y1 => "$y1_brand%", x2 => "$x_brand%", y2 => (($y1_brand+$y2_brand)/2)."%", stroke=>"black","stroke-width"=>3*$linewidth."%");
			my $cx = $width*($x_brand+0.3)/100;
			my $cy = $height*$y1_brand/100;
			$svg->text(x => ($x_brand+0.3)."%", y => "$y1_brand%", transform =>"rotate(90, $cx, $cy)", "font-family"=>"Times New Roman","font-size"=>$brfontsize, "-cdata" => "not express");

			$svg->line(x1 => "$x_brand%", y1 => (($y1_brand+$y2_brand)/2)."%", x2 => "$x_brand%", y2 => "$y2_brand%", stroke=>"red","stroke-width"=>3*$linewidth."%");
			$cy = $height*(($y1_brand + $y2_brand)/2)/100;
			$svg->text(x => ($x_brand+0.3)."%", y => (($y1_brand + $y2_brand)/2)."%", transform =>"rotate(90, $cx, $cy)", "font-family"=>"Times New Roman","font-size"=>$brfontsize, "-cdata" => "express");
		}else{
			my $i = 1;
			foreach (@linecolors){
				$svg -> line(x1 => "$x_brand%", y1 => ($y1_brand+(($i-1)/10)*($y2_brand-$y1_brand))."%", x2 => "$x_brand%", y2 => ($y1_brand+($i/10)*(int($y2_brand-$y1_brand)))."%", stroke => $linecolors[$i-1], "stroke-width"=>3*$linewidth."%");
				$i++;
			}
			my $cx = $width*($x_brand+0.3)/100;
			my $cy = $height*(($y2_brand-$y1_brand)/2-9)/100;
			$svg->text(x => ($x_brand+0.3)."%", y => (($y2_brand-$y1_brand)/2-9)."%", transform =>"rotate(90, $cx, $cy)", "font-family"=>"Times New Roman","font-size"=>$brfontsize, "-cdata" => "expression level");
		}
	}


	## Save svg
	my $out = $svg->xmlify;
	open SVGFILE, ">$outpath";
	print SVGFILE $out;

}


sub set_linecolor {
	my ($exp,$linecolors) = @_;
	
	if ($exp >= 5000){
		return $$linecolors[9];
	}elsif($exp <= 0){
		return $$linecolors[0];
	}else{
		return $$linecolors[int($exp*10/5000)];
	}
}


