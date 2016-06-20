use strict;
use warnings;
use feature 'say';
use DBI;
use Getopt::Long;

#
# Import2db_sqlite.pl
# Reqirements: Perl5
# Author: Xiao-Tai Huang
# Email: xthuang226@gmail.com
# Version: 1.0
# Date: 2016/06/19
#

my $file        = "Import2db_sqlite.pl";
my $program     = uc($&) if $file =~ /^\w+/;
my $version     = "1.0";
my $date        = "2016/06/19";
my $author      = "Xiao-Tai Huang";

### program identification
my $title = "$program $version, $date $author\n";


### option variables
$::opt_help         = 0;
$::opt_inpath       = '../data/single_cell_data/';
$::opt_outpath      = '../data/single_cell_data/cdfiles_nhr25.db3';
$::opt_expinfopath  = '../data/single_cell_data/experiment_info.csv';


my $usage = <<"END_OF_USAGE";
${title}Syntax:   \L$program\E [options]
Function: Import cell lineage data into SQLite.
Options:                                                             (defaults:)
  -h, --help            print usage
  -i, --inpath          set input cell lineage data .csv file path   ($::opt_inpath)
  -o, --outpath         set output SQLite .db3 file path             ($::opt_outpath)
  -x, --expinfopath     set input experiment_info.csv file path      ($::opt_expinfopath)
Examples:
  \L$program\E
  \L$program\E -i=../data/single_cell_data/ -o=../data/single_cell_data/cdfiles_nhr25.db3 -x=../data/single_cell_data/experiment_info.csv
END_OF_USAGE


### process options
my @OrgArgv = @ARGV;
GetOptions(
  "help|h!",
  "inpath|i=s",
  "outpath|o=s",
  "expinfopath|x=s"
) or die $usage;
!$::opt_help or die $usage;


### Set parameters
my ($inpath, $outpath, $expinfopath) = ($::opt_inpath, $::opt_outpath, $::opt_expinfopath);


### Load experiment information.
my %exper_info;
open(INPUT, $expinfopath)||die $!;
# open(INPUT, $expinfopath)||die $!;
foreach (<INPUT>) {
	chomp;
	unless (/^exp/) {
		my @row = split(/,/, $_);
		$exper_info{"CD$row[0]"} = \@row;
	}
}
close INPUT;


my $output = $outpath;
$output =~ s/\/([^\/]+)$//;

my $output_tmp = $output;
my $output_whole;

while ($output_tmp =~ s/^([\.\/]*[^\.\/]+\/?)//) {
	$output_whole .= $&;
	mkdir $output_whole;
}


### Connect to database.
my $conn = DBI->connect("dbi:SQLite:dbname=$outpath","","",{ RaiseError => 1 }) or die $DBI::errstr;
my %tables;
foreach (split (/[\s\n]+/, `sqlite3 $outpath .tables`)) {
	$tables{$_} = $_;
}


foreach (qw/w m/) {			### 'w' for wild-type '.csv' files import, while 'm' for mutant.

	my $etype = $_;

	### Build a table named 'table_info' for mapping the original data table to the processed data table, as wel as other experiment information.
	my ($sql,$query);
	unless (defined $tables{'table_info'}) {
		$sql = "
			CREATE TABLE table_info(
				curr_name TEXT NOT NULL,
				orig_name TEXT PRIMARY KEY  NOT NULL,
				lineaging_strain TEXT,
				marker_gene TEXT,
				ko_gene TEXT,
				imaging_time_point INT,
				edited_cell_number INT,
				edited_time_point INT
			);
		";
		$query = $conn->prepare($sql);
		$query->execute;
		undef $query;
		$tables{'table_info'} = 'table_info';
	}


	### Import .csv files.
	opendir (DIR, "$inpath$etype") || die "Cannot open dir : $!";
	my @files = readdir(DIR);
	close DIR;


	my $i = 0;
	foreach (@files){
		if (/^((.+)\.csv)$/){
			unless (defined $tables{$etype.'_'.$2}) {
				$i++;
				say "\n$1";
				&insert1cdf ("$etype/$1","$etype\_$2", $inpath, $outpath);
				say "$etype\_t$i";
				&info_table ($etype,$2,$i);
				say "Adding cell supplementary information...";
				&add_supp_info("$etype\_t$i");
			}
		}
	}
}


map{say}sort keys %tables;



$conn->disconnect;
undef $conn;






### Load data in .csv file into tables of database.
### Table name starting with 'w_' represents wild-type, while 'm_' for mutant, e.g. 'w_CD120911NHR25p4' and 'm_CD120927NHR25fkb3ip1'.
sub insert1cdf {
	my ($file, $table, $inpath) = @_;

	my $sql = "
		CREATE TABLE \'[$table]\' (
			cellTime	TEXT PRIMARY KEY NOT NULL,
			cell	TEXT,
			time	INT,
			none	INT,
			global	INT,
			local	INT,
			blot	INT,
			cross	INT,
			z	REAL,
			x	INT,
			y	INT,
			size	INT,
			gweight	INT
		);
	";
	my $query = $conn->prepare($sql);
	$query->execute;
	undef $query;

	open(PIPE, "|sqlite3 $outpath") or die "Open pipe error: $!";
	say PIPE ".mode csv\n.import $inpath$file [$table]\n.exit";
	close PIPE;
	$tables{$table} = $table;

	$sql = "
		DROP TABLE \'$table\';
	";
	$query = $conn->prepare($sql);
	$query->execute;
	undef $query;

	$sql = "
		ALTER TABLE \'[$table]\' RENAME TO \'$table\';
	";
	$query = $conn->prepare($sql);
	$query->execute;
	undef $query;

}


### Create tables with useful information from original tables.
### Table name starting with 'w_' represents wild-type, while 'm_' for mutant, e.g. 'w_t1' and 'm_t1'.
sub info_table {
	my ($etype,$exper,$i) = @_;

	unless (defined $tables{"$etype\_t$i"}) {
		my $sql = "
			CREATE TABLE $etype\_t$i
			(
			  cell TEXT PRIMARY KEY,
			  time INT,
			  len INT,
			  exp VARYING TEXT,
			  avgexp REAL,
			  cellpath TEXT
			);
		";
		my $query = $conn->prepare($sql);
		$query->execute;
		undef $query;
		$tables{"$etype\_t$i"} = "$etype\_t$i";

		$sql = "
			INSERT INTO $etype\_t$i (cell, time, len, exp, avgexp)
			SELECT cell, min(time), count(*), group_concat(blot), avg(blot)
			FROM \'$etype\_$exper\'
			WHERE cell GLOB \'[ACDEMPZ]*\'
			GROUP BY cell;
		";
		$query = $conn->prepare($sql);
		$query->execute;
		undef $query;

		my $row = $exper_info{$exper};
		$sql = "
			INSERT INTO table_info VALUES (\'$etype\_t$i\',\'$etype\_$exper\',\'$$row[1]\',\'$$row[2]\',\'$$row[3]\',$$row[4],$$row[5],$$row[6]);
		";
		$query = $conn->prepare($sql);
		$query->execute;
		undef $query;
	}
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


### Add some supplementary information about cells into tables, such as cell path, left or right child.
sub add_supp_info {
	my ($table) = @_;
	my $sql = "
		SELECT cell from $table
	";
	my $query = $conn->prepare($sql);
	$query->execute;

	while (my @row = $query->fetchrow_array) {
		my $sql = "
			UPDATE $table
			SET cellpath = ?
			WHERE cell IS ?;
		";
		my $query = $conn->prepare($sql);
		$query->execute(&get_cellpath($row[0]), $row[0]);
		undef $query;		
	}
	undef $query;
}

