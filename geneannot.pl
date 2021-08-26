#!/usr/local/bin/perl -w

# BINF2010 Assignment Term 3 2020
# by Samuel Thorley (z5257239)

use JSON;
use LWP::Simple;
use HTTP::Request;

system 'source ~binftools/setup-env.sh';

# File handling:
if ($#ARGV < 0) {
    die "./geneannot.pl: one file must be supplied\n";
} elsif ($#ARGV > 0) {
    die "./geneannot.pl: only one file can be supplied\n";
}

my $file = $ARGV[0];

if (!($file =~ /^\S*\.fasta$/)) {
    die "./geneannot.pl: file supplied must be of .fasta format\n";
}

# Prints error message and exits if getorf fails
system "getorf -sequence $file -minsize 150 -outseq temp1.fasta 2>/dev/null" 
and die "./geneannot.pl: $file invalid file format\n";

# Extract lengths from first line of every ORF
# Sizes stored in %sizes hash

my %sizes = ();
open(FILE1, '<', "temp1.fasta") or die "./geneannot.pl: could not open temporary file\n";
while (my $line = <FILE1>) {
    if ($line =~ /^>.*/) {
        my $index = $line;
        $index =~ /\.(\S*).*/;
        $index = $1;
        my $start = $line;
        $start =~ /\[(\d*).*/;
        $start = $1;
        my $end = $line;
        $end =~ /(\d*)\].*/;
        $end = $1;
        $length = abs($end - $start);
        $sizes{"$index"} = $length;
    }
}
close(FILE1);

# @sizes sorted with 3 highest lengths stored in array @highest
# if less than 4 ORFs, all ORFs stored in @highest
my @lengths = sort { $sizes{$a} <=> $sizes{$b} } keys %sizes;
my @highest;
my $i = 1;
if ($#lengths >= 3) {
    while ($i < 4) {
        push(@highest, $lengths[-$i]);
        $i++;
    }
} else {
    @highest = @lengths;
}

# Loops through first temp file (temp1.fasta)
# 3 highest ORFs stored in new temp file (temp2.fasta) based on values in @highest
# All lines copied if less than 3
my $output = 0;
open(FILE1, '<', "temp1.fasta") or die "./geneannot.pl: could not open temporary file\n";
open(FILE2, '>', "temp2.fasta") or die "./geneannot.pl: could not open new temporary file\n";
while (my $line = <FILE1>) {
    if ($#highest < 2) {
        print FILE2 "$line";
    } else {
        if ($line =~ /^>.*/) {
            $output = 0;
            my $index = $line;
            $index =~ /\.(\S*).*/;
            $index = $1;
            if ("$index" eq "$highest[0]" or "$index" eq "$highest[1]" or "$index" eq "$highest[2]") {
                $output = 1;
                print FILE2 "$line";
            }
        } elsif ($output == 1) {
            print FILE2 "$line";
        }
    }
}
close(FILE1);
close(FILE2);

# Loop through highest ORFs (temp2.fasta) 
# Responses from RCSB PDB Search API stored in new file (temp3.fasta) to be used
# Corresponding start index, end index, and strand direction also stored
# Start arrays of ORFs also stored in new array @start_indexes
open(FILE2, '<', "temp2.fasta") or die "./geneannot.pl: could not open temporary file\n";
open(FILE3, '>', "temp3.fasta") or die "./geneannot.pl: could not open new temporary file\n";
my $sequence = "";
my @start_indexes;
while (my $line = <FILE2>) {
    chomp($line);
    if ($line =~ /^>.*/) {
        if ($sequence ne "") {
            my $result = create_json($sequence);
            print FILE3 "$result\n";
        }
        $sequence = "";
        my $start = $line;
        $start =~ /\[(\d*).*/;
        $start = $1;
        push @start_indexes, $start;
        print FILE3 "\nstart: $start\n";
        my $end = $line;
        $end =~ /(\d*)\].*/;
        $end = $1;
        print FILE3 "end: $end\n";
        my $strand = "FORWARD";
        if ($line =~ /\(REVERSE SENSE\)/) {
            $strand = "REVERSE";
        }
        print FILE3 "strand: $strand\n";
    } elsif ($line =~ /\w*/) {
        $sequence .= $line;
    }
}
my $result = create_json($sequence);
print FILE3 "$result\n";
close(FILE2);
close(FILE3);

# Array of start indexes sorted to be in order numerically
my @sorted_start_indexes = sort { $a <=> $b } @start_indexes;

# Resulting temp file (temp3.fasta) looped with results listed in .csv format
# Order based on @sorted_start_indexes
print "Start,End,Strand,PDB_ID,E-value\n";
foreach my $start (@sorted_start_indexes) {
    my $found = 0;
    my $ID = "";
    my $evalue = "";
    open(FILE3, '<', "temp3.fasta") or die "./geneannot.pl: could not open new temporary file\n";
    print "$start,";
    while (my $line = <FILE3>) {
        if ($line =~ /^start: $start/) {
            $found = 1;
        } elsif ($line =~ /^start:.*/) {
            $found = 0;
        }
        if ($found == 1) {
            if ($line =~ /^end: .*/) {
                my $end = $line;
                $end =~ /(\d+)/;
                $end = $1;
                print "$end,";
            } elsif ($line =~ /^strand: .*/) {
                my $strand = $line;
                $strand =~ / (\w+)$/;
                $strand = $1;
                print "$strand,";
            } elsif ($line =~ /.*"identifier".*/ and $ID eq "") {
                $ID = $line;
                $ID =~ /: "(\S+)",/;
                $ID = $1;
                print "$ID,";
            } elsif ($line =~ /.*"evalue".*/ and $evalue eq "") {
                $evalue = $line;
                $evalue =~ / : (\S+),/;
                $evalue = $1;
                print "$evalue";
            }
        }
    }
    if ($ID eq "" and $evalue eq "") {
        print "-,-";
    }
    print "\n";
    close(FILE3);
}

# Delete temporary files
unlink("temp1.fasta") or die "./geneannot.pl: temp1.fasta could not be deleted\n";
unlink("temp2.fasta") or die "./geneannot.pl: temp2.fasta could not be deleted\n";
unlink("temp3.fasta") or die "./geneannot.pl: temp3.fasta could not be deleted\n";


# Helper function given sequence creates json to be parsed to RCSB PDB Search API
# Returns results from search
sub create_json {
    my $sequence = $_[0];
    my %parameters = (
        "evalue_cutoff" => 1,
        "identity_cutoff" => 0.4,
        "target" => "pdb_protein_sequence",
        "value" => $sequence
    );
    my %query = (
        "type" => "terminal", 
        "service" => "sequence",
        "parameters" => \%parameters
    );
    my %request_options = (
        "scoring_strategy" => "sequence"
    );
    my %json = (
        "query" => \%query,
        "request_options" => \%request_options,
        "return_type" => "polymer_entity"
    );
    my $data = encode_json \%json;
    my $url = "https://search.rcsb.org/rcsbsearch/v1/query";
    my $header = ['Content-Type' => 'application/json; charset=UTF-8'];
    my $ua = LWP::UserAgent->new;
    my $resp = HTTP::Request->new('POST', $url, $header, $data);
    die "./geneannot.pl: no results returned from RCSB PDB Search API\n" unless defined $resp;
    $r = $ua->request($resp);
    $results = $r->content;
    return $results;
}