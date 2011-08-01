#!/usr/bin/perl
use strict;
use warnings;

my $USAGE = "\n\tUSAGE: perl $0 [annotation_file] [probe_id]\n\n";
my $file = $ARGV[0];
my $id = $ARGV[1];
die $USAGE unless $file;
die $USAGE unless -e $file;
die $USAGE unless $id;

my $tmp = 'probe.ann.tmp';
my $command = "grep \'$id\' $file \> $tmp";
# Add check for field (ID or MODEL?)

system('clear');
system('clear');
system('clear');
system($command);

die "\n\nProbe id does not exists!\n\n" if -z $tmp;

open(IN,$tmp);
my $row = <IN>;
chomp $row;
my @field = split(/\t/,$row);

# probe_id model_id syno orto blast go design model_seq

my $model_seq = $field[7];
$model_seq =~ s/\,/\n\n/g;
$model_seq =~ s/\:/\n/g;

print "
----------------------------------------------------------------------------
Probe ID
---------
$field[0]

Matching Models
----------------
$field[1]

Syno
-----
$field[2]

Orto
-----
$field[3]

Blast
------
$field[4]

Gene Ontology
--------------
$field[5]

Design Sequence
----------------
$field[6]

Model Sequences
----------------

$model_seq
--------------------------------------------------------

";

system("rm $tmp");
