#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use DBI;
use DBD::mysql;
use Bio::SeqIO;

my $design = 'CINT06a520380F_sif.fa';

# DATABASE SPECIFIC SETTINGS
my $DB = '';
my $USR = 'mysql_dev';
my $PWD = 'dEvEl0pEr';
my $HOST;
my $DBH = connect_to_db($DB,$USR,$PWD,$HOST);

$DBH->do('CREATE DATABASE ciona_affy_annotation');
$DBH->do('USE ciona_affy_annotation');

my $create_association = q!
CREATE TABLE association (
  probe_id VARCHAR(256),
  model_id VARCHAR(256),
  syno VARCHAR(256),
  orto VARCHAR(256),
  blast VARCHAR(256),
  KEY probe_id(probe_id),
  KEY model_id(model_id)
)
!;
$DBH->do($create_association);
$DBH->do("LOAD DATA LOCAL INFILE \'PROBE_ANNOTATION.xls\' INTO TABLE association IGNORE 1 LINES");

my $create_seq = q!
CREATE TABLE seq (
  model_id VARCHAR(256) PRIMARY KEY,
  seq TEXT
)
!;
$DBH->do($create_seq);
$DBH->do("LOAD DATA LOCAL INFILE \'aniseed_models_seq.xls\' INTO TABLE seq IGNORE 1 LINES");

my $create_model_go = q!
CREATE TABLE model_go (
  model_id VARCHAR(256),
  go_id VARCHAR(30),
  KEY model_id(model_id),
  KEY go_id(go_ID)
)
!;
$DBH->do($create_model_go);
$DBH->do("LOAD DATA LOCAL INFILE \'aniseed_models_go.xls\' INTO TABLE model_go IGNORE 1 LINES");

my $create_go = q!
CREATE TABLE go (
  go_id VARCHAR(30) PRIMARY KEY,
  secondary VARCHAR(256),
  definition VARCHAR(256),
  class VARCHAR(1),
  KEY definition(definition)
)
!;
$DBH->do($create_go);
system("grep -v \'\\\!\' GO.terms_alt_ids \> GO.terms_alt_ids.mod");
$DBH->do("LOAD DATA LOCAL INFILE \'GO.terms_alt_ids.mod\' INTO TABLE go");

open(OUT,">$design\.xls");
my $seqio = Bio::SeqIO->new(-file => $design,
                            -format => 'fasta');
while(my $seq = $seqio->next_seq) {
  print OUT join("\t",$seq->id,$seq->seq);
  print OUT "\n";
}
close(OUT);
my $create_design = q!
CREATE TABLE design (
  probe_id VARCHAR(256) PRIMARY KEY,
  seq TEXT
)
!;
$DBH->do($create_design);
$DBH->do("LOAD DATA LOCAL INFILE \'$design\.xls\' INTO TABLE design IGNORE 1 LINES");

sub connect_to_db {
  my $db = shift;
  my $usr = shift;
  my $pwd = shift;
  my $host = shift;
  my $dsn = 'dbi:mysql:'.$db;
  $dsn .= ':'.$host if $host; # IN THE CURRENT DBI POD VERSION THERE IS THE '@' IN THE PLACE OF ':'
  my $dbh = DBI->connect($dsn,$usr,$pwd,{PrintWarn=>1,PrintError=>1,RaiseError=>1}) or die $DBI::errstr;
  return $dbh;
}
