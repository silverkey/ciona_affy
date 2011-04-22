#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use DBI;

my $id = $ARGV[0];
my $dbh = DBI->connect("dbi:SQLite:dbname=ciona_affy_annotation.sqlite3");

my $anno = "select distinct probe_id, model_id, syno, orto, blast from association where probe_id=\'$id\'";
$dbh->do($anno);





__END__
select distinct probe_id, model_id, syno, orto, blast from association where probe_id='ENSCINT00000021990_at';
select model_id, seq from seq where model_id in (select model_id from association where probe_id='ENSCINT00000021990_at');
select probe_id, seq from design where probe_id = 'ENSCINT00000021990_at';
select distinct definition from model_go m, go g where m.go_id=g.go and model_id in (select model_id from association where probe_id='ENSCINT00000021990_at');
