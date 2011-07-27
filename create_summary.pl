#!/usr/bin/perl;
use strict;
use warnings;
use DBI;
use Data::Dumper;

# WHEN THE DATABASE IS READY THIS SCRIPT CREATES 2 TABLES WITH THE INFORMATION COLLECTED
# BASED ON THE ANISEED ANNOTATIONS OF THE MODELS WRITING ONE SINGLE ROW FOR EACH PROBE
# THE DIFFERENCE BETWEEN THE 2 TABLES IS THAT THE SECOND DOES NOT CONTAINS THE SEQUENCES
# OF THE MODEL

my $dbh = connect_to_db('ciona_affy_annotation','mysql_dev','dEvEl0pEr','');

my $probes = get_all_probe($dbh);

open(OUT,">CIONA_AFFY_ANNOTATION_SUMMARY.xls");
print OUT join("\t",qw(probe_id model_id syno orto blast go design model_seq));
print OUT "\n";

open(SMALL,">CIONA_AFFY_ANNOTATION_SUMMARY_SMALL.xls");
print SMALL join("\t",qw(probe_id model_id syno orto blast go design));
print SMALL "\n";

foreach my $probe(keys %$probes) {
  my $models = get_all_from_association($dbh,$probe,'model_id');
  my $synos = get_all_from_association($dbh,$probe,'syno');
  my $ortos = get_all_from_association($dbh,$probe,'orto');
  my $blasts = get_all_from_association($dbh,$probe,'blast');
  my $gos = get_all_go($dbh,$models);
  my $seqs = get_all_seq($dbh,$models);
  $models =~ s/\'//g;
  $synos =~ s/\'//g;
  $ortos =~ s/\'//g;
  $blasts =~ s/\'//g;
  print OUT join("\t",$probe,$models,$synos,$ortos,$blasts,$gos,$probes->{$probe},$seqs);
  print OUT "\n";
  print SMALL join("\t",$probe,$models,$synos,$ortos,$blasts,$gos,$probes->{$probe});
  print SMALL "\n";
}

sub get_all_go {
  my $dbh = shift;
  my $models = shift;
  return 'NA' if $models eq 'NA';
  my $res;
  my $sth = $dbh->prepare("SELECT DISTINCT definition FROM model_go m, go_definition g WHERE m.go_id=g.go_id AND model_id IN($models)");
  $sth->execute;
  while(my $row = $sth->fetchrow_hashref) {
    next unless $row->{definition};
#   $res .= "'".$row->{definition}."'".',';
    $res .= $row->{definition}.',';
  }
  $res =~ s/\,$// if $res;
  return $res || 'NA';
}

sub get_all_seq {
  my $dbh = shift;
  my $models = shift;
  return 'NA' if $models eq 'NA';
  my $res;
  my $sth = $dbh->prepare("SELECT DISTINCT * FROM seq WHERE model_id IN($models)");
  $sth->execute;
  while(my $row = $sth->fetchrow_hashref) {
    next unless $row->{model_id};
#   $res .= "'".$row->{model_id}.':'.$row->{seq}."'".',';
    $res .= $row->{model_id}.':'.$row->{seq}.',';
  }
  $res =~ s/\,$// if $res;
  return $res || 'NA';
}

sub get_all_from_association {
  my $dbh = shift;
  my $probe = shift;
  my $field = shift;
  my $models;
  my $sth = $dbh->prepare('SELECT DISTINCT '.$field.' FROM association WHERE probe_id = '.$dbh->quote($probe));
  $sth->execute;
  while(my $row = $sth->fetchrow_hashref) {
    next unless $row->{$field};
    $models .= "'".$row->{$field}."'".',';
#   $models .= $row->{$field}.',';
  }
  $models =~ s/\,$// if $models;
  return $models || 'NA';
}

sub get_all_probe {
  my $dbh = shift;
  my $href = {};
  my $sth = $dbh->prepare('SELECT DISTINCT * FROM design');
  $sth->execute;
  while(my $row = $sth->fetchrow_hashref) {
    next unless $row->{probe_id};
    $href->{$row->{probe_id}} = $row->{seq};
  }
  return $href || 'NA';
}

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

__END__
select distinct probe_id, model_id, syno, orto, blast from association where probe_id='ENSCINT00000021990_at';
select model_id, seq from seq where model_id in (select model_id from association where probe_id='ENSCINT00000021990_at');
select probe_id, seq from design where probe_id = 'ENSCINT00000021990_at';
select distinct definition from model_go m, go g where m.go_id=g.go and model_id in (select model_id from association where probe_id='ENSCINT00000021990_at');
