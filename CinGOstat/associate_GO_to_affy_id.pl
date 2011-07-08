#!/usr/bin/perl;
use strict;
use warnings;
use DBI;
use Data::Dumper;

# WHEN THE DATABASE IS READY THIS SCRIPT ASSOCIATE TO EACH PROBE THE CORRESPONDING GO ID
# BASED ON THE ANISEED ANNOTATIONS OF THE MODELS
# REDUNDANCY IS ELIMINATED EVERY PROBE HAS ONLY ONE OCCURRENCE FOR EACH TERM ASSOCIATED

my $dbh = connect_to_db('ciona_affy_annotation','mysql_dev','dEvEl0pEr','');

my $probes = get_all_probe_id($dbh);

foreach my $probe(keys %$probes) {
  my $models = get_all_model_id($dbh,$probe);
  my $gos = get_all_go_id($dbh,$models);
  foreach my $go(keys %$gos) {
    print "$probe\t$go\n";
  }
}

sub get_all_go_id {
  my $dbh = shift;
  my $models = shift;
  my $href = {};
  my $sth = $dbh->prepare("SELECT DISTINCT go_id FROM model_go WHERE model_id IN($models)");
  $sth->execute;
  while(my $row = $sth->fetchrow_hashref) {
    $href->{$row->{go_id}} ++;
  }
  return $href;
}

sub get_all_model_id {
  my $dbh = shift;
  my $probe = shift;
  my $models;
  my $sth = $dbh->prepare('SELECT DISTINCT model_id FROM association WHERE probe_id = '.$dbh->quote($probe));
  $sth->execute;
  while(my $row = $sth->fetchrow_hashref) {
    $models .= "'".$row->{model_id}."'".',';
  }
  $models =~ s/\,$//;
  return $models;
}

sub get_all_probe_id {
  my $dbh = shift;
  my $href = {};
  my $sth = $dbh->prepare('SELECT DISTINCT probe_id FROM association');
  $sth->execute;
  while(my $row = $sth->fetchrow_hashref) {
    $href->{$row->{probe_id}} ++;
  }
  return $href;
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
