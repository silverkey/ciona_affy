#!/usr/bin/perl
use strict;
use warnings;

use Bio::SearchIO;

my $out = 'ASSOCIATIONS_PROBE_ANISEED.xls';
open(OUT,">$out");

print OUT "PROBE_ID\tTRANSCRIPT_ID\tPERCENTAGE\tCOVERAGE\n";

my $blast_out = Bio::SearchIO->new(-file => 'MEGABLAST_DEFAULT.out',
                                   -format => 'blast');

while(my $result = $blast_out->next_result) {

  my $id = $result->query_name;
  my $db = $result->database_name;
  my $querylength = $result->query_length;

  while(my $hit = $result->next_hit) {
    my $hname = $hit->name;

    while(my $hsp = $hit->next_hsp) {
      my $percentage = $hsp->percent_identity;
      my $hsplength = $hsp->length;
      my $qcoverage = $hsplength/$querylength;
      my @info = ($id,$hname,$percentage,$qcoverage);
      next unless $percentage>=90 && $qcoverage>=0.75;
      print OUT join("\t",@info);
      print OUT "\n";
    }
  }
}

