#!/usr/bin/perl
use strict;
use warnings;

# This script is to parse the official GO.terms_alt_ids file
# in order to take both the primary and the secondary ID for
# the annotations

open(IN,"GO.terms_alt_ids");
while(<IN>) {
  next unless $_ =~ /^GO/;
  my @f = split(/\t/);
  my @secondary = split(/ /,$f[1]);
  print "$f[0]\t$f[2]\n";
  foreach my $s(@secondary){
    print "$s\t$f[2]\n";
  }
}
