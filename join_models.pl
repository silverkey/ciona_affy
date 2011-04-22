#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;

# aniseed link: http://www.aniseed.cnrs.fr/molecule-gene.php?name=ci0100130007

# Takes the gene models annotations ans sequences downloaded from ANISEED
# using the script 'download_models.sh' and produce two files joining the
# annotations with the sequences. The files are one table and one fasta
# where some annotations are in the description.

# The folder in which there are the sequences and the annotation downloaded
# Basically the folder in which 'download_models.sh' has been run

# Modified to add GO and to build 3 tables instead of only 1:
# 1) GO
# 2) Annotations
# 3) Sequences

my @fasta = glob('*.fa');
my @ann = glob('*.ann');

#print "\nGoing to parse files:\n\n";
#print join("\n",@fasta)."\n".join("\n",@ann)."\n\n";

my $href = {};
open(ANNO,">aniseed_models.xls");
print ANNO "id\tsyno\torto\tblast\n";
open(SEQ,">aniseed_models_seq.xls");
print SEQ "id\tseq\n";
open(GO,">aniseed_models_go.xls");
print GO "id\tgo\n";

my $outfa = Bio::SeqIO->new(-file => ">aniseed_models.fa",
                            -format => 'fasta');

foreach my $file(@fasta) {
  next if $file eq 'CINT06a520380F_sif.fa';
  my $fasta = Bio::SeqIO->new(-file => $file,
                              -format => 'fasta');
  while(my $seq = $fasta->next_seq) {
    $href->{$seq->id}->{seq} = $seq->seq;
  }
}

foreach my $file(@ann) {
  open(IN,$file);
  my $id = 0;
  my $seen = {};
  while(my $row = <IN>) {
    chomp($row);
    if($row =~ /^ID\t(.+)$/) {
      $id = "$1";
      print "\nERROR: $id seen twice!!!\n" if exists $seen->{$id};
      $seen->{$id}++;
      $href->{$id}->{syno} = '\N';
      $href->{$id}->{orto} = '\N';
      $href->{$id}->{blast} = '\N';
      $href->{$id}->{seq} = '\N';
    }
    if($row =~ /^SYNO\t(.+)\tORIGIN\t.+$/) {
      my $syno = "$1";
      if(exists $href->{$id}->{syno} && $href->{$id}->{syno} eq '\N') {
        $href->{$id}->{syno} = $syno unless($syno eq $id);
      }
    }
    if($row =~ /^ORTHOLOGY\t.+\tFAMILLY\t(.+)\tSCORE\t.+$/) {
      my $orto = "$1";
      if(exists $href->{$id}->{orto} && $href->{$id}->{orto} eq '\N') {
        $href->{$id}->{orto} = $orto;
      }
    }
    if($row =~ /^BLAST\t.+\tDESCRIPTION\t(.+)\tPROTEIN.+$/) {
      my $blast = "$1";
      if(exists $href->{$id}->{blast} && $href->{$id}->{blast} eq '\N') {
        $href->{$id}->{blast} = $blast;
      }
    }
    if($row =~ /^GO\tGO_ID\t(.+)\tORIGIN.+$/) {
      my $go = "$1";
      print GO "$id\t$go\n";
    }
    # Needed because some ID from the annotation files are not present in the fasta
    if($row =~ /^SEQ$/) {
      my $seq;
      my $end = 0;
      while(!$end) {
        my $newrow = <IN>;
        chomp($newrow);
        $end ++ if $newrow =~ /\/\//;
        $seq .= $newrow;
      }
      if(exists $href->{$id}->{seq} && $href->{$id}->{seq} eq '\N') {
        $href->{$id}->{seq} = $seq;
      }
    }
  }
}

foreach my $id(keys %$href) {
  my $info = $href->{$id};
  print ANNO "$id\t\\N\t\\N\t\\N\n" unless $info->{seq};
  my $string = $info->{seq};
  $string =~ s/\/\///;
  $string =~ s/\_//g;
  my $seqobj = Bio::Seq->new(-id => $id,
                             -seq => $string);
  $outfa->write_seq($seqobj);
  print SEQ "$id\t$string\n";
  my @info = ($id,$info->{syno},$info->{orto},$info->{blast});
  print ANNO join("\t",@info);
  print ANNO "\n";
}
