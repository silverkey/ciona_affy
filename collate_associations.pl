#/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

my $aass = 'aniseed_models.xls'; # transcript ID are in 1st column
my $bass = 'ASSOCIATIONS_PROBE_ANISEED.xls'; # transcript ID are in 2th column, probe ID are in the 1st
my $out = 'PROBE_ANNOTATION.xls';
open(OUT,">$out");

my ($ahref,$ahead) = read_table($aass,1);

my @hstring = qw(probe_id transcript_id);
push(@hstring,@$ahead);
print OUT join("\t",@hstring);
print OUT "\n";

my $associated = collate_files($bass,2,1,$ahref);

sub collate_files {
  my $file = shift;
  my $key_col = shift;
  my $pid_col = shift;
  my $ahref = shift;

  open(IN,$file);
  my $h = <IN>;
  # my $header = get_head($h);
  while(my $row = <IN>) {
    chomp($row);
    my @field = split(/\t/,$row);
    my $pid = $field[$pid_col-1];
    my $tid = $field[$key_col-1];
    my @string = ($pid,$tid);
    if (exists $ahref->{$tid}) {
      my $el = $ahref->{$tid};
      foreach my $head(@$ahead) {
        if (exists $el->{$head}) {
	  push(@string,$el->{$head});
        }
	else {
	    push(@string,'NA');
        }
      }
      print OUT join("\t",@string);
      print OUT "\n";
    }
  }
}    

# make an hashref of hashref from a table using as a key the id in the column specified
# and using as keys of the internal hashref the headers of the table avoiding the id column
sub read_table {
  my $file = shift;
  my $col = shift;
  my $href = {};
  open(IN,$file);
  my $h = <IN>;
  my ($head,$titles) = get_head($h,$col);
  while(my $row = <IN>) {
    chomp($row);
    my @field = split(/\t/,$row);
    my $id = $field[$col-1];
    for (my $c=0;$c<=$#field;$c++) {
      $href->{$id}->{$head->[$c]} = $field[$c] unless $c == $col-1;
    }
  }
  return ($href,$titles);
}

sub get_head {
  my $row = shift;
  # strip out the primary key field only for the hader!
  my $strip = shift;
  chomp($row);
  my @field = split(/\t/,$row);
  my @new_field;
  for (my $c=0;$c<=$#field;$c++) {
    push(@new_field,$field[$c]) unless $c == $strip-1;
  }
  return (\@field,\@new_field);
}
