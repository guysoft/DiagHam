#!/usr/bin/perl -w

use strict 'vars';


if (!(defined($ARGV[0])))
  {
    die ("usage: MagneticTranslation2Translation filename\n");
  }

unless (open (INFILE, $ARGV[0]))
  {
    die ("can't open ".$ARGV[0]."\n");
  }
my $TmpLine;
my %Spectrum;
foreach $TmpLine (<INFILE>)
  {
    chomp ($TmpLine);
    my @TmpArray = split (/ /, $TmpLine);
    my $TmpKey = $TmpArray[1];
    if (!(defined($Spectrum{$TmpKey})))
      {
	my @Array = ($TmpArray[2]);
	$Spectrum{$TmpKey} = \@Array;
      }
    else
      {
	my $Array = $Spectrum{$TmpKey};
	push (@$Array, $TmpArray[2]);
      }
  }
close (INFILE);

foreach $TmpLine (sort {$a <=> $b} (keys(%Spectrum)))
  {
    my $Array = $Spectrum{$TmpLine};
    my $Value;
    foreach $Value (sort {$a <=> $b} (@$Array))
      {
	print $TmpLine." ".$Value."\n";
      }
  }
