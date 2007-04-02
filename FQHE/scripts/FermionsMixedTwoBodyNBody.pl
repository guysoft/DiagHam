#!/usr/bin/perl -w

use strict 'vars';

my $NbrPoints = 11;
my $NbrParticles= 10;
my $SValue = 17;

my $Index = 0;
my $Step = 1.0 / ($NbrPoints - 1);
my $Value = 0.0;
#my $ThreeBodyWeight = 0.306582478478;
#my $ThreeBodyWeight = 0.644340439407;
my $ThreeBodyWeight = 1.0;
#my $InteractionName = "mixed_coulomb_1_dv1_0.05_nbody_3";
my $InteractionName = "mixed_v1_nbody_3";
my $PseudopotentialFileName = "pseudopotential_v1_2s_15.dat";
#my $PseudopotentialFileName = "pseudopotential_coulomb_l_1_dv1_0.05_2s_15.dat";
#my $PseudopotentialFileName = "pseudopotential_coulomb_l_1_2s_15.dat";
my $GapFileName = "fermions_".$InteractionName."_n_".$NbrParticles."_2s_".$SValue.".gap.dat";

my @Pseudopotentials;

unless (open (INFILE, $PseudopotentialFileName))
  {
    die ("can't open ".$PseudopotentialFileName."\n");
  }
my $TmpLine;
while (defined($TmpLine = <INFILE>))
  {
    if ($TmpLine =~ /^\s*Pseudopotentials\s*\=/)
      {
	chomp($TmpLine);
	$TmpLine =~ s/^\s*Pseudopotentials\s*\=\s*//;
	$TmpLine =~ s/\s*$//;
	@Pseudopotentials = split (/\s+/, $TmpLine);
      }
  }
close (INFILE);

unless (open (OUTFILE, ">".$GapFileName))
  {
    die ("can't open ".$GapFileName."\n");
  }
close (OUTFILE);

while ($Index < $NbrPoints)
  {
    unless (open (OUTFILE, ">pseudopotential_mixed_".$Value."_2s_".$SValue.".dat"))
      {
	die ("can't create file pseudopotential_mixed_".$Value."_2s_".$SValue.".dat\n");
      }
    print OUTFILE "# pseudopotentials on the sphere for coulomb interaction 
# in the Landau level N=1 for 2S=".$SValue." flux quanta
#
# Pseudopotentials = V_0 V_1 ...

Pseudopotentials =";
    my $TmpPseudopotential;
    foreach $TmpPseudopotential (@Pseudopotentials)
      {
	print OUTFILE " ".((1.0 - $Value) * $TmpPseudopotential);
      }
    print OUTFILE "\n\nNbrNBody=3

Weights=0 0 0 ".($ThreeBodyWeight * $Value)."\n";
    close (OUTFILE);
    my $Command = "/home/regnault/development/Physics/DiagHam/build/FQHE/src/Programs/FQHEOnSphere/QHEFermionsNBodyHardCore -p ".$NbrParticles." -l ".$SValue." --nbr-lz 2 --eigenstate --memory 2000 --force-reorthogonalize --nbr-lz 2 -n 1 --nbody-file pseudopotential_mixed_".$Value."_2s_".$SValue.".dat";
    system ($Command);
    my $OutputName = "fermions_".$InteractionName."_".$Value."_n_".$NbrParticles."_2s_".$SValue."_lz";
    rename ("fermions_hardcore_nbody_3_n_".$NbrParticles."_2s_".$SValue."_lz_0.0.vec", $OutputName."_0.0.vec");    
    rename ("fermions_hardcore_nbody_3_n_".$NbrParticles."_2s_".$SValue."_lz_0.0.vec", $OutputName."_2.0.vec");    
    $OutputName .= ".dat";
    rename ("fermions_hardcore_nbody_3_n_".$NbrParticles."_2s_".$SValue."_lz.dat", $OutputName);
    unless (open (INFILE, $OutputName))
      {
	die ("can't open ".$OutputName."\n");
      }
    $TmpLine = <INFILE>;
    chomp ($TmpLine);
    my @TmpArray = split (/ /, $TmpLine);
    my $Gap = -$TmpArray[1];
    while ((defined($TmpLine = <INFILE>)) && ($TmpLine =~/^0 /))
      {
      }
    chomp ($TmpLine);
    @TmpArray = split (/ /, $TmpLine);
    $Gap += $TmpArray[1];
    close (INFILE);
    unless (open (OUTFILE, ">>".$GapFileName))
      {
	die ("can't open ".$GapFileName."\n");
      }
    print OUTFILE $Value." ".$Gap."\n";
    close (OUTFILE);
    $Value += $Step;
    $Index++;
  }
