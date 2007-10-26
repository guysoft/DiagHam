#!/usr/bin/perl -w
# Batch for calculations of the combined V0-V1 interaction for fermions with spin
use strict 'vars';
use Math::Trig;

# global settings
my $DiagonalizationProgram="/home/gunnar/DiagHam/build1/FQHE/src/Programs/FQHEOnSphere/FQHESphereFermionsWithSpin";
my $memory = 0;

if (!(defined($ARGV[2])))
  {
    print ("Run FQHESphereFermionsWithSpin with pseudopotentials V0=cos(phi), V1=sin(phi)\n");
    die "usage: ClockModelRun.pl nbr_fermions nbr_flux nbr_points [RescaleFirst] [print-command]\n";
  }

my $NbrFermions = $ARGV[0];
my $NbrFlux = $ARGV[1];
my $NbrPoints = $ARGV[2];

my $R0=($NbrFlux+1)*($NbrFlux+1)/(4*pi*(2*$NbrFlux+1));
my $R1=($NbrFlux+1)*($NbrFlux+1)/($NbrFlux*(2*$NbrFlux+1)/2.0);

my $RawRescale;
if (!defined($ARGV[3]))
  {
   $RawRescale=1.0;
  }
else
  {
    print ("Rescaling pseudopotentials with finite size scaling\n");
    if ( $NbrFlux == 2*$NbrFermions-4 )
      {
	# shift of the Haldane-Rezayi state -> filling factor 1/2
	$RawRescale=sqrt(0.5*$NbrFlux/$NbrFermions);
	print ("1/2 state\n");
      }
    else
      {
	if ( $NbrFlux == 3*$NbrFermions/2-3 )
	  {
	    # shift of the 2/3 state -> filling factor 2/3
	    $RawRescale=sqrt(2.0*$NbrFlux/$NbrFermions/3.0);
	    print ("2/3 state\n");
	  }
	else
	  {
	    $RawRescale=1.0;
	    print("Attention: No state is known at this shift! Scaling factor set to one!");
	  }
      }
    }

my $R20=$RawRescale*$RawRescale;
my $R21=$RawRescale*$RawRescale*$RawRescale*$RawRescale;
  
my $V0eff;
my $V1eff;
my $Phi;

for ( my $point=0; $point<$NbrPoints; $point++)
  {
    $Phi=$point*2*pi/($NbrPoints-1);
    $V0eff=$R0*$R20*cos($Phi);
    $V1eff=$R1*$R21*sin($Phi);
    # write pseudopotentials to file
    my $PseudopotentialFile = "pseudopotentials_2s_".$NbrFlux."_phi_".$Phi/pi."pi.dat";
    unless (open (OUTFILE, ">".$PseudopotentialFile))
      {
	die ("can't create pseudopotential file\n");
      }
    print OUTFILE ("# Pseudopotentials for clock model without finite size correction for 2s=".$NbrFlux." at angle Phi=".$Phi/pi." pi\n");
    print OUTFILE ("Pseudopotentials =");
    print OUTFILE ($V0eff." ".$V1eff);
    for ( my $i=2; $i<=$NbrFlux; $i++)
      {
	print OUTFILE (" 0");
      }
    close (OUTFILE);
    # run calculation
    my $Command = $DiagonalizationProgram." -p ".$NbrFermions." -l ".$NbrFlux." -s 0 --show-itertime --memory ".$memory." --interaction-file ".$PseudopotentialFile." --interaction-name clock_phi_".$Phi/pi;
    #$Command = $Command." --szsymmetrized-basis";
    if (!defined($ARGV[4]))
      {
	system ($Command);
      }
    else
      {
	print ("To run for Phi=".$Phi/pi."pi, type: \n".$Command."\n");
      }
  }
