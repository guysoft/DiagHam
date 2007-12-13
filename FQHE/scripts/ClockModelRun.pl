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
    die "usage: ClockModelRun.pl nbr_fermions nbr_flux nbr_points filling [GS] [print-command]\n";
  }

my $NbrFermions = $ARGV[0];
my $NbrFlux = $ARGV[1];
my $NbrPoints = $ARGV[2];
my $Filling = $ARGV[3];

# previously utilized version:
#my $R0=($NbrFlux+1)*($NbrFlux+1)/(4*pi*(2*$NbrFlux+1));
#my $R1=($NbrFlux+1)*($NbrFlux+1)/($NbrFlux*(2*$NbrFlux+1)/2.0);

# version provided by Thierry
my $R0=($NbrFlux+1)*($NbrFlux+1)/($NbrFlux*(2*$NbrFlux+1)/2.0);
my $R1=($NbrFlux+1)*($NbrFlux+1)/($NbrFlux*(2*$NbrFlux-1)/2.0);

my $RawRescale;
if (1==2) # always use scaling factor now...
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
	if ($Filling =~ "1/2" )
	  {
	    print ("1/2 state\n");
	  }
	else
	  {
	    die ("wrong flux for nu=1/2.\n");
	  }
      }
    else
      {
	if ($Filling =~ "2/3" )
	  {
	    $RawRescale=sqrt(2.0*$NbrFlux/$NbrFermions/3.0);
	    if ( $NbrFlux == 3*$NbrFermions/2-3 )
	      {
		# shift of the 2/3 state -> filling factor 2/3
		print ("2/3 state at 3/2 N - 3.\n");
	      }
	    else
	      {
		if ( $NbrFlux == 3*$NbrFermions/2-1 )
		  {
		    # shift of the 2/3 state -> filling factor 2/3
		    print ("2/3 state at 3/2 N - 1.\n");
		  }
		else
		  {
		    die ("wrong flux for nu=2/3.\n");
		  }
	      }
	  }
	else
	  {
	    if ($Filling =~ "1" )
	      {
		$RawRescale=sqrt($NbrFlux/$NbrFermions);
		if ( $NbrFlux == $NbrFermions-1 )
		  {
		    # shift of the p-wave paired state
		    print ("nu=1 ferromagnetic state at N - 1.\n");
		  }
		else
		  {
		    if ( $NbrFlux == $NbrFermions-2 )
		      {
			# shift of the s-wave paired state
			print ("nu=1 s-wave paired state at N - 2.\n");
		      }
		    else
		      {
			die ("wrong flux for nu=1.\n");
		      }
		  }
	      }
	    else
	      {
		$RawRescale=1.0;
		print("Attention: Unknown state! Scaling factor set to one!");
	      }
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
    $Phi=$point*2*pi/($NbrPoints);
    $V0eff=$R0*$R20*cos($Phi);
    $V1eff=$R1*$R21*sin($Phi);
    # write pseudopotentials to file
    $Filling =~ s/\//\_/g;
    print("Filling:".$Filling);
    my $PseudopotentialFile = "pseudopotentials_2s_".$NbrFlux."_phi_".$Phi/pi."pi_nu_".$Filling.".dat";
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
    if (defined($ARGV[4]))
      {
	$Command = $Command." --nbr-lz 1 --eigenstate -n2 --force-reorthogonalize";
      }
    if (!defined($ARGV[5]))
      {
	system ($Command);
      }
    else
      {
	print ("To run for Phi=".$Phi/pi."pi, type: \n".$Command."\n");
      }
  }
