#!/usr/bin/perl -w
#
# script for optimizing overlaps for vectors associated with given protocol files
#
use strict 'vars';
use File::stat;
use Math::Complex;


# hardwire which state to look at
my $Program_32="FQHELatticeBosonsGeneric";
my $Program_64="FQHELatticeBosonsGeneric_64";
my $MatrixProgram="MatrixElement";

my $DiceLattice="NbrSites = 6\n"
  ."Dimension = 2\n"
  ."LatticeVector0 = 1.73205080756888,0\n"
  ."LatticeVector1 = 0,3\n"
  ."SubLatticeVector0 = 0,0\n"
  ."SubLatticeVector1 = 0,1\n"
  ."SubLatticeVector2 = 0,2\n"
  ."SubLatticeVector3 = 0.866025403784439,0.5\n"
  ."SubLatticeVector4 = 0.866025403784439,1.5\n"
  ."SubLatticeVector5 = 0.866025403784439,2.5\n"
  ."NeighborsInCell = 0,1 | 1,2 | 1,3 | 1,4 | 2,5 | 4,5\n"
  ."NeighborCells = 0,1 | 1,0 | 0,-1 | -1,0 | 1,1 | -1,-1\n"
  ."NeighborsAcrossBoundary0_1 = 5,0 | 5,3\n"
  ."NeighborsAcrossBoundary0_-1 = 0,5 | 3,5\n"
  ."NeighborsAcrossBoundary1_1 = 5,0\n"
  ."NeighborsAcrossBoundary-1_-1 = 0,5\n"
  ."NeighborsAcrossBoundary1_0 = 3,1 | 4,1 | 5,2\n"
  ."NeighborsAcrossBoundary-1_0 = 1,3 | 1,4 | 2,5\n"
  ."UseGauge = yes\n"
  ."GaugeAyx = 1.15470053837925\n";

my $EffectiveTriangularLattice="Descriptor = eff_triang\n"
  ."NbrSites = 2\n"
  ."Dimension = 2\n"
  ."LatticeVector0 = 1.73205080756888,0\n"
  ."LatticeVector1 = 0,3\n"
  ."SubLatticeVector0 = 0,0\n"
  ."SubLatticeVector1 = 0.866025403784439,1.5\n"
  ."NeighborsInCell = 0,1\n"
  ."NeighborCells = 0,1 | 1,0 | 0,-1 | -1,0 | 1,1 | -1,-1\n"
  ."NeighborsAcrossBoundary0_1 = 1,0\n"
  ."NeighborsAcrossBoundary0_-1 = 0,1\n"
  ."NeighborsAcrossBoundary1_1 = 1,0\n"
  ."NeighborsAcrossBoundary-1_-1 = 0,1\n"
  ."NeighborsAcrossBoundary1_0 =  0,0 | 1,0 | 1,1\n"
  ."NeighborsAcrossBoundary-1_0 = 0,0 | 0,1 | 1,1\n"
  ."UseGauge = yes\n"
  ."GaugeAyx = 1.15470053837925\n";

my $Descriptor = "dice_doubled";
my $Trapping = 0.001;
my $RatioU = 1.0;

my $Directory="";
my @UnitCells;
$UnitCells[0]=3;
$UnitCells[1]=2;
my $DisplayHelp=0;

while( (defined($ARGV[0])&&$ARGV[0] =~ /^-/ ))
  {
    if ( $ARGV[0] =~ /-d/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $Directory = $ARGV[0];
	    $Directory =~ s/-d//;
	  }
	else
	  {
	    shift(@ARGV);
	    $Directory = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-h/ )
      {
	$DisplayHelp=1;
      }
    if ( $ARGV[0] =~ /-u/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $RatioU = $ARGV[0];
	    $RatioU =~ s/-u//;
	  }
	else
	  {
	    shift(@ARGV);
	    $RatioU = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-p/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $Trapping = $ARGV[0];
	    $Trapping =~ s/-p//;
	  }
	else
	  {
	    shift(@ARGV);
	    $Trapping = $ARGV[0];
	  }
	print ("set trapping to ".$Trapping."\n");
      }
    if ( $ARGV[0] =~ /-C/ )
      {
	my $TmpStr;
	if (length($ARGV[0])>2)
	  {
	    $TmpStr = $ARGV[0];
	    $TmpStr =~ s/-C//;
	  }
	else
	  {
	    shift(@ARGV);
	    $TmpStr = $ARGV[0];
	  }
	@UnitCells = split(/,/,$TmpStr);	
	my $TmpLength = $#UnitCells+1;
	if ( $TmpLength != 2)
	  {
	    die ("need a lattice dimension of two: -C Lx,Ly\n");
	  }
	print ("Evaluating matrix elements for $UnitCells[0]x$UnitCells[1] unit cells\n");
      }
    shift(@ARGV);
  }

if ($DisplayHelp)
  {
    print("usage DiceLatticeModel.pl -C Lx,Ly [-p V_p] [-d dir]\n");
    print("option -C: indicate length in x- and y-directions\n");
    print("       -p: strength of pinning potential (default: 1e-4)\n");
    print("       -d: name of directory to generate files in (default: ./dice_Lx_Ly)\n");
    print("       -u: ratio of U6 to U3 on 6-fold and 3-fold connected sites of dice lattice\n");
    exit(1);
  }


# set default directory name, if not given
if ($Directory eq "")
  {
    $Directory = "./dice_".($UnitCells[0])."_".($UnitCells[1]);
  }

my $NbrCells = $UnitCells[0]*$UnitCells[1];
my $NbrSites = 2*$NbrCells;

# append lattice geometry to lattice definition
$DiceLattice.="PeriodicRepeat = ".($UnitCells[0]).",".($UnitCells[1])."\n";
$EffectiveTriangularLattice.="PeriodicRepeat = ".($UnitCells[0]).",".($UnitCells[1])."\n";
$EffectiveTriangularLattice.="NbrFlux = ".($NbrSites/2)."\n";
my $Program;
my $Have64Bits=0;
my $tmp = "";
$tmp = `status`;
if ( $tmp =~ /x86_64/ )
  {
    $Program = $Program_64;
  }
else
  {
    $Program = $Program_32;	
  }

# create directory if not existent
if ( ! -e $Directory )
  {
    system ("mkdir -p $Directory");
  }

chdir($Directory);

#calculate localized single-particle states:
for (my $i=0; $i<$NbrSites; ++$i)
  {
    # generate lattice definition
    my $LatticeFile = "DiceDoubledPhases_S".$i."_V_-".abs($Trapping)."_on_$UnitCells[0]x$UnitCells[1].dat";
    open (DEFINITION, ">$LatticeFile");
    my $CurrentDescriptor=$Descriptor."_S".$i."_V_-".abs($Trapping);
    print DEFINITION ("Descriptor = ".$CurrentDescriptor."\n");
    print DEFINITION ("LocalPotentials = ".GetSiteIndex($i).",-".abs($Trapping)."\n");
    print DEFINITION ($DiceLattice);
    close(DEFINITION);
    # run single particle calculation
    my $Command = "$Program -p 1 -L $LatticeFile -q ".(3*$NbrCells)." -c --eigenstate -n 1";
    system($Command);
  }

#write definition of effective lattice
my $LatticeFile = "DiceDoubledEffective_$UnitCells[0]x$UnitCells[1].dat";
open (DEFINITION, ">$LatticeFile");
print DEFINITION ($EffectiveTriangularLattice);
close(DEFINITION);

#generate description of interaction, here: delta interaction
my $Interaction="HilbertSpaceDimension = ".(6*$NbrCells)."\n"
  ."NbrMatrixElements = 2\n"
  ."# Matrix element 0 -> U3, Matrix element 1 -> U6\n"
  ."MatrixElements = 1 1\n"
  ."Sparse = false\n"
  ."HaveOffDiagonal = false\n"
  ."DiagonalEntries = ";
for (my $i=0; $i<$NbrCells; ++$i)
  {
    $Interaction.="0 1 0 0 0 1 ";
  }
$Interaction.="\n";
my $InteractionFile = "DeltaInteraction_$UnitCells[0]x$UnitCells[1].dat";
open (DEFINITION, ">$InteractionFile");
print DEFINITION ($Interaction);
close (DEFINITION);

my $MatrixFile = "MatrixElements_Delta_u_".$RatioU."_".$UnitCells[0]."x".$UnitCells[1].".dat";
open (MATRIX, ">$MatrixFile");

# calculate all matrix elements
for (my $Index1=0; $Index1<$NbrSites; ++$Index1)
  {
    for (my $Index2=0; $Index2<=$Index1; ++$Index2)
      {
	for (my $Index3=0; $Index3<$NbrSites; ++$Index3)
	  {
	    for (my $Index4=0; $Index4<=$Index3; ++$Index4)
	      {
		# call MatrixElement evaluator: 4,3: creation operators, 2,1 annihilation operators
		my $Command = "$MatrixProgram -c --gauge --quiet --interaction $InteractionFile ".GetLocalWavefunction($Index4)." "
		  .GetLocalWavefunction($Index3)." ".GetLocalWavefunction($Index2)." ".GetLocalWavefunction($Index1);
		my $Output = `$Command`;
		my @OutputLines = split(/\n/,$Output);
		chomp (@OutputLines);
		my $TmpLength = $#OutputLines+1;

		#print ("Output for element $Index1 $Index2 $Index3 $Index4:\n");
		#for (my $i=0; $i<$TmpLength; ++$i)
		#  {
		#    print ($OutputLines[$i]."\n");
		#  }
		
		#get matrix element for U3:
		my @Complex3 = split(/ /,$OutputLines[0]);
		if ((abs($Complex3[1])>1e-5)&&(abs(abs($Complex3[0]/$Complex3[1])-1.0)<1e-5))
		  {
		    $Complex3[0]=$Complex3[1];
		  }
		my $AbsVal = sqrt($Complex3[0]*$Complex3[0]+$Complex3[1]*$Complex3[1]);
		my $RoundedMultiples3 = sprintf("%.0f", $AbsVal*144.0);
		my $Factor = $RoundedMultiples3/144.0/$AbsVal;
		$Complex3[0]*=$Factor;
		$Complex3[1]*=$Factor;
		#get matrix element for U6:
		my @Complex6 = split(/ /,$OutputLines[1]);
		if ((abs($Complex6[1])>1e-5)&&(abs(abs($Complex6[0]/$Complex6[1])-1.0)<1e-5))
		  {
		    $Complex6[0]=$Complex6[1];
		  }
		$AbsVal = sqrt($Complex6[0]*$Complex6[0]+$Complex6[1]*$Complex6[1]);
		my $RoundedMultiples6 = sprintf("%.0f", $AbsVal*144.0);
		$Factor = $RoundedMultiples6/144.0/$AbsVal;
		$Complex6[0]*=$Factor;
		$Complex6[1]*=$Factor;

		if (($RoundedMultiples3 != 0)||($RoundedMultiples6!=0))
		  {
		    my $Real = $Complex3[0]+$RatioU*$Complex6[0];
		    if (abs($Real)<1e-8)
		      {
			$Real=0.0;
		      }
		    my $Imag = $Complex3[1]+$RatioU*$Complex6[1];
		    if (abs($Imag)<1e-8)
		      {
			$Imag=0.0;
		      }
		    printf("$Index1 $Index2 $Index3 $Index4 ".$Real." ".$Imag." "
			   .$RoundedMultiples3." ".$RoundedMultiples6."\n");
		    printf MATRIX ("$Index1 $Index2 $Index3 $Index4 ".$Real." ".$Imag."\n"); #." ".$RoundedMultiples3." ".$RoundedMultiples6."\n");
		  }
	      }
	  }
      }
  }
close (MATRIX);



# get coordinate of trapping site
#
# $_[0] = index of trapping site

sub GetSiteIndex
  {
    my $SiteIndex = $_[0];
    if ( $SiteIndex % 2 == 0)
      {
	return (6*($SiteIndex/2)+1);
      }
    else
      {
	return (6*(($SiteIndex-1)/2)+5)
      }
  }


# get filename of wavefunction at given trapping site
#
# $_[0] = index of trapping site

sub GetLocalWavefunction
  {
    my $SiteIndex = $_[0];
    return "bosons_lattice_dice_doubled_S".$SiteIndex."_V_-".abs($Trapping)."_".$UnitCells[0]."x".$UnitCells[1]."_n_1_hardcore_q_"
      .(3*$NbrCells).".0.vec";
  }


