#!/usr/bin/perl -w
#
# script for optimizing overlaps for vectors associated with given protocol files
#
use strict 'vars';
use File::stat;
use Math::Complex;


# hardwire which state to look at
my $Program_32="FQHELatticeBosons";
my $Program_64="FQHELatticeBosons_64";
my $OverlapExe="GenericOverlap";

my $CalculateVectors=0;
my $NbrGrid=20;
my $GridString="";
my $ReferenceString="0,0,1,1";
my $Degeneracy=1;
my $Memory=1000;
my $NbrCalculate=1;

while( (defined($ARGV[0])&&$ARGV[0] =~ /^-/ ))
  {
    if ( $ARGV[0] =~ /-d/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $Degeneracy = $ARGV[0];
	    $Degeneracy =~ s/-d//;
	  }
	else
	  {
	    shift(@ARGV);
	    $Degeneracy = $ARGV[0];
	  }
      }
      if ( $ARGV[0] =~ /-g/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $GridString = $ARGV[0];
	    $GridString =~ s/-g//;
	  }
	else
	  {
	    shift(@ARGV);
	    $GridString = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-m/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $Memory = $ARGV[0];
	    $Memory =~ s/-m//;
	  }
	else
	  {
	    shift(@ARGV);
	    $Memory = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-n/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $NbrGrid = $ARGV[0];
	    $NbrGrid =~ s/-n//;
	  }
	else
	  {
	    shift(@ARGV);
	    $NbrGrid = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-r/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $ReferenceString = $ARGV[0];
	    $ReferenceString =~ s/-r//;
	  }
	else
	  {
	    shift(@ARGV);
	    $ReferenceString = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-s/ )
      {
	if (length($ARGV[0])>2)
	  {
	    $NbrCalculate = $ARGV[0];
	    $NbrCalculate =~ s/-s//;
	  }
	else
	  {
	    shift(@ARGV);
	    $NbrCalculate = $ARGV[0];
	  }
      }
    if ( $ARGV[0] =~ /-c/ )
      {
	$CalculateVectors=1;
	print("Will calculate missing vectors!\n");
      }
    shift(@ARGV);
  }

if (!defined($ARGV[0]))
  {
    print("usage ChernNumberLattice.pl [-g GRIDPOINTS] [-n NbrPoints] [-c] [-r s1x,s1y,s2x,s2y] basename_*SX_SY*\n");
    print("option -g: list of discrete points used for both x- and y- directions\n");
    print("       -n: number of gridpoints to be used\n");    
    print("       -r: reference points A,B as theta1A,theta2A,theta1B,theta2B (default 0,0,1,1)\n");
    print("       -c: optionally calculate missing vector files\n");
    print("       -m: memory for precalculations when calculating vectors\n");
    print("       -s: number of states to be calculated at each point\n");
    exit(1);
  }


my @GridPoints;

if (length($GridString)>0)
  {
    @GridPoints = split (/,/, $GridString);
    $NbrGrid = $#GridPoints + 1;
  }
else
  {
    my $Sep=2.0/$NbrGrid;
    for (my $i=0; $i<$NbrGrid; ++$i)
      {
	$GridPoints[$i]=-1.0+($i+1)*$Sep;
      }
  }
if ($NbrCalculate<$Degeneracy)
  {
    $NbrCalculate=$Degeneracy;
  }

my $Program;
my $Have64Bits=0;
my $tmp = "";
# $tmp = `status`;
if ( $tmp =~ /x86_64/ )
  {
    $Program = $Program_64;
  }
else
  {
    $Program = $Program_32;	
  }

my $TmpFile;
foreach $TmpFile (@ARGV)
  {
    if ($TmpFile =~ m/bosons\_lattice/)
      {
	print ("Analyzing states with base name ".$TmpFile."\n");
	&AnalyzeChern($TmpFile);
      }
  }


# calculate gauge and prepare output for plotting
#
# $_[0] = base file name

sub AnalyzeChern
  {
    my $BaseName = $_[0];
    my $HardCore;
    my $N;
    my $x;
    my $y;
    my $u;
    my $q;
    if ($BaseName =~ m/hardcore/)
      {

	$BaseName =~ /n\_(\d+)\_x\_(\d*)\_y\_(\d*)\_.*\_q\_(\d*)/;
	$N = $1;
	$x = $2;
	$y = $3;
	$u = 0;
	$HardCore=1;
	$q = $4;
      }
    else
      {
	$BaseName =~ /n\_(\d+)\_x\_(\d*)\_y\_(\d*)\_u\_(-*\d*[\.]*\d*)\_.*\_q\_(\d*)/;
	$N = $1;
	$x = $2;
	$y = $3;
	$u = $4;
	$HardCore=0;
	$q = $5;
      }
    my $TotalSolenoid="";
    my $SolenoidX=0.0;
    my $SolenoidY=0.0;
#    if ($BaseName =~ m/\_s\_/)
#      {	
#	$BaseName =~ /\_s\_(-*\d*[\.]*\d*e*-*\d*)\_(-*\d*[\.]*\d*e*-*\d*)/;
#	$SolenoidX = $1;
#	$SolenoidY = $2;
#	$TotalSolenoid="--solenoid-flux $SolenoidX,$SolenoidY";
#      }
    my $Interaction;
    if ( $HardCore == 1)
      {
	print "n = ".$N."  x = ".$x."  y = ".$y."  q = ".$q."  (hardcore bosons)\n";
	$Interaction ="-c";
      }
    else
      {
	print "n = ".$N."  x = ".$x."  y = ".$y."  q = ".$q."  u = ".$u."\n";
	$Interaction ="";
      }
    
    my @ReferenceVals = split (/,/, $ReferenceString);

    my $RefS1x = $ReferenceVals[0];
    my $RefS1y = $ReferenceVals[1];
    my $RefS2x = $ReferenceVals[2];
    my $RefS2y = $ReferenceVals[3];

    my $CommandLine = "$Program -p $N -x $x -y $y $Interaction -q $q -n $NbrCalculate -m $Memory";

    TestVectors ($BaseName, $RefS1x, $RefS1y, $Degeneracy, $CalculateVectors, $CommandLine);
    TestVectors ($BaseName, $RefS2x, $RefS2y, $Degeneracy, $CalculateVectors, $CommandLine);

    my $LogFileName = GetVectorName($BaseName, 0, 0, 0);
    $LogFileName =~ s/0\.vec/cn/;
    open (LOGFILE, ">$LogFileName");

    print LOGFILE ("# Evaluation of Chern-Number for state at N = $N, x = $x, y = $y, q = $q  ($Interaction)\n");
    print LOGFILE ("# Theta_x\tTheta_y\tOmega_x\tOmega_y\n");

    my $GridX;
    my $GridY;
    for ($GridX=0; $GridX<$NbrGrid;++$GridX)
      {
	for ($GridY=0; $GridY<$NbrGrid;++$GridY)
	  {
	    my $SolenoidX = $GridPoints[$GridX];
	    my $SolenoidY = $GridPoints[$GridY];
	    print ("Analysing point $SolenoidX,$SolenoidY \n");

	    # make sure we have all necessary data
	    TestVectors ($BaseName, $SolenoidX, $SolenoidY, $Degeneracy, $CalculateVectors, $CommandLine);
	    my $Alpha;
	    my $Beta;
	    my $OvlCommand="GenericOverlap --quiet -s -c ";
	    my @Matrix1;
	    my @Matrix2;
	    my @Overlaps;
	    for ($Alpha=0; $Alpha<$Degeneracy;++$Alpha)
	      {
		my @Column1;
		my @Column2;
		my $VectorA1=GetVectorName($BaseName, $RefS1x, $RefS1y, $Alpha);
		my $VectorA2=GetVectorName($BaseName, $SolenoidX, $SolenoidY, $Alpha);
		for ($Beta=0; $Beta<$Degeneracy;++$Beta)
		  {
		    my $VectorB1=GetVectorName($BaseName, $SolenoidX, $SolenoidY, $Beta);
		    my $VectorB2=GetVectorName($BaseName, $RefS2x, $RefS2y, $Beta);
		    my $TmpCmd = $OvlCommand." ".$VectorA1." ".$VectorB1;
		    my $OvlString= `$TmpCmd`;
		    @Overlaps = split(/ /,$OvlString);
		    my $z = cplx($Overlaps[0],$Overlaps[1]);
		    push(@Column1,$z);
		    $TmpCmd = $OvlCommand." ".$VectorA2." ".$VectorB2;
		    $OvlString= `$TmpCmd`;
		    @Overlaps = split(/ /,$OvlString);
		    $z = cplx($Overlaps[0],$Overlaps[1]);
		    push(@Column2,$z);
		  }
		push(@Matrix1,\@Column1);
		push(@Matrix2,\@Column2);
	      }
	    my $Det1 = det (\@Matrix1);
	    my $Det2 = det (\@Matrix2);
	    my $Arg1 = arg($Det1);
	    my $Arg2 = arg($Det2);
	    print ("det1 = ".$Det1."\n");
	    print ("det2 = ".$Det2."\n");
	    my $TotalX = sin($Arg1+$Arg2);
	    my $TotalY = cos($Arg1+$Arg2);
	    print LOGFILE ("$SolenoidX\t$SolenoidY\t$TotalX\t$TotalY\n");
	  }
      }
    close(LOGFILE);
  } # end of AnalyzeChern



# test if all required vectors at a given point are present
# if not, calculate
sub TestVectors {
  my $BaseName = $_[0];
  my $SolenoidX = $_[1];
  my $SolenoidY = $_[2];
  my $Degeneracy = $_[3];
  my $Calculate = $_[4];
  my $Command = $_[5];

  my $HavePoint=1;
  for (my $i=0; $i<$Degeneracy; ++$i)
    {
      my $VectorName = GetVectorName($BaseName, $SolenoidX, $SolenoidY, $i);
      if ( ! -e $VectorName )
	{
	  print ("State $VectorName not found!\n");
	  $HavePoint=0;
	}
      else
	{
	  print ("State $VectorName found!\n");
	}
    }
  if ( $HavePoint == 0)
    {
      if ( $Calculate == 1 )
	{
	  print ("Missing vectors at $SolenoidX,$SolenoidY ... recalculating\n");
	  print ("Command executed: ".$Command." --eigenstate --show-itertime --solenoid-flux $SolenoidX,$SolenoidY\n");
	  system($Command." --eigenstate --show-itertime --solenoid-flux $SolenoidX,$SolenoidY");
	}
      else
	{
	  die("Missing vectors at $SolenoidX,$SolenoidY ...\nAborting - please supply the missing files manually!\n");
	}
    }
}

# generate the name of a vector-file
#
sub GetVectorName
  {
    my $BaseName=$_[0];
    my $SolenoidX = $_[1];
    my $SolenoidY = $_[2];
    my $ID = $_[3];

    my $VectorName = $BaseName.".".$ID.".vec";
    $VectorName =~ s/XX/$SolenoidX/;
    $VectorName =~ s/YY/$SolenoidY/;
    if (($SolenoidX==0)&&($SolenoidY==0))
      {
	$VectorName =~ s/\_s\_0\_0//;
      }
    return $VectorName;
  }


# calculate the determinant of a matrix given as a double list
# call this function by reference: det (\@matrix)
sub det {
    my $matrix = shift;
    my $size   = $#{ $matrix } + 1;

    foreach (@$matrix) {
#        print ("Line ".$_.", size ".($#{ $_ } + 1)."\n");
#	foreach (@$_)
#	  {
#	    print("Entry".$_."\n");
#	  }
        die "det(Matrix) requires n x n matrix!" if @$_ != $size;
        }

    return $matrix->[0][0] if $size == 1;
    return $matrix->[0][0] * $matrix->[1][1] - $matrix->[1][0] * $matrix->[0][1]
      if $size == 2;
    return _det_helper( $matrix, $size );
}

sub _det_helper {
    my $matrix = shift;
    my $size   = shift;

    return $matrix->[0][0] * $matrix->[1][1] * $matrix->[2][2] + $matrix->[1][0]
      * $matrix->[2][1] * $matrix->[0][2] + $matrix->[2][0] * $matrix->[0][1] *
      $matrix->[1][2] - $matrix->[0][2] * $matrix->[1][1] * $matrix->[2][0] -
      $matrix->[1][2] * $matrix->[2][1] * $matrix->[0][0] - $matrix->[2][2] *
      $matrix->[0][1] * $matrix->[1][0]
      if $size == 3;

    my $det;
    foreach ( 0 .. $size - 1 ) {
        if ( $_ % 2 ) {
            $det -=
              $matrix->[0][$_] *
              _det_helper( _matrix_slice( $matrix, 0, $_ ), $size - 1 );
        }
        else {
            $det +=
              $matrix->[0][$_] *
              _det_helper( _matrix_slice( $matrix, 0, $_ ), $size - 1 );
        }
    }
    return $det;
}

sub _matrix_slice {
    my $matrix = shift;
    my $x      = shift;
    my $y      = shift;

    return [ map { [ @{$_}[ 0 .. $y - 1, $y + 1 ... $#$_ ] ] }
          @{$matrix}[ 0 .. $x - 1, $x + 1 .. $#$matrix ] ];
}
