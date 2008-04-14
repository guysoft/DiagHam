#!/usr/bin/perl -w

use strict 'vars';

my $CommandFile="/rscratch/gm360/latticeQHE/ParametersToLaunch";
my $LaunchedFile="/rscratch/gm360/latticeQHE/ParametersLaunched";
my $FinishedFile="/rscratch/gm360/latticeQHE/ParametersFinished";

my $Program_32="/rscratch/gm360/bin/FQHELatticeBosons";
my $Program_64="/rscratch/gm360/bin/FQHELatticeBosons_64";
my $Program_S5="/scratch/gm360/DiagHam/buildSMP/FQHE/src/Programs/FQHEOnLattice/FQHELatticeBosons";

if (!defined($ARGV[1]))
  {
    print("usage ManageBosonsOnLattice.pl Machine# #Processors [PrecalculationMemory=0]\n");
    exit(1);
  }

my $Machine;
my $Program;
my $MaxProcessor=1;
if ($ARGV[0] =~ /^[0-9]/)
  {
    $Machine="tcmpc".$ARGV[0];
    my $tmp = `ssh ${Machine} status`;
    if ( $tmp =~ /x86_64/ )
      {
	$Program = $Program_64;
	if ( $tmp =~ /Core2/ )
	  {
	    $MaxProcessor=2;
	  }
	print ("found 64-bit ".$MaxProcessor."-processor machine\n");
      }
    else
      {
	$Program = $Program_32;	
      }
  }
else
  {
    if ($ARGV[0] =~ /s5/)
      {
	$Machine="s5.tcm.phy.private";
	$MaxProcessor=8;
      }
    else
      {
	print ("Machine not recognised!");
	exit(1);
      }
  }
print ("Running on machine ".$Machine."\n");
my $Processors=" ";
if ( $ARGV[1] > 1 )
  {
    if ($ARGV[1]<=$MaxProcessor)
      {
	$Processors = " -S --processors ".$ARGV[1]." ";
      }
    else
      {
	$Processors = " -S --processors ".$MaxProcessor." ";
      }
  }

# read options file
system ("cp ".$CommandFile." ".$CommandFile.".save");
open(MYCOMMANDS, $CommandFile) or die("Error: cannot open file '$CommandFile'\n");
open(LEFTCOMMANDS, ">${CommandFile}.remain") or die("Error: cannot open file '$CommandFile.remain'\n");
my @Lines = <MYCOMMANDS>;
close (MYCOMMANDS);
my $TmpLine;
my $Parameters="";
foreach $TmpLine (@Lines)
  {
    if (( $TmpLine =~ /^#/) || ( $Parameters ne "" ))
      {
	print LEFTCOMMANDS ($TmpLine);
      }
    else
      {
	chomp($TmpLine);
	$Parameters=$TmpLine;
      }
  }
if ( $Parameters eq "" )
  {
    print ("Error: No parameters found\n");
    exit(1);
  }
# extract individual parameters
my @AllParam=split(/\t\s*/,$Parameters);
my $paramR = $AllParam[0];
my $paramT = $AllParam[1];
my $paramLx = $AllParam[2];
my $paramLy = $AllParam[3];
my $paramQ = $AllParam[4];
my $paramU = $AllParam[5];
my $paramN1 = $AllParam[6];
my $paramN2 = $AllParam[7];

my $StatesDir = "/rscratch/gm360/latticeQHE/states/n_${paramR}_${paramT}/";
my $SpecDir = "/rscratch/gm360/latticeQHE/spectra/n_${paramR}_${paramT}/";
my $WorkDir="";
my $EigenVectors="";
my $QString="q";
if ($paramQ>0)
  {
    $QString="q_${paramQ}";
  }
my $NbrBosons = $paramLx*$paramLy*$paramR/$paramT;
my $OutputName = "bosons_lattice_n_${NbrBosons}_x_${paramLx}_y_${paramLy}_u_${paramU}_${QString}.dat";

if ( ! -e $SpecDir )
  {
    system ("mkdir ${SpecDir}");
  }
if ( ! -e $StatesDir )
  {
    system ("mkdir ${StatesDir}");
  }
if ( $paramN1 == $paramN2 ) # need to calculate only states
  {
    $WorkDir = $StatesDir;
    $EigenVectors= "--eigenstate -n ${paramN1}";
  }
else
  {
    if ( -e $SpecDir.$OutputName )
      {
	$WorkDir = $StatesDir;
	$EigenVectors= "--eigenstate -n ${paramN2}";
      }
    else
      {
	$WorkDir = $SpecDir;
	$EigenVectors= "-n ${paramN1}";
	system ("touch ${SpecDir}${OutputName}");
	print LEFTCOMMANDS ($Parameters."\tLS\n");
      }
  }
close(LEFTCOMMANDS);
system ("cp ".$CommandFile.".remain ".$CommandFile);

my $Command = "ssh $Machine \"cd ${WorkDir}; nohup nice -n15 $Program -p ${NbrBosons} -x ${paramLx} -y ${paramLy} -q $paramQ -u $paramU ${Processors} ${EigenVectors} --show-itertime > log_p_${NbrBosons}_u_${paramU} \" &";

open(LAUNCHEDCOMMANDS, ">>${LaunchedFile}") or die("Error: cannot open file '$LaunchedFile'\n");
print LAUNCHEDCOMMANDS ($Parameters." ".$Command."\n");
close(LAUNCHEDCOMMANDS);

print $Command."\n"; # launch here after testing!

open(FINISHEDCOMMANDS, ">>${FinishedFile}") or die("Error: cannot open file '$FinishedFile'\n");
print FINISHEDCOMMANDS ($Parameters." ".$Command."\n");
close (FINISHEDCOMMANDS);


