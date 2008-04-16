#!/usr/bin/perl -w

use strict 'vars';

my $CommandFile="/rscratch/gm360/latticeQHE/Parameters";
my $LaunchedFile="/rscratch/gm360/latticeQHE/Parameters.launch";
my $FinishedFile="/rscratch/gm360/latticeQHE/Parameters.finish";

my $Program_32="/rscratch/gm360/bin/FQHELatticeBosons";
my $Program_64="/rscratch/gm360/bin/FQHELatticeBosons_64";
my $Program_S5="/scratch/gm360/DiagHam/buildSMP/FQHE/src/Programs/FQHEOnLattice/FQHELatticeBosons";

if (!defined($ARGV[1]))
  {
    print("usage ManageBosonsOnLattice.pl Machine# #Processors [PrecalculationMemory=0] [ParameterFile]\n");
    exit(1);
  }

my $Machine;
my $Program;
my $Memory=0;
my $MaxProcessor=1;
my $Have64Bits=0;
my $RunLocal=0;
if ($ARGV[0] =~ /^[0-9]/)
  {
    $Machine="ssh tcmpc".$ARGV[0];
    my $tmp = `${Machine} status`;
    if ( $tmp =~ /x86_64/ )
      {
	$Program = $Program_64;
	if ( $tmp =~ /Core2/ )
	  {
	    $MaxProcessor=2;
	  }
	print ("found 64-bit ".$MaxProcessor."-processor machine\n");
	$Have64Bits=1;
      }
    else
      {
	$Program = $Program_32;	
      }
  }
else
  {
    if ($ARGV[0] =~ /^s5/)
      {
	$Machine="ssh s5.tcm.phy.private";
	$Program = $Program_S5;	
	$MaxProcessor=8;
	$Have64Bits=1;
      }
    else
      {
	if ($ARGV[0] =~ /locals5/)
	  {
	    $Machine="";
	    $Program = $Program_S5;	
	    $MaxProcessor=8;
	    $Have64Bits=1;
	    $RunLocal=1;
	  }
	else
	  {
	    print ("Machine not recognised!");
	    exit(1);
	  }
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
if (defined($ARGV[2]))
  {
    $Memory=$ARGV[2];
  }
if (defined($ARGV[3]))
  {
    $CommandFile=$ARGV[3];
    $LaunchedFile="${CommandFile}.launch";
    $FinishedFile="${CommandFile}.finish";
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
my $NbrBosons = $paramLx*$paramLy*$paramR/$paramT;

if ($Have64Bits==0)
  {
    if ($paramLx*$paramLy+$NbrBosons > 32)
      {
	print ("Need 64 bit processor for next command\n");
	system ("cp ".${CommandFile}.".save ".$CommandFile);
	exit(1);
      }	
  }
else
  {
    if ($paramLx*$paramLy+$NbrBosons > 64)
      {
	print ("Hilbert space cannot be coded in a word of 64 bits!\n");
	exit(1);
      }	
  }

my $StatesDir = "/rscratch/gm360/latticeQHE/states/n_${paramR}_${paramT}/";
my $SpecDir = "/rscratch/gm360/latticeQHE/spectra/n_${paramR}_${paramT}/";
if ($RunLocal==1)
  {
    $StatesDir = "states/n_${paramR}_${paramT}/";
    $SpecDir = "spectra/n_${paramR}_${paramT}/";
  }
my $WorkDir="";
my $EigenVectors="";
my $QString="q";
if ($paramQ>0)
  {
    $QString="q_${paramQ}";
  }
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

my $Command = "$Machine \"cd ${WorkDir}; touch log_p_${NbrBosons}_u_${paramU}; nohup nice -n15 $Program -p ${NbrBosons} -x ${paramLx} -y ${paramLy} -q $paramQ -u $paramU ${Processors} ${EigenVectors} --show-itertime >> log_p_${NbrBosons}_u_${paramU} \" &";

open(LAUNCHEDCOMMANDS, ">>${LaunchedFile}") or die("Error: cannot open file '$LaunchedFile'\n");
print LAUNCHEDCOMMANDS ($Parameters." ".$Command."\n");
close(LAUNCHEDCOMMANDS);

print ("running: ".$Command."\n"); # launch here after testing!
system ($Command);

open(FINISHEDCOMMANDS, ">>${FinishedFile}") or die("Error: cannot open file '$FinishedFile'\n");
print FINISHEDCOMMANDS ($Parameters." ".$Command."\n");
close (FINISHEDCOMMANDS);


