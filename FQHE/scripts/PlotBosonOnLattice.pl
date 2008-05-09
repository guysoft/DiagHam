#!/usr/bin/perl -w
#
# script for analysis of eigenvector files designed to run on TCM group PC's
#
use strict 'vars';


if (!defined($ARGV[0]))
  {
    print("usage PlotBosonOnLattice.pl protocol1.gs protocol2.gs ...\n");
    exit(1);
  }


my @ListFiles;
my $TmpFile;

foreach $TmpFile (@ARGV)
  {
    if ($TmpFile =~ /bosons\_lattice.*\.gs/)
      {
	push (@ListFiles, $TmpFile);
      }
    else
      {
	print ("$TmpFile not recognized as a protocol file with ending .gs\n");
      }
  }

my $TmpFileName = "tmp".time().".xm";

open (OUTFILE, ">$TmpFileName");

print OUTFILE ("with g0\nview 0.15, 0.535, 0.60, 0.85\nlegend 0.4, 0.85\n");
print OUTFILE ("with g1\nview 0.7, 0.535, 1.15, 0.85\nyaxes scale Logarithmic\nlegend off\ng1 on\n");
print OUTFILE ("with g2\nview 0.15, 0.325, 0.60, 0.5\nxaxis ticklabel off\nlegend off\ng2 on\n");
print OUTFILE ("with g3\nview 0.15, 0.15, 0.60, 0.325\nlegend off\ng3 on\n");
print OUTFILE ("xaxis label \"n\\s\\xf\"\n");
print OUTFILE ("with g4\nview 0.7, 0.325, 1.15, 0.5\nxaxis ticklabel off\nlegend off\ng4 on\n");
print OUTFILE ("with g5\nview 0.7, 0.15, 1.15, 0.325\nlegend off\ng5 on\n");
print OUTFILE ("xaxis label \"n\\s\\xf\"\n");

print OUTFILE ("with string\nstring on\nstring 0.175, 0.81\nstring def \"Gap \\xD\"\n");
print OUTFILE ("with string\nstring on\nstring 0.71, 0.81\nstring def \"\\xr\\s0\\N/r\\s1\"\n");
print OUTFILE ("with string\nstring on\nstring 0.175, 0.46\nstring def \"\\xr\\s\\f{}ave\"\n");
print OUTFILE ("with string\nstring on\nstring 0.175, 0.29\nstring def \"\\xN(r\\s\\f{}ave\\N)/N\\se\"\n");
print OUTFILE ("with string\nstring on\nstring 0.71, 0.46\nstring def \"d(GS)\"\n");
print OUTFILE ("with string\nstring on\nstring 0.71, 0.29\nstring def \"|K|\\sGS\"\n");

foreach $TmpFile (@ListFiles)
  {
    print ("Adding ".$TmpFile."\n");
    $TmpFile =~ /n\_(\d+)\_x\_(\d*)\_y\_(\d*)\_u\_(-*\d*[\.]*\d*)\_/;
    my $N = $1;
    my $x = $2;
    my $y = $3;
    my $u = $4;
    my $Factor=1.0/$x/$y;
    print OUTFILE ("READ BLOCK \"${TmpFile}\"\n");
    print OUTFILE ("with g0\n");
    print OUTFILE ("BLOCK xy \"1:4\"\n");
    print OUTFILE ("S_ legend \"N=${N}, u=${u}\"\n");
    print OUTFILE ("S_ comment \"${TmpFile} :: 1:4\"\n");
    print OUTFILE ("S_.x = S_.x*${Factor}\n");

    print OUTFILE ("with g1\n");
    print OUTFILE ("BLOCK xy \"1:5\"\n");
    print OUTFILE ("S_ legend \"N=${N}, u=${u}\"\n");
    print OUTFILE ("S_ comment \"${TmpFile} :: 1:5\"\n");
    print OUTFILE ("S_.x = S_.x*${Factor}\n");

    print OUTFILE ("with g2\n");
    print OUTFILE ("BLOCK xy \"1:6\"\n");
    print OUTFILE ("S_ legend \"N=${N}, u=${u}\"\n");
    print OUTFILE ("S_ comment \"${TmpFile} :: 1:6\"\n");
    print OUTFILE ("S_.x = S_.x*${Factor}\n");

    print OUTFILE ("with g3\n");
    print OUTFILE ("BLOCK xy \"1:7\"\n");
    print OUTFILE ("S_ legend \"N=${N}, u=${u}\"\n");
    print OUTFILE ("S_ comment \"${TmpFile} :: 1:7\"\n");
    print OUTFILE ("S_.x = S_.x*${Factor}\n");
    print OUTFILE ("S_.y = S_.y/${N}\n");

    print OUTFILE ("with g4\n");
    print OUTFILE ("BLOCK xy \"1:9\"\n");
    print OUTFILE ("S_ legend \"N=${N}, u=${u}\"\n");
    print OUTFILE ("S_ comment \"${TmpFile} :: 1:9\"\n");
    print OUTFILE ("S_.x = S_.x*${Factor}\n");

    print OUTFILE ("with g5\n");
    print OUTFILE ("BLOCK xy \"1:8\"\n");
    print OUTFILE ("S_ legend \"N=${N}, u=${u}\"\n");
    print OUTFILE ("S_ comment \"${TmpFile} :: 1:8\"\n");
    print OUTFILE ("S_.x = S_.x*${Factor}\n");
  }

for (my $i=0; $i<6; ++$i)
  {
    print OUTFILE ("with G${i}\n");
    print OUTFILE ("G${i} AUTOSCALE TYPE SPEC\n");
    print OUTFILE ("WORLD XMIN 0\n");
    print OUTFILE ("WORLD XMAX 0.5\n");
    print OUTFILE ("AUTOSCALE YAXES\n");
    print OUTFILE ("AUTOTICKS\n");
  }

close(OUTFILE);

system("xmgrace -batch ${TmpFileName} -nosafe -noask &");
system("sleep 2");

system ("rm $TmpFileName");


