#!/usr/bin/perl -w

use strict 'vars';

use XML::Simple;


unless (open (INFILE, "docs/built_in_progams/header.tex"))
  {
    die ("can't open ".$ARGV[0]."\n");
  }
my $TmpLine;
foreach $TmpLine (<INFILE>)
  {
    print $TmpLine;
  }
close (INFILE);


my $XMLContent = XMLin($ARGV[0]);

my $Key;
my $Value;

if (defined($$XMLContent{'name'}))
  {
    print "{\\bf{\\underline{".$$XMLContent{'name'}."}:}}\\\\\n";
  }
if (defined($$XMLContent{'shortdesc'}))
  {
    print $$XMLContent{'shortdesc'}."\\\\\n";
  }
if (defined($$XMLContent{'location'}))
  {
    print "{\\bf{\\underline{Location}:}} ".$$XMLContent{'location'}."\\\\\n";
  }
if (defined($$XMLContent{'option'}))
  {
    print "{\\bf{\\underline{Options}:}}\\\\

\\begin{description}\n";
    my $Key;
    my $Value;
    my $OptionRef = $$XMLContent{'option'};
    while (($Key, $Value) = each (%$OptionRef))
      {
	print "\\item[";
	if (defined($$Value{'short'}))
	  {
	    print "-".$$Value{'short'}.", ";
	  }
	print "--".$Key."] ".$$Value{'description'};
	if (defined($$Value{'default'}))
	  {
	    print " (default value is set to ".$$Value{'default'}.")";
	  }
	print "\n";
      }
    print "\\end{description}\n";
  }
if (defined($$XMLContent{'longdesc'}))
  {
    print "{\\bf{\\underline{Description}:}}\\\\\n";
    print $$XMLContent{'longdesc'}."\\\\\n";
  }
if (defined($$XMLContent{'accuracy'}))
  {
    print "{\\bf{\\underline{Accuracy}:}}\\\\\n";
    print $$XMLContent{'accuracy'}."\\\\\n";
  }
if (defined($$XMLContent{'remarks'}))
  {
    print "{\\bf{\\underline{Remarks}:}}\\\\\n";
    print $$XMLContent{'remarks'}."\\\\\n";
  }
if (defined($$XMLContent{'author'}))
  {
    my $Key;
    my $Value;
    my $OptionRef = $$XMLContent{'author'};
    while (($Key, $Value) = each (%$OptionRef))
      {
	print $Key."\\\\\n";
      }
  }


unless (open (INFILE, "docs/built_in_progams/footer.tex"))
  {
    die ("can't open ".$ARGV[0]."\n");
  }
foreach $TmpLine (<INFILE>)
  {
    print $TmpLine;
  }
close (INFILE);

