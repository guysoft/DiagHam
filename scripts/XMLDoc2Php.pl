#!/usr/bin/perl -w

use strict 'vars';

use XML::Simple;
use Data::Dumper;
  
my $XMLContent = XMLin($ARGV[0]);

print Dumper($XMLContent);
