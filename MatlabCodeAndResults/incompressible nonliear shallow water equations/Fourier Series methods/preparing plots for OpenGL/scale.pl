#!/usr/bin/perl -w
use strict;
use File::Slurp;

# A crufty little program which will, for at least some EPS images,
# scale them by a given factor (e.g., 0.5 for half size).
# Note that it doesn't give up if it can't understand the EPS.  It
# just carries on regardless and hopes for the best.
#
# I think there is a bug in the scaling.  I'll fix it sometime.
#
# This program is public domain.  Do what you want with it.  Don't sue me.
#
# --Qef, 2004-03-04.

die "Usage: $0 scale-factor input-filename output-filename\n"
  unless @ARGV == 3;

my ($scale, $input_filename, $output_filename) = @ARGV;
$input_filename = \*STDIN if $input_filename eq '-';
$output_filename = \*STDOUT if $output_filename eq '-';

local $_ = read_file($input_filename);

# Adjust bounding boxes.
s[^(%%(?:Page)?BoundingBox:)\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)]
 [join ' ', $1, map { int($_ * $scale) } $2, $3, $4, $5 ]meg;

# Insert scaling command.
s[^(%%Page:.*\n(?:%%.*\n)*)][$1\ngsave $scale $scale scale\n]m
or
s[^(%!.*\n(?:%%.*\n)*)][$1\ngsave $scale $scale scale\n]
or
s[^][gsave $scale $scale scale\n];

# Restore from the scaling.
s[^(%%PageTrailer)][grestore\n$1]m
or
s[$][\ngrestore\n];

write_file($output_filename, $_);

# vi:ts=4 sw=4 expandtab
