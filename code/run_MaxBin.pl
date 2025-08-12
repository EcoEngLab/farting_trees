#!/usr/bin/env perl

# Wrapper script for MaxBin2 to provide the expected run_MaxBin.pl interface
# This script acts as a bridge between the expected Perl script interface
# and the actual MaxBin binary provided by the conda package

use strict;
use warnings;
use File::Basename;

# Get the conda environment path
my $conda_prefix = $ENV{'CONDA_PREFIX'} || die "Error: CONDA_PREFIX not set. Make sure you're in the maxbin2_env environment.\n";

# Path to MaxBin binary
my $maxbin_binary = "$conda_prefix/bin/MaxBin";

# Check if MaxBin binary exists
unless (-x $maxbin_binary) {
    die "Error: MaxBin binary not found at $maxbin_binary\n";
}

# Pass all arguments to the MaxBin binary
exec($maxbin_binary, @ARGV) or die "Error: Could not execute MaxBin: $!\n";
