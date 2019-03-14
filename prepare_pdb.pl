#!/usr/bin/perl

# Matthias Schmidt
# 12/5/11
# matthias.rene.schmidt@gmail.com

# Renumbers a file correctly
# Adds chains, if missing

# Todo: Renumber independent of input numbers
# Todo: Don't subtract if numbering starts at 1,
#       this need manual script changing now.

#use warnings;
use strict;

unless ($ARGV[2] )
{
	die "Usage:script.pl <pdbfile> <truncated chain-length> <n-truncation> <ligand> eg.
	script.pl run2-1.pdb 335 34 PIP or
	script.pl run2-1.pdb 297 42 POP \n";
	
}
unless ($ARGV[3])
{
	$ARGV[3] = "none";
}

my $pdbfile = $ARGV[0];
my $length = $ARGV[1];
my $trunc = $ARGV[2];
my $ligand = $ARGV[3];
my $chain = "@"; # the letter that comes before "A".
my $atom_prev = "";
my $atom_endofchain = "";

open (PDB, "$pdbfile") || die "Cannot open $pdbfile\n";

while ( my $line = <PDB>)
{
	unless ($line =~ m/^ATOM/ and not $line =~ m/$ligand/ 
				  and not $line =~ m/NA\+/ 
				  and not $line =~ m/K\+/ 
				  and not $line =~ m/SOL/ 
				  and not $line =~ m/POP/ )
	{
		print $line;	
		next;
	}
	
	else
	{
		my $aa = substr $line, 22, 4;
		
		# if numbering starts at 1
		$aa = $aa % $length;
		
		# if numbering does not start at 1
		# $aa = ($aa - $trunc) % $length;
			
		if ($aa == 0) 
		{
			$aa = $length;
			$atom_endofchain = substr $line, 6, 5;
		}
		if ($aa == 1 and $atom_prev eq $atom_endofchain)
		{
			$chain = chr(ord($chain)+1);
		}
		
		$aa += $trunc;
				
		substr( $line, 22, 4 ) = sprintf( "%4s", $aa);
		substr( $line, 21, 1 ) = $chain;
		print "$line";
		
		$atom_prev = substr $line, 6, 5;
				
	}

}
close (PDB);
