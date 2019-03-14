#!/usr/bin/perl

#use warnings;
#use strict;

unless ($ARGV[3] )
{
	die "Usage:script.pl <pdbfile> <bfactorfile> <indextoread> <outpdb>\n";
}

my $pdbfile = $ARGV[0];
my $bfact = $ARGV[1];
my $col = $ARGV[2];
my $outpdb = $ARGV[3];

my @bfactdata = ();
my $temp = ();

`gmx editconf -f $pdbfile -o labelZ.gro -label Z`;
`gmx editconf -f labelZ.gro -o labelZ.pdb`;

open LABZ, "<", 'labelZ.pdb' or die;
open NOLAB, ">", 'nolabel.pdb' or die;

my $nochain = 'nolabel.pdb';

my @z = <LABZ>;
s/ Z/  /g for @z;
print NOLAB @z;

close (LABZ);
close (NOLAB);

open (BFACT, "$bfact") || die "Cannot open $bfact\n";

my $i = 0;
while ( my $line = <BFACT> )
{
	next if ( $line =~ m/^#/ );
	my @temp = split ( /,/, $line );
	$bfactdata[$temp[0]] = $temp[$col];
	$i++;		
}
close (BFACT);


open (PDB, "$nochain") || die "Cannot open $nochain\n";
open (OUT, "> $outpdb") || die "Cannot open $outpdb\n";

while ( my $line = <PDB>)
{
	unless ($line =~ m/^ATOM/ )
	{
		print OUT $line;
		next;
	}
	if ($line =~ m/DPP|POP|PPP|DLP|PLP|DDM|BOG|PIP|CHO|DHP|DMP|DOP|W/ )
	{
		print OUT $line;
	}
	else
	{
		my @temp = split ( /\s+/, $line );
		my $string= sprintf( "%.3f", $bfactdata[$temp[4]]);
		substr( $line, 61, 5 ) = "$string";
		print OUT $line;
	}

}

system("rm -r nolabel*");
system("rm -r label*");
close (PDB);
close (OUT);
