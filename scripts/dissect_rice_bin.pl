use strict;
use warnings;

my $binfile = shift;

my $winsize = 100 * 1000;

open(IN,"<$binfile") or die $!;

my %hash_bin;
my %hash_count;
my $binsum = 0;

while(<IN>){
	chomp;
	next if($_ =~ /^#/);
	my($chr,$bin,$locinum,$aus,$ind,$tej,$trj) = split/\t/;
	my $start = $bin * $winsize + 1;
	my $end = ($bin + 1) * $winsize;
	my $type = "NA";
	if($locinum < 5){
		goto OUTPUT;
	}
	if(($aus + $ind) > 0.75){
		$type = "ss_INDx";
		if($ind > ($aus + $ind) * 0.75){
			$type = "sp_ind";
		}elsif($aus > ($aus + $ind) * 0.75){
			$type = "sp_aus";
		}
	}elsif(($tej + $trj) > 0.75){
		$type = "ss_JAPx";
		if($tej > ($tej + $trj) * 0.75){
			$type = "sp_tej";
		}elsif($trj > ($tej + $trj) * 0.75){
			$type = "sp_trj";
		}
	}else{
		$type = "admixed";
	}
	OUTPUT:
	$binsum ++;
	$hash_count{$type}{count} ++;
	$hash_bin{$chr}{$start}{end} = $end;
	$hash_bin{$chr}{$start}{type} = $type;
	print "$chr\t$start\t$end\t$type\n";
}
close IN;

=cut
my @types = sort {$hash_count{$b}{count} <=> $hash_count{$a}{count}} keys %hash_count;
my $smptype = "admixed";

if($types[0] eq "INDICA"){
	if($types[1] eq "ind"){
		$smptype = "ind";
	}elsif($types[1] eq "aus"){
		$smptype = "aus";
	}
}elsif($types[0] eq "JAPONICA"){
	if($types[1] eq "tej"){
		$smptype = "tej";
	}elsif($types[1] eq "trj"){
		$smptype = "trj";
	}
}
