use strict;
use warnings;

my($ratiofile,$bedfile,$outfile) = @ARGV;
my $usage = "USAGE:\nperl $0 <sample haplotype stdratio file> <gene bed file> <outfile>\n";

die $usage unless(@ARGV == 3);

=ratiofile
#TAG    RANK    aus     bas     ind     int     jap     tej     trj
LOC_Os01g01010  0       0.040   0.161   0.137   0.151   0.167   0.175   0.169
LOC_Os01g01019  0       0.038   0.160   0.136   0.153   0.168   0.175   0.169
LOC_Os01g01030  0       0.040   0.160   0.137   0.153   0.167   0.175   0.168
=cut

=bedfile
Chr1    2903    10817   LOC_Os01g01010
Chr1    11218   12435   LOC_Os01g01019
Chr1    12648   15915   LOC_Os01g01030
Chr1    16292   20323   LOC_Os01g01040
=cut

my $window = 100 * 1000;
my %hash_bed;
my %hash_window;
my %hash_group;
my @groupnames;
my @null_groups;
my @chrs;

open(IN,"<$bedfile") or die $!;
while(<IN>){
	chomp;
	my($chr,$start,$end,$tag) = split/\t/;
	next if(exists $hash_bed{$tag});
	$hash_bed{$tag}{chr} = $chr;
	$hash_bed{$tag}{pos} = $start;
}
close IN;

open(IN,"<$ratiofile") or die $!;
while(<IN>){
	chomp;
	my($tag,$rank,@groups) = split/\t/;
	if($tag =~ /^#/){
		@groupnames = @groups;
		for(my $i = 0; $i < @groups; $i++){
			$null_groups[$i] = 0;
		}
		next;
	}
	next if($rank eq "-");
	my $chr = $hash_bed{$tag}{chr};
	my $pos = $hash_bed{$tag}{pos};
	my $windowrank = int(($pos - 1)/$window);
	unless(exists $hash_window{$chr}){
		push @chrs, $chr;
	}
	unless(exists $hash_window{$chr}{$windowrank}){
		@{$hash_window{$chr}{$windowrank}{groups}} = @null_groups;
	}
	push @{$hash_window{$chr}{$windowrank}{tags}}, $tag;
	for(my $i = 0; $i < @groups; $i++){
		${$hash_window{$chr}{$windowrank}{groups}}[$i] += $groups[$i];
	}
}
close IN;

open(OUT,">$outfile");
print OUT "#CHROM\tRANK\tLOCINUM\t".join("\t",@groupnames)."\n";

foreach my $chr(@chrs){
	my @ranks = sort {$a <=> $b} keys %{$hash_window{$chr}};
	for(my $rank = 0; $rank <=$ranks[-1]; $rank++){
		my $lociNum = 0;
		unless(exists $hash_window{$chr}{$rank}){
			@{$hash_window{$chr}{$rank}{groups}} = @null_groups;
		}else{
			$lociNum = @{$hash_window{$chr}{$rank}{tags}};
		}
		my @outs = ();
		my $sum = 0;
		foreach(@{$hash_window{$chr}{$rank}{groups}}){
			$sum += $_;
		}
		if($sum > 0){
			foreach(@{$hash_window{$chr}{$rank}{groups}}){
				my $newratio = $_/$sum;
				$newratio = sprintf("%.3f",$newratio);
				push @outs, $newratio;
			}
		}else{
			@outs = @null_groups;
		}
		print OUT "$chr\t$rank\t$lociNum\t".join("\t",@outs)."\n";
	}
}
close OUT;
