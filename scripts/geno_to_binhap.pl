use strict;
use warnings;
use Getopt::Long;

my($in,$out,$windowSize);
my $usage = "A specified haplotype caller for homozygous bi-allelic SNP data.\n";
$usage .= "USAGE:\nperl $0 --in <input geno file> --out <output path>\n";
$usage .= "--window <window size>. Default=10000\n";

#defaults
$windowSize = 10 * 1000;
my $nameLen = 5; # Chr12_00001 Chr1_02347


GetOptions(
	"in=s" => \$in,
	"out=s" => \$out,
	"window=s" => \$windowSize,
) or die $usage;

die $usage unless(defined $in and defined $out);

open(IN,"<$in") or die $!;

# data for each bin
my $chrTmp = "";
my $binTmp = "";
my %hash_sample = ();
my %hash_bin = ();
# data for the whole file
my @allSamples = ();

while(<IN>){
	chomp;
	my($chr,$pos,$ref,$alt,@datas) = split/\t/;
	if($_ =~ /^#/){
		@allSamples = @datas;
		next;
	}
	my $binRank = int(($pos-1)/$windowSize);
	
	# output data when entering a new bin
	if($chr ne $chrTmp or $binRank ne $binTmp){
		goto RESET if($binTmp eq "");
		
		# make haplotype and output datas
		foreach my $sample(sort keys %hash_sample){
			my $haplotype = join("|",@{$hash_sample{$sample}{genos}});
			push @{$hash_bin{$haplotype}{samples}}, $sample;
			$hash_bin{$haplotype}{count}++;
		}
		#rank the haplotypes
		my @sorted_hap = sort {$hash_bin{$b}{count} <=> $hash_bin{$a}{count}} keys %hash_bin;
		for(my $i = 0; $i < @sorted_hap; $i++){
			my $rank = $i + 1;
			my $haplotype = $sorted_hap[$i];
			if($haplotype !~ /[1-9]/ and $haplotype !~ /\-/){
				$rank = 0;
			}
			$hash_bin{$haplotype}{rank} = $rank;
		}
		#output datas
		my $binName = (sprintf "%0$nameLen"."d", $binTmp);
		$binName = "$chrTmp\_$binName";
		unless(-e "$out/$chrTmp/"){
			system("mkdir -p $out/$chrTmp/");
		}
		open(OUT,">$out/$chrTmp/$binName.haplotype");
		@sorted_hap = sort {$hash_bin{$a}{rank} <=> $hash_bin{$b}{rank}} keys %hash_bin;
		foreach my $haplotype(@sorted_hap){
			my $rank = $hash_bin{$haplotype}{rank};
			my $samples_join = join(",",@{$hash_bin{$haplotype}{samples}});
			my $count = $hash_bin{$haplotype}{count};
			print OUT "$binName\t$rank\t$haplotype\t-\t$count\t$samples_join\n";
		}
		close OUT;
		
		# reset data
		RESET:
		$chrTmp = $chr;
		$binTmp = $binRank;
		%hash_sample = ();
		%hash_bin = ();
	}
	
	# data input
	for(my $i = 0; $i < @datas; $i++){
		my $sample = $allSamples[$i];
		push @{$hash_sample{$sample}{genos}}, $datas[$i];
	}
}

# output the last bin
if($binTmp ne ""){
	# make haplotype and output datas
	foreach my $sample(sort keys %hash_sample){
		my $haplotype = join("|",@{$hash_sample{$sample}{genos}});
		push @{$hash_bin{$haplotype}{samples}}, $sample;
		$hash_bin{$haplotype}{count}++;
	}
	#rank the haplotypes
	my @sorted_hap = sort {$hash_bin{$b}{count} <=> $hash_bin{$a}{count}} keys %hash_bin;
	for(my $i = 0; $i < @sorted_hap; $i++){
		my $rank = $i + 1;
		my $haplotype = $sorted_hap[$i];
		if($haplotype !~ /[1-9]/ and $haplotype !~ /\-/){
			$rank = 0;
		}
		$hash_bin{$haplotype}{rank} = $rank;
	}
	#output datas
	my $binName = (sprintf "%0$nameLen"."d", $binTmp);
	$binName = "$chrTmp\_$binName";
	unless(-e "$out/$chrTmp/"){
		system("mkdir -p $out/$chrTmp/");
	}
	open(OUT,">$out/$chrTmp/$binName.haplotype");
	@sorted_hap = sort {$hash_bin{$a}{rank} <=> $hash_bin{$b}{rank}} keys %hash_bin;
	foreach my $haplotype(@sorted_hap){
		my $rank = $hash_bin{$haplotype}{rank};
		my $samples_join = join(",",@{$hash_bin{$haplotype}{samples}});
		my $count = $hash_bin{$haplotype}{count};
		print OUT "$binName\t$rank\t$haplotype\t-\t$count\t$samples_join\n";
	}
	close OUT;
}

close IN;
