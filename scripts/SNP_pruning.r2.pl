use strict;
use warnings;
use Getopt::Long;

my($in,$out,$windowSize,$min_r2,$minCover,$minMaf,$keepBed,$keepList,$pass,$pass1,$varList);

# default values
$windowSize = 10 * 1000;
$min_r2 = 0.9;
$minCover = 0.9;
$minMaf = 0.2;
# default values

my $usage = "USAGE:\nperl $0 --in <input VCF> --out <output file> \n";
$usage .= "--window <window size>. Default=10000 \n";
$usage .= "--r2 <minor geno-r2>: minor geno-r2 between two SNPs to be merged. Default=0.9 \n";
$usage .= "--cov <minor coverage>: minor coverage of a SNP to be kept. Default=0.9 \n";
$usage .= "--maf <minor maf>: minor maf of a SNP to be kept. Default=0.2 \n";
$usage .= "--bed <selected region BED>. Selected region to be analyzed in BED format. [Optional] \n";
$usage .= "--keep <keep sample list>. Selected sample list to be analyzed. [Optional] \n";
$usage .= "--pass to keep SNPs with filter TAG [PASS] or [SnpCluster].\n";
$usage .= "--pass1 to keep SNPs with filter TAG [PASS].\n";
$usage .= "--varlist <var list>\n\n";

GetOptions(
	"in=s" => \$in,
	"out=s" => \$out,
	"window=s" => \$windowSize,
	"r2=s" => \$min_r2,
	"cov=s" => \$minCover,
	"maf=s" => \$minMaf,
	"bed=s" => \$keepBed,
	"keep=s" => \$keepList,
	"pass!" => \$pass,
	"pass1!" => \$pass1,
	"varlist=s" => \$varList,
) or die $usage;

die $usage unless(defined $in and defined $out);

my %hash_keep;
if(defined $keepList){
	open(IN,"<$keepList") or die $!;
	while(<IN>){
		chomp;
		my($sample,$other) = split/\t/;
		$hash_keep{$sample}{exist} = 0;
	}
	close IN;
}

my %hash_bed;
my %hash_chr;
if(defined $keepBed){
	open(IN,"<$keepBed") or die $!;
	while(<IN>){
		chomp;
		my($chr,$start,$end,$other) = split/\t/;
		next if(exists $hash_bed{$chr}{$start} and $hash_bed{$chr}{$start} > $end);
		$hash_bed{$chr}{$start} = $end;
	}
	close IN;
	# merge regions
	foreach my $chr(keys %hash_bed){
		$hash_chr{$chr} = "";
		my @starts = sort {$a <=> $b} keys %{$hash_bed{$chr}};
		for(my $i = 0; $i < @starts; $i++){
			my $starti = $starts[$i];
			my $endi = $hash_bed{$chr}{$starti};
			for(my $j = $i+1; $j < @starts; $j++){
				my $startj = $starts[$j];
				my $endj = $hash_bed{$chr}{$startj};
				last if($startj > $endi);
				if($endj > $endi){
					$endi = $endj;
					$hash_bed{$chr}{$starti} = $endj;
				}
				delete($hash_bed{$chr}{$startj});
				splice(@starts,$j,1);
				$j--;
			}
		}
	}
}

my %hash_var;
if(defined $varList){
	open(IN,"<$varList") or die $!;
	while(<IN>){
		chomp;
		next if($_ =~ /^#/);
		my($chr,$pos,$ref,$alt) = split/\t/;
		$hash_var{$chr}{$pos}{$ref}{$alt} = 0;
	}
	close IN;
}

open(IN,"<$in") or die $!;
open(OUT,">$out");
# fundamental data
my @pickedRanks = ();
my %hash_bin;
my $binTmp = -1;
my $chrTmp = "NA";
my $snpRank = -1;
# data for region filtering
my @remained_chrs = ();
my $chr_tmp = "NA";
my @regions = ();

while(<IN>){
	chomp;
	next if($_ =~ /^##/);
	my($chr,$pos,$id,$ref,$alts_join,$qual,$filter,$info,$format,$datas_join) = split/\t/,$_,10;
	
	# dissect samples
	if($chr =~ /#CHROM/){
		my @samples = split/\t/,$datas_join;
		my $tail = @samples - 1;
		my @pickedSamples = ();
		
		print @pickedRanks."\n";
		unless(defined $keepList){
			@pickedSamples = @samples;
			@pickedRanks = (0..$tail);
			goto NOKEEP;
		}
		for(my $i = 0; $i < @samples; $i++){
			next unless(exists $hash_keep{$samples[$i]});
			push @pickedRanks, $i;
			push @pickedSamples, $samples[$i];
		}
		NOKEEP:
		print "#Keeping ".@pickedRanks." samples from ".@samples." samples.\n";
		print OUT "#CHROM\tPOS\tREF\tALT\t".join("\t",@pickedSamples)."\n";
		next;
	}
	
	if(@pickedRanks == 0){
		die "#ERROR: No sample picked in the vcf file, please check your inputs!\n";
	}
	
	my $bin = int(($pos-1)/$windowSize);
	
	# output the previous bin when entering a new bin
	if($chr ne $chrTmp or $bin ne $binTmp){
		my @snpRanks = sort {$a <=> $b} keys %hash_bin;
		my $count = @snpRanks;
		for(my $i = 0; $i < @snpRanks; $i++){
			for(my $j = $i+1; $j < @snpRanks; $j++){
				my $snpRanki = $snpRanks[$i];
				my $snpRankj = $snpRanks[$j];
				my @genosi = @{$hash_bin{$snpRanki}{genos}};
				my @genosj = @{$hash_bin{$snpRankj}{genos}};
				#my($sum,$match) = compareSNP(\@genosi,\@genosj);
				my($indv,$r2) = geno_r2(\@genosi,\@genosj);
				next if($r2 < $min_r2);
				my $covi = $hash_bin{$snpRanki}{cov};
				my $covj = $hash_bin{$snpRankj}{cov};
				if($covi > $covj){
					splice(@snpRanks,$j,1);
					delete($hash_bin{$snpRankj});
					$j--;
				}else{
					splice(@snpRanks,$i,1);
					$i--;
					delete($hash_bin{$snpRanki});
					last;
				}
			}
		}
		my $remained = @snpRanks;
		print "# CHROM:$chrTmp BIN:$binTmp INPUT:$count OUTPUT:$remained\n";
		for(my $i = 0; $i < @snpRanks; $i++){
			my $snpRanki = $snpRanks[$i];
			my $posi = $hash_bin{$snpRanki}{pos};
			my $refi = $hash_bin{$snpRanki}{ref};
			my $alti = $hash_bin{$snpRanki}{alt};
			print OUT "$chrTmp\t$posi\t$refi\t$alti\t".join("\t",@{$hash_bin{$snpRanki}{genos}})."\n";
		}
		$snpRank = -1;
		$binTmp = $bin;
		$chrTmp = $chr;
		%hash_bin = ();
	}
	
	# Check if the position locates in the candidate regions
	goto NOBED unless(defined $keepBed);
	
	next unless(exists $hash_bed{$chr});
	if($chr ne $chr_tmp){
		@regions = ();
		foreach my $start(sort {$a <=> $b} keys %{$hash_bed{$chr}}){
			my $end = $hash_bed{$chr}{$start};	
			push @regions, "$start,$end";
		}
		$chr_tmp = $chr;
		if(exists $hash_chr{$chr}){
			delete($hash_chr{$chr});
		}
		@remained_chrs = keys %hash_chr;
		print "\t# Reading chr:$chr\n";
		if(@regions == 0){
			print "\t# Skipping chr:$chr\n";
		}
	}
	
	if(@remained_chrs == 0 and @regions == 0){
		last;
	}

	if(@regions == 0){
		next;
	}
	
	my $pickIt = 0;
	for(my $i = 0; $i < @regions; $i++){
		my($start,$end) = split/,/,$regions[$i];
		if($pos < $start){
			last;
		}elsif($pos > $end){
			splice(@regions, $i, 1);
			$i--;
			next;
		}else{
			$pickIt = 1;
			last;
		}
	}

	next unless($pickIt == 1);
	NOBED:
	
	# Check if the variant is in the varList
	if(defined $varList){
		next unless(defined $hash_var{$chr}{$pos}{$ref}{$alts_join});
	}
	
	# Check if the variant is a bi-allelic SNP
	next unless(length($ref) == 1 and length($alts_join) == 1);
	
	# Check if the SNP
	if(defined $pass){
		next unless($filter eq "PASS" or $filter eq "SnpCluster");
	}

	if(defined $pass1){
		next unless($filter eq "PASS");
	}
	
	# Check maf and coverage
	my @datas = split/\t/, $datas_join;

	my @pickedDatas = ();
	foreach my $rank(@pickedRanks){
		push @pickedDatas, $datas[$rank];
	}
	my($cov_,$maf_) = get_cov_maf(@pickedDatas);
	next unless($cov_ >= $minCover and $maf_ >= $minMaf);
	
	# push the geno data of the SNP to the hash_bed
	$snpRank++;
	$hash_bin{$snpRank}{cov} = $cov_;
	$hash_bin{$snpRank}{maf} = $maf_;
	$hash_bin{$snpRank}{pos} = $pos;
	$hash_bin{$snpRank}{ref} = $ref;
	$hash_bin{$snpRank}{alt} = $alts_join;
	foreach my $rank(@pickedRanks){
		my($geno) = vcf2geno($datas[$rank]);
		push @{$hash_bin{$snpRank}{genos}}, $geno;
	}
}

if($snpRank >= 0){
	my @snpRanks = sort {$a <=> $b} keys %hash_bin;
	my $count = @snpRanks;
	for(my $i = 0; $i < @snpRanks; $i++){
		for(my $j = $i+1; $j < @snpRanks; $j++){
			my $snpRanki = $snpRanks[$i];
			my $snpRankj = $snpRanks[$j];
			my @genosi = @{$hash_bin{$snpRanki}{genos}};
			my @genosj = @{$hash_bin{$snpRankj}{genos}};
			my($indv,$r2) = geno_r2(\@genosi,\@genosj);
			next if($r2 < $min_r2);
			my $covi = $hash_bin{$snpRanki}{cov};
			my $covj = $hash_bin{$snpRankj}{cov};
			if($covi > $covj){
				splice(@snpRanks,$j,1);
				delete($hash_bin{$snpRankj});
				$j--;
			}else{
				splice(@snpRanks,$i,1);
				$i--;
				delete($hash_bin{$snpRanki});
				last;
			}
		}
	}
	my $remained = @snpRanks;
	print "# CHROM:$chrTmp BIN:$binTmp INPUT:$count OUTPUT:$remained\n";
	for(my $i = 0; $i < @snpRanks; $i++){
		my $snpRanki = $snpRanks[$i];
		my $posi = $hash_bin{$snpRanki}{pos};
		my $refi = $hash_bin{$snpRanki}{ref};
		my $alti = $hash_bin{$snpRanki}{alt};
		print OUT "$chrTmp\t$posi\t$refi\t$alti\t".join("\t",@{$hash_bin{$snpRanki}{genos}})."\n";
	}
}

close IN;

sub geno_r2{
	my($arr1,$arr2) = @_;
	unless(@{$arr1} == @{$arr2}){
		die "#ERROR: the two arrys must be the same length";
	}
	my $sum = 0;
	my $r2 = 0;
	my $pa = 0;
	my $pb = 0;
	my $pab = 0;
	my $ga;
	my $gb;
	for(my $i = 0; $i < @{$arr1}; $i++){
		my $va1 = ${$arr1}[$i];
		my $va2 = ${$arr2}[$i];
		next if($va1 eq "NA" || $va1 eq "n/a" || $va1 eq "-" || $va1 eq "null");
		next if($va2 eq "NA" || $va2 eq "n/a" || $va2 eq "-" || $va2 eq "null");
		$sum++;
		unless(defined $ga){
			$ga = $va1;
			$gb = $va2;
		}
		if($va1 eq $ga){
			$pa++;
		}
		if($va2 eq $gb){
			$pb++;
		}
		if($va1 eq $ga and $va2 eq $gb){
			$pab++;
		}
	}
	if($sum > 0 and $pa < $sum and $pb < $sum){
		$pa = $pa/$sum;
		$pb = $pb/$sum;
		$pab = $pab/$sum;
		$r2 = (($pab - $pa * $pb) ** 2) / ($pa * (1-$pa) * $pb * (1-$pb));
	}
	return($sum,$r2);
}

sub get_cov_maf{
	my @genos = @_;
	my $altcount = 0;
	my $covcount = 0;
	my $lostcount = 0;
	
	foreach my $geno(@genos){
		my $gt = "-";
		if($geno =~ /^(\d+)\/(\d+)/ and $1 eq $2){
			$gt = $1;
			if($gt != 0){
				$altcount ++;
			}
		}
		if($gt ne '-'){
			$covcount ++;
		}else{
			$lostcount ++;
		}
	}
	
	my $cov = 1 - $lostcount/@genos;
	my $maf = $altcount/@genos;
	my $maf2 = $cov - $maf;
	if($maf > $maf2){
		$maf = $maf2;
	}
	return($cov,$maf);
}

sub compareSNP{
	my($in1,$in2) = @_;
	my @arr1 = @$in1;
	my @arr2 = @$in2;
	if(@arr1 != @arr2){
		die "#ERROR: wrong array length.\n";
	}
	my $sum = 0;
	my $value1 = 0;
	my $value2 = 0;
	for(my $i = 0; $i < @arr1; $i++){
		next if($arr1[$i] eq "-" or $arr2[$i] eq "-");
		$sum++;
		if($arr1[$i] eq $arr2[$i]){
			$value1++;
		}else{
			$value2++;
		}
	}
	if($value1 < $value2){
		$value1 = $value2;
	}
	return($sum,$value1);
}

sub vcf2geno{
	my($data) = @_;
	my $geno = "-";
	if($data =~ /(\d+)\/(\d+)/){
		if($1 == $2 and $1 == 0){
			$geno = 0;
		}elsif($1 == $2 and $1 != 0){
			$geno = 1;
		}
	}
	return($geno);
}
