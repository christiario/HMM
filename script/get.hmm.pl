use strict;
use Statistics::Basic qw(:all);
use File::Basename;

my $dir_work = `pwd`; chomp $dir_work; 
$dir_work =~ s/\/[^\/]+$//;
require "$dir_work/script/get.hmm.sub_functions.pl";
use Array::Utils qw(:all);

#### Parameter Setting
my $minCnum  = 10;
my $species  = "hg19";
my $mCH_dist = 1000000;
my $bsize = 180;

#### Directory Setting
if(@ARGV != 2){die "ERROR!:perl $0 [input_file (full path)] [output_dir]\n";}
my $infile = $ARGV[0];
my $sample = basename($infile);
my $outdir = $ARGV[1]."/";
my $dir_gblock = "$outdir/gblock";
my $dir_emit = "$outdir/emit";
my $dir_TRN = "$outdir/TRN";
my $dir_initProb = "$outdir/initProb";
my $dir_viterbi = "$outdir/viterbi";
my $dir_statistics = "$outdir/statistics";
system("mkdir -p $dir_gblock")unless(-s $dir_gblock);
system("mkdir -p $dir_emit")unless(-s $dir_emit);
system("mkdir -p $dir_TRN")unless(-s $dir_TRN);
system("mkdir -p $dir_initProb")unless(-s $dir_initProb);
system("mkdir -p $dir_viterbi")unless(-s $dir_viterbi);
system("mkdir -p $dir_statistics")unless(-s $dir_statistics);
print "Input file:\n$infile\nOutput directory created under:\n$outdir\nReference data loaded from:\n$dir_work\nSample Name:$sample\n";

### Check target chromosomes
my %vChr;
open DAT, "gunzip -c $infile |" or die "cant open Input file:$!\n";
<DAT>;
while(<DAT>){
	if(m/^(\S+)\s/){ $vChr{$1}=1; }
}close(DAT);
my @chr = keys %vChr;

print "Target chromosomes: @chr \n";

foreach my $chr (@chr){
	#### Average mCG/mCH, number of CG/CH in each block
	my (%refCGn, %refCHn);
	my $ref_file = "$dir_work/Data/hg19_b180/chr$chr.txt.gz";
	open REF,"gunzip -c $ref_file |" or die "ERROR!:cant open REF:$!\nREF FILE:$ref_file\n";
	while(<REF>){
		if(m/chr$chr\_(\S+)\s(\d+)\s(\d+)/){	
	 		if(($2>0) | ($3>0)){ $refCGn{$1}=$2; $refCHn{$1}=$3; }
		}
	}close(REF);
	my $tsize = keys %refCGn;
	print "hg19 chr$chr Cytoine # read\n";
	
	#### Load mCG/mCH read count at each block 
	my (%mCG,%CG,%mCH,%CH);
	open DAT, "gunzip -c $infile |" or die "cant open Input file:$!\n";
	while(<DAT>){
		if(m/$chr\s(\d+)\s(\S+)\s(\d+)\s(\d+)/){
			my $pos = int($1/$bsize)*$bsize;
			if($2 eq "CG"){	$mCG{$pos} += $3; $CG{$pos} +=$4;}
			else{	$mCH{$pos} += $3; $CH{$pos} +=$4;}
		}
	}close(DAT);

	#### Average and Median
	my ($null_mcg, $ave_mcg, $sd_mcg) = &get_statistic(\%mCG,\%CG);
	my ($ave_mch, $null_mch, $sd_mch) = &get_statistic(\%mCH,\%CH);
	print "chr$chr statics: mCG Median: $ave_mcg\tmCG sd: $sd_mcg\tmCH ave: $ave_mch\tmCH sd: $sd_mch\n";

	##=============== select valid blocks
	my @cg_key = sort{$a<=>$b} keys %CG;
	my @ch_key = sort{$a<=>$b} keys %CH;
	my @p = intersect(@cg_key, @ch_key);
	my (@pos, @mCG, @CG, @mCH, @CH, @refCGn, @refCHn);
	my $c=0;
	foreach my $pos(@p){
		if(($CG{$pos}>=$minCnum) & ($CH{$pos}>=$minCnum)){
			($pos[$c], $mCG[$c], $CG[$c], $mCH[$c], $CH[$c], $refCGn[$c], $refCHn[$c]) = ($pos, $mCG{$pos}, $CG{$pos}, $mCH{$pos}, $CH{$pos}, $refCGn{$pos}, $refCHn{$pos});
			$c++;
		}
	}
	print "chr$chr validBlock#: $c ===\n";

	##============= dividing chromosome based on distance between mCs (>1M)
	my @group;
	my $g=0; $group[0]=0;
	my $bf = $pos[0];
	for(my $i=1; $i<@pos; $i++){
		if(($pos[$i]-$bf) > $mCH_dist){$g++;}
		$group[$i] = $g;
		$bf = $pos[$i];
	}
	my $g_num = $g+1; 
	print "chr$chr HMM Block#: $g_num ===\n";
	
	##============ save valid and group-divided blocks
	for(my $gr=0; $gr<$g_num; $gr++){
		my $incount = 0;
		open GBLOCK,">$dir_gblock/chr$chr.$gr.txt" or die "ERROR!:cant open OUT $dir_gblock/chr$chr.$gr.txt\n$!\n";
		print GBLOCK "pos\tmCGn\tCGn\tmCHn\tCHn\trefCGn\trefCHn\n";
		for(my $i=0; $i<@pos; $i++){
			if($group[$i] == $gr){
				print GBLOCK "$pos[$i]\t$mCG[$i]\t$CG[$i]\t$mCH[$i]\t$CH[$i]\t$refCGn[$i]\t$refCHn[$i]\n";
			}
		}
		close(GBLOCK);
		
		##============= get emission probability
		&get_emit($dir_gblock,$dir_emit,"chr$chr",$gr,$ave_mcg,$sd_mcg,$ave_mch,$sd_mch);
	
		##============= transision probability learning
		if(-e "$dir_emit/chr$chr.$gr.txt"){
			&get_em($dir_emit,$dir_TRN,$dir_initProb,"chr$chr",$gr);
		}else{ die "ERROR!:EMISSION FILE MISSING: $dir_emit/chr$chr.$gr.txt"; }
	
		##============= viterbi
		if((-e "$dir_TRN/chr$chr.$gr.txt") || (-e "$dir_initProb/chr$chr.$gr.txt")){
			&get_viterbi($dir_emit, $dir_TRN, $dir_initProb, $dir_viterbi, "chr$chr", $gr);
		}else{ 	die "ERROR!:TRN or initProb MISSING: \n$dir_TRN/chr$chr.$gr.txt\n$dir_initProb/chr$chr.$gr.txt\n"; }
	
		##=========== write statistics
		my (%emit,%viterbi);
		($emit{0},$emit{1},$emit{2})=(0,0,0);	
		($viterbi{0},$viterbi{1},$viterbi{2})=(0,0,0);	
		my ($e_t,$v_t) = (0,0);
		open IN,"$dir_emit/chr$chr.$gr.txt" or die "ERROR!:cant open IN for simple test:$!\n";
		while(<IN>){
			if(m/(\d)$/){$emit{$1}++;$e_t++;}
		}close(IN);
		open IN,"$dir_viterbi/chr$chr.$gr.txt" or die "ERROR!:cant open IN for simple test:$!\n";
		while(<IN>){
			if(m/(\d)$/){$viterbi{$1}++;$v_t++;}
		}close(IN);

		open OUT,">$dir_statistics/chr$chr.$gr.txt" or die "cant open statistics:$!\n$dir_statistics/chr$chr.$gr.txt\n";
		print OUT "Emission_total\tEmission_P\tEmission_N\tEmission_I\tViterbi_total\tViterbi_P\tViterbi_N\tViterbi_I\n";
		print OUT "$e_t\t$emit{0}\t$emit{1}\t$emit{2}\t$v_t\t$viterbi{0}\t$viterbi{1}\t$viterbi{2}\n";
		close(OUT);
	}
}

sub get_statistic(){
	my %mc    = %{$_[0]};
	my %c     = %{$_[1]};
	my ($mc_sum,$c_sum)=(0,0);
	my @key = keys %c;
	my @ave;
	foreach my $key(@key){
		$mc_sum += $mc{$key};
		$c_sum += $c{$key};
		push @ave, $mc{$key}/$c{$key};
	}
	my $ave = ($c_sum!=0)? ($mc_sum/$c_sum) : "NA";
	my $med = median(@ave);
	my $sd = stdev(\@ave);
	return ($ave, $med, $sd);
}
sub stdev{
	my $data = $_[0];
	if(@$data == 1){ return 0; }
	my $average = &average($data);
	my $sqtotal = 0;
	foreach(@$data) { $sqtotal += ($average-$_) ** 2; }
	my $std = ($sqtotal / (@$data-1)) ** 0.5;
	return $std;
}

sub average{
	my($data) = @_;
	if (not @$data) { die("Empty array"); }
	my $total = 0;
	foreach (@$data) { $total += $_; }
	my $average = $total / @$data;
	return $average;
}

