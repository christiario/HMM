use strict;
use List::Util 'max';
use Math::Counting ':student';
use Statistics::ROC;

my $MAX_DIFF = 30;
my $NEG_INF = -2e20;

sub get_emit(){
	my ($dir_gblock,$dir_emit,$chr,$gr,$ave_mcg,$sd_mcg,$ave_mch,$sd_mch) = @_;
	my (@mCGlv,@CG,@mCH,@CH,@refCGn,@refCHn); my $i=0;

	open EMIT,">$dir_emit/$chr.$gr.txt" or die "ERROR!:cant open OUT(EMIT):$dir_emit/$chr.$gr.txt\n$!\n";
	print EMIT "pos\tmCGlv\trefCGnum\tmCHlv\trefCHnum\te_pcor\te_ncor\te_indi\tmax_stat\n";

	open GBLOCK,"$dir_gblock/$chr.$gr.txt" or die "ERROR!:cant open IN(GBLOCK):$dir_gblock/$chr.$gr.txt\n$!\n";
	while(<GBLOCK>){
		if(m/^(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)/){
			my ($pos,$mCG,$CG,$mCH,$CH,$refCGn,$refCHn) = ($1,$2,$3,$4,$5,$6,$7);
			my $mCGlv=$mCG/$CG;
			my $x= $sd_mch*($mCGlv-$ave_mcg)/$sd_mcg;
			my ($p_mu, $n_mu) = ($ave_mch+$x, $ave_mch-$x);

			if($p_mu<=0){$p_mu = 0.00000000000001;}
			elsif($p_mu>1){$p_mu = 1;}
			if($n_mu<=0){$n_mu = 0.00000000000001;}
			elsif($n_mu>1){$n_mu = 1;}
		
			my $refn = $CG;	
			my ($a_p,$b_p) = ($refn*$p_mu, $refn*(1-$p_mu));
			my ($a_n,$b_n) = ($refn*$n_mu, $refn*(1-$n_mu));

			my $p_pcor = ($ave_mch*$CH + $a_p)/($CH + $a_p + $b_p);
			my $p_ncor = ($ave_mch*$CH + $a_n)/($CH + $a_n + $b_n);
			
			if($p_pcor<=0){$p_pcor = 0.00000000000001;}
			elsif($p_pcor>1){$p_pcor = 1;}
			if($p_ncor<=0){$p_ncor = 0.00000000000001;}
			elsif($p_ncor>1){$p_ncor = 1;}

			my $comb_num = loggamma($CH+1) - loggamma($mCH+1) - loggamma($CH-$mCH+1);	
			my $e_indi = ($comb_num) + $mCH*log($ave_mch) + ($CH-$mCH)*log(1-$ave_mch);
			my $e_pcor = ($comb_num) + $mCH*log($p_pcor) + ($CH-$mCH)*log(1-$p_pcor);
			my $e_ncor = ($comb_num) + $mCH*log($p_ncor) + ($CH-$mCH)*log(1-$p_ncor);

			$e_indi = ($e_indi > $NEG_INF)? $e_indi : $NEG_INF;
			$e_pcor = ($e_pcor > $NEG_INF)? $e_pcor : $NEG_INF;
			$e_ncor = ($e_ncor > $NEG_INF)? $e_ncor : $NEG_INF;
			my $max_stat=2;
			if(($e_indi>=$e_pcor) & ($e_indi>=$e_ncor)){ $max_stat = 2;}
			elsif(($e_pcor>$e_indi) & ($e_pcor>$e_ncor)){ $max_stat = 0;} 
			elsif(($e_ncor>$e_indi) & ($e_ncor>$e_pcor)){ $max_stat = 1;} 

			my $mCHlv = sprintf "%.3f", ($mCH/$CH);
			print EMIT "$pos\t$mCGlv\t$refCGn\t$mCHlv\t$refCHn\t$e_pcor\t$e_ncor\t$e_indi\t$max_stat\n";
		}
	}close(GBLOCK);
	close(EMIT);
	return 1;
}

sub get_em(){
	my ($dir_emit,$dir_TRN,$dir_initProb,$chr,$gr) = @_;

	#========= parameters
	my $statN=3;
	my $iterateN=1;
	my $tolerance=0.0001;
	my (@FW,@BW);

	#==== get emission probabiliy
	my @emit; my $i=0;
	open EMIT,"$dir_emit/$chr.$gr.txt" or die "ERROR!:cant open IN:$dir_emit/$chr.$gr.txt\n$!\n";
	while(<EMIT>){
		if(m/^\d+\s\S+\s\d+\s\S+\s\d+\s(\S+)\s(\S+)\s(\S+)/){
			$emit[$i][0]=$1; $emit[$i][1]=$2; $emit[$i][2]=$3; $i++;
		}
	}close(EMIT);
	my $emit_len=$i;	
	#==== set inition probability
	my @initProb = (log(0.4),log(0.3),log(0.3));
	
	#==== set transition probability
	my @TRN;
	for(my $i=0; $i<$statN; $i++){
		for(my $j=0; $j<$statN; $j++){
			$TRN[$i][$j] = ($i==$j)? log(0.4) : log(0.3);
		}
	}

	&printLOG("EM: group:$gr| iter_0",$statN,\@TRN,\@initProb);
	#==== run em
	my $it=1;
	while(){
		my $converge=1;
		#=== foward
		for(my $j=0; $j<$statN; $j++){ 	$FW[0][$j]= $initProb[$j] + $emit[0][$j]; }
		for(my $i=1; $i<$emit_len; $i++){
			for(my $l=0; $l<$statN; $l++){
				my $val = $NEG_INF;
				for(my $k=0; $k<$statN; $k++){
					my $x = $FW[$i-1][$k] + $TRN[$k][$l] + $emit[$i][$l];
					if($x>$NEG_INF){
						if($k==0){$val = $x
						}else{$val = &logSum($val,$x);}
					}
				}
				$FW[$i][$l] = $val
			}
		}
		
		#=== backward
		for(my $k=0; $k<$statN; $k++){ $BW[$emit_len-1][$k]=log(1.0); }
		for(my $i=1; $i<$emit_len; $i++){
			for(my $j=0; $j<$statN; $j++){
				my $val = $NEG_INF;
				for(my $k=0; $k<$statN; $k++){
					my $x = $TRN[$j][$k] + $BW[$emit_len-$i][$k] + $emit[$emit_len-$i][$k];
					if($x >$NEG_INF){
						if($k==0){$val = $x
						}else{$val = &logSum($val,$x);}
					}
				}
				$BW[$emit_len-$i-1][$j] = $val;
			}
		}
		
		#=== total probability
		my ($t_f,$t_b) = ($NEG_INF,$NEG_INF);
		for(my $l=0; $l<$statN; $l++){
			my $x = $FW[$emit_len-1][$l];
			if($x>$NEG_INF){
				if($l==0){$t_f = $x
				}else{$t_f = &logSum($t_f,$x);}
			}
		}
		for(my $l=0; $l<$statN; $l++){
			my $x = $initProb[$l] + $emit[0][$l] + $BW[0][$l];
			if($x>$NEG_INF){
				if($l==0){$t_b = $x
				}else{$t_b = &logSum($t_b,$x);}
			}
		}
		my $total = ($t_f+$t_b)/2;

		#=== get expectation values
		my @ex_init;
		for(my $j=0; $j<$statN; $j++){
			my $x = $initProb[$j] + $emit[0][$j] + $BW[0][$j] - $total;
			$ex_init[$j] = ($x>$NEG_INF)? exp($x) : 0;
		}

		my @ex_trans;
		for(my $j=0; $j<$statN; $j++){
			for(my $k=0; $k<$statN; $k++){
				my $val=$NEG_INF;
				for(my $i=1; $i<$emit_len; $i++){
					my $x = $FW[$i-1][$j] + $TRN[$j][$k] + $emit[$i][$k] + $BW[$i][$k] - $total;
					if($x>$NEG_INF){
						if($i==1){$val = $x
						}else{$val = &logSum($val,$x);}
					}
				}
				$ex_trans[$j][$k] = ($val>$NEG_INF)? exp($val): 0;
			}
		}
		#=== parameter update
		my $sum=0;
		for(my $j=0; $j<$statN; $j++){ $sum += $ex_init[$j];}
		for(my $j=0; $j<$statN; $j++){
			if($ex_init[$j] > 0){
				my $val = $ex_init[$j]/$sum;
				if(abs(exp($initProb[$j]) - $val) > $tolerance){
					$converge = 0;
					$initProb[$j] = log($val);
				}
			}else{
				$initProb[$j] = $NEG_INF; 
			}
		}
		for(my $j=0; $j<$statN; $j++){
			my $sum_trans = 0;
			for(my $k=0; $k<$statN; $k++){$sum_trans += $ex_trans[$j][$k];}
			for(my $k=0; $k<$statN; $k++){
				if($ex_trans[$j][$k]>0){
					my $val = $ex_trans[$j][$k]/$sum_trans;
					if(abs(exp($TRN[$j][$k]) - $val) > $tolerance){ 
						$converge = 0;
						$TRN[$j][$k] = ($val == 0)? $NEG_INF : log($val);
					}
				}else{
					if(exp($TRN[$j][$k]) > $tolerance){ 
						$converge = 0;
						$TRN[$j][$k] = $NEG_INF; print "TRN unstable\n";
					}
				}
			}
		}
			&printLOG("EM: group:$gr| iter_$it",$statN,\@TRN,\@initProb);
		$it++;
		if($converge == 1){
			"conversion success! \n"; last;
		}elsif($it>10000){die "ERROR!:conversion fail!: $chr $gr $it\n";}
		
	}
	open TRN,">$dir_TRN/$chr.$gr.txt" or die "ERROR!:cant open TRN OUT:$dir_TRN/$chr.$gr.txt\n$!\n";
	for(my $j=0; $j<$statN; $j++){
			print TRN "$TRN[$j][0]\t$TRN[$j][1]\t$TRN[$j][2]\n";
	}close(TRN);

	open INIT_PROB,">$dir_initProb/$chr.$gr.txt" or die "ERROR!:cant open INIT_PROB:$dir_initProb/$chr.$gr.txt\n:$!\n";
	print INIT_PROB "$initProb[0]\t$initProb[1]\t$initProb[2]\n";
	close(INIT_PROB);

	print "$chr $gr conversion success!\n";
	return 1;
}

sub get_viterbi(){
	my ($dir_emit, $dir_TRN, $dir_initProb, $dir_viterbi, $chr, $gr) = @_;
	my $statN=3;

	#=== get Transmission probability
	my @TRN;
	my $i=0;
	open IN,"$dir_TRN/$chr.$gr.txt" or die "ERROR!:viterbi step:cant open TRN_IN:$!\n";
	while(<IN>){
		if(m/(\S+)\s(\S+)\s(\S+)/){
			$TRN[$i][0] = $1; $TRN[$i][1] = $2; $TRN[$i][2] = $3; $i++;
		}
	}close(IN);

	#=== get inition probability
	my @initProb;
	$i=0;
	open IN,"$dir_initProb/$chr.$gr.txt" or die "ERROR!:viterbi step: cant open initProb_IN:$!\n";
	while(<IN>){
		if(m/(\S+)\s(\S+)\s(\S+)/){ @initProb = ($1,$2,$3); }
	}close(IN);

	#=== get emission probability
	my (@pos,@mCG,@CGn,@mCH,@CHn);
	my @emit;
	$i=0;
	open IN,"$dir_emit/$chr.$gr.txt" or die "ERROR!:viterbi step: cant open EMIT_IN:$!\n";
	while(<IN>){
		if(m/^(\d+)\s(\S+)\s(\d+)\s(\S+)\s(\d+)\s(\S+)\s(\S+)\s(\S+)/){
			$pos[$i]=$1; $mCG[$i]=$2; $CGn[$i]=$3; $mCH[$i]=$4; $CHn[$i]=$5;
			$emit[0][$i]=$6; $emit[1][$i]=$7; $emit[2][$i]=$8; 
			$i++;
		}
	}close(IN);
	my $emit_len = $i;
	print "viterbi: emit length: $emit_len\n";
	&printLOG("viterbi: TRN: ",$statN,\@TRN,\@initProb);

	#=== viterbi
	my (@vtb,@trc_vtb);
	for(my $j=0; $j<$statN; $j++){ 	$vtb[$j][0] = $emit[$j][0] + $initProb[$j]; }
	for(my $i=1; $i<$emit_len; $i++){
		for(my $j=0; $j<$statN; $j++){
			my $max  = $NEG_INF;
			my $idx = 0;
			for(my $k=0; $k<$statN; $k++){
				my $v = $vtb[$k][$i-1] + $TRN[$k][$j] + $emit[$j][$i];
				if($v > $max){$max = $v; $idx=$k;}
			}
			$vtb[$j][$i] = $max;
			$trc_vtb[$j][$i] = $idx;
		}
	}
	
	#=== viterbi trace back
	my @path;
	my $max = $NEG_INF;
	my $idx;
	for(my $j=0; $j<$statN; $j++){
		if($vtb[$j][-1] > $max){
			$max = $vtb[$j][-1];
			$idx = $j;
		}
	}
	$path[$emit_len-1]=$idx;
	for(my $i=1; $i<$emit_len; $i++){
		$path[$emit_len-1-$i] = $trc_vtb[$path[$emit_len-$i]][$emit_len-$i];
	}

	#write results
	open OUT,">$dir_viterbi/$chr.$gr.txt" or die "ERROR!:cant open viterbi out:$dir_viterbi/$chr.$gr.txt\n$!\n";
	for(my $i=0; $i<@pos; $i++){
		print OUT "$pos[$i]\t$mCG[$i]\t$CGn[$i]\t$mCH[$i]\t$CHn[$i]\t$path[$i]\n";
	}close(OUT);
	print "viterbi done: $chr.$gr\n";
	return 1;
}

sub logSum(){
	my ($a,$b) = ($_[0],$_[1]);
	if($a<$b){($a,$b) = ($b,$a);}
	if(($a-$b)<$MAX_DIFF){
		my $out = log(exp($a-$b)+1) + $b;
		if($out == 0){print "exceed float: $a:\t$b:\t$out:\n";}
		elsif($out eq "inf"){die "ERROR!:exceed float:$a:\t$b:\t$out:\n";}
		return $out;
	}else{
		return $a;
	}
}

sub printLOG(){
	my ($pre,$statN,$TRN_pt,$initProb_pt) = @_;
	my @TRN = @$TRN_pt;
	my @initProb = @$initProb_pt;
	my $out = $pre;
	for(my $i=0; $i<$statN; $i++){
		for(my $j=0; $j<$statN; $j++){
			my $val = sprintf "%.4f", exp($TRN[$i][$j]);
			$out = $out."\t$val";
		}
	}
	for(my $i=0; $i<$statN; $i++){
			my $val = sprintf "%.4f", exp($initProb[$i]);
			$out = $out."\t$val";
	}
	print LOG "$out\n";
	print "$out\n";
}

