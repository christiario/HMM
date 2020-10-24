# HMM
HMM for detection of CpG-CpH differentially methylated regions
details are described in https://openlooper.hgc.jp/cphHMM/

perl script/get.hmm.pl [path/input file] [output_dir]

path/input file: full path to the input file (e.x. example/br81y.chr19.50M.txt.gz).
output_dir: full path to output directory

Input: compressed (gzip) file with tap-seperated 4 columns, CHR BP CXX mC_read and totalC_read 
CHR: chromosome number (1..19,X,Y)
BP: position of the cytosine
CXX: cytosine context. CG, CHH, or CHG (both CHH and CHG are recognized as CH)
mC_read: methylated read count mapped at this position
totalC_read: total read count mapped at this position

Example:
CHR BP  CXX mC_read totalC_read
19  60001 CHH 0 10
19  60006 CHG 0 17
19  60120 CG  2 2
...


Output: 
output_dir/gblock/  
  pos: starting position of 180bp bin 
  mCGn: sum of methylated read count aligned at CG in the bin
  CGn: sum of total read count aligned at CG in the bin
  mCHn: sum of methylated read count aligned at CH in the bin
  CHn: sum of total read count aligned at CH in the bin
  refCGn: number of CG in the bin in hg19
  refCHn: number of CH in the bin in hg19
  
output_dir/emit/
  pos: starting position of 180bp bin 
  mCGlv: average methylation level at CG
  refCGn: number of CG in the bin in hg19
  mCHlv: average methylation level at CH
  refCHn: number of CH in the bin in hg19
  e_p: probability that the bin is belong to P-state
  e_n: probability that the bin is belong to N-state
  e_i: probability that the bin is belong to I-state
  max_stat: state of top probability (0:P, 1:N, 2:I)
  
output_dir/initProb/
  randomly set initial probability of P-, N-, and I-state

output_dir/TRN/
  log-scaled transition rates between states. column order: P, N, I, row order: P, N, I
 
output_dir/viterbi/
  pos: starting position of 180bp bin 
  mCGlv: average methylation level at CG
  refCGn: number of CG in the bin in hg19
  mCHlv: average methylation level at CH
  refCHn: number of CH in the bin in hg19
  state: state designated by Viterbi decoding (0:P-, 1:N-, 2:I-state)
  
output_dir/statistics/  Number of bins deteced as P, N, and I-state by top emission probability and Viterbi decoding. 
  
 
  
  

