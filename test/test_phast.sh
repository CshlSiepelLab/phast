******************** phyloP ********************

phyloFit hmrc.ss --tree "(human, (mouse,rat), cow)" -i SS --quiet
@phyloP --null 10 phyloFit.mod
msa_view -i SS -o SS --end 100 hmrc.ss > hmrc_short.ss
@phyloP -i SS phyloFit.mod hmrc_short.ss
@phyloP -i SS --method LRT --base-by-base phyloFit.mod hmrc.ss
@phyloP -i SS --method LRT --mode CONACC --wig-scores phyloFit.mod hmrc.ss
@phyloP -i SS --method SCORE --mode NNEUT --base-by-base phyloFit.mod hmrc.ss
@phyloP -i SS --method GERP --mode ACC --wig-scores phyloFit.mod hmrc.ss
@phyloP -i SS --method GERP --base-by-base phyloFit.mod hmrc.ss
@phyloP -i SS --method SCORE --wig-scores phyloFit.mod hmrc.ss
@phyloP -i SS --method LRT --wig-scores phyloFit.mod hmrc.ss
@phyloP -i SS --method LRT --base-by-base phyloFit.mod hmrc.ss
@phyloP -i SS --method SCORE --wig-scores --refidx 2 phyloFit.mod hmrc.ss
echo -e "chr1\t0\t10\nchr1\t50\t100\nchr1\t200\t300" > temp.bed
@phyloP -i SS --method LRT --mode CONACC --features temp.bed phyloFit.mod hmrc.ss
@phyloP -i SS --method SCORE --features temp.bed -g phyloFit.mod hmrc.ss
tree_doctor --name-ancestors phyloFit.mod > phyloFit-named.mod
@phyloP -i SS --method LRT --mode CONACC --subtree mouse-rat -w phyloFit-named.mod hmrc.ss
@phyloP -i SS --method LRT --mode ACC --branch mouse-rat -w phyloFit-named.mod hmrc.ss
@phyloP --posterior -i SS phyloFit.mod hmrc_short.ss
@phyloP --fit-model -w -i SS --subtree mouse-rat phyloFit-named.mod hmrc_short.ss
@phyloP --epsilon 0.0001 -i SS phyloFit-named.mod hmrc_short.ss
@phyloP --confidence-interval 0.05 -i SS phyloFit-named.mod hmrc_short.ss
@phyloP --quantiles --null 100 phyloFit-named.mod
rm -f hmrc_short.ss phyloFit.mod phyloFit-named.mod temp.bed



******************** phastCons ********************

# these are the two examples in the original test Makefile
!elements.bed @phastCons hpmrc.ss hpmrc-rev-dg-global.mod --nrates 20 --transitions .08,.008  --viterbi elements.bed --seqname chr22
tree_doctor hpmrc-rev-dg-global.mod --prune galGal2 > hpmr.mod
!elements-4way.bed @phastCons hpmrc.ss hpmr.mod --nrates 20 --transitions .08,.008  --viterbi elements-4way.bed --seqname chr22


# let's go through the examples given in phastCons --help
tree_doctor --scale 3.0 hpmr.mod > hpmr_fast.mod
@phastCons  hpmrc.ss hpmr.mod,hpmr_fast.mod
@phastCons  hpmrc.ss --rho  0.5 hpmr.mod
msa_view -i SS --end 1000 -o SS hpmrc.ss > hpmrc_short.ss
!tempTree.cons.mod !tempTree.noncons.mod @phastCons  hpmrc_short.ss hpmr.mod --estimate-trees tempTree
@phastCons --target-coverage 0.25 --expected-length 12 hpmrc.ss hpmr.mod,hpmr_fast.mod
@phastCons --transitions 0.01,0.02 hpmrc.ss hpmr.mod,hpmr_fast.mod
!tempRho.cons.mod !tempRho.noncons.mod @phastCons --target-coverage 0.25 --expected-length 12 --estimate-rho tempRho --no-post-probs hpmrc.ss hpmr.mod
!tempTree.cons.mod !tempTree.noncons.mod @phastCons --target-coverage 0.25 --expected-length 12 --estimate-trees tempTree --no-post-probs hpmrc_short.ss hpmr.mod
@phastCons hpmrc.ss hpmr.mod,hpmr_fast.mod 
!tempTree.cons.mod !tempTree.noncons.mod @phastCons --target-coverage 0.25 --estimate-rho tempTree hpmrc.ss hpmr.mod

# now do at least minimal testing on each option.  (or at least the ones that seem useful)
#--gc,-G
!tempRho.cons.mod !tempRho.noncons.mod @phastCons --gc 0.8 --estimate-rho tempRho hpmrc_short.ss hpmr.mod
!tempTree.cons.mod !tempTree.noncons.mod @phastCons --gc 0.7 --estimate-trees tempTree hpmrc_short.ss hpmr.mod
#--nrates,-K
!tempTree.cons.mod !tempTree.noncons.mod @phastCons -k 12 --estimate-trees tempTree hpmrc_short.ss hpmr.mod
!tempTree.cons.mod !tempTree.noncons.mod @phastCons -k 2,3 --estimate-trees tempTree hpmrc_short.ss hpmr.mod 
#--transitions,-t
@phastCons -t 0.01,0.02 hpmrc_short.ss hpmr.mod
@phastCons -t ~0.02,0.03 hpmrc_short.ss hpmr.mod
#--target-coverage is tested in examples above
#--expected-length is tested in examples above
#--msa-format,-i
msa_view -i SS -o FASTA hpmrc.ss > hpmrc.fa
@phastCons -i FASTA hpmrc.fa hpmr.mod
rm hpmrc.fa
#--viterbi,--most-conserved
!elements.bed @phastCons --most-conserved elements.bed hpmrc.ss hpmr.mod
!elements.gff @phastCons --viterbi elements.gff hpmrc.ss hpmr.mod
#--score and --most-conserved
!elements.bed @phastCons --most-conserved elements.bed --score hpmrc.ss hpmr.mod
#--lnl
!likeFile.txt @phastCons --lnl likeFile.txt hpmrc.ss hpmr.mod
#--no-post-probs
!likeFile.txt @phastCons --lnl likeFile.txt --no-post-probs hpmrc.ss hpmr.mod
!elements.bed @phastCons --most-conserved elements.bed --no-post-probs hpmrc.ss hpmr.mod
#--log.  But don't compare the log files because they include runtime information.
!tempTree.cons.mod !tempTree.noncons.mod  @phastCons --estimate-trees tempTree --log log.txt hpmrc_short.ss hpmr.mod
rm log.txt
#--refidx
@phastCons --refidx 0 hpmrc_short.ss hpmr.mod
@phastCons --refidx 2 hpmrc_short.ss hpmr.mod
#--seqname and --idpref
!elements.gff @phastCons --seqname testSeq --idpref idprefstr --viterbi elements.gff hpmrc_short.ss hpmr.mod

# experimental options
#--indels
@phastCons --indels hpmrc_short.ss hpmr.mod
#--max-micro-indel
@phastCons --max-micro-indel 2 --indels hpmrc_short.ss hpmr.mod
#--indel-params
@phastCons --indels --indel-params 0.02,0.01,0.3,0.01,0.015,0.4 hpmrc_short.ss hpmr.mod
@phastCons --indels --indel-params ~0.02,0.01,0.3,0.01,0.015,0.4 hpmrc_short.ss hpmr.mod
#--indels-only
!likeFile.txt @phastCons --indels-only --lnl likeFile.txt hpmrc_short.ss hpmr.mod

# Felsenstein/Churchill model
#--FC
@phastCons --FC hpmrc_short.ss hpmr.mod
#--lambda
@phastCons --FC --lambda 0.1 hpmrc_short.ss hpmr.mod
@phastCons --FC --lambda ~0.1 hpmrc_short.ss hpmr.mod

#--coding-potential
msa_view -i FASTA -o SS --tuple-size 3 hmrc_correct.fa > hmrc_correct_tuple3.ss
# don't compare stderr because it contains references to the phast home dir
-stderr @phastCons --coding-potential hmrc_correct_tuple3.ss
#@phastCons --coding-potential -i FASTA hmrc_correct.fa
#--extrapolate and alias
tree_doctor --prune mm3 hpmr.mod --rename "hg16 -> human; panTro1 -> chimp; rn3 -> rat" > hpmr_pruned.mod
-stderr @phastCons --alias hg16=human,panTro1=chimp,rn3=rat,mm3=mouse,galGal2=chicken --extrapolate default hpmrc_short.ss hpmr_pruned.mod

#--hmm
tree_doctor --scale 0.1 hpmr.mod > hpmr_slow.mod
@phastCons --hmm ../data/phastCons/simple-coding.hmm hpmrc_short.ss hpmr.mod,hpmr_fast.mod,hpmr_slow.mod,hpmr_fast.mod,hpmr.mod
# TODO: --catmap tests (not sure how to use)
#--states 
@phastCons --hmm ../data/phastCons/simple-coding.hmm --states 2,3,4 hpmrc_short.ss hpmr.mod,hpmr_fast.mod,hpmr_slow.mod,hpmr_fast.mod,hpmr.mod
@phastCons --hmm ../data/phastCons/simple-coding.hmm --states 0 hpmrc_short.ss hpmr.mod,hpmr_fast.mod,hpmr_slow.mod,hpmr_fast.mod,hpmr.mod
#--reflect-strand
@phastCons --hmm ../data/phastCons/simple-coding.hmm --reflect-strand 2,3 hpmrc_short.ss hpmr.mod,hpmr_fast.mod,hpmr_slow.mod,hpmr_fast.mod,hpmr.mod

#missing data (rarely used options)
@phastCons --require-informative 0 --not-informative panTro1,hg16,mm3 hpmrc_short.ss hpmr.mod
@phastCons --ignore-missing hpmrc.ss hpmr.mod

#remove temporary files
rm -f hpmr.mod hpmr_fast.mod hpmrc_short.ss hmrc_correct_tuple3.ss hpmr_pruned.mod hpmr_slow.mod



******************** phyloFit ********************

# These are the same tests implemented in $PHAST/test/Makefile
# Should add more thorough testing at some point

!phyloFit.mod @phyloFit hmrc.ss --subst-mod JC69 --tree "(human, (mouse,rat), cow)" -i SS
!phyloFit.mod @phyloFit hmrc.ss --subst-mod JC69 --tree "((((human,chimp), (mouse,rat)), cow), chicken)" -i SS
!phyloFit.mod @phyloFit hmrc.ss --subst-mod F81 --tree "(human, (mouse,rat), cow)" -i SS
!phyloFit.mod @phyloFit hmrc.ss --subst-mod HKY85 --tree "(human, (mouse,rat), cow)" -i SS
!phyloFit.mod @phyloFit hmrc.ss --subst-mod REV --tree "(human, (mouse,rat), cow)" -i SS
!phyloFit.mod @phyloFit hmrc.ss --subst-mod UNREST --tree "(human, (mouse,rat), cow)" -i SS
!phyloFit.mod @phyloFit hmrc.ss --subst-mod HKY85 --tree "(human, (mouse,rat), cow)" -i SS -k 4
!phyloFit.mod @phyloFit hmrc.ss --subst-mod REV --tree "(human, (mouse,rat), cow)" -i SS -k 4
!phyloFit.mod @phyloFit hmrc.ss --subst-mod HKY85 --tree "(human, (mouse,rat), cow)" -i SS --EM
!phyloFit.mod @phyloFit hmrc.ss --subst-mod REV --tree "(human, (mouse,rat), cow)" -i SS --EM
!phyloFit.mod @phyloFit hpmrc.ss --subst-mod REV --tree "(hg16, (mm3,rn3), galGal2)" -i SS --gaps-as-bases
!phyloFit.mod !phyloFit.postprob @phyloFit hmrc.ss --subst-mod REV -i SS --init-model rev.mod --post-probs --lnl
!phyloFit.mod @phyloFit hmrc.ss --subst-mod REV --tree "(human, (mouse,rat))" -i SS
msa_view hmrc.ss -i SS --seqs human,mouse,rat --unordered -o SS > hmr.ss
!phyloFit.mod @phyloFit hmr.ss -i SS
msa_view hmrc.ss -i SS --seqs human,mouse --unordered -o SS > hm.ss
!phyloFit.mod @phyloFit hm.ss -i SS
!phyloFit.mod @phyloFit hmrc.ss --subst-mod UNREST --tree "((human, mouse), cow)" -i SS --ancestor cow
rm -f phyloFit.mod phyloFit.postprob hmr.ss hm.ss

******************** msa_view ********************

# These are the same tests implemented in $PHAST/test/Makefile
# Should add more thorough testing at some point

# first make some files to use (all these commands are tested later)
msa_view hmrc.ss -i SS --end 10000 > hmrc.fa
msa_view hmrc.fa -o PHYLIP > hmrc.ph
msa_view hmrc.ph --out-format MPM --in-format PHYLIP > hmrc.mpm
msa_view hmrc.fa -o SS > hmrc_short_a.ss



@msa_view hmrc.ss -i SS --end 10000
@msa_view hmrc.fa -o PHYLIP
@msa_view hmrc.ph -i PHYLIP
@msa_view hmrc.ph --out-format MPM --in-format PHYLIP
@msa_view hmrc.mpm --in-format MPM
@msa_view hmrc.fa -o SS
@msa_view hmrc.ss --end 10000 -i SS -o SS
@msa_view --seqs human,cow hmrc_short_a.ss -i SS -o SS | msa_view - -i SS
@msa_view --seqs human,cow hmrc.fa
@msa_view hmrc.ss --gap-strip ANY -i SS -o SS | msa_view - -i SS
@msa_view hmrc.ss -i SS | msa_view - --gap-strip ANY
@msa_view hmrc.fa --seqs human,cow --gap-strip ALL
@msa_view hmrc_short_a.ss --seqs human,cow -i SS --gap-strip ALL
@msa_view hmrc.ss --start 10000 --end 20000 -i SS
@msa_view hmrc.ss -i SS | msa_view - --start 10000 --end 20000
@msa_view hmrc.ss --summary -i SS

rm -f hmrc.fa hmrc.ph hmrc.mpm hmrc_short_a.ss
