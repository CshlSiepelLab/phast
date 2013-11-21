# Basic syntax:
# Commands starting with @ mean "test this command", ie, run it in
# two different versions of PHAST and compare results.
# commands starting with @ can be preceded with files starting with !,
# which indicates to compare these files as well as stdout/stderr.
# You can also preceded a command with "-stdout" or "-stderr" to avoid
# comparing stdout/stderr, otherwise they are compared by default.
# Commands that do not start with @ are run only once (using the default
# $PATH).  Lines starting with "*" are comments.

******************* phastOdds ******************

phyloFit --tree "((hg16, panTro1), (mm3, rn3))" -o hpmr hpmrc.ss --quiet

# let's go through the examples given in phastOdds--help
tree_doctor --scale 3.0 hpmr.mod > hpmr_fast.mod
@phastOdds --background-mods hpmr.mod --feature-mods hpmr_fast.mod --features elemcfa.bed hpmrc.fa
@phastOdds --background-mods hpmr.mod --feature-mods hpmr_fast.mod --window 100 hpmrc.fa -v
@phastOdds --background-mods hpmr.mod --feature-mods hpmr_fast.mod -y hpmrc.fa
#@phastOdds --background-mods hpmr.mod --feature-mods hpmr_fast.mod --background-hmm coding.hmm --feature-hmm noncoding.hmm --features elemcfa.bed --output-bed hpmrc.fa

#@phastCons  hpmrc.ss hpmr.mod,hpmr_fast.mod
#@phastCons  hpmrc.ss --rho  0.5 hpmr.mod
#msa_view --end 1000 -o SS hpmrc.ss > hpmrc_short.ss
#!tempTree.cons.mod !tempTree.noncons.mod @phastCons  hpmrc.ss hpmr.mod --estimate-trees tempTree
#@phastCons --target-coverage 0.25 --expected-length 12 hpmrc.ss hpmr.mod,hpmr_fast.mod
#@phastCons --transitions 0.01,0.02 hpmrc.ss hpmr.mod,hpmr_fast.mod
#!tempRho.cons.mod !tempRho.noncons.mod @phastCons --target-coverage 0.25 --expected-length 12 --estimate-rho tempRho --no-post-probs hpmrc.ss hpmr.mod
#!tempTree.cons.mod !tempTree.noncons.mod @phastCons --target-coverage 0.25 --expected-length 12 --estimate-trees tempTree --no-post-probs hpmrc_short.ss $
#@phastCons hpmrc.ss hpmr.mod,hpmr_fast.mod 
#!tempTree.cons.mod !tempTree.noncons.mod @phastCons --target-coverage 0.25 --estimate-rho tempTree hpmrc.ss hpmr.mod

# now do at least minimal testing on each option.  (or at least the ones that seem useful)
#--gc,-G
#!tempRho.cons.mod !tempRho.noncons.mod @phastCons --gc 0.8 --estimate-rho tempRho hpmrc_short.ss hpmr.mod

******************** phyloP ********************

phyloFit hmrc.ss --tree "(human, (mouse,rat), cow)" --quiet
@phyloP --null 10 phyloFit.mod
msa_view -o SS --end 100 hmrc.ss > hmrc_short.ss
@phyloP phyloFit.mod hmrc_short.ss
@phyloP --method LRT --base-by-base phyloFit.mod hmrc.ss
@phyloP --method LRT --mode CONACC --wig-scores phyloFit.mod hmrc.ss
@phyloP -d 12345 --method SCORE --mode NNEUT --base-by-base phyloFit.mod hmrc.ss
@phyloP --method GERP --mode ACC --wig-scores phyloFit.mod hmrc.ss
@phyloP --method GERP --base-by-base phyloFit.mod hmrc.ss
@phyloP -d 12345 --method SCORE --wig-scores phyloFit.mod hmrc.ss
@phyloP --method LRT --wig-scores phyloFit.mod hmrc.ss
@phyloP --method LRT --base-by-base phyloFit.mod hmrc.ss
@phyloP -d 12345 --method SCORE --wig-scores --refidx 2 phyloFit.mod hmrc.ss
echo -e "chr1\t0\t10\nchr1\t50\t100\nchr1\t200\t300" > temp.bed
@phyloP --method LRT --mode CONACC --features temp.bed phyloFit.mod hmrc.ss
@phyloP --method SCORE --features temp.bed -g phyloFit.mod hmrc.ss
tree_doctor --name-ancestors phyloFit.mod > phyloFit-named.mod
@phyloP --method LRT --mode CONACC --subtree mouse-rat --base-by-base phyloFit-named.mod hmrc.ss
@phyloP --method LRT --mode ACC --branch mouse-rat -w phyloFit-named.mod hmrc.ss
@phyloP --posterior phyloFit.mod hmrc_short.ss
@phyloP --fit-model -w --subtree mouse-rat phyloFit-named.mod hmrc_short.ss
@phyloP --epsilon 0.0001 phyloFit-named.mod hmrc_short.ss
@phyloP --confidence-interval 0.05 phyloFit-named.mod hmrc_short.ss
@phyloP --quantiles --null 100 phyloFit-named.mod
rm -f hmrc_short.ss phyloFit.mod phyloFit-named.mod temp.bed



******************** phastCons ********************

# these are the two examples in the original test Makefile
!elements.bed @phastCons hpmrc.ss hpmrc-rev-dg-global.mod --nrates 20 --transitions .08,.008  --viterbi elements.bed --seqname chr22
phyloFit --tree "((hg16, panTro1), (mm3, rn3))" -o hpmr hpmrc.ss --quiet
!elements-4way.bed @phastCons hpmrc.ss hpmr.mod --nrates 20 --transitions .08,.008  --viterbi elements-4way.bed --seqname chr22


# let's go through the examples given in phastCons --help
tree_doctor --scale 3.0 hpmr.mod > hpmr_fast.mod
@phastCons  hpmrc.ss hpmr.mod,hpmr_fast.mod
@phastCons  hpmrc.ss --rho  0.5 hpmr.mod
msa_view --end 1000 -o SS hpmrc.ss > hpmrc_short.ss
!tempTree.cons.mod !tempTree.noncons.mod @phastCons  hpmrc.ss hpmr.mod --estimate-trees tempTree
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
phyloFit --quiet --tree "((hg16, panTro1), (mm3, rn3))" -o hpmr-ratevar --alpha 2.0 --nrates 3 hpmrc.ss
!tempTree.cons.mod !tempTree.noncons.mod @phastCons -k 4 --estimate-trees tempTree hpmrc_short.ss hpmr-ratevar.mod
!tempTree.cons.mod !tempTree.noncons.mod @phastCons -k 2,3 --estimate-trees tempTree hpmrc_short.ss hpmr-ratevar.mod 
#--transitions,-t
@phastCons -t 0.01,0.02 hpmrc_short.ss hpmr.mod
@phastCons -t ~0.02,0.03 hpmrc_short.ss hpmr.mod
#--target-coverage is tested in examples above
#--expected-length is tested in examples above
#--msa-format,-i
@phastCons hpmrc.fa hpmr.mod
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
rm -f log.txt
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

@phastCons --FC hpmrc_short.ss hpmr-ratevar.mod
#--lambda
@phastCons --FC --lambda 0.1 hpmrc_short.ss hpmr-ratevar.mod
@phastCons --FC --lambda ~0.1 hpmrc_short.ss hpmr-ratevar.mod

#--coding-potential
msa_view -o SS --tuple-size 3 hmrc_correct.fa > hmrc_correct_tuple3.ss
# don't compare stderr because it contains references to the phast home dir
-stderr @phastCons --coding-potential hmrc_correct_tuple3.ss
#@phastCons --coding-potential hmrc_correct.fa
#--extrapolate and alias
tree_doctor --prune mm3 hpmr.mod --rename "hg16 -> human; panTro1 -> chimp; rn3 -> rat" > hpmr_pruned.mod
-stderr @phastCons --alias "hg16=human; panTro1=chimp; rn3=rat; mm3=mouse; galGal2=chicken" --extrapolate default hpmrc_short.ss hpmr_pruned.mod

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
rm -f hpmr.mod hpmr-ratevar.mod hpmr_fast.mod hpmrc_short.ss hmrc_correct_tuple3.ss hpmr_pruned.mod hpmr_slow.mod



******************** phyloFit ********************

# These are the same tests implemented in $PHAST/test/Makefile
# Should add more thorough testing at some point

seed=124
!phyloFit.mod @phyloFit hmrc.ss --subst-mod JC69 --tree "(human, (mouse,rat), cow)"
# try at least one other kind of input
!phyloFit.mod @phyloFit hmrc.ss --subst-mod JC69 --tree "((((human,chimp), (mouse,rat)), cow), chicken)"
!phyloFit.mod @phyloFit hmrc.ss --subst-mod F81 --tree "(human, (mouse,rat), cow)"
!phyloFit.mod @phyloFit hmrc.ss --subst-mod HKY85 --tree "(human, (mouse,rat), cow)"
!phyloFit.mod @phyloFit -D 12345 hmrc.ss --subst-mod REV --tree "(human, (mouse,rat), cow)"
!phyloFit.mod @phyloFit -D 12345 hmrc.ss --subst-mod UNREST --tree "(human, (mouse,rat), cow)"
!phyloFit.mod @phyloFit hmrc.ss --subst-mod HKY85 --tree "(human, (mouse,rat), cow)" -k 4
!phyloFit.mod @phyloFit -D 12345 hmrc.ss --subst-mod REV --tree "(human, (mouse,rat), cow)" -k 4
!phyloFit.mod @phyloFit -D 12345 hpmrc.ss --subst-mod REV --tree "(hg16, (mm3,rn3), galGal2)" --gaps-as-bases
!phyloFit.mod !phyloFit.postprob @phyloFit hmrc.ss --subst-mod REV --init-model rev.mod --post-probs --lnl
!phyloFit.mod @phyloFit -D 12345 hmrc.ss --subst-mod REV --tree "(human, (mouse,rat))"
!phyloFit.mod @phyloFit -D 12345 hpmrc.fa --tree "(((hg16,panTro2),(rn3,mm3)),galGal2)"

# try most of the above again with --EM
!phyloFit.mod @phyloFit hmrc.ss --EM --subst-mod JC69 --tree "(human, (mouse,rat), cow)"
!phyloFit.mod @phyloFit hmrc.ss --EM --subst-mod JC69 --tree "((((human,chimp), (mouse,rat)), cow), chicken)"
!phyloFit.mod @phyloFit hmrc.ss --EM --subst-mod F81 --tree "(human, (mouse,rat), cow)"
!phyloFit.mod @phyloFit hmrc.ss --EM --subst-mod HKY85 --tree "(human, (mouse,rat), cow)"
!phyloFit.mod @phyloFit hmrc.ss -D 12345 --EM --subst-mod REV --tree "(human, (mouse,rat), cow)"
!phyloFit.mod @phyloFit hmrc.ss -D 12345 --EM --subst-mod UNREST --tree "(human, (mouse,rat), cow)"
!phyloFit.mod @phyloFit hmrc.ss --EM --subst-mod HKY85 --tree "(human, (mouse,rat), cow)" -k 4
!phyloFit.mod @phyloFit hmrc.ss -D 12345 --EM --subst-mod REV --tree "(human, (mouse,rat), cow)" -k 4
!phyloFit.mod @phyloFit hpmrc.ss -D 12345 --EM --subst-mod REV --tree "(hg16, (mm3,rn3), galGal2)" --gaps-as-bases
!phyloFit.mod !phyloFit.postprob @phyloFit hmrc.ss -D 12345 --subst-mod REV --EM --init-model rev.mod --post-probs --lnl
!phyloFit.mod @phyloFit hmrc.ss -D 12345 --subst-mod REV --EM --tree "(human, (mouse,rat))"
!phyloFit.mod @phyloFit -D 12345 hpmrc.fa --EM --tree "(((hg16,panTro2),(rn3,mm3)),galGal2)"
!phyloFit.mod @phyloFit -D 12345 hpmrc.fa --EM --nrates 3 --tree "(((hg16,panTro2),(rn3,mm3)),galGal2)"

# test some of the higher order models (they are slow so use small simulated data set)
base_evolve --nsites 100 rev.mod > simulated.fa
# move all these commands to another section so they aren't run automatically with phastCons check
#!phyloFit.mod @phyloFit --subst-mod R2 --tree "((human,(mouse,rat)),cow)" simulated.fa
#!phyloFit.mod @phyloFit --subst-mod R2S --tree "((human,(mouse,rat)),cow)" simulated.fa
#!phyloFit.mod @phyloFit --subst-mod U2 --tree "((human,(mouse,rat)),cow)" simulated.fa
#!phyloFit.mod @phyloFit --subst-mod U2S --tree "((human,(mouse,rat)),cow)" simulated.fa
#!phyloFit.mod @phyloFit --subst-mod R3 --tree "((human,(mouse,rat)),cow)" simulated.fa
#!phyloFit.mod @phyloFit --subst-mod R3S --tree "((human,(mouse,rat)),cow)" simulated.fa
#!phyloFit.mod @phyloFit --subst-mod U3 --tree "((human,(mouse,rat)),cow)" simulated.fa
#!phyloFit.mod @phyloFit --subst-mod U3S --tree "((human,(mouse,rat)),cow)" simulated.fa

msa_view hmrc.ss --seqs human,mouse,rat --unordered -o SS > hmr.ss
!phyloFit.mod @phyloFit hmr.ss -D 12345
msa_view hmrc.ss --seqs human,mouse --unordered -o SS > hm.ss
!phyloFit.mod @phyloFit hm.ss -D 12345
!phyloFit.mod @phyloFit hmrc.ss -D 12345 --subst-mod UNREST --tree "((human, mouse), cow)" --ancestor cow
!phyloFit.mod @phyloFit hmrc.ss -D 12345 --subst-mod UNREST --tree "((human, (mouse, rat)mouse-rat), cow)" --ignore-branches mouse-rat
!phyloFit.mod @phyloFit hmrc.ss -D 12345 --precision LOW --tree "((human, (mouse, rat)), cow)"
!phyloFit.mod @phyloFit hmrc.ss --init-model rev-em.mod
#!phyloFit.mod @phyloFit hmrc.ss --init-random --tree "((human, (mouse, rat)), cow)"
#parsimony
!phyloFit.mod @phyloFit hmrc.ss -D 12345 --init-parsimony --tree "((human, (mouse, rat)), cow)"
!parsimony.txt @phyloFit hmrc.ss --print-parsimony parsimony.txt --tree "((human, (mouse, rat)), cow)"
#clock
!phyloFit.mod @phyloFit hmrc.ss -D 12345 --tree "((human, (mouse,rat)), cow)" --clock
#scale-only
!phyloFit.mod @phyloFit hmrc.ss --init-mod rev-em.mod --scale-only 
#scale-subtree
tree_doctor --name-ancestors rev-em.mod --scale 2.0 > rev-em-scaled-named.mod
!phyloFit.mod @phyloFit hmrc.ss --init-mod rev-em-scaled-named.mod --scale-only --scale-subtree mouse-rat
!phyloFit.mod @phyloFit hmrc.ss --init-mod rev-em-scaled-named.mod --scale-only --scale-subtree mouse-rat:gain
!phyloFit.mod @phyloFit hmrc.ss --init-mod rev-em-scaled-named.mod --scale-only --scale-subtree mouse-rat:loss
# estimate-freqs
!phyloFit.mod @phyloFit hmrc.ss --init-mod rev-em.mod --estimate-freqs
# sym-freqs
!phyloFit.mod @phyloFit hmrc.ss --sym-freqs --init-mod rev-em.mod
# no-freqs
!phyloFit.mod @phyloFit hmrc.ss --no-freqs --init-mod rev-em.mod
# no-rates
!phyloFit.mod @phyloFit hmrc.ss --no-rates --init-mod rev-em.mod
# ancestor tested above
# error
!errors.txt !phyloFit.mod @phyloFit hmrc.ss --no-rates --init-mod rev-em.mod --error errors.txt
# no-opt
!phyloFit.mod @phyloFit hmrc.ss --init-mod rev-em.mod --no-opt branches
!phyloFit.mod @phyloFit hmrc.ss --init-mod rev-em.mod --no-opt backgd
!phyloFit.mod @phyloFit hmrc.ss --init-mod rev-em.mod --no-opt ratematrix
# bound
!phyloFit.mod @phyloFit hmrc.ss --init-mod rev-em.mod --bound "branches[0.1,0.2]"
!phyloFit.mod @phyloFit hmrc.ss --init-mod rev-em.mod --bound "branches[0.1,]"
#--nrates
!phyloFit.mod @phyloFit hmrc.ss --init-mod rev-em.mod --nrates 4
#--alpha
!phyloFit.mod @phyloFit hmrc.ss --init-mod rev-em.mod --nrates 4 --alpha 5.2
#--rate-constants
!phyloFit.mod @phyloFit hmrc.ss --init-mod rev-em.mod --nrates 4 --rate-constants 10.0,6.0,1.0,0.1

#--features
!phyloFit.bed_feature.mod !phyloFit.background.mod @phyloFit hpmrc.ss -D 12345 --tree "(((hg16,panTro2),(rn3,mm3)),galGal2)" --features elements_correct.bed

#--markov
!phyloFit.mod @phyloFit --markov --subst-mod R2 --tree "((human,(mouse,rat)),cow)" simulated.fa
#--non-overlapping, --min-informative
!phyloFit.mod @phyloFit --non-overlapping --subst-mod R2 --tree "((human,(mouse,rat)),cow)" --min-informative 25 simulated.fa

#--alt-mod
!phyloFit.mod @phyloFit -D 12345 --tree "((human,(mouse#MR,rat#MR)mouse-rat),cow)" --alt-mod MR:HKY85 --subst-mod REV simulated.fa
!phyloFit.mod @phyloFit -D 12345 --tree "((human,(mouse,rat)mouse-rat),cow)" --alt-mod MR:HKY85 --label-subtree mouse-rat:MR --subst-mod REV simulated.fa
!phyloFit.mod @phyloFit -D 12345 --tree "((human,(mouse#MR,rat#MR)mouse-rat),cow)" --alt-mod MR:ratematrix,backgd --subst-mod REV simulated.fa
!phyloFit.mod @phyloFit -D 12345 --tree "((human,(mouse#MR,rat#MR)mouse-rat#MR),cow)" --alt-mod MR:backgd --subst-mod REV simulated.fa --estimate-freqs
!phyloFit.mod @phyloFit -D 12345 --tree "((human,(mouse,rat)mouse-rat),cow)" --label-subtree "mouse-rat+:MR" --alt-mod MR:bgc,sel --subst-mod REV simulated.fa --estimate-freqs
!phyloFit.mod @phyloFit -D 12345 --tree "((human,(mouse,rat#MR)mouse-rat#MR),cow)" --alt-mod MR:backgd --subst-mod REV simulated.fa --estimate-freqs
!phyloFit.mod @phyloFit -D 12345 --tree "((human,(mouse,rat)mouse-rat),cow)" --label-branches rat,mouse-rat:MR --alt-mod MR:backgd --subst-mod REV simulated.fa --estimate-freqs


#--post-probs 
!phyloFit.mod !phyloFit.postprob @phyloFit -D 12345 --tree "((human,(mouse,rat)mouse-rat),cow)" --post-probs simulated.fa
!phyloFit.mod !phyloFit.postprob @phyloFit -D 12345 --tree "((human,(mouse,rat)mouse-rat),cow)" --post-probs --alt-mod mouse:HKY85 simulated.fa

#--expected-subs, --expected-total-subs
!phyloFit.mod !phyloFit.expsub !phyloFit.exptotsub @phyloFit -D 12345 --tree "((human,(mouse,rat)mouse-rat),cow)" --expected-subs --expected-total-subs simulated.fa
#--column-probs
!phyloFit.mod !phyloFit.colprobs @phyloFit --init-mod rev.mod hmrc.ss --column-probs

#--windows
!phyloFit.win-1.mod !phyloFit.win-2.mod !phyloFit.win-3.mod !phyloFit.win-4.mod !phyloFit.win-5.mod !phyloFit.win-6.mod @phyloFit --tree "((human,(mouse,rat)mouse-rat),cow)" --windows 20,15 simulated.fa --min-informative 15 -D 12345
#--windows-explicit
!phyloFit.win-1.mod !phyloFit.win-2.mod @phyloFit  --tree "((human,(mouse,rat)mouse-rat),cow)" --windows-explicit 1,20,25,45 simulated.fa --min-informative 15 -D 12345
echo -e "1\t20\n25\t45" > windows.txt
!phyloFit.win-1.mod !phyloFit.win-2.mod @phyloFit  --tree "((human,(mouse,rat)mouse-rat),cow)" --windows-explicit '*windows.txt' simulated.fa --min-informative 15 -D 12345
rm -f windows.txt


rm -f phyloFit.mod phyloFit.postprob hmr.ss hm.ss rev-em-scaled-named.mod simulated.fa


# TODO: phyloFit options not currently tested above:
# --min-informative,--log,--catmap,--do-cats,--reverse-groups

******************** phyloFit-slow ********************
# These models are so slow only test them if we really need to
base_evolve --nsites 50 rev.mod > simulated.fa
!phyloFit.mod @phyloFit --subst-mod R2 --tree "((human,(mouse,rat)),cow)" simulated.fa
!phyloFit.mod @phyloFit --subst-mod R2S --tree "((human,(mouse,rat)),cow)" simulated.fa
!phyloFit.mod @phyloFit --subst-mod U2 --tree "((human,(mouse,rat)),cow)" simulated.fa -D 12345
!phyloFit.mod @phyloFit --subst-mod U2S --tree "((human,(mouse,rat)),cow)" simulated.fa -D 12345
#!phyloFit.mod @phyloFit --subst-mod R3 --tree "((human,(mouse,rat)),cow)" simulated.fa
#!phyloFit.mod @phyloFit --subst-mod R3S --tree "((human,(mouse,rat)),cow)" simulated.fa
#!phyloFit.mod @phyloFit --subst-mod U3 --tree "((human,(mouse,rat)),cow)" simulated.fa
#!phyloFit.mod @phyloFit --subst-mod U3S --tree "((human,(mouse,rat)),cow)" simulated.fa
rm -f simulated.fa

******************** msa_view ********************

# These are the same tests implemented in $PHAST/test/Makefile
# Should add more thorough testing at some point

# first make some files to use (all these commands are tested later)
msa_view hmrc.ss --end 10000 > hmrc.fa
msa_view hmrc.fa -o PHYLIP > hmrc.ph
msa_view hmrc.ph --out-format MPM --in-format PHYLIP > hmrc.mpm
msa_view hmrc.fa -o SS > hmrc_short_a.ss



@msa_view hmrc.ss --end 10000
@msa_view hmrc.fa -o PHYLIP
@msa_view hmrc.ph
@msa_view hmrc.ph --out-format MPM
@msa_view hmrc.mpm
@msa_view hmrc.fa -o SS
@msa_view hmrc.ss --end 10000 -o SS
@msa_view --seqs human,cow hmrc_short_a.ss -o SS | msa_view -
@msa_view --seqs human,cow hmrc.fa
@msa_view hmrc.ss --gap-strip ANY -o SS | msa_view -
@msa_view hmrc.ss | msa_view - --gap-strip ANY
@msa_view hmrc.fa --seqs human,cow --gap-strip ALL
@msa_view hmrc_short_a.ss --seqs human,cow --gap-strip ALL
@msa_view hmrc.ss --start 10000 --end 20000
@msa_view hmrc.ss | msa_view - --start 10000 --end 20000
@msa_view hmrc.ss --summary
@msa_view -o SS chr22.14500000-15500000.maf
@msa_view -o SS --refseq chr22.14500000-15500000.fa chr22.14500000-15500000.maf
@msa_view -o SS --gap-strip 1 chr22.14500000-15500000.maf

refeature chr22.14500000-15500000.gp | awk -v OFS="\t" '{start=$4-14500000; end=$5-14500000; print "hg17."$1,$2,$3,start,end,$6,$7,$8,$9}' > temp.gff
@msa_view -o SS --features temp.gff chr22.14500000-15500000.maf
@msa_view -o SS --features temp.gff --4d chr22.14500000-15500000.maf

rm -f hmrc.fa hmrc.ph hmrc.mpm hmrc_short_a.ss temp.gff


******************** tree_doctor ********************

# this is just a start for some recently added options; many more tests could/should
# be added
echo '(human, (mouse, rat), cow)' > tree.nh
@tree_doctor --newick tree.nh
@tree_doctor --newick --label-branches mouse,rat:MR tree.nh
@tree_doctor --newick --name-ancestors --label-branches mouse-rat,mouse,rat:MR tree.nh
@tree_doctor --newick --name-ancestors --label-subtree mouse-rat:MR tree.nh
@tree_doctor --newick --name-ancestors --label-subtree mouse-rat+:MR tree.nh
@tree_doctor --newick --name-ancestors --label-subtree mouse-rat+:MR --label-branches mouse:mouse tree.nh
phyloFit hmrc.ss --tree "(human, (mouse,rat), cow)" --quiet
@tree_doctor phyloFit.mod
@tree_doctor --label-branches mouse,rat:MR phyloFit.mod
@tree_doctor --name-ancestors --label-branches mouse-rat,mouse,rat:MR phyloFit.mod
@tree_doctor --name-ancestors --label-subtree mouse-rat:MR phyloFit.mod
@tree_doctor --name-ancestors --label-subtree mouse-rat+:MR phyloFit.mod

rm -f phyloFit.mod tree.nh
