# Basic syntax:
# commands starting with @ mean "test this command", ie, run it in
# two different versions of PHAST and compare results.
# commands starting with @ can be preceded with files starting with !,
# which indicates to compare these files as well as stdout/stderr.

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
@phyloP -i SS --method LRT --mode CONACC --subtree mouse-rat --base-by-base phyloFit-named.mod hmrc.ss
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
phyloFit -i SS --tree "((hg16, panTro1), (mm3, rn3))" -o hpmr hpmrc.ss --quiet
!elements-4way.bed @phastCons hpmrc.ss hpmr.mod --nrates 20 --transitions .08,.008  --viterbi elements-4way.bed --seqname chr22


# let's go through the examples given in phastCons --help
tree_doctor --scale 3.0 hpmr.mod > hpmr_fast.mod
@phastCons  hpmrc.ss hpmr.mod,hpmr_fast.mod
@phastCons  hpmrc.ss --rho  0.5 hpmr.mod
msa_view -i SS --end 1000 -o SS hpmrc.ss > hpmrc_short.ss
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
phyloFit --quiet -i SS --tree "((hg16, panTro1), (mm3, rn3))" -o hpmr-ratevar --alpha 2.0 --nrates 3 hpmrc.ss
!tempTree.cons.mod !tempTree.noncons.mod @phastCons -k 12 --estimate-trees tempTree hpmrc_short.ss hpmr-ratevar.mod
!tempTree.cons.mod !tempTree.noncons.mod @phastCons -k 2,3 --estimate-trees tempTree hpmrc_short.ss hpmr-ratevar.mod 
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
msa_view -i FASTA -o SS --tuple-size 3 hmrc_correct.fa > hmrc_correct_tuple3.ss
# don't compare stderr because it contains references to the phast home dir
-stderr @phastCons --coding-potential hmrc_correct_tuple3.ss
#@phastCons --coding-potential -i FASTA hmrc_correct.fa
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

!phyloFit.mod @phyloFit hmrc.ss --subst-mod JC69 --tree "(human, (mouse,rat), cow)" -i SS
# try at least one other kind of input
!phyloFit.mod @phyloFit -i FASTA --subst-mod REV+GC --tree "(human, (mouse, rat), cow)" hmrc_correct.fa
!phyloFit.mod @phyloFit hmrc.ss --subst-mod JC69 --tree "((((human,chimp), (mouse,rat)), cow), chicken)" -i SS
!phyloFit.mod @phyloFit hmrc.ss --subst-mod F81 --tree "(human, (mouse,rat), cow)" -i SS
!phyloFit.mod @phyloFit hmrc.ss --subst-mod HKY85 --tree "(human, (mouse,rat), cow)" -i SS
!phyloFit.mod @phyloFit hmrc.ss --subst-mod REV --tree "(human, (mouse,rat), cow)" -i SS
!phyloFit.mod @phyloFit hmrc.ss --subst-mod UNREST --tree "(human, (mouse,rat), cow)" -i SS
!phyloFit.mod @phyloFit hmrc.ss --subst-mod HKY85 --tree "(human, (mouse,rat), cow)" -i SS -k 4
!phyloFit.mod @phyloFit hmrc.ss --subst-mod REV --tree "(human, (mouse,rat), cow)" -i SS -k 4
!phyloFit.mod @phyloFit hpmrc.ss --subst-mod REV --tree "(hg16, (mm3,rn3), galGal2)" -i SS --gaps-as-bases
!phyloFit.mod !phyloFit.postprob @phyloFit hmrc.ss --subst-mod REV -i SS --init-model rev.mod --post-probs --lnl
!phyloFit.mod @phyloFit hmrc.ss --subst-mod REV --tree "(human, (mouse,rat))" -i SS
msa_view  -i SS -o FASTA hpmrc.ss > hpmrc.fa
!phyloFit.mod @phyloFit -i FASTA hpmrc.fa --tree "(((hg16,panTro2),(rn3,mm3)),galGal2)"

# try most of the above again with --EM
!phyloFit.mod @phyloFit hmrc.ss --EM --subst-mod JC69 --tree "(human, (mouse,rat), cow)" -i SS
!phyloFit.mod @phyloFit hmrc.ss --EM --subst-mod JC69 --tree "((((human,chimp), (mouse,rat)), cow), chicken)" -i SS
!phyloFit.mod @phyloFit hmrc.ss --EM --subst-mod F81 --tree "(human, (mouse,rat), cow)" -i SS
!phyloFit.mod @phyloFit hmrc.ss --EM --subst-mod HKY85 --tree "(human, (mouse,rat), cow)" -i SS
!phyloFit.mod @phyloFit hmrc.ss --EM --subst-mod REV --tree "(human, (mouse,rat), cow)" -i SS
!phyloFit.mod @phyloFit hmrc.ss --EM --subst-mod UNREST --tree "(human, (mouse,rat), cow)" -i SS
!phyloFit.mod @phyloFit hmrc.ss --EM --subst-mod HKY85 --tree "(human, (mouse,rat), cow)" -i SS -k 4
!phyloFit.mod @phyloFit hmrc.ss --EM --subst-mod REV --tree "(human, (mouse,rat), cow)" -i SS -k 4
!phyloFit.mod @phyloFit hpmrc.ss --EM --subst-mod REV --tree "(hg16, (mm3,rn3), galGal2)" -i SS --gaps-as-bases
!phyloFit.mod !phyloFit.postprob @phyloFit hmrc.ss --subst-mod REV --EM -i SS --init-model rev.mod --post-probs --lnl
!phyloFit.mod @phyloFit hmrc.ss --subst-mod REV --EM --tree "(human, (mouse,rat))" -i SS
!phyloFit.mod @phyloFit -i FASTA hpmrc.fa --tree "(((hg16,panTro2),(rn3,mm3)),galGal2)"

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

msa_view hmrc.ss -i SS --seqs human,mouse,rat --unordered -o SS > hmr.ss
!phyloFit.mod @phyloFit hmr.ss -i SS
msa_view hmrc.ss -i SS --seqs human,mouse --unordered -o SS > hm.ss
!phyloFit.mod @phyloFit hm.ss -i SS
!phyloFit.mod @phyloFit hmrc.ss --subst-mod UNREST --tree "((human, mouse), cow)" -i SS --ancestor cow
!phyloFit.mod @phyloFit hmrc.ss --subst-mod UNREST --tree "((human, (mouse, rat)mouse-rat), cow)" -i SS --ignore-branches mouse-rat
!phyloFit.mod @phyloFit hmrc.ss --precision LOW --tree "((human, (mouse, rat)), cow)" -i SS
!phyloFit.mod @phyloFit hmrc.ss --init-model rev-em.mod -i SS
#!phyloFit.mod @phyloFit hmrc.ss --init-random --tree "((human, (mouse, rat)), cow)" -i SS
#parsimony
!phyloFit.mod @phyloFit hmrc.ss --init-parsimony --tree "((human, (mouse, rat)), cow)" -i SS
!parsimony.txt @phyloFit hmrc.ss --print-parsimony parsimony.txt --tree "((human, (mouse, rat)), cow)" -i SS
#clock
!phyloFit.mod @phyloFit hmrc.ss --tree "((human, (mouse,rat)), cow)" -i SS --clock
#scale-only
!phyloFit.mod @phyloFit hmrc.ss --init-mod rev-em.mod --scale-only -iSS 
#scale-subtree
tree_doctor --name-ancestors rev-em.mod --scale 2.0 > rev-em-scaled-named.mod
!phyloFit.mod @phyloFit hmrc.ss --init-mod rev-em-scaled-named.mod --scale-only --scale-subtree mouse-rat -i SS
!phyloFit.mod @phyloFit hmrc.ss --init-mod rev-em-scaled-named.mod --scale-only --scale-subtree mouse-rat:gain -i SS
!phyloFit.mod @phyloFit hmrc.ss --init-mod rev-em-scaled-named.mod --scale-only --scale-subtree mouse-rat:loss -i SS
# estimate-freqs
!phyloFit.mod @phyloFit hmrc.ss --init-mod rev-em.mod --estimate-freqs -i SS
# sym-freqs
!phyloFit.mod @phyloFit hmrc.ss -i SS --sym-freqs --init-mod rev-em.mod
# no-freqs
!phyloFit.mod @phyloFit hmrc.ss -i SS --no-freqs --init-mod rev-em.mod
# no-rates
!phyloFit.mod @phyloFit hmrc.ss -i SS --no-rates --init-mod rev-em.mod
# ancestor tested above
# error
!errors.txt !phyloFit.mod @phyloFit hmrc.ss -i SS --no-rates --init-mod rev-em.mod --error errors.txt
# no-opt
!phyloFit.mod @phyloFit hmrc.ss -i SS --init-mod rev-em.mod --no-opt branches
!phyloFit.mod @phyloFit hmrc.ss -i SS --init-mod rev-em.mod --no-opt backgd
!phyloFit.mod @phyloFit hmrc.ss -i SS --init-mod rev-em.mod --no-opt ratematrix
# bound
!phyloFit.mod @phyloFit hmrc.ss -i SS --init-mod rev-em.mod --bound "branches[0.1,0.2]"
!phyloFit.mod @phyloFit hmrc.ss -i SS --init-mod rev-em.mod --bound "branches[0.1,]"
#--nrates
!phyloFit.mod @phyloFit hmrc.ss -i SS --init-mod rev-em.mod --nrates 4
#--alpha
!phyloFit.mod @phyloFit hmrc.ss -i SS --init-mod rev-em.mod --nrates 4 --alpha 5.2
#--rate-constants
!phyloFit.mod @phyloFit hmrc.ss -i SS --init-mod rev-em.mod --nrates 4 --rate-constants 10.0,6.0,1.0,0.1

#--features
!phyloFit.bed_feature.mod !phyloFit.background.mod @phyloFit hpmrc.ss --tree "(((hg16,panTro2),(rn3,mm3)),galGal2)" --features elements_correct.bed -i SS

#--markov
!phyloFit.mod @phyloFit --markov --subst-mod R2 --tree "((human,(mouse,rat)),cow)" simulated.fa
#--non-overlapping, --min-informative
!phyloFit.mod @phyloFit --non-overlapping --subst-mod R2 --tree "((human,(mouse,rat)),cow)" --min-informative 25 simulated.fa

#--alt-mod
!phyloFit.mod @phyloFit --tree "((human,(mouse,rat)mouse-rat),cow)" --alt-mod mouse-rat:HKY85 --subst-mod REV simulated.fa
!phyloFit.mod @phyloFit --tree "((human,(mouse,rat)mouse-rat),cow)" --alt-mod mouse-rat:ratematrix --subst-mod REV simulated.fa
!phyloFit.mod @phyloFit --tree "((human,(mouse,rat)mouse-rat),cow)" --alt-mod mouse-rat+:backgd --subst-mod REV simulated.fa --estimate-freqs
!phyloFit.mod @phyloFit --tree "((human,(mouse,rat)mouse-rat),cow)" --alt-mod mouse-rat+:backgd --subst-mod REV --alt-mod mouse+ simulated.fa --estimate-freqs

#--post-probs 
!phyloFit.mod !phyloFit.postprob @phyloFit --tree "((human,(mouse,rat)mouse-rat),cow)" --post-probs simulated.fa
!phyloFit.mod !phyloFit.postprob @phyloFit --tree "((human,(mouse,rat)mouse-rat),cow)" --post-probs --alt-mod mouse+:HKY85 simulated.fa

#--expected-subs, --expected-total-subs
!phyloFit.mod !phyloFit.expsub !phyloFit.exptotsub @phyloFit --tree "((human,(mouse,rat)mouse-rat),cow)" --expected-subs --expected-total-subs simulated.fa
#--column-probs
!phyloFit.mod !phyloFit.colprobs @phyloFit --init-mod rev.mod -i SS hmrc.ss --column-probs

#--windows
!phyloFit.win-1.mod !phyloFit.win-2.mod !phyloFit.win-3.mod !phyloFit.win-4.mod !phyloFit.win-5.mod !phyloFit.win-6.mod @phyloFit --tree "((human,(mouse,rat)mouse-rat),cow)" --windows 20,15 simulated.fa --min-informative 15
#--windows-explicit
!phyloFit.win-1.mod !phyloFit.win-2.mod @phyloFit  --tree "((human,(mouse,rat)mouse-rat),cow)" --windows-explicit 1,20,25,45 simulated.fa --min-informative 15
echo "1\t20\n25\t45" > windows.txt
!phyloFit.win-1.mod !phyloFit.win-2.mod @phyloFit  --tree "((human,(mouse,rat)mouse-rat),cow)" --windows-explicit '*windows.txt' simulated.fa --min-informative 15
rm -f windows.txt


rm -f phyloFit.mod phyloFit.postprob hmr.ss hm.ss rev-em-scaled-named.mod simulated.fa


# TODO: phyloFit options not currently tested above:
# --min-informative,--log,--catmap,--do-cats,--reverse-groups

******************** phyloFit-slow ********************
# These models are so slow only test them if we really need to
base_evolve --nsites 50 rev.mod > simulated.fa
!phyloFit.mod @phyloFit --subst-mod R2 --tree "((human,(mouse,rat)),cow)" simulated.fa
!phyloFit.mod @phyloFit --subst-mod R2S --tree "((human,(mouse,rat)),cow)" simulated.fa
!phyloFit.mod @phyloFit --subst-mod U2 --tree "((human,(mouse,rat)),cow)" simulated.fa
!phyloFit.mod @phyloFit --subst-mod U2S --tree "((human,(mouse,rat)),cow)" simulated.fa
#!phyloFit.mod @phyloFit --subst-mod R3 --tree "((human,(mouse,rat)),cow)" simulated.fa
#!phyloFit.mod @phyloFit --subst-mod R3S --tree "((human,(mouse,rat)),cow)" simulated.fa
#!phyloFit.mod @phyloFit --subst-mod U3 --tree "((human,(mouse,rat)),cow)" simulated.fa
#!phyloFit.mod @phyloFit --subst-mod U3S --tree "((human,(mouse,rat)),cow)" simulated.fa


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
@msa_view -i MAF -o SS chr22.14500000-15500000.maf
@msa_view -i MAF -o SS --refseq chr22.14500000-15500000.fa chr22.14500000-15500000.maf
@msa_view -i MAF -o SS --gap-strip 1 chr22.14500000-15500000.maf

refeature chr22.14500000-15500000.gp | awk -v OFS="\t" '{start=$4-14500000; end=$5-14500000; print $1,$2,$3,start,end,$6,$7,$8,$9}' > temp.gff
@msa_view -i MAF -o SS --features temp.gff chr22.14500000-15500000.maf
@msa_view -i MAF -o SS --features temp.gff --4d chr22.14500000-15500000.maf

rm -f hmrc.fa hmrc.ph hmrc.mpm hmrc_short_a.ss temp.gff
