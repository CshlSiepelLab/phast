# acs 8/4/04
# simple test cases, designed to catch obvious errors
# add cases as needed

all: msa_view phyloFit phastCons

msa_view:
	@echo "*** Testing msa_view ***"
	msa_view hmrc.ss -i SS --end 10000 > hmrc.fa
	@if [[ -n `diff --brief hmrc.fa hmrc_correct.fa` ]] ; then echo "ERROR" ; exit 1 ; fi
	msa_view hmrc.fa -o PHYLIP > hmrc.ph
	msa_view hmrc.ph -i PHYLIP > hmrc.fa
	@if [[ -n `diff --brief hmrc.fa hmrc_correct.fa` ]] ; then echo "ERROR" ; exit 1 ; fi
	msa_view hmrc.ph --out-format MPM --in-format PHYLIP > hmrc.mpm
	msa_view hmrc.mpm --in-format MPM > hmrc.fa
	@if [[ -n `diff --brief hmrc.fa hmrc_correct.fa` ]] ; then echo "ERROR" ; exit 1 ; fi
	msa_view hmrc.fa -o SS > hmrc_short_a.ss
	msa_view hmrc.ss --end 10000 -i SS -o SS > hmrc_short_b.ss
	@if [[ -n `diff --brief hmrc_short_[ab].ss` ]] ; then echo "ERROR" ; exit 1 ; fi
	msa_view --seqs human,cow hmrc_short_a.ss -i SS -o SS | msa_view - -i SS > hm_a.fa
	msa_view --seqs human,cow hmrc.fa > hm_b.fa
	@if [[ -n `diff --brief hm_[ab].fa` ]] ; then echo "ERROR" ; exit 1 ; fi
	msa_view hmrc.ss --gap-strip ANY -i SS -o SS | msa_view - -i SS > hmrc_nogaps_a.fa
	msa_view hmrc.ss -i SS | msa_view - --gap-strip ANY > hmrc_nogaps_b.fa
	@if [[ -n `diff --brief hmrc_nogaps_[ab].fa` ]] ; then echo "ERROR" ; exit 1 ; fi
	msa_view hmrc.fa --seqs human,cow --gap-strip ALL > hm_a.fa
	msa_view hmrc_short_a.ss --seqs human,cow -i SS --gap-strip ALL > hm_b.fa
	@if [[ -n `diff --brief hm_[ab].fa` ]] ; then echo "ERROR" ; exit 1 ; fi
	msa_view hmrc.ss --start 10000 --end 20000 -i SS > hmrc_sub_a.fa
	msa_view hmrc.ss -i SS | msa_view - --start 10000 --end 20000 > hmrc_sub_b.fa
	@if [[ -n `diff --brief hmrc_sub_[ab].fa` ]] ; then echo "ERROR" ; exit 1 ; fi
	msa_view hmrc.ss --summary -i SS > hmrc_summary
	@if [[ -n `diff --brief hmrc_summary hmrc_summary_correct` ]] ; then echo "ERROR" ; exit 1 ; fi
	@echo -e "Passed all tests.\n"
	@rm -f hmrc.ph hmrc.fa hmrc.mpm hmrc_short_[ab].ss hm_[ab].fa hmrc_nogaps_[ab].fa hmrc_sub_[ab].fa hmrc_summary


# still need test for MAFs, reverse complement, aggregate, dinucs,
# category-specific stats, missing data, reorder rows.  Also, maybe
# should choose a better short test case (first 10000 bases all gaps
# in mouse and rat)

phyloFit:
	@echo "*** Testing phyloFit ***"
	phyloFit hmrc.ss --subst-mod JC69 --tree "(human, (mouse,rat), cow)" -i SS --quiet
	@if [[ -n `diff --brief phyloFit.mod jc.mod` ]] ; then echo "ERROR" ; exit 1 ; fi
	phyloFit hmrc.ss --subst-mod JC69 --tree "((((human,chimp), (mouse,rat)), cow), chicken)" -i SS --quiet
	@if [[ -n `diff --brief phyloFit.mod jc.mod` ]] ; then echo "ERROR" ; exit 1 ; fi
	phyloFit hmrc.ss --subst-mod F81 --tree "(human, (mouse,rat), cow)" -i SS --quiet
	@if [[ -n `diff --brief phyloFit.mod f81.mod` ]] ; then echo "ERROR" ; exit 1 ; fi
	phyloFit hmrc.ss --subst-mod HKY85 --tree "(human, (mouse,rat), cow)" -i SS --quiet
	@if [[ -n `diff --brief phyloFit.mod hky.mod` ]] ; then echo "ERROR" ; exit 1 ; fi
	phyloFit hmrc.ss --subst-mod REV --tree "(human, (mouse,rat), cow)" -i SS --quiet
	@if [[ -n `diff --brief phyloFit.mod rev.mod` ]] ; then echo "ERROR" ; exit 1 ; fi
	phyloFit hmrc.ss --subst-mod UNREST --tree "(human, (mouse,rat), cow)" -i SS --quiet
	@if [[ -n `diff --brief phyloFit.mod unrest.mod` ]] ; then echo "ERROR" ; exit 1 ; fi
	phyloFit hmrc.ss --subst-mod HKY85 --tree "(human, (mouse,rat), cow)" -i SS -k 4 --quiet
	@if [[ -n `diff --brief phyloFit.mod hky-dg.mod` ]] ; then echo "ERROR" ; exit 1 ; fi
	phyloFit hmrc.ss --subst-mod REV --tree "(human, (mouse,rat), cow)" -i SS -k 4 --quiet 
	@if [[ -n `diff --brief phyloFit.mod rev-dg.mod` ]] ; then echo "ERROR" ; exit 1 ; fi
	phyloFit hmrc.ss --subst-mod HKY85 --tree "(human, (mouse,rat), cow)" -i SS --EM --quiet
	@if [[ -n `diff --brief phyloFit.mod hky-em.mod` ]] ; then echo "ERROR" ; exit 1 ; fi
	phyloFit hmrc.ss --subst-mod REV --tree "(human, (mouse,rat), cow)" -i SS --EM --quiet
	@if [[ -n `diff --brief phyloFit.mod rev-em.mod` ]] ; then echo "ERROR" ; exit 1 ; fi
	phyloFit hmrc.ss --subst-mod REV --tree "(human, (mouse,rat), cow)" -i SS --gaps-as-bases --quiet  # could use better test case (missing data)
	@if [[ -n `diff --brief phyloFit.mod rev-gaps.mod` ]] ; then echo "ERROR" ; exit 1 ; fi
	phyloFit hmrc.ss --subst-mod REV -i SS --init-model rev.mod --post-probs --lnl --quiet
	@if [[ -n `diff --brief phyloFit.mod rev-lnl.mod` ]] ; then echo "ERROR" ; exit 1 ; fi
	@if [[ -n `diff --brief phyloFit.postprob rev.postprob` ]] ; then echo "ERROR" ; exit 1 ; fi
	phyloFit hmrc.ss --subst-mod REV --tree "(human, (mouse,rat))" -i SS --quiet
	@if [[ -n `diff --brief phyloFit.mod rev-hmr.mod` ]] ; then echo "ERROR" ; exit 1 ; fi
	msa_view hmrc.ss -i SS --seqs human,mouse,rat --unordered -o SS > hmr.ss
	phyloFit hmr.ss -i SS --quiet
	@if [[ -n `diff --brief phyloFit.mod rev-hmr2.mod` ]] ; then echo "ERROR" ; exit 1 ; fi  # FIXME: rev-hmr.mod and rev-hmr2.mod should be equal (eq freqs being obtained from all seqs)
	msa_view hmrc.ss -i SS --seqs human,mouse --unordered -o SS > hm.ss
	phyloFit hm.ss -i SS --quiet
	@if [[ -n `diff --brief phyloFit.mod rev-hm.mod` ]] ; then echo "ERROR" ; exit 1 ; fi  
	phyloFit hmrc.ss --subst-mod UNREST --tree "((human, mouse), cow)" -i SS --ancestor cow --quiet
	@if [[ -n `diff --brief phyloFit.mod unrest-cow-anc.mod` ]] ; then echo "ERROR" ; exit 1 ; fi
	@echo -e "Passed all tests.\n"
	@rm -f phyloFit.mod phyloFit.postprob hmr.ss hm.ss

# still need tests for dinucs, functional categories, scale-only,
# estimate-freqs, empirical rate variation, reverse-groups,
# expected subs, column-probs, windows

phastCons:
	@echo "*** Testing phastCons ***"
	# FIXME: should get a test case with missing data and with more rate variation
	phyloFit -s HKY85 -k 4 hmrc.ss --tree "(human, (mouse,rat), cow)" -i SS --quiet
	phastCons hmrc.ss phyloFit.mod --nrates 10 --transitions .08,.008 --quiet > cons.dat
	@if [[ -n `diff --brief cons.dat cons_correct.dat` ]] ; then echo "ERROR" ; exit 1 ; fi  
	phyloFit -s HKY85 -k 4 hmrc.ss --tree "(human, mouse, cow)" -i SS --quiet
	phastCons hmrc.ss phyloFit.mod --nrates 10 --transitions .08,.008 --quiet > cons-hmr.dat
	@if [[ -n `diff --brief cons-hmr.dat cons-hmr_correct.dat` ]] ; then echo "ERROR" ; exit 1 ; fi  
	@echo -e "Passed all tests.\n"
	@rm -f phyloFit.mod cons.dat cons-hmr.dat

# still need to test esimation of MLE for transition probs, coding potential, felsenstein/churchill model, prediction of conserved elements

# msa_split

# refeature

# exoniphy


#...