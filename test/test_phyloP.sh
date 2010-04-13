# phyloP test commands, to be run with the testPhast.pl program
# commands beginning with @ are executed once, and output is not captured
# by the program.
# all other commands are run in two different versions of phast, and
# output is captured in each, and compared.

@phyloFit hmrc.ss --tree "(human, (mouse,rat), cow)" -i SS --quiet
phyloP --null 10 phyloFit.mod
@msa_view -i SS -o SS --end 100 hmrc.ss > hmrc_short.ss
phyloP -i SS phyloFit.mod hmrc_short.ss
phyloP -i SS --method LRT --base-by-base phyloFit.mod hmrc.ss
phyloP -i SS --method LRT --mode CONACC --wig-scores phyloFit.mod hmrc.ss
phyloP -i SS --method SCORE --mode NNEUT --base-by-base phyloFit.mod hmrc.ss
phyloP -i SS --method GERP --mode ACC --wig-scores phyloFit.mod hmrc.ss
phyloP -i SS --method GERP --base-by-base phyloFit.mod hmrc.ss
phyloP -i SS --method SCORE --wig-scores phyloFit.mod hmrc.ss
phyloP -i SS --method LRT --wig-scores phyloFit.mod hmrc.ss
phyloP -i SS --method LRT --base-by-base phyloFit.mod hmrc.ss
phyloP -i SS --method SCORE --wig-scores --refidx 2 phyloFit.mod hmrc.ss
@echo -e "chr1\t0\t10\nchr1\t50\t100\nchr1\t200\t300" > temp.bed
phyloP -i SS --method LRT --mode CONACC --features temp.bed phyloFit.mod hmrc.ss
phyloP -i SS --method SCORE --features temp.bed -g phyloFit.mod hmrc.ss
@tree_doctor --name-ancestors phyloFit.mod > phyloFit-named.mod
phyloP -i SS --method LRT --mode CONACC --subtree mouse-rat -w phyloFit-named.mod hmrc.ss
phyloP -i SS --method LRT --mode ACC --branch mouse-rat -w phyloFit-named.mod hmrc.ss
phyloP --posterior -i SS phyloFit.mod hmrc_short.ss
phyloP --fit-model -w -i SS --subtree mouse-rat phyloFit-named.mod hmrc_short.ss
phyloP --epsilon 0.0001 -i SS phyloFit-named.mod hmrc_short.ss
phyloP --confidence-interval 0.05 -i SS phyloFit-named.mod hmrc_short.ss
phyloP --quantiles --null 100 phyloFit-named.mod
@rm -f hmrc_short.ss phyloFit.mod phyloFit-named.mod temp.bed
