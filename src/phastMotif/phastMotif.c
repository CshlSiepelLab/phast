/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include <trees.h>
#include <tree_model.h>
#include <msa.h>
#include <motif.h>
#include <ctype.h>
#include <sufficient_stats.h>
#include <bed.h>

#define DEFAULT_SIZE 10
#define DEFAULT_NUMBER 3
#define PRIOR 0.7

void usage(char *prog) {
  printf("\n\
PROGRAM:      %s\n\
\n\
DESCRIPTION:  Predicts motifs from a set of multiple alignments.  Uses\n\
              an EM algorithm similar to that of MEME, but a motif is\n\
              defined by phylogenetic models rather than multinomial\n\
              distributions.  The specified multiple alignments may\n\
              actually be single sequences (see -m).  Various parameters\n\
              control the strategy for initialization (see below).\n\
              Currently, the F81 substitution model is assumed.\n\
\n\
USAGE:        %s [-t <treefile>] [OPTIONS] <msa_list>\n\
\n\
OPTIONS:\n\
    -t <file> (Required unless -m or -p) Use specified tree topology for\n\
              all phylogenetic models (Newick format).\n\
\n\
    -i <fmt>  Input format for alignment.  May be FASTA, PHYLIP, MPM, SS,\n\
              or MAF (default FASTA).\n\
\n\
    -b <file> Read background model from specified file (.mod format).\n\
              By default, the background model is estimated\n\
              in a preprocessing step, by pooling all data.\n\
\n\
    -s        Estimate a separate background model for each multiple alignment.\n\
              (Not yet implemented.)\n\
\n\
    -k <size> Learn motifs of the specified size (default is %d).\n\
\n\
    -B <n>    Report best <n> motifs (default %d).\n\
\n\
    -m        MEME mode.  Use multinomial rather than phylogenetic\n\
              models.  Causes multiple alignments to be ignored -- any\n\
              gaps are discarded and all sequences are assumed\n\
              independent.\n\
\n\
    -d <+lst> Use the discriminative training method of Segal et\n\
              al. (RECOMB'02), rather than EM.  The specified list\n\
              should contain the filenames from msa_list that are to\n\
              be considered *positive* examples (containing the\n\
              desired motif); all others will be considered negative\n\
              examples.  Can be used with or without -m.\n\
\n\
    -p        Use \"profile\" models rather than phylogenetic models\n\
              (characters in each alignment column assumed\n\
              independent).  The resulting model is a hybrid of the\n\
              full model and MEME's model.  Essentially, it uses the\n\
              multiple alignments but not the phylogeny.  NOT YET IMPLEMENTED.\n\
\n\
    -n <n>    Perform <n> random restarts and report the motif with highest\n\
              likelihood.  Default number is 10.  Ignored with -I, -P, and\n\
              -R unless -S is specified (see below).\n\
\n\
    -I <mlst> Run the algorithm after a \"soft\" initialization with\n\
              each of the consensus sequences in the specified list.\n\
              At each position, <pc> pseudocounts (see -c) are given\n\
              to the consensus base and 1 pseudocount to all other\n\
              bases.  Each string must have length at most equal to\n\
              the size of the motif.  If shorter, it is used as a\n\
              \"seed\" for a motif, with flanking positions treated as\n\
              wildcards.\n\
\n\
    -P <x,y>  Initialize with the x most prevalent y-tuples.  A soft\n\
              initialization is performed, as above.  If y is less\n\
              than the motif size, y-tuples are used as a \"seed\" for\n\
              a motif, as above.\n\
\n\
    -R <x,y>  Initialize with a random sample of x y-tuples.  A soft\n\
              initialization is performed, as above.  If y is less\n\
              than the motif size, y-tuples are used as a \"seed\" for\n\
              a motif, as above.\n\
\n\
    -w <n>    (for use with -I, -P, -R) Winnow initialization sequences\n\
              to the top <n> based on the unmaximized likelihood.\n\
\n\
    -c <pc>   (for use with -I, -P, -R) Number of pseudocounts for\n\
              consensus bases (default 5).\n\
\n\
    -S        (for use with -I, -P, -R) Instead of doing a deterministic\n\
              initialization based on a consensus sequence, sample\n\
              parameters from a Dirichlet distribution defined by the\n\
              pseudocounts (see -c).  In this case, random restarts\n\
              are performed, as specified by -n.\n\
\n\
    -o <pref> Use the specified prefix for all output files (dflt. \"phastm\").\n\
    -H        Produce HTML formatted output, in addition to ordinary output.\n\
              One file is produced per predicted motif, as well as a \n\
              single HTML-formatted summary file.\n\
\n\
    -D        Produce a BED file with predicted motifs, for use in the \n\
              UCSC browser.  Currently, sequence names must be\n\
              formatted such as \"chr10:102553847-102554897+\", with\n\
              the final '+' or '-' indicating strand.\n\
\n\
    -x        (For use with -H or -D) Suppress ordinary output to stdout.\n\
\n\
    -h        Print this help message.\n\n", prog, prog, DEFAULT_SIZE, 
         DEFAULT_NUMBER);
  exit(0);
}

int main(int argc, char *argv[]) {
  TreeNode *tree = NULL;
  TreeModel *backgd_mod = NULL;
  int i, j, separate_backgd = 0, 
    size = DEFAULT_SIZE, meme_mode = 0, profile_mode = 0, 
    nrestarts = 10, npseudocounts = 5, nsamples = -1, 
    nmostprevalent = -1, tuple_size = -1, nbest = -1, sample_parms = 0,
    nmotifs = DEFAULT_NUMBER, nseqs = -1, do_html = 0, do_bed = 0, 
    suppress_stdout = 0;
  List *msa_name_list = NULL, *pos_examples = NULL, *init_list = NULL, *tmpl;
  List *msas, *motifs;
  SeqSet *seqset = NULL;
  PooledMSA *pmsa = NULL;
  msa_format_type msa_format = UNKNOWN_FORMAT;
  Vector *backgd_mnmod = NULL;
  Hashtable *hash=NULL;
  String *output_prefix = str_new_charstr("phastm.");
  double *has_motif = NULL;
  double prior = PRIOR;
  char c;
  GFF_Set *bedfeats = NULL;

  while ((c = getopt(argc, argv, "t:i:b:sk:md:pn:I:R:P:w:c:SB:o:HDxh")) != -1) {
    switch (c) {
    case 't':
      tree = tr_new_from_file(phast_fopen(optarg, "r"));
      break;
    case 'i':
      msa_format = msa_str_to_format(optarg);
      if (msa_format == UNKNOWN_FORMAT) 
	die("ERROR: bad input format.\n");
      break;
    case 'b':
      backgd_mod = tm_new_from_file(phast_fopen(optarg, "r"), 1);
      break;
    case 's':
      separate_backgd = 1;
      break;
    case 'k':
      size = get_arg_int(optarg);
      break;
    case 'm':
      meme_mode = 1;
      break;
    case 'd':
      pos_examples = get_arg_list(optarg);
      break;
    case 'p':
      profile_mode = 1;
      break;
    case 'n':
      nrestarts = get_arg_int(optarg);
      break;
    case 'I':
      init_list = get_arg_list(optarg);
      break;
    case 'P':
      tmpl = str_list_as_int(get_arg_list(optarg));
      if (lst_size(tmpl) != 2) die("ERROR: bad argument to -P.\n");
      nmostprevalent = lst_get_int(tmpl, 0);
      tuple_size = lst_get_int(tmpl, 1);
      if (!(nmostprevalent > 0 && tuple_size > 0))
	die("ERROR: bad argument nmostprevalent=%i tuple_size=%i\n", 
	    nmostprevalent, tuple_size);
      lst_free(tmpl);
      break;
    case 'R':
      tmpl = str_list_as_int(get_arg_list(optarg));
      if (lst_size(tmpl) != 2) die("ERROR: bad argument to -R.\n");
      nsamples = lst_get_int(tmpl, 0);
      tuple_size = lst_get_int(tmpl, 1);
      if (!(nsamples > 0 && tuple_size > 0))
	die("ERROR nsamples=%i tuple_sizse=%i\n", nsamples, tuple_size);
      lst_free(tmpl);
      break;
    case 'c':
      npseudocounts = get_arg_int(optarg);
      break;
    case 'w':
      nbest = get_arg_int(optarg);
      break;
    case 'S':
      sample_parms = 1;
      break;
    case 'B':
      nmotifs = get_arg_int(optarg);
      break;
    case 'o': 
      str_free(output_prefix);
      output_prefix = str_new_charstr(optarg);
      str_append_char(output_prefix, '.'); 
      break;
    case 'H': 
      do_html = 1;
      break;
    case 'D': 
      do_bed = 1;
      break;
    case 'x':
      suppress_stdout = 1;
      break;
    case 'h':
      usage(argv[0]);
    case '?':
      die("Bad argument.  Try '%s -h'.\n", argv[0]);
    }
  }

  if (optind != argc - 1) 
    die("ERROR: List of alignment files required.  Try '%s -h'.\n", argv[0]);

  if ((nsamples > 0 && nmostprevalent > 0) || 
      (nsamples > 0 && init_list != NULL) || 
      (nmostprevalent > 0 && init_list != NULL)) 
    die("ERROR: -I, -P, and -R are mutually exclusive.");

  set_seed(-1);
    
  msa_name_list = get_arg_list(argv[optind]);

  if (backgd_mod != NULL && tree == NULL) tree = backgd_mod->tree;

  if (tree == NULL && !meme_mode && !profile_mode) 
    die("ERROR: Must specify -t, -m, or -p.\n");

  if ((init_list != NULL || nsamples > 0 || nmostprevalent > 0) && 
      !sample_parms)
    nrestarts = 1;

  if (pos_examples != NULL) {
    hash = hsh_new(lst_size(pos_examples));
    for (i = 0; i < lst_size(pos_examples); i++)
      hsh_put_int(hash, ((String*)lst_get_ptr(pos_examples, i))->chars, 1);
    has_motif = smalloc(lst_size(msa_name_list) * sizeof(double));
  }

  /* open all MSAs */
  msas = lst_new_ptr(lst_size(msa_name_list));
  fprintf(stderr, "Reading alignment(s) ...\n");
  for (i = 0, j = 0; i < lst_size(msa_name_list); i++) {
    String *name = lst_get_ptr(msa_name_list, i);
    FILE *mfile = phast_fopen(name->chars, "r");
    msa_format_type temp_format;
    MSA *msa;
    if (msa_format == UNKNOWN_FORMAT)
      temp_format = msa_format_for_content(mfile, 1);
    else temp_format = msa_format;
    msa = msa_new_from_file_define_format(mfile, temp_format, NULL);
    phast_fclose(mfile);
    if (nseqs == -1) nseqs = msa->nseqs;
    if (!meme_mode &&
        (msa->length - msa_num_gapped_cols(msa, STRIP_ANY_GAPS, -1, -1) < 300 ||
        msa->nseqs != nseqs)) {
      fprintf(stderr, "WARNING: ignoring alignment '%s' -- too few informative sites.\n", name->chars);
      msa_free(msa);
      continue;
    }

    if (msa_alph_has_lowercase(msa)) msa_toupper(msa); 
    msa_remove_N_from_alph(msa); /* Ns can be a problem */
    lst_push_ptr(msas, msa);
    if (has_motif != NULL) {
      int k, hm = (hsh_get_int(hash, name->chars) == 1);
      if (meme_mode) {          /* here need to record at individ seq level */
        has_motif = srealloc(has_motif, 
                             (j + msa->nseqs + 1) * sizeof(double)); /* FIXME */
        for (k = 0; k < msa->nseqs; k++) has_motif[j++] = hm;
      }
      else has_motif[j++] = hm;
    }
  }
  if (!meme_mode) {
    fprintf(stderr, "Extracting and pooling sufficient statistics ...\n");
    pmsa = ss_pooled_from_msas(msas, 1, size, NULL, 0);
    msa_remove_N_from_alph(pmsa->pooled_msa);
  }

  /* obtain individual sequences, if necessary */
  if (nmostprevalent > 0 || nsamples > 0 || meme_mode) {
    if (meme_mode) fprintf(stderr, "Converting to individual sequences ...\n");
    else fprintf(stderr, "Obtaining reference sequences for pre-processing ...\n");
    seqset = mtf_get_seqset(msas, meme_mode ? -1 : 1, 10 * size);
                                /* for now, assume 1st seq is reference */
    msa_remove_N_from_alph(seqset->set); 
  }

  if (nmostprevalent > 0) {
    fprintf(stderr, "Obtaining %d most prevalent %d-tuples ...\n", 
            nmostprevalent, tuple_size);
    init_list = lst_new_ptr(nmostprevalent);
    mtf_get_common_ntuples(seqset, init_list, tuple_size, nmostprevalent);
  }
  else if (nsamples > 0) {
    fprintf(stderr, "Sampling %d %d-tuples ...\n", nsamples, tuple_size);
    init_list = lst_new_ptr(nsamples);
    mtf_sample_ntuples(seqset, init_list, tuple_size, nsamples);
  }

  /* in meme_mode, backgd model can be specified as eq freqs in a .mod file */
  if (meme_mode && backgd_mod != NULL && has_motif == NULL)
    backgd_mnmod = backgd_mod->backgd_freqs;

  /* estimate background model, if necessary */
  else if (backgd_mod == NULL && (!meme_mode || has_motif == NULL)) {
    fprintf(stderr, "Fitting background model%s ...\n", 
            has_motif == NULL ? "" : " (for use in initialization)");
                                /* if discriminative, be clear
                                   backgd isn't really part of the
                                   estimation procedure */
    if (meme_mode) {
      backgd_mnmod = vec_new(strlen(seqset->set->alphabet));
      mtf_estim_backgd_mn(seqset, backgd_mnmod);
    }
    else {
      backgd_mod = tm_new(tr_create_copy(tree), NULL, NULL, F81, 
                          pmsa->pooled_msa->alphabet, 1, 0, NULL, -1);
      tm_fit(backgd_mod, pmsa->pooled_msa, 
             tm_params_init(backgd_mod, .1, 5, 0), 
             -1, OPT_MED_PREC, NULL, 0, NULL);
    }
  }

  /* select subset of init strings, if necessary */
  if (nbest > 0 && init_list != NULL) {
    fprintf(stderr, "Winnowing candidate start strings ...\n");
    tmpl = lst_new_ptr(nbest);
    mtf_winnow_starts(meme_mode ? (void*)seqset : (void*)pmsa,
                      init_list, nbest, tmpl, !meme_mode, size, tree,
                      meme_mode ? (void*)backgd_mnmod : (void*)backgd_mod, 
                      has_motif);
    lst_free(init_list);
    init_list = tmpl;
  }

  /* Now find motifs */
  motifs = mtf_find(meme_mode ? (void*)seqset : (void*)pmsa, 
                    !meme_mode, size, nmotifs, tree,
                    meme_mode ? (void*)backgd_mnmod : (void*)backgd_mod, 
                    has_motif, prior, nrestarts, init_list, sample_parms, 
                    npseudocounts);
     
  fprintf(stderr, "\n\n");
  if (do_bed)
    bedfeats = gff_new_set_init("phast_motif", "0.1b");

  /* generate output */
  for (i = 0; i < lst_size(motifs); i++) {
    Motif *m = lst_get_ptr(motifs, i);

    if (!suppress_stdout) {
      if (lst_size(motifs) > 1) 
        printf("\n**********\nMOTIF #%d\n**********\n\n", i+1);

      mtf_print(stdout, m);
    }

    if (do_html) {
      String *fname = str_dup(output_prefix);
      str_append_int(fname, i+1);
      str_append_charstr(fname, ".html");
      mtf_print_html(phast_fopen(fname->chars, "w+"), m);
      str_free(fname);
    }

    if (do_bed) 
      mtf_add_features(m, bedfeats);
  }
  if (do_html) {
    String *fname = str_dup(output_prefix);
    str_append_charstr(fname, "index.html");
    mtf_print_summary_html(phast_fopen(fname->chars, "w+"), 
                           motifs, output_prefix);
    str_free(fname);
  }
  if (do_bed) {
    String *fname = str_dup(output_prefix);
    str_append_charstr(fname, "bed");
    gff_print_bed(phast_fopen(fname->chars, "w+"),
                  bedfeats, FALSE);
    str_free(fname);
  }

  return 0;
}
