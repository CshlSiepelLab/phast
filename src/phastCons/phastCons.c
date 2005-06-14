#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include <phylo_hmm.h>
#include <em.h>
#include <sufficient_stats.h>
#include <bed.h>
#include <dgamma.h>
#include <tree_likelihoods.h>
#include "phastCons.help"

#define DEFAULT_RHO 0.3

/* functions implemented below and used internally */
void setup_two_state(HMM **hmm, CategoryMap **cm, double mu, double nu);

double fit_two_state(PhyloHmm *phmm, MSA *msa, int estim_func, int estim_indels,
                     int estim_trees, int estim_rho, double *mu, double *nu, 
                     double *alpha_0, double *beta_0, double *tau_0, 
                     double *alpha_1, double *beta_1, double *tau_1, 
                     double *rho, double gamma, FILE *logf);

void reestimate_trees(void **models, int nmodels, void *data, 
                      double **E, int nobs, FILE *logf);

void reestimate_rho(void **models, int nmodels, void *data, 
		    double **E, int nobs, FILE *logf);

void phmm_estim_trans_em_coverage(HMM *hmm, void *data, double **A);

void collapse_cats(CategoryMap *cm, List *cats_to_merge);

void init_eqfreqs(TreeModel *mod, MSA *msa, double gc);


int main(int argc, char *argv[]) {

  /* arguments and defaults */
  int post_probs = TRUE, score = FALSE, quiet = FALSE, 
    gff = FALSE, FC = FALSE, estim_lambda = TRUE, 
    estim_transitions = TRUE, two_state = TRUE, indels = FALSE,
    coding_potential = FALSE, indels_only = FALSE, estim_indels = TRUE,
    estim_trees = FALSE, ignore_missing = FALSE, estim_rho = FALSE,
    set_transitions = FALSE;
  int nrates = -1, nrates2 = -1, refidx = 1, max_micro_indel = 20;
  double lambda = 0.9, mu = 0.01, nu = 0.01, alpha_0 = 0.05, beta_0 = 0.05, 
    tau_0 = 0.45, alpha_1 = 0.05, beta_1 = 0.05, tau_1 = 0.2, gc = -1,
    gamma = -1, rho = DEFAULT_RHO, omega = -1;
  msa_format_type msa_format = SS;
  FILE *viterbi_f = NULL, *lnl_f = NULL, *log_f = NULL;
  List *states = NULL, *pivot_states = NULL, *inform_reqd = NULL, 
    *mod_fname_list = NULL, *not_informative = NULL;
  char *seqname = NULL, *idpref = NULL, *estim_trees_fname_root = NULL,
    *extrapolate_tree_fname = NULL;
  HMM *hmm = NULL;
  Hashtable *alias_hash = NULL;
  TreeNode *extrapolate_tree = NULL;

  struct option long_opts[] = {
    {"states", 1, 0, 'S'},
    {"hmm", 1, 0, 'H'},
    {"viterbi", 1, 0, 'V'},
    {"most-conserved", 1, 0, 'V'}, /* same as --viterbi */
    {"no-post-probs", 0, 0, 'n'},
    {"msa-format", 1, 0, 'i'},
    {"FC", 0, 0, 'X'},
    {"lambda", 1, 0, 'l'},
    {"target-coverage", 1, 0, 'C'},
    {"transitions", 1, 0, 't'},
    {"expected-length", 1, 0, 'E'},
    {"expected-lengths", 1, 0, 'E'}, /* for backward compatibility */
    {"estimate-trees", 1, 0, 'T'},
    {"estimate-rho", 1, 0, 'O'},
    {"rho", 1, 0, 'R'},
    {"gc", 1, 0, 'G'},
    {"ignore-missing", 0, 0, 'z'},
    {"nrates", 1, 0, 'k'},
    {"log", 1, 0, 'g'},
    {"refidx", 1, 0, 'r'},
    {"suppress-missing", 0, 0, 'x'}, /* for backward compatibility */
    {"reflect-strand", 1, 0, 'U'},
    {"catmap", 1, 0, 'c'},
    {"extrapolate", 1, 0, 'e'},
    {"indels", 0, 0, 'I'},
    {"max-micro-indel", 1, 0, 'Y'},
    {"indel-params", 1, 0, 'D'},
    {"min-informative-types", 1, 0, 'M'}, /* for backward compatibility */
    {"require-informative", 1, 0, 'M'},
    {"not-informative", 1, 0, 'F'},
    {"lnl", 1, 0, 'L'},
    {"seqname", 1, 0, 'N'},
    {"idpref", 1, 0, 'P'},
    {"score", 0, 0, 's'},
    {"coding-potential", 0, 0, 'p'},
    {"indels-only", 0, 0, 'J'},
    {"alias", 1, 0, 'A'},
    {"quiet", 0, 0, 'q'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  /* other vars */
  char c;
  int opt_idx, i, j, last;
  List *tmpl = NULL;
  MSA *msa = NULL;
  double lnl = INFTY;
  String *tmpstr;
  TreeModel **mod;
  PhyloHmm *phmm;
  CategoryMap *cm = NULL;
  char *mods_fname = NULL, *newname;
  indel_mode_type indel_mode;

  while ((c = getopt_long(argc, argv, 
			  "S:H:V:ni:k:l:C:G:zt:E:R:T:O:r:xL:s:N:P:g:U:c:e:IY:D:JM:F:pA:Xqh", 
                          long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'S':
      states = get_arg_list(optarg);
      break;
    case 'H':
      hmm = hmm_new_from_file(fopen_fname(optarg, "r"));
      two_state = FALSE;
      break;
    case 'V':
      viterbi_f = fopen_fname(optarg, "w+");
      tmpstr = str_new_charstr(optarg);
      if (str_ends_with_charstr(tmpstr, ".gff")) gff = TRUE;
      str_free(tmpstr);
      break;
    case 'n':
      post_probs = FALSE;
      break;
    case 'i':
      msa_format = msa_str_to_format(optarg);
      if (msa_format == -1) die("ERROR: bad argument to --msa-format\n");
      break;
    case 'X':
      FC = TRUE;
      two_state = FALSE;
      break;
    case 'l':
      if (optarg[0] != '~') estim_lambda = FALSE;
      else optarg = &optarg[1];
      lambda = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'C':
      gamma = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'G':
      gc = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 't':
      set_transitions = TRUE;
      if (optarg[0] != '~') estim_transitions = FALSE;
      else optarg = &optarg[1];
      tmpl = get_arg_list_dbl(optarg);
      if (lst_size(tmpl) != 2) 
        die("ERROR: bad argument to --transitions.\n");
      mu = lst_get_dbl(tmpl, 0);
      nu = lst_get_dbl(tmpl, 1);
      if (mu <= 0 || mu >= 1 || nu <= 0 || nu >= 1)
        die("ERROR: bad argument to --transitions.\n");
      lst_free(tmpl);
      break;
    case 'E':
      if (optarg[0] != '~') estim_transitions = FALSE;
      else optarg = &optarg[1];
      omega = get_arg_dbl_bounds(optarg, 1, INFTY);
      mu = 1/omega;
      break;
    case 'T':
      estim_trees = TRUE;
      estim_trees_fname_root = optarg;
      break;
    case 'O':
      estim_rho = TRUE;
      estim_trees_fname_root = optarg;
      break;
    case 'z':
      ignore_missing = TRUE;
      break;
    case 'k':
      tmpl = get_arg_list_int(optarg);
      if (lst_size(tmpl) > 2) 
        die("ERROR: too many arguments with --nrates.\n");
      nrates = lst_get_int(tmpl, 0);
      if (nrates <= 0) 
        die("ERROR: bad argument to --nrates (%d).\n", nrates);
      if (lst_size(tmpl) == 2) {
        nrates2 = lst_get_int(tmpl, 1);
        if (nrates2 <= 0) 
          die("ERROR: bad argument to --nrates (%d).\n", nrates2);
      }
      lst_free(tmpl);
      break;
    case 'R':
      rho = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'g':
      if (!strcmp(optarg, "-")) log_f = stderr;
      else log_f = fopen_fname(optarg, "w+");
      break;
    case 'r':
      refidx = get_arg_int_bounds(optarg, 0, INFTY);
      break;
    case 'x':
      /* do nothing; left in for backward compatibility */
      break;
    case 'U':
      pivot_states = get_arg_list(optarg); /* we want strings not ints
                                             for phmm_new */
      break;
    case 'e':
      extrapolate_tree_fname = optarg;
      break;
    case 'I':
      indels = TRUE;
      break;
    case 'Y':
      max_micro_indel = get_arg_int_bounds(optarg, 1, INFTY);
      break;
    case 'D':
      if (optarg[0] != '~') estim_indels = FALSE;
      else optarg = &optarg[1];
      tmpl = get_arg_list_dbl(optarg);
      if (lst_size(tmpl) != 6) die("ERROR: bad argument to --indel-params.\n");
      alpha_0 = lst_get_dbl(tmpl, 0);
      beta_0 = lst_get_dbl(tmpl, 1);
      tau_0 = lst_get_dbl(tmpl, 2);
      alpha_1 = lst_get_dbl(tmpl, 3);
      beta_1 = lst_get_dbl(tmpl, 4);
      tau_1 = lst_get_dbl(tmpl, 5);
      if (alpha_0 < 0 || beta_0 < 0 || tau_0 < 0 || 
          alpha_1 < 0 || beta_1 < 0 || tau_1 < 0)
        die("ERROR: bad argument to --indel-params.\n");
      lst_free(tmpl);
      break;
    case 'J':
      indels_only = TRUE;
      two_state = FALSE;
      indels = TRUE;
      post_probs = FALSE;
      break;
    case 'M':
      inform_reqd = get_arg_list(optarg);
      break;
    case 'F':
      not_informative = get_arg_list(optarg);
      break;
    case 'c':
      cm = cm_new_string_or_file(optarg);
      break;
    case 'L':
      lnl_f = fopen_fname(optarg, "w+");
      break;
    case 'N':
      seqname = optarg;
      break;
    case 'P':
      idpref = optarg;
      break;
    case 's':
      score = TRUE;
      break;
    case 'p':
      coding_potential = TRUE;
      break;
    case 'A':
      alias_hash = make_name_hash(optarg);
      break;
    case 'q':
      quiet = TRUE;
      break;
    case 'h':
      printf(HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try '%s -h'.\n", argv[0]);
    }
  }

  /* enforce usage rules */
  if ((hmm != NULL && FC))
    die("ERROR: --hmm and --FC are mutually exclusive.\n");

  if (indels_only && (hmm != NULL || FC))
    die("ERROR: --indels-only cannot be used with --hmm or --FC.\n");

  if ((estim_trees || gamma != -1 || estim_rho || 
       omega != -1 || set_transitions) && !two_state)
    die("ERROR: --estimate-trees, --target-coverage, --expected-length, --transitions,\nand --estimate-rho can only be used with default two-state HMM.\n");

  if (set_transitions && (gamma != -1 || omega != -1))
    die("ERROR: --transitions and --target-coverage/--expected-length cannot be used together.\n");

  if (omega != -1 && gamma == -1) 
    die("ERROR: --expected-length requires --target-coverage.\n");

  if (cm != NULL && hmm == NULL) 
    die("ERROR: --catmap can only be used with --hmm.\n");
  
  if (indels == TRUE && FC)
    die("ERROR: --indels cannot be used with --FC.\n");

  if (nrates != -1 && hmm != NULL)
    die("ERROR: --nrates currently can't be used with --hmm.\n");

  if ((!coding_potential && optind != argc - 2) ||
      (coding_potential && optind != argc - 2 && optind != argc - 1))
    die("ERROR: extra or missing arguments.  Try '%s -h'.\n", argv[0]);

  if (!indels) estim_indels = FALSE;

  if (extrapolate_tree_fname != NULL &&
      !strcmp(extrapolate_tree_fname, "default")) {
    extrapolate_tree_fname = smalloc(1000 * sizeof(char));
    sprintf(extrapolate_tree_fname, 
            "%s/data/exoniphy/mammals/cftr25_hybrid.nh", PHAST_HOME);
  }
  if (extrapolate_tree_fname != NULL)
    extrapolate_tree = tr_new_from_file(fopen_fname(extrapolate_tree_fname, "r"));

  mods_fname = (optind == argc - 2 ? argv[argc - 1] : NULL);
  /* if there are two args, mods are the second one; otherwise will
     use default mods for coding potential (see below) */
  
  /* set defaults for coding-potential mode */
  if (coding_potential) {
    char tmp[5000];
    two_state = FALSE;
    if (cm == NULL) cm = cm_new_string_or_file("NCATS=4; CNS 1; CDS 2-4");
    if (hmm == NULL) {
      sprintf(tmp, "%s/data/phastCons/%s", PHAST_HOME,
              indels ? "simple-coding-indels.hmm" : "simple-coding.hmm");
      if (!quiet) fprintf(stderr, "Reading HMM from %s...\n", tmp);
      hmm = hmm_new_from_file(fopen_fname(tmp, "r"));
    }
    if (mods_fname == NULL) {
      sprintf(tmp, "\
%s/data/exoniphy/mammals/r3.ncns.mod,\
%s/data/exoniphy/mammals/r3.cns.mod,\
%s/data/exoniphy/mammals/r3.cds-1.mod,\
%s/data/exoniphy/mammals/r3.cds-2.mod,\
%s/data/exoniphy/mammals/r3.cds-3.mod", 
              PHAST_HOME, PHAST_HOME, PHAST_HOME, PHAST_HOME, PHAST_HOME);
      mods_fname = tmp;
    }
    if (states == NULL) states = get_arg_list("CDS");
    if (pivot_states == NULL) pivot_states = get_arg_list("background,CNS");
    if (inform_reqd == NULL) inform_reqd = get_arg_list("CDS");
  }
  
  /* read alignment */
  if (!quiet)
    fprintf(stderr, "Reading alignment from %s...\n", argv[optind]);
  msa = msa_new_from_file(fopen_fname(argv[optind], "r"), msa_format, NULL);
  msa_remove_N_from_alph(msa);  /* for backward compatibility */
  if (msa_format == SS && msa->ss->tuple_idx == NULL) 
    die("ERROR: Ordered representation of alignment required.\n");
  if (msa->ss == NULL) ss_from_msas(msa, 1, TRUE, NULL, NULL, NULL, -1);
                                /* SS assumed below */

  /* rename if aliases are defined */
  if (alias_hash != NULL) {
    for (i = 0; i < msa->nseqs; i++) {
      if ((newname = hsh_get(alias_hash, msa->names[i])) != (char*)-1) {
        free(msa->names[i]);
        msa->names[i] = strdup(newname);
      }
    }
  }

  /* mask out macro-indels, if necessary */
  if (indels) {
    /* this little hack allows gaps in refseq to be restored before
       output (needed for proper coord conversion) */
    if (msa->seqs == NULL) { ss_to_msa(msa); ss_free(msa->ss); msa->ss = NULL; }
    assert(strlen(msa->missing) >= 2);
    for (i = 0; i < msa->length; i++) 
      if (msa->is_missing[(int)msa->seqs[0][i]]) msa->seqs[0][i] = msa->missing[1];
                                /* msa->missing[0] is used in msa_mask_macro_indels */

    msa_mask_macro_indels(msa, max_micro_indel);
  }

  /* Set up array indicating which seqs are informative, if necessary */
  if (not_informative != NULL)
    msa_set_informative(msa, not_informative);

  /* strip missing columns, if necessary */
  if (ignore_missing)
    ss_strip_missing(msa, refidx);

  /* read tree models */
  mod_fname_list = get_arg_list(mods_fname);

  if ((FC || indels_only) && lst_size(mod_fname_list) != 1)
    die("ERROR: only one tree model allowed with --FC and --indels-only.\n");

  if (two_state && lst_size(mod_fname_list) > 2)
    die("ERROR: must specify either one or two tree models with default two-state model.\n");
    
  mod = (TreeModel**)smalloc(sizeof(TreeModel*) * lst_size(mod_fname_list));
  for (i = 0; i < lst_size(mod_fname_list); i++) {
    String *fname = lst_get_ptr(mod_fname_list, i);
    int old_nnodes, found;
    List *pruned_names = lst_new_ptr(msa->nseqs);

    if (!quiet)
      fprintf(stderr, "Reading tree model from %s...\n", fname->chars);
    mod[i] = tm_new_from_file(fopen_fname(fname->chars, "r"));
    mod[i]->use_conditionals = 1;     
    old_nnodes = mod[i]->tree->nnodes;

    /* extrapolate tree and/or prune away extra species */
    if (extrapolate_tree != NULL) {
      double scale = tm_extrapolate_and_prune(mod[i], extrapolate_tree, 
                                              msa, pruned_names);
      if (!quiet) 
        fprintf(stderr, "Extrapolating based on %s (scale=%f)...\n", 
                extrapolate_tree_fname, scale);
    }
    else
      tm_prune(mod[i], msa, pruned_names);

    if (lst_size(pruned_names) == (old_nnodes + 1) / 2)
      die("ERROR: no match for leaves of tree in alignment (leaf names must match alignment names).\n");
    if (!quiet && lst_size(pruned_names) > 0) {
      fprintf(stderr, "WARNING: pruned away leaves of tree with no match in alignment (");
      for (j = 0; j < lst_size(pruned_names); j++)
        fprintf(stderr, "%s%s", ((String*)lst_get_ptr(pruned_names, j))->chars, 
                j < lst_size(pruned_names) - 1 ? ", " : ").\n");
    }

    /* also make sure match for reference sequence in tree */
    if (refidx > 0) {
      for (j = 0, found = FALSE; !found && j < mod[i]->tree->nnodes; j++) {
	TreeNode *n = lst_get_ptr(mod[i]->tree->nodes, j);
	if (n->lchild == NULL && n->rchild == NULL && 
	    n->name != NULL && !strcmp(n->name, msa->names[refidx-1]))
	  found = TRUE;
      }
      if (!found) die("ERROR: no match for reference sequence in tree.\n");
    }

    lst_free_strings(pruned_names);
    lst_free(pruned_names);
  }

  /* initial checks and setup of tree models for two-state and FC */
  if (two_state) {
    if (lst_size(mod_fname_list) == 2 && (estim_trees || estim_rho))
      die("ERROR: If re-estimating tree models, pass in only one model for initialization.\n");
    if (mod[0]->empirical_rates || 
        (lst_size(mod_fname_list) == 2 && mod[1]->empirical_rates))
      die("ERROR: nonparameteric rate variation not allowed with default two-state HMM.\n");

    /* set equilibrium frequencies if estimating tree models */
    if (estim_trees || (gc != -1 && estim_rho)) 
      init_eqfreqs(mod[0], msa, gc);

    if (lst_size(mod_fname_list) == 1) { /* create 2nd tree model &
                                            rescale first */
      mod = srealloc(mod, 2 * sizeof(void*));
      mod[1] = tm_create_copy(mod[0]);
      if (!estim_rho) tm_scale(mod[0], rho, TRUE);
    }
    if (nrates != -1 && nrates != mod[0]->nratecats) 
      tm_reinit(mod[0], mod[0]->subst_mod, nrates, mod[0]->alpha, NULL, NULL);
    if (nrates2 != -1 && nrates2 != mod[1]->nratecats) 
      tm_reinit(mod[1], mod[1]->subst_mod, nrates2, mod[1]->alpha, NULL, NULL);
  }
  else if (FC) {
    if (mod[0]->nratecats <= 1)
      die("ERROR: a tree model allowing for rate variation is required.\n");
    if (nrates != -1 && mod[0]->empirical_rates)
      die("ERROR: can't use --nrates with nonparameteric rate model.\n");
    if (nrates == -1) nrates = mod[0]->nratecats;
  }

  /* use file name root for default seqname */
  if (viterbi_f != NULL && (seqname == NULL || idpref == NULL)) {
    String *tmp = str_new_charstr(argv[optind]);
    if (!str_equals_charstr(tmp, "-")) {
      str_remove_path(tmp);
      str_root(tmp, '.');
      if (idpref == NULL) idpref = strdup(tmp->chars);
      str_root(tmp, '.');         /* apply one more time for double suffix */
      if (seqname == NULL) seqname = tmp->chars;    
    }
    else if (seqname == NULL) seqname = "refseq";
  }

  /* set up states */
  if (states == NULL) {
    states = lst_new_ptr(1);
    lst_push_ptr(states, str_new_charstr("0"));
  }

  if (two_state) {
    if (!quiet) 
      fprintf(stderr, "Creating 'conserved' and 'nonconserved' states in HMM...\n");
    if (gamma != -1) {
      nu = gamma/(1-gamma) * mu;
      if (nu >= 1) 
        die("ERROR: mu=%f and gamma=%f imply nu >= 1.\n", mu, gamma);
    }
    setup_two_state(&hmm, &cm, mu, nu);
  }
  else if (cm == NULL)
    cm = cm_create_trivial(lst_size(mod_fname_list)-1, NULL);

  /* set up PhyloHmm */
  if (!indels) indel_mode = MISSING_DATA;
  else if (hmm == NULL || hmm->nstates == cm->ncats + 1)
    indel_mode = PARAMETERIC;
  else indel_mode = NONPARAMETERIC;

  phmm = phmm_new(hmm, mod, cm, pivot_states, indel_mode);

  if (FC) {
    if (!quiet) 
      fprintf(stderr, "Creating %d scaled versions of tree model...\n", nrates);
    phmm_rates_cross(phmm, nrates, lambda, TRUE);
  }

  /* set inform_reqd, if necessary.  This has to be done
     *after* the set of models is expanded (two-state or FC) */
  if (inform_reqd != NULL) {
    List *l = cm_get_category_list(cm, inform_reqd, 0);
    for (i = 0; i < lst_size(l); i++) {
      int modno = lst_get_int(l, i);
      if (modno < 0 || modno >= phmm->nmods) 
        die("ERROR: illegal argument to --require-informative.\n");
      phmm->mods[modno]->inform_reqd = TRUE;
    }
    lst_free(l);
  }        

  /* compute emissions */
  phmm_compute_emissions(phmm, msa, quiet);

  /* estimate lambda, if necessary */
  if (FC && estim_lambda) {
    if (!quiet) fprintf(stderr, "Finding MLE for lambda...");
    lnl = phmm_fit_lambda(phmm, &lambda, log_f);
    if (!quiet) fprintf(stderr, " (lambda = %f)\n", lambda);
    phmm_update_cross_prod(phmm, lambda);
  }

  /* estimate mu and nu and indel params, if necessary */
  else if (two_state && 
	   (estim_transitions || estim_indels || estim_trees || estim_rho)) {
    char cons_fname[STR_MED_LEN], noncons_fname[STR_MED_LEN];
    if (!quiet) {
      fprintf(stderr, "Finding MLE for (");
      if (estim_transitions) 
        fprintf(stderr, "mu, nu%s", estim_indels || estim_trees || estim_rho 
		? ", " : "");
      if (estim_indels) 
        fprintf(stderr, "alpha_0, beta_0, tau_0, alpha_1, beta_1, tau_1%s",
                estim_trees || estim_rho ? ", " : "");
      if (estim_trees)
        fprintf(stderr, "[tree models]");
      else if (estim_rho) 
        fprintf(stderr, "rho");
      fprintf(stderr, ")...\n");
    }
    lnl = fit_two_state(phmm, msa, estim_transitions, estim_indels, 
			estim_trees, estim_rho,
                        &mu, &nu, &alpha_0, &beta_0, &tau_0, 
                        &alpha_1, &beta_1, &tau_1, &rho,
                        gamma, log_f);
    if (!quiet && (estim_transitions || estim_indels || estim_rho)) {      
      fprintf(stderr, "(");
      if (estim_transitions)
        fprintf(stderr, "mu = %f. nu = %f%s", mu, nu, 
		estim_indels || estim_rho ? ", " : "");
      if (estim_indels)
        fprintf(stderr, 
		"alpha_0 = %f, beta_0 = %f, tau_0 = %f, alpha_1 = %f, beta_1 = %f, tau_1 = %f%s", 
		alpha_0, beta_0, tau_0, alpha_1, beta_1, tau_1, 
		estim_rho ? ", " : "");
      if (estim_rho) 
        fprintf(stderr, "rho = %f", rho);
      fprintf(stderr, ")\n");
    }

    if (estim_trees || estim_rho) {
      sprintf(cons_fname, "%s.cons.mod", estim_trees_fname_root);
      sprintf(noncons_fname, "%s.noncons.mod", estim_trees_fname_root);
      if (!quiet)
        fprintf(stderr, "Writing re-estimated tree models to %s and %s...\n", 
                cons_fname, noncons_fname);
      tm_print(fopen_fname(cons_fname, "w+"), phmm->mods[0]);
      tm_print(fopen_fname(noncons_fname, "w+"), phmm->mods[1]);
    }
  }

  /* estimate indel parameters only, if necessary */
  else if (indels_only) {
    if (!quiet) fprintf(stderr, "Estimating parameters for indel model...");
    lnl = phmm_fit_em(phmm, msa, TRUE, FALSE, log_f);
    if (!quiet) fprintf(stderr, "...\n");
  }

  /* still have to set indel params if not estimating */
  else if (indel_mode == PARAMETERIC) {
    phmm->alpha[0] = alpha_0; phmm->beta[0] = beta_0; phmm->tau[0] = tau_0;
    phmm->alpha[1] = alpha_1; phmm->beta[1] = beta_1; phmm->tau[1] = tau_1;
    phmm_reset(phmm);
  }

  /* before output, have to restore gaps in reference sequence, for
     proper coord conversion */
  if (indels && (post_probs || viterbi_f != NULL)) {
    ss_free(msa->ss); msa->ss = NULL; /* msa->seqs must already exist */
    for (i = 0; i < msa->length; i++) 
      if (msa->seqs[0][i] == msa->missing[0]) msa->seqs[0][i] = GAP_CHAR;
  }
    
  /* Viterbi */
  if (viterbi_f != NULL) {
    GFF_Set *predictions;

    if (lst_size(states) > 1) 
      /* if possible, we want to merge categories of specified states,
         so that the "union" of states is automatically considered by
         phmm_viterbi_features.  Reduces potential for generation of
         huge numbers of features (could be as many as one per site) */
      collapse_cats(phmm->cm, states);

    if (!quiet) fprintf(stderr, "Running Viterbi algorithm...\n");
    predictions = phmm_predict_viterbi_cats(phmm, states, seqname, NULL,
                                            idpref, NULL, "phastCons_predicted");
    /* note that selected state numbers are also cat numbers  */
   
    /* score predictions, if necessary */
    if (score) { 
      if (!quiet) fprintf(stderr, "Scoring predictions...\n");            
      phmm_score_predictions(phmm, predictions, states, NULL, NULL, FALSE);
    }

    /* convert GFF to coord frame of reference sequence and adjust
       coords by idx_offset, if necessary  */
    if (refidx != 0 || msa->idx_offset != 0)
      msa_map_gff_coords(msa, predictions, 0, refidx, msa->idx_offset, NULL);

    if (refidx != 0) 
      gff_flatten(predictions);	
    /* necessary because coord conversion might create overlapping
       features (can happen in deletions in reference sequence) */

    /* now output predictions */
    if (gff)
      gff_print_set(viterbi_f, predictions);
    else                        /* BED format */
      gff_print_bed(viterbi_f, predictions, FALSE); 
  }

  /* posterior probs */
  if (post_probs) {
    int j, k;
    double *postprobs;

    if (!quiet) fprintf(stderr, "Computing posterior probabilities...\n");

    postprobs = phmm_postprobs_cats(phmm, states, &lnl);

    /* print to stdout */
    last = -INFTY;
    for (j = 0, k = 0; j < msa->length; j++) {
      if (refidx == 0 || msa_get_char(msa, refidx-1, j) != GAP_CHAR) {
        if (!msa_missing_col(msa, refidx, j)) {
	  if (k > last + 1) 
	    printf("fixedStep chrom=%s start=%d step=1\n", seqname, 
		   k + msa->idx_offset + 1);
          printf("%.3f\n", postprobs[j]);
	  last = k;
	}
        k++;
      }
    }
  }

  /* likelihood */
  if (lnl_f != NULL) {
    if (lnl > 0) {              /* may have already been computed */
      if (!quiet) fprintf(stderr, "Computing total log likelihood...\n");
      lnl = phmm_lnl(phmm); 
    }
    fprintf(lnl_f, "lnL = %.4f\n", lnl); 
    if (FC) fprintf(lnl_f, "(lambda = %f)\n", lambda);
    else if (two_state && (estim_transitions || estim_indels)) {
      fprintf(lnl_f, "(");
      if (estim_transitions)
        fprintf(lnl_f, "mu = %f, nu = %f%s", mu, nu, estim_indels ? ", " : "");
      if (estim_indels)
        fprintf(lnl_f, "alpha_0 = %f, beta_0 = %f, tau_0 = %f, alpha_1 = %f, beta_1 = %f, tau_1 = %f", alpha_0, beta_0, tau_0, alpha_1, beta_1, tau_1);
      fprintf(lnl_f, ")\n");
    }
  }

  if (!quiet)
    fprintf(stderr, "Done.\n");

  return 0;
}

/* Set up HMM and category map for two-state case */
void setup_two_state(HMM **hmm, CategoryMap **cm, double mu, double nu) {

  *hmm = hmm_new_nstates(2, TRUE, FALSE);

  /* set HMM transitions according to mu and nu */
  mm_set((*hmm)->transition_matrix, 0, 0, 1-mu);
  mm_set((*hmm)->transition_matrix, 0, 1, mu);
  mm_set((*hmm)->transition_matrix, 1, 0, nu);
  mm_set((*hmm)->transition_matrix, 1, 1, 1-nu);

  /* just use stationary distribution for begin transitions */
  gsl_vector_set((*hmm)->begin_transitions, 0, nu/(mu+nu));
  gsl_vector_set((*hmm)->begin_transitions, 1, mu/(mu+nu));

  hmm_reset(*hmm);

  /* define two-category category map */
  *cm = cm_create_trivial(1, "cons_");
}

/* Version of compute_emissions for use when estimating rho only (see
   fit_two_state, below); makes use of fact that emissions for
   nonconserved state need not be recomputed */
void compute_emissions_estim_rho(double **emissions, void **models, 
				 int nmodels, void *data, int sample, 
				 int length) {
  PhyloHmm *phmm = (PhyloHmm*)data;
  tl_compute_log_likelihood(phmm->mods[0], phmm->em_data->msa, 
			    phmm->emissions[0], -1, NULL);
}


/* Estimate parameters for the two-state model using an EM algorithm.
   Any or all of the parameters 'mu' and 'nu', the indel parameters, and
   the tree models themselves may be estimated.  Returns ln
   likelihood. */
double fit_two_state(PhyloHmm *phmm, MSA *msa, int estim_func, int estim_indels,
                     int estim_trees, int estim_rho, double *mu, double *nu, 
                     double *alpha_0, double *beta_0, double *tau_0, 
                     double *alpha_1, double *beta_1, double *tau_1, 
                     double *rho, double gamma, FILE *logf) {
  double retval;

  mm_set(phmm->functional_hmm->transition_matrix, 0, 0, 1-*mu);
  mm_set(phmm->functional_hmm->transition_matrix, 0, 1, *mu);
  mm_set(phmm->functional_hmm->transition_matrix, 1, 0, *nu);
  mm_set(phmm->functional_hmm->transition_matrix, 1, 1, 1-*nu);
                                /* note that phmm->functional_hmm ==
                                   phmm->hmm if no indel model */

  phmm->em_data = smalloc(sizeof(EmData));
  phmm->em_data->msa = msa;
  phmm->em_data->fix_functional = !estim_func;
  phmm->em_data->fix_indel = !estim_indels;
  phmm->em_data->rho = *rho;
  phmm->em_data->gamma = gamma;
  phmm->em_data->H = NULL;      /* will be defined as needed */

  if (phmm->indel_mode == PARAMETERIC) {
    phmm->alpha[0] = *alpha_0;
    phmm->beta[0] = *beta_0;
    phmm->tau[0] = *tau_0;
    phmm->alpha[1] = *alpha_1;
    phmm->beta[1] = *beta_1;
    phmm->tau[1] = *tau_1;
  }

  phmm_reset(phmm); 

  if (estim_trees || estim_rho) {
    msa->ncats = phmm->nmods - 1;   /* ?? */
    if (msa->ss == NULL) 
      ss_from_msas(msa, phmm->mods[0]->order+1, TRUE, NULL, NULL, NULL, -1);
    else if (msa->ss->cat_counts == NULL)
      ss_realloc(msa, msa->ss->tuple_size, msa->ss->ntuples, TRUE, TRUE);
  }

  if (estim_trees) {
    /* force re-initialization of tree models with rate variation;
       this is a hack that helps keep the parameterization simple */
    if (phmm->mods[0]->nratecats == 1) {
      tm_reinit(phmm->mods[0], phmm->mods[0]->subst_mod, 2, 
                phmm->mods[0]->alpha, NULL, NULL);
      phmm->mods[0]->nratecats = 1; phmm->mods[0]->alpha = -1; /* ignore rate variation */
      phmm->mods[0]->freqK[0] = phmm->mods[0]->rK[0] = 1;
    }
    if (phmm->mods[1]->nratecats == 1) {
      tm_reinit(phmm->mods[1], phmm->mods[1]->subst_mod, 2, 
                phmm->mods[1]->alpha, NULL, NULL);
      phmm->mods[1]->nratecats = 1; phmm->mods[1]->alpha = -1;
      phmm->mods[1]->freqK[0] = phmm->mods[1]->rK[0] = 1;
    }

    retval = hmm_train_by_em(phmm->hmm, phmm->mods, phmm, 1, &phmm->alloc_len, NULL, 
                             phmm_compute_emissions_em, reestimate_trees,
                             gamma > 0 ? 
                             phmm_estim_trans_em_coverage : phmm_estim_trans_em, 
                             phmm_get_obs_idx_em, 
                             phmm_log_em, phmm->emissions, logf) * log(2);

    /* have to do final rescaling of tree models to get units of subst/site */
    if (phmm->mods[0]->subst_mod != JC69 && phmm->mods[0]->subst_mod != F81) {   
                                /* JC69 and F81 are exceptions */
      double scale = tm_scale_rate_matrix(phmm->mods[0]); 
      tm_scale_rate_matrix(phmm->mods[1]); /* will be the same */
      tm_scale(phmm->mods[0], scale, 0); 
      tm_scale(phmm->mods[1], scale, 0);
    }

    phmm->mods[0]->lnL = phmm->mods[1]->lnL = retval;
  }

  else if (estim_rho) {
    phmm->mods[0]->estimate_branchlens = TM_SCALE_ONLY;
    phmm->mods[0]->scale = phmm->em_data->rho;
    tm_set_subst_matrices(phmm->mods[0]);

    retval = hmm_train_by_em(phmm->hmm, phmm->mods, phmm, 1, &phmm->alloc_len, NULL, 
                             compute_emissions_estim_rho, reestimate_rho,
                             gamma > 0 ? 
                             phmm_estim_trans_em_coverage : phmm_estim_trans_em, 
                             phmm_get_obs_idx_em, 
                             phmm_log_em, phmm->emissions, logf) * log(2);

    /* do final rescaling of conserved tree */
    tm_scale(phmm->mods[0], phmm->em_data->rho, FALSE);  
    phmm->mods[0]->scale = 1;

    phmm->mods[0]->lnL = phmm->mods[1]->lnL = retval;
  }

  else {                        /* not estimating tree models */
    retval = hmm_train_by_em(phmm->hmm, phmm->mods, phmm, 1, 
                             &phmm->alloc_len, NULL, NULL, NULL,
                             gamma > 0 ? 
                             phmm_estim_trans_em_coverage : phmm_estim_trans_em, 
                             NULL, phmm_log_em, phmm->emissions, logf) * log(2);
  }

  *mu = mm_get(phmm->functional_hmm->transition_matrix, 0, 1);
  *nu = mm_get(phmm->functional_hmm->transition_matrix, 1, 0);
  *rho = phmm->em_data->rho;

  if (phmm->indel_mode == PARAMETERIC) {
    *alpha_0 = phmm->alpha[0];
    *beta_0 = phmm->beta[0];
    *tau_0 = phmm->tau[0];
    *alpha_1 = phmm->alpha[1];
    *beta_1 = phmm->beta[1];
    *tau_1 = phmm->tau[1];
  }

  return retval;
}

/* Special-purpose unpack function, adapted from tm_unpack_params */
void unpack_params_mod(TreeModel *mod, gsl_vector *params) {
  TreeNode *n;
  int assigned = 0, nodeidx, i;
  List *traversal;

  assert(!mod->estimate_backgd && !mod->empirical_rates);

  /* check parameter values */
  for (i = 0; i < params->size; i++) {
    double mu = gsl_vector_get(params, i);
    if (mu < 0 && abs(mu) < TM_IMAG_EPS) /* consider close enough to 0 */
      gsl_vector_set(params, i, mu=0);
    if (mu < 0) die("ERROR: parameter %d has become negative (%f).\n", i, mu);
    if (!finite(mu)) die("ERROR: parameter %d is no longer finite (%f).\n", i, mu);
  }
  i = 0;

  if (mod->estimate_branchlens == TM_SCALE_ONLY) 
    mod->scale = gsl_vector_get(params, i++);
  else if (mod->estimate_branchlens == TM_BRANCHLENS_ALL) {
    traversal = tr_preorder(mod->tree);
    for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
      n = lst_get_ptr(traversal, nodeidx);

      if (n->parent != NULL) {
        if ((n == mod->tree->lchild || n == mod->tree->rchild) && 
            tm_is_reversible(mod->subst_mod)) {
          n->dparent = gsl_vector_get(params, 0)/2;
          if (!assigned) {
            i++;     /* only increment the first time */
            assigned = 1;
          }
        }
        else if (n->id == mod->root_leaf_id) 
          n->dparent = 0;
        else 
          n->dparent = gsl_vector_get(params, i++);
      }
    }
  }

  /* next parameters are for rate variation */
  if (mod->nratecats > 1) 
    mod->alpha = gsl_vector_get(params, i++);
  else i++;                     /* there's always a placeholder
                                   here */

  tm_set_rate_matrix(mod, params, i);

  /* diagonalize, if necessary */
  if (mod->subst_mod != JC69 && mod->subst_mod != F81)
    mm_diagonalize(mod->rate_matrix);
}

/* Unpack all params for two-state HMM.  Used during M step of EM */
void unpack_params_phmm(PhyloHmm *phmm, gsl_vector *params) {
  unpack_params_mod(phmm->mods[0], params);
  unpack_params_mod(phmm->mods[1], params);
  phmm->em_data->rho = gsl_vector_get(params, params->size - 1);
  tm_scale(phmm->mods[0], gsl_vector_get(params, params->size-1), FALSE);
  
  if (phmm->mods[0]->nratecats > 1) 
    DiscreteGamma(phmm->mods[0]->freqK, phmm->mods[0]->rK, phmm->mods[0]->alpha, 
                  phmm->mods[0]->alpha, phmm->mods[0]->nratecats, 0); 
                                /* mods[0]->alpha will already be set */
  if (phmm->mods[1]->nratecats > 1) {
    phmm->mods[1]->alpha = gsl_vector_get(params, params->size - 2);
    DiscreteGamma(phmm->mods[1]->freqK, phmm->mods[1]->rK, phmm->mods[1]->alpha, 
                  phmm->mods[1]->alpha, phmm->mods[1]->nratecats, 0); 
  }
  tm_set_subst_matrices(phmm->mods[0]);
  tm_set_subst_matrices(phmm->mods[1]);
}
 
/* Wrapper for computation of likelihood, for use by reestimate_trees (below) */
double likelihood_wrapper(gsl_vector *params, void *data) {
  PhyloHmm *phmm = (PhyloHmm*)data;
  double retval0, retval1;

  unpack_params_phmm(phmm, params);

  retval0 = -tl_compute_log_likelihood(phmm->mods[0], phmm->em_data->msa, NULL, 0, NULL);
  retval1 = -tl_compute_log_likelihood(phmm->mods[1], phmm->em_data->msa, NULL, 1, NULL);
  
  return retval0 + retval1;
                                /* FIXME: what happens when not one to
                                   one cats and mods? */
}

/* Re-estimate phylogenetic model based on expected counts (M step of EM) */
void reestimate_trees(void **models, int nmodels, void *data, 
                      double **E, int nobs, FILE *logf) {

  PhyloHmm *phmm = (PhyloHmm*)data;
  int k, obsidx;
  gsl_vector *params, *lower_bounds, *upper_bounds;
  double ll;

  /* FIXME: what about when multiple states per model?  Need to
     collapse sufficient stats.  Could probably be done generally...
     need to use state_to_cat, etc. in deciding which categories to
     use */

  for (k = 0; k < phmm->nmods; k++) 
    for (obsidx = 0; obsidx < nobs; obsidx++) 
      phmm->em_data->msa->ss->cat_counts[k][obsidx] = E[k][obsidx];

  params = gsl_vector_alloc(tm_get_nparams(phmm->mods[1]) + 2);
  tm_params_init_from_model(phmm->mods[1], params, 0); /* unscaled branch lens */

  /* special initialization of rate-variation parameters and rho */
  gsl_vector_set(params, tm_get_nbranchlenparams(phmm->mods[0]), 
                 phmm->mods[0]->nratecats > 1 ? phmm->mods[0]->alpha : 0);
  gsl_vector_set(params, params->size - 2,  
                 phmm->mods[1]->nratecats > 1 ? phmm->mods[1]->alpha : 0);
  gsl_vector_set(params, params->size - 1, phmm->em_data->rho);

  lower_bounds = gsl_vector_calloc(params->size);
  upper_bounds = gsl_vector_alloc(params->size);
  gsl_vector_set_all(upper_bounds, INFTY);
  gsl_vector_set(upper_bounds, params->size - 1, 1); /* 0 < rho < 1 */

  if (logf != NULL)
    fprintf(logf, "\nRE-ESTIMATION OF TREE MODEL:\n");

  /* keep Hessian arround so it can be used from one iteration to the
     next */
  if (phmm->em_data->H == NULL) {
    phmm->em_data->H = gsl_matrix_alloc(params->size, params->size);
    gsl_matrix_set_identity(phmm->em_data->H);
  }

  if (opt_bfgs(likelihood_wrapper, params, phmm, &ll, lower_bounds, 
               NULL, logf, NULL, OPT_MED_PREC, phmm->em_data->H) != 0)
    die("ERROR returned by opt_bfgs.\n");

  if (logf != NULL) 
    fprintf(logf, "END RE-ESTIMATION OF TREE MODEL\n\n");

  unpack_params_phmm(phmm, params);

  if (phmm->indel_mode == PARAMETERIC)
    phmm_set_branch_len_factors(phmm);

  gsl_vector_free(params); 
  gsl_vector_free(lower_bounds);
  gsl_vector_free(upper_bounds);
}

/* Wrapper for computation of likelihood, for use by reestimate_rho (below) */
double likelihood_wrapper_rho(double rho, void *data) {
  PhyloHmm *phmm = (PhyloHmm*)data;
  phmm->mods[0]->scale = rho;
  tm_set_subst_matrices(phmm->mods[0]);
  return -tl_compute_log_likelihood(phmm->mods[0], phmm->em_data->msa, 
				    NULL, 0, NULL);
}

/* Similar to reestimate_trees, but re-estimate only scale parameter
   rho; only the first model (for the conserved state) needs to be
   considered */
void reestimate_rho(void **models, int nmodels, void *data, 
		    double **E, int nobs, FILE *logf) {

  PhyloHmm *phmm = (PhyloHmm*)data;
  int obsidx;
  double ll, ax, bx, cx, fa, fb, fc;

  for (obsidx = 0; obsidx < nobs; obsidx++) 
    phmm->em_data->msa->ss->cat_counts[0][obsidx] = E[0][obsidx];

  if (logf != NULL)
    fprintf(logf, "\nRE-ESTIMATION OF RHO (BRENT'S METHOD):\n");

  bx = phmm->em_data->rho;
  ax = max(0.1, phmm->em_data->rho - .05);
  mnbrak(&ax, &bx, &cx, &fa, &fb, &fc, likelihood_wrapper_rho, phmm, logf);
  ll = opt_brent(ax, bx, cx, likelihood_wrapper_rho, 5e-3, 
		 &phmm->em_data->rho, phmm, logf);

  if (logf != NULL) 
    fprintf(logf, "END RE-ESTIMATION OF RHO\n\n");

  assert(phmm->indel_mode != PARAMETERIC);
  /* FIXME: to make work with parameteric indel model, will have to
     propagate scale parameter through phmm_set_branch_len_factors */
}

/* Maximize HMM transition parameters subject to constrain implied by
   target coverage (M step of EM).  For use with two-state HMM.  This
   function is passed to hmm_train_by_em in phmm_fit_em */
void phmm_estim_trans_em_coverage(HMM *hmm, void *data, double **A) {

  PhyloHmm *phmm = data;
  IndelEstimData *ied = NULL;
  double **C;

  if (phmm->em_data->fix_functional && phmm->em_data->fix_indel) return;

  if (phmm->indel_mode == PARAMETERIC) {
    ied = phmm_new_ied(phmm, A);
    C = ied->fcounts;
  }
  else C = A;

  /* estimate transition probs for functional cats subject to
     constraint on coverage */
  if (!phmm->em_data->fix_functional) {
    double a, b, c, mu, nu, nu1, nu2, z, tmp;
    /* if you take the first derivative wrt nu of the expression inside
       the argmax and set it to zero, you get a quadratic eqn which
       can be solved using the quadratic formula */
    z = (1-phmm->em_data->gamma)/phmm->em_data->gamma;
    a = z * (C[0][0] + C[0][1] + C[1][0] + C[1][1]);
    b = -C[0][1] - C[1][0] - C[1][1] - z * (C[0][0] + C[0][1] + C[1][0]);
    c = C[0][1] + C[1][0];

    tmp = b*b - 4*a*c;
    assert (tmp >= 0);
    tmp = sqrt(tmp);
    nu1 = (-b + tmp) / (2*a);
    nu2 = (-b - tmp) / (2*a);
    /* only one root can be valid */
    if (nu1 < 1e-10 || z * nu1 > 1 - 1e-10)                                 
      nu = nu2;                   /* (allow for rounding errors) */
    else nu = nu1;

    /* double check that derivative is really zero */
    assert(fabs(-z*C[0][0]/(1-z*nu) + (C[0][1] + C[1][0])/nu - C[1][1]/(1-nu)) < 1e-6);

    mu = z * nu;
    assert(nu >= 0 && nu <= 1 && mu >= 0 && mu <= 1);

    mm_set(phmm->functional_hmm->transition_matrix, 0, 0, 1-mu);
    mm_set(phmm->functional_hmm->transition_matrix, 0, 1, mu);
    mm_set(phmm->functional_hmm->transition_matrix, 1, 0, nu);
    mm_set(phmm->functional_hmm->transition_matrix, 1, 1, 1-nu);

    /* use stationary distribution for begin transitions */
    gsl_vector_set(phmm->functional_hmm->begin_transitions, 0, nu/(mu+nu));
    gsl_vector_set(phmm->functional_hmm->begin_transitions, 1, mu/(mu+nu));
  }

  if (phmm->indel_mode == PARAMETERIC) {
    if (!phmm->em_data->fix_indel) phmm_em_estim_indels(phmm, ied);
    phmm_free_ied(ied);
  }

  phmm_reset(phmm);
}

/* attempts to merge cats by incorporating together into ranges */
void collapse_cats(CategoryMap *cm, List *cats_to_merge) {
  /* assumes all ranges are of size one initially */
  int i, beg, end = -INFTY;
  lst_qsort_int(cats_to_merge, ASCENDING);
  for (i = 0; i < lst_size(cats_to_merge); i++) {
    int cat = lst_get_int(cats_to_merge, i);
    if (cat == end + 1) {
      cm_free_category_range(cm->ranges[cat]);
      cm->ranges[cat] = cm->ranges[beg];
      end = cat;
      cm->ranges[beg]->end_cat_no = end;
    }
    else beg = end = cat;
  }
}

/* initialize equilibrium freqs for tree model; either make consistent
   with given G+C content or estimate from alignment */
void init_eqfreqs(TreeModel *mod, MSA *msa, double gc) {
  if (gc != -1) {               /* gc specified */
    if (strlen(mod->rate_matrix->states) != 4 || 
        mod->rate_matrix->inv_states[(int)'A'] < 0 ||
        mod->rate_matrix->inv_states[(int)'C'] < 0 ||
        mod->rate_matrix->inv_states[(int)'G'] < 0 ||
        mod->rate_matrix->inv_states[(int)'T'] < 0)
      die("ERROR: Four-character DNA alphabet required with --gc.\n");
    assert(gc > 0 && gc < 1);
    gsl_vector_set(mod->backgd_freqs, 
                   mod->rate_matrix->inv_states[(int)'G'], gc/2);
    gsl_vector_set(mod->backgd_freqs, 
                   mod->rate_matrix->inv_states[(int)'C'], gc/2);
    gsl_vector_set(mod->backgd_freqs, 
                   mod->rate_matrix->inv_states[(int)'A'], (1-gc)/2);
    gsl_vector_set(mod->backgd_freqs, 
                   mod->rate_matrix->inv_states[(int)'T'], (1-gc)/2);
  }
  else {                        /* estimate from alignment */
    if (mod->subst_mod == JC69 || mod->subst_mod == K80)
      gsl_vector_set_all(mod->backgd_freqs, 1.0/mod->backgd_freqs->size);
    else
      msa_get_base_freqs_tuples(msa, mod->backgd_freqs, mod->order+1, -1);
  }
}
