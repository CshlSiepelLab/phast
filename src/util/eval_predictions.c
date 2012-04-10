/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: eval_predictions.c,v 1.14 2008-11-12 02:07:58 acs Exp $ */

#include <stdlib.h>
#include <stdio.h>
#include "gff.h"
#include <getopt.h>
#include <math.h>
#include <misc.h>

void print_usage() {
  printf("USAGE: eval_predictions -r <real_fname_list> -p <pred_fname_list>\n\
                        -l <seq_len_list> [OPTIONS] \n\
\nDESCRIPTION:\n\
Compares predicted genes with \"real\" (or annotated) genes.\n\
Reports standard measures of prediction quality.  The following\n\
measures are reported:\n\
\n\
    - nucleotide sensitivity (Sn)\n\
    - nucleotide specificity (Sp)\n\
    - approximate correlation (AC)\n\
    - correlation coefficient (CC)\n\
    - exon sensitivity (ESn)\n\
    - exon specificity (ESp)\n\
    - proportion of real exons correctly predicted (CRa)\n\
    - proportion of real exons partially predicted (PCa)\n\
    - proportion of real exons with overlapping predictions (OLa)\n\
    - missed exons (ME)\n\
    - proportion of predicted exons that are correct (CRp)\n\
    - proportion of predicted exons that are partially correct (PCp)\n\
    - proportion of predicted exons that overlap real ones (OLp)\n\
    - wrong exons (WE)\n\
\n\
All quantities are computed as described in \"Evaluation of Gene-Finding\n\
Programs on Mammalian Sequences,\" by Rogic et al. (Genome Research\n\
11:817-832).  Note that CRa + PCa + OLa + ME = 1 and CRp + PCp + OLp +\n\
WE = 1.  Note also that each set (predicted and real) should consist of\n\
non-overlapping groups of features (see 'refeature').\n\
\n\
OPTIONS:\n\
    -r <real_fname_list> \n\
        (required) List of names of files defining real genes (GFF). \n\
\n\
    -p <pred_fname_list> \n\
        (required) List of names of files defining predicted genes\n\
        (GFF).  Must correspond in order to <real_fname_list>.\n\
\n\
    -l <seq_len_list>\n\
        (required) List of lengths of sequences.  Needed to compute\n\
        certain nucleotide-level statistics.\n\
\n\
    -f <feat_list>  \n\
        List of names of all features denoting exon regions.  By\n\
        default, equal to the single name \"CDS\".\n\
\n\
    -d <fname_prefix>  \n\
        Dump full coords of correct, partially correct, wrong, missed, \n\
        and overlapping exons to a set of files having the specified \n\
        file name prefix.\n\
\n\
    -n <nbases>\n\
        Also report stats on \"nearly correct\" exons, that is, incorrect\n\
        exons whose boundaries are within <nbases> of being correct.\n\
        Columns will be labeled \"NCa\" and \"NCp\".\n\
\n\
    -h\n\
        Print this help message.\n\
\n\
\n\
NOTE: be sure stop codons are included in CDSs in both the predicted\n\
and real sets, or in neither set.\n\n");  
} 

int is_exon(GFF_Feature *feat, List *l) {
  int i;
  for (i = 0; i < lst_size(l); i++) 
    if (str_equals_nocase(feat->feature, (String*)lst_get_ptr(l, i)))
      return 1;
  return 0;
}

/* types of predictions */
typedef enum {CR, PC, WE, ME, OL, NC} pred_type;
static FILE *SUMF = NULL;       /* hack to allow total number of bases
                                   in (real) exons to be written to
                                   summary file */

/* dump record of predictions to files */
void dump(char *prefix, GFF_Feature *real, GFF_Feature *pred, pred_type t,
          double pct_correct) {
  char *type=NULL;
  static FILE *CRF = NULL, *PCF = NULL, *NCF = NULL, *WEF = NULL, 
    *MEF = NULL, *OLF = NULL;

  /* set up files on first call */
  if (CRF == NULL) {
    String *fname = str_new(STR_SHORT_LEN);
    str_cpy_charstr(fname, prefix);
    str_append_charstr(fname, ".log.cr");
    CRF = phast_fopen(fname->chars, "w+");
    str_cpy_charstr(fname, prefix);
    str_append_charstr(fname, ".log.pc");
    PCF = phast_fopen(fname->chars, "w+");
    str_cpy_charstr(fname, prefix);
    str_append_charstr(fname, ".log.nc");
    NCF = phast_fopen(fname->chars, "w+");
    str_cpy_charstr(fname, prefix);
    str_append_charstr(fname, ".log.we");
    WEF = phast_fopen(fname->chars, "w+");
    str_cpy_charstr(fname, prefix);
    str_append_charstr(fname, ".log.me");
    MEF = phast_fopen(fname->chars, "w+");
    str_cpy_charstr(fname, prefix);
    str_append_charstr(fname, ".log.ol");
    OLF = phast_fopen(fname->chars, "w+");
    str_cpy_charstr(fname, prefix);
    str_append_charstr(fname, ".summary");
    SUMF = phast_fopen(fname->chars, "w+");
    str_free(fname);
  }

  switch (t) {
  case CR:
    gff_print_feat(CRF, real);
    gff_print_feat(CRF, pred);
    type = "CR";
    break;
  case PC:
    gff_print_feat(PCF, real);
    gff_print_feat(PCF, pred);
    type = "PC";
    break;
  case NC:
    gff_print_feat(NCF, real);
    gff_print_feat(NCF, pred);
    type = "NC";
    break;
  case WE:
    gff_print_feat(WEF, pred);
    type = "WE";
    break;
  case ME:
    gff_print_feat(MEF, real);
    break;
  case OL:
    gff_print_feat(OLF, real);
    gff_print_feat(OLF, pred);
    type = "OL";
    break;
  }

  /* also write per-prediction summary, with scores */
  if (t != ME) 
    fprintf(SUMF, "%-12s %12d %12d %12.3f %10s %10.4f\n", pred->feature->chars, 
            pred->start, pred->end, pred->score_is_null ? -1 : pred->score, 
            type, pct_correct);
}

void compute_and_print_stats(FILE *F, String *real_name, String *pred_name, 
                             int tp, int fp, int nreal_pos, int npred_pos, 
                             int seqlen, int ncr, int npca, int nola, 
                             int nme, int npcp, int nolp, int nwe, 
                             int nexons_real, int nexons_pred, int nnc) {

    double Sn, Sp, CC, ACP, AC, ESn, ESp, CRa, PCa, OLa, ME, CRp, PCp, 
      OLp, WE, NCa, NCp;
    double fn, tn;              /* use doubles to avoid overflow */

    Sn = (double)tp/nreal_pos;
    Sp = (double)tp/npred_pos;
    fn = nreal_pos - tp;
    tn = seqlen - nreal_pos;
    if (tp + fn > 0 && tn + fp > 0 && tp + fp > 0 && tn+fn > 0) {
      CC = ((tp * tn) - (fn * fp)) / 
        sqrt((tp+fn) * (tn+fp) * (tp+fp) * (tn+fn)); 
      ACP = 0.25 * (tp/(tp+fn) + (double)tp/(tp+fp) + 
                    tn/(tn+fp) + tn/(tn+fn));
      AC = (ACP - 0.5) * 2;
    }
    else {
      CC = AC = 0;
    }

    ESn = (double)ncr/nexons_real;
    ESp = (double)ncr/nexons_pred;

    CRa = (double)ncr/nexons_real;
    PCa = (double)npca/nexons_real;
    OLa = (double)nola/nexons_real;
    ME = (double)nme/nexons_real;

    CRp = (double)ncr/nexons_pred;
    PCp = (double)npcp/nexons_pred;
    OLp = (double)nolp/nexons_pred;
    WE = (double)nwe/nexons_pred;

    if (nnc != -1) {
      NCa = (double)nnc/nexons_real;
      NCp = (double)nnc/nexons_pred;
    }

    fprintf(F, "%-25s %-25s %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f", real_name->chars, pred_name->chars, Sn, Sp, 
           AC, CC, ESn, ESp, CRa, PCa, OLa, ME, CRp, PCp, OLp, WE);
    if (nnc != -1) fprintf(F, " %7.4f %7.4f %7.4f %7.4f", NCa, NCp, NCa+CRa, NCp+CRp);
    fprintf(F, "\n");


}

int main(int argc, char* argv[]) {
  FILE* F;
  GFF_Set *gff_real=NULL, *gff_pred=NULL;
  char c;
  List *real_fname_list = NULL, *pred_fname_list = NULL, 
    *feat_list = NULL, *seq_len_list = NULL, *l = NULL;
  int nfile, i, j;
  char *prefix = NULL;
  int tot_tp = 0, tot_fp = 0, tot_nreal_pos = 0, tot_npred_pos = 0, 
    tot_seqlen = 0, tot_ncr = 0, tot_npca = 0, tot_nola = 0, tot_nme = 0, 
    tot_npcp = 0, tot_nolp = 0, tot_nwe = 0, tot_nexons_real = 0, 
    tot_nexons_pred = 0, dump_exons = 0, nnc = -1, tot_nnc = -1, 
    nc_threshold = 0;

  while ((c = getopt(argc, argv, "r:p:f:l:d:n:h")) != -1) {
    switch(c) {
    case 'r':
      real_fname_list = get_arg_list(optarg);
      break;
    case 'p':
      pred_fname_list = get_arg_list(optarg);
      break;
    case 'l':
      l = get_arg_list(optarg);
      /* convert to ints */
      seq_len_list = lst_new_int(lst_size(l));
      for (i = 0; i < lst_size(l); i++) {
        int tmp;
        if (str_as_int((String*)lst_get_ptr(l, i), 
                       &tmp) != 0) {
          die("ERROR: Bad integer in <seq_len_list>.\n"); 
        }
        lst_push_int(seq_len_list, tmp);
      }
      break;
    case 'f':
      feat_list = get_arg_list(optarg);
      break;
    case 'd':
      dump_exons = 1;
      prefix = optarg;
      break;
    case 'n':
      nnc = tot_nnc = 0;
      nc_threshold = get_arg_int(optarg);
      break;
    case 'h':
      print_usage();
      exit(0);
    case '?':
      die("Unrecognized option.  Try \"eval_predictions -h\" for help.\n");
    }
  }

  set_seed(-1);

  if (feat_list == NULL) {
    feat_list = lst_new_ptr(1);
    lst_push_ptr(feat_list, str_new_charstr(GFF_CDS_TYPE));
  }
  
  if (real_fname_list == NULL || pred_fname_list == NULL || 
      seq_len_list == NULL) {
    die("ERROR: Must specify -r, -p, and -l.  Try \"eval_predictions -h\" for help.\n");
  }

  if (lst_size(real_fname_list) != lst_size(pred_fname_list)) {
    die("ERROR: Must specify lists of equal length for real and predicted filenames.\n\n.");
  }

  if (lst_size(seq_len_list) == 1 && lst_size(real_fname_list) > 1)
    for (i = 1; i < lst_size(real_fname_list); i++)
      lst_push_int(seq_len_list, lst_get_int(seq_len_list, 0));
  else if (lst_size(seq_len_list) != lst_size(real_fname_list))
    die("ERROR: List of sequence lengths does not match lists of real and predicted filenames.\n");

  /* print header */
  printf("%-25s %-25s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s %7s", "Real_fname", "Pred_fname", "Sn", "Sp", "AC", "CC", "ESn", "ESp", "CRa", "PCa", "OLa", "ME", "CRp", "PCp", "OLp", "WE");
  if (nnc != -1) printf(" %7s %7s %7s %7s", "NCa", "NCp", "CR+NCa", "CR+NCp");
  printf("\n");

  for (nfile = 0; nfile < lst_size(real_fname_list); nfile++) {
    int tp, fp, nexons_real, nexons_pred, nwe, nme, ncr, npca, 
      npcp, nola, nolp, nreal_pos, npred_pos, len_real, len_pred, seqlen,
      already_counted_real;
    String *real_fname, *pred_fname;
    GFF_Feature *feat_real, *feat_pred=NULL;

    real_fname = (String*)lst_get_ptr(real_fname_list, nfile);
    F = phast_fopen(real_fname->chars, "r");
    if ((gff_real = gff_read_set(F)) == NULL) {
      die("ERROR: Unable to read file \"%s\".\n", 
	  real_fname->chars);
    }
    phast_fclose(F);

    pred_fname = (String*)lst_get_ptr(pred_fname_list, nfile);
    F = phast_fopen(pred_fname->chars, "r");
    if ((gff_pred = gff_read_set(F)) == NULL) {
      die("ERROR: Unable to read file \"%s\".\n", 
	  pred_fname->chars);
    }
    phast_fclose(F);

    seqlen = lst_get_int(seq_len_list, nfile);

    /* sort ungrouped -- only cds exons will be considered, and each
       one will be considered individually */
    gff_ungroup(gff_real); 
    gff_ungroup(gff_pred);
    gff_sort(gff_real);
    gff_sort(gff_pred);

    nexons_real = nexons_pred = nwe = nme = ncr = npca = npcp = nola = 
      nolp = tp = fp = nreal_pos = npred_pos = 0;
    if (nnc != -1) nnc = 0;
    i = j = 0;
    already_counted_real = 0;
    while (i < lst_size(gff_real->features)) {
      feat_real = (GFF_Feature*)lst_get_ptr(gff_real->features, i);
      if (!is_exon(feat_real, feat_list)) { i++; continue; }

      len_real = feat_real->end - feat_real->start + 1;

      if (!already_counted_real) {
        nexons_real++;
        nreal_pos += len_real;
      }

      /* look at all predicted exons up to and overlapping this real exon */
      while (j < lst_size(gff_pred->features)) {
        feat_pred = (GFF_Feature*)lst_get_ptr(gff_pred->features, j);
        if (!is_exon(feat_pred, feat_list)) {
          j++;
          continue;
        }
        else if (feat_pred->start > feat_real->end) {
          if (!already_counted_real) {
            nme++;
            if (dump_exons) dump(prefix, feat_real, NULL, ME, -1);
          }
          break;
        }

        /* otherwise we have a predicted exon to count (start of pred
           <= end of real) */
        nexons_pred++;
        len_pred = feat_pred->end - feat_pred->start + 1;
        npred_pos += len_pred;
        j++;                    /* we'll be done with this prediction
                                   one way or another; next time
                                   through look at a new one */

        if (feat_pred->end < feat_real->start) { /* WE */
          nwe++;
          fp += len_pred;
          if (dump_exons) dump(prefix, NULL, feat_pred, WE, 0);
        }
        else if (feat_pred->start == feat_real->start && /* CR */
                 feat_pred->end == feat_real->end) {
          ncr++;
          tp += len_pred;
          if (dump_exons) dump(prefix, feat_real, feat_pred, CR, 1);
          break;
        }
        else if (feat_pred->start == feat_real->start || /* PC */
                 feat_pred->end == feat_real->end) {
          pred_type type;
          npca++;
          npcp++;
          if (nnc != -1 && 
              max(abs(feat_pred->start - feat_real->start), 
                  abs(feat_pred->end - feat_real->end)) <= nc_threshold) {
            nnc++; 
            type = NC;
          }
          else type = PC;
          if (len_pred < len_real) 
            tp += len_pred;
          else {
            tp += len_real;
            fp += (len_pred - len_real);
          }
          if (dump_exons) dump(prefix, feat_real, feat_pred, type, 
                               min(1, (double)len_real/len_pred));
          break;
        }
        else {                  /* OL */
          int overlap_size;
          pred_type type;
          nola++;
          nolp++;
          if (nnc != -1 && 
              max(abs(feat_pred->start - feat_real->start), 
                  abs(feat_pred->end - feat_real->end)) <= nc_threshold) {
            nnc++; 
            type = NC;
          }
          else type = PC;

          overlap_size = min(feat_pred->end, feat_real->end) - 
            max(feat_pred->start, feat_real->start) + 1;
          tp += overlap_size;
          fp += len_pred - overlap_size;
          if (dump_exons) dump(prefix, feat_real, feat_pred, type,
                               (double)overlap_size/len_pred);
          break;
        }
        /* NOTE: I'm ignoring the possibility that a predicted exon
           could be a PC and/or OL with respect to multiple real
           exons.  The effect on the exon-level stats will be fairly
           minor (at worst a predicted exon is scored as an OL when it
           should be scored as an PC, and a real exon is erroneously
           counted as a ME), but the effect on the nucleotide-level Sn
           and Sp could conceivably be significant.  */
      }

      /* if we have counted at least one prediction (and thus failed
         to reach the end of the list), but the last prediction did
         not extend as far as the end of the real exon, then delay
         moving on to the next real exon */
      if (j < lst_size(gff_pred->features) && feat_pred->end < feat_real->end) 
          already_counted_real = 1;
      else {
        /* if we reached the end of the list of predictions, then it
           must not have contained any exons, and the real exon in
           question is a ME (if it hasn't already been counted) */
        if (j == lst_size(gff_pred->features) && !already_counted_real) 
          nme++; 

        i++;
        already_counted_real = 0;
      }
    }
    
    /* any remaining predictions must be wrong */
    for (; j < lst_size(gff_pred->features); j++) {
      if (is_exon((GFF_Feature*)lst_get_ptr(gff_pred->features, j), 
                  feat_list)) {
        nexons_pred++;
        nwe++;
      }
    }

    compute_and_print_stats(stdout, real_fname, pred_fname, 
                            tp, fp, nreal_pos, npred_pos, seqlen, ncr, 
                            npca, nola, nme, npcp, nolp, nwe, 
                            nexons_real, nexons_pred, nnc);

    tot_tp += tp;
    tot_fp += fp;
    tot_nreal_pos += nreal_pos;
    tot_npred_pos += npred_pos;
    tot_seqlen += seqlen;
    tot_ncr += ncr;
    tot_npca += npca;
    tot_nola += nola;
    tot_nme += nme;
    tot_npcp += npcp;
    tot_nolp += nolp;
    tot_nwe += nwe;
    tot_nexons_real += nexons_real;
    tot_nexons_pred += nexons_pred;
    if (nnc != -1) tot_nnc += nnc;

    if (dump_exons && SUMF != NULL)
      fprintf(SUMF, "# Total number of bases in real exons: %d\n", nreal_pos);

    gff_free_set(gff_real);
    gff_free_set(gff_pred);
  }

  if (lst_size(real_fname_list) > 1)
    compute_and_print_stats(stdout, str_new_charstr("TOTAL"), str_new_charstr(""), 
                            tot_tp, tot_fp, tot_nreal_pos, tot_npred_pos, 
                            tot_seqlen, tot_ncr, tot_npca, tot_nola, tot_nme, 
                            tot_npcp, tot_nolp, tot_nwe, tot_nexons_real, 
                            tot_nexons_pred, tot_nnc);

  return 0;
}

