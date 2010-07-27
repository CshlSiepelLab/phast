#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <getopt.h>
#include <misc.h>
#include <maf.h>
#include <sufficient_stats.h>
#include <dmotif_phmm.h>
#include "dmTranslateToRealigned.help"

#define DEFAULT_REFSPEC "hg18";

int main(int argc, char *argv[]) {
  int i, j, k, opt_idx, msa_len_orig, msa_len_trans, refseq_len,
    **idx_map, refstart_orig, refend_orig, start_trans, end_trans,
    refidx_orig, refidx_trans, end_coord, md_offset, inf_len, algn_end,
    inf_end;
  char c, *msa_path_orig, *msa_path_trans, *refseq_path, *msa_suffix,
    *refseq_suffix, c_msa, c_ref, *gff_trans_fname, *gff_unmapped_fname,
    seq1_end[11], seq2_end[11];
  List *seqlist;
  MSA *msa_orig, *msa_trans, *refseq, *tmp_msa;
  GFF_Set *gff_orig, *gff_trans, *gff_unmapped;
  GFF_Feature *f_orig;
  String *f_orig_seqname, *msa_fname_orig, *msa_fname_trans, *refseq_fname;
  FILE *tmp_f;

  struct option long_opts[] = {
    {"msa-format", 1, 0, 'i'},
    {"refseq-format", 1, 0, 'f'},
    {"refspecies", 1, 0, 'r'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };
  
  /* arguments and defaults for options */
  char *refspec = DEFAULT_REFSPEC;
  msa_format_type msa_format = FASTA,
    refseq_format = FASTA;

#ifdef RPHAST
  GetRNGstate(); //seed R's random number generator
#endif

  while ((c = getopt_long(argc, argv,
			  "i:f:r:h", 
			  long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'i':
      msa_format = msa_str_to_format(optarg);
      if (msa_format == -1)
        die("ERROR: unrecognized alignment format.\n");
      break;
    case 'f':
      refseq_format = msa_str_to_format(optarg);
      if (refseq_format == -1)
        die("ERROR: unrecognized refseq format.\n");
      break;
    case 'r':
      refspec = optarg;
      break;
    case 'h':
      printf(HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try 'dmTranslateToRealigned -h'.\n");
    }
  }

  if (optind != argc - 6) 
    die("Six arguments required.  Try 'dmTranslateToRealigned -h'.\n");

  /* To-Do:
     3) Step through GFF features to cross-map
         -- Pull out hg18 sequence and attempt to map within hg18 in msa_2
         -- If don't find match in hg18, try next seq...

	 -- OR --
	 -- Create mapping of columns from msa_orig to msa_trans
	 -- Use mapping to relocate features from gff_orig to gff_trans

	 -- Create output feature and add to gff

     4) Write output files
     5) Clean up and exit
  */

  /* Read in the original GFF */
  tmp_f = fopen_fname(argv[optind], "r");
  gff_orig = gff_read_set(tmp_f);
  fclose(tmp_f);

  /* Get the orig and trans paths for the msa files */
  msa_path_orig = argv[optind+1];
  msa_path_trans = argv[optind+2];

  /* Get the path to the reference sequence files */
  refseq_path = argv[optind+3];

  /* Get the names of the output files */
  gff_trans_fname = argv[optind+4];
  gff_unmapped_fname = argv[optind+5];

  /* Get the suffix for MSA and refseq files */
  msa_suffix = msa_suffix_for_format(msa_format);
  refseq_suffix = msa_suffix_for_format(refseq_format);

  /* Allocate strings */
  msa_fname_orig = str_new(STR_MED_LEN);
  msa_fname_trans = str_new(STR_MED_LEN);
  refseq_fname = str_new(STR_MED_LEN);

  /* Allocate space for column to refseq mappings */
  idx_map = smalloc(2 * sizeof(int*));

  /* Create GFF's for output and unmapped features */
  gff_trans = gff_new_set();
  gff_unmapped = gff_new_set();


  /* Step through the original gff features, opening and closing corresponding
     msa files as needed */
  for (i = 0; i < lst_size(gff_orig->features); i++) {
    f_orig = lst_get_ptr(gff_orig->features, i);
    
    /* Read in the MSA's if we're working with a new sequence and free the
       previous set, if needed. */
    if (i == 0 || !str_equals(f_orig->seqname, f_orig_seqname)) {
      if (i > 0) {
	str_free(f_orig_seqname);
	msa_free(msa_orig);
	msa_free(msa_trans);
	msa_free(refseq);
	free(idx_map[0]);
	free(idx_map[1]);
      }
      f_orig_seqname = str_dup(f_orig->seqname);

      str_clear(msa_fname_orig);
      str_append_charstr(msa_fname_orig, msa_path_orig);
      if (msa_fname_orig->chars[msa_fname_orig->length-1] != '/')
	str_append_charstr(msa_fname_orig, "/");
      str_append_charstr(msa_fname_orig, f_orig_seqname->chars);
      str_append_charstr(msa_fname_orig, ".");
      str_append_charstr(msa_fname_orig, msa_suffix);

      str_clear(msa_fname_trans);
      str_append_charstr(msa_fname_trans, msa_path_trans);
      if (msa_fname_trans->chars[msa_fname_trans->length-1] != '/')
	str_append_charstr(msa_fname_trans, "/");
      str_append_charstr(msa_fname_trans, f_orig_seqname->chars);
      str_append_charstr(msa_fname_trans, ".");
      str_append_charstr(msa_fname_trans, msa_suffix);

      str_clear(refseq_fname);
      str_append_charstr(refseq_fname, refseq_path);
      if (msa_fname_trans->chars[refseq_fname->length-1] != '/')
	str_append_charstr(refseq_fname, "/");
      str_append_charstr(refseq_fname, f_orig_seqname->chars);
      str_append_charstr(refseq_fname, ".");
      str_append_charstr(refseq_fname, refspec);
      str_append_charstr(refseq_fname, ".");
      str_append_charstr(refseq_fname, refseq_suffix);

/*       fprintf(stderr, "msa_fname_orig %s, msa_fname_trans %s,  refseq_fname %s\n", */
/* 	      msa_fname_orig->chars, msa_fname_trans->chars, refseq_fname->chars); */

      tmp_f = fopen_fname(msa_fname_orig->chars, "r");
      msa_orig = msa_new_from_file(tmp_f, msa_format, NULL);
      fclose(tmp_f);

      tmp_f = fopen_fname(msa_fname_trans->chars, "r");
      msa_trans = msa_new_from_file(tmp_f, msa_format, NULL);
      fclose(tmp_f);

      tmp_f = fopen_fname(refseq_fname->chars, "r");
      refseq = msa_new_from_file(tmp_f, refseq_format, NULL);
      fclose(tmp_f);
      
      msa_len_orig = msa_orig->length;
      msa_len_trans = msa_trans->length;
      refseq_len = refseq->length;

/*       fprintf(stderr, "msa_len_orig %d, msa_len_trans %d, refseq_len %d\n", */
/* 	      msa_len_orig, msa_len_trans, refseq_len); */

      idx_map[0] = smalloc(msa_len_orig * sizeof(int));
      idx_map[1] = smalloc(refseq_len * sizeof(int));
      for (j = 0; j < refseq_len; j++)
	idx_map[1][j] = -1;

      /* Look up the indices of the reference species */
      refidx_orig = msa_get_seq_idx(msa_orig, refspec);
      refidx_trans = msa_get_seq_idx(msa_trans, refspec);
      
      /* Map columns in msa_orig to the reference sequence */
      k = 0; /* For index in reference sequence */
      for (j = 0; j < msa_len_orig; j++) {
	c_msa = msa_get_char(msa_orig, refidx_orig, j);
	if (c_msa == GAP_CHAR || msa_orig->is_missing[(int)c_msa]) {
	  idx_map[0][j] = -1;
	} else {
	  idx_map[0][j] = k;
	  k++;
	}
/* 	fprintf(stderr, "%d ", idx_map[0][j]); */
      }
/*       fprintf(stderr, "\n\n"); */

      /* Some sequences will have missing and/or gap chars at one or both ends
	 in the realigned sequence. The missing columns from the reference 
	 sequence need to be accounted for. */
      md_offset = 0;
      seqlist = lst_new_int(2);
      lst_push_int(seqlist, refidx_trans);
      tmp_msa = msa_sub_alignment(msa_trans, seqlist, 1, 0,
				  msa_trans->length);
      msa_strip_gaps(tmp_msa, STRIP_ALL_GAPS);
      inf_len = tmp_msa->length;
      for (j = 0; j < tmp_msa->length; j++) {
	c = msa_get_char(tmp_msa, 0, j);
	if (tmp_msa->is_missing[(int)c])
	  inf_len--;
	/* stop when we encounter an informative character */
/* 	else */
/* 	  break; */
      }
      
      md_offset = (refseq->length - inf_len);
/*       fprintf(stderr, "refseq->length %d, msa_len_trans %d, inf_len %d, md_offset %d\n", */
/* 	      refseq->length, msa_len_trans, inf_len, md_offset); */
      /* Check the end of the realigned sequence for missing data characters */
      inf_end = tmp_msa->length;
      for (j = tmp_msa->length-1; j >= 0; j--) {
	c = msa_get_char(tmp_msa, 0, j);
	if (tmp_msa->is_missing[(int)c])
	  inf_end--;
	else
	  break;
      }

      /* Check that the ends of the reference and realigned seqs are the same
	 and adjust the offset if characters have been dropped in the realigned
	 version */
      k = algn_end = 0;
      seq1_end[10] = seq2_end[10] = '\0';
      for (j = 0; j < 10; j++) {
	c = msa_get_char(tmp_msa, 0, ((inf_end-10) + j));
	seq1_end[j] = c;
      }
/*       fprintf(stderr, "seq1_end %s\n", seq1_end); */
      while (algn_end == 0) {
	for (j = 0; j < 10; j++) {
	  c = msa_get_char(refseq, 0, ((refseq->length - (10 + k)) + j));
	  seq2_end[j] = c;
	}
/* 	fprintf(stderr, "k %d, seq2_end %s\n", k, seq2_end); */
	if (!strcmp(seq1_end, seq2_end))
	  algn_end = (msa_trans->length - k);
	else
	  k++;
      }
      md_offset -= k;
/*       fprintf(stderr, "refseq->length %d, tmp_msa->length %d, inf_len %d, md_offset %d\n", */
/* 	      refseq->length, tmp_msa->length, inf_len, md_offset); */
      
      lst_free(seqlist);
      msa_free(tmp_msa);


      /* Map columns from the reference sequence to msa_trans */
      k = md_offset; /* For index in refseq */
      for (j = 0; j < msa_len_trans; j++) {
	c_msa = msa_get_char(msa_trans, refidx_trans, j);
	if (c_msa == GAP_CHAR || msa_trans->is_missing[(int)c_msa]) {
	  continue;
	} else {
	  c_ref = msa_get_char(refseq, 0, k);
/* 	  fprintf(stderr, "j %d, k %d, c_msa %c, c_ref %c\n", */
/* 		  j, k, c_msa, c_ref); */
	  if (c_msa == c_ref) {
	    idx_map[1][k] = j;
	  } else {
	    idx_map[1][k] = -1;
	  }
/* 	  fprintf(stderr, "idx_map[1][%d] %d\n", */
/* 		  k, idx_map[1][k]); */
	  k++;
	}
      }
/*       fprintf(stderr, "\n\n");       */
    }

    /* Remap the start and end coordinates */ 
    
    /* this hack suppresses errors due to features that overlap the sequence
       end */
    if (f_orig->end == msa_len_orig)
      end_coord = f_orig->end - 1;
    else
      end_coord = f_orig->end;

    refstart_orig = idx_map[0][f_orig->start];
    refend_orig = idx_map[0][end_coord];
    start_trans = end_trans = -1;

/*     fprintf(stderr, "f_orig->start %d, f_orig->end %d, refstart_orig %d, refend_orig %d\n", */
/* 	    f_orig->start, f_orig->end, refstart_orig, refend_orig); */

    if (refstart_orig == -1 && refend_orig == -1) {
      /* Scan the motif window for any mappable (non-gap) character */
/*       offset = 1; */
/*       for (j = (f_orig->start + 1); j < f_orig->end; j++) { */
/* 	if (idx_map[0][j] != -1) { /\* Found a mappable character *\/ */
/* 	  start_trans = idx_map[0][j] - offset; */
/* 	  end_trans = start_trans + (f_orig->end - f_orig->start); */
/* 	  break; */
/* 	} else { */
/* 	  offset++; */
/* 	} */
/*       } */
/*       if (start_trans == -1) { */
	lst_push_ptr(gff_unmapped->features, gff_new_feature_copy(f_orig));
	continue;
/*       } */
    } else if (refstart_orig == -1) {
      end_trans = idx_map[1][refend_orig];
      start_trans = end_trans - (f_orig->end - f_orig->start);

    } else if (refend_orig == -1) {
      start_trans = idx_map[1][refstart_orig];
      end_trans = start_trans + (f_orig->end - f_orig->start);

    } else {
      start_trans = idx_map[1][refstart_orig];
      end_trans = idx_map[1][refend_orig];
    }

    if (start_trans == -1 || end_trans == -1) {
      lst_push_ptr(gff_unmapped->features, gff_new_feature_copy(f_orig));
      continue;
    }
    
/*     fprintf(stderr, "start %d, end %d, refstart_orig %d, refend_orig %d, start_trans %d, end_trans %d\n", */
/* 	    f_orig->start, f_orig->end, refstart_orig, refend_orig, */
/* 	    start_trans, end_trans); */

    lst_push_ptr(gff_trans->features, 
		 gff_new_feature(str_dup(f_orig->seqname), 
				 str_dup(f_orig->source),
				 str_dup(f_orig->feature),
				 start_trans, end_trans, f_orig->score, 
				 f_orig->strand, f_orig->frame,
				 str_dup(f_orig->attribute), TRUE));
  }
  
  /* Print the GFF's to their respective output files */
  tmp_f = fopen_fname(gff_trans_fname, "w");
  gff_print_set(tmp_f, gff_trans);
  fclose(tmp_f);
  
  tmp_f = fopen_fname(gff_unmapped_fname, "w");
  gff_print_set(tmp_f, gff_unmapped);
  fclose(tmp_f);

  gff_free_set(gff_orig);
  str_free(msa_fname_orig);
  str_free(msa_fname_trans);
  str_free(refseq_fname);
  gff_free_set(gff_trans);
  gff_free_set(gff_unmapped);
  msa_free(msa_orig);
  msa_free(msa_trans);
  str_free(f_orig_seqname);
  free(idx_map);
  
  return 0;
}

