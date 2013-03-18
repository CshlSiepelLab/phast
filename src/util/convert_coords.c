/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: convert_coords.c,v 1.4 2008-11-12 02:07:58 acs Exp $ */

#include <stdlib.h>
#include <stdio.h>
#include <msa.h>
#include <gff.h>
#include <getopt.h>
#include <local_alignment.h>

void print_usage() {
  fprintf(stderr, "USAGE: convert_coords -m <msa_fname> -f <feature_fname> [-s <src_frame>] [-d <dest_frame>] [-p] [-n] [-i PHYLIP|FASTA|MPM]\n\
\n\
DESCRIPTION:\n\
Converts coordinates of features in a GFF file according to a multiple\n\
alignment.  Will map from the coordinate system of any sequence to any\n\
other sequence; can also map to or from the coordinate system of the\n\
entire alignment.  In addition, supports translation of coordinates by\n\
specified offset.\n\
\n\
OPTIONS:\n\
    -m <msa_fname>  (required) Name of file in which alignment is defined.\n\
    -f <gff_fname>  (required) Name of file in which features are defined (GFF). \n\
    -s <src_frame>  Index of frame of reference for feature coordinates, as \n\
                    defined in the GFF file.  Use an integer 1-N (if N seqs) \n\
                    or 0 to indicate the coordinate system of the alignment \n\
                    as a whole.  Default behavior is to match features with \n\
                    alignment sequences by name (feature by feature).\n\
    -d <dest_frame> Index of destination frame of reference.  Default is 0\n\
                    (whole MSA).\n\
    -p <coord_off>  Positive coordinate offset.  This value will be\n\
                    added to all coordinates.  Useful when \n\
                    the alignment (or sequence) for which the coordinates \n\
                    are specified is a sub-alignment of yours. \n\
    -n <coord_off>  Negative coordinate offset.  This value will be\n\
                    subtracted from all coordinates.  Useful when your\n\
                    alignment is a sub-alignment of the alignment (or \n\
                    sequence) for which the coordinates are specified.\n\
    -i FASTA|PHYLIP|MPM|SS\n\
                    Alignment format.  Default is to guess format from file\n\
                    contents\n\n");  
} 

int main(int argc, char* argv[]) {
  FILE* F;
  MSA *msa;
  msa_format_type format = UNKNOWN_FORMAT;
  int src_ref = -1, dest_ref = 0, offset = 0;
  char *msa_fname = NULL, *feat_fname = NULL;
  GFF_Set *gff;
  char c;

  while ((c = getopt(argc, argv, "hm:f:s:d:i:p:n:")) != -1) {
    switch(c) {
    case 'm':
      msa_fname = optarg;
      break;
    case 'f':
      feat_fname = optarg;
      break;
    case 's':
      src_ref = get_arg_int(optarg);
      break;
    case 'd':
      dest_ref = get_arg_int(optarg);
      break;
    case 'i':
      format = msa_str_to_format(optarg);
      if (format == UNKNOWN_FORMAT) die("ERROR: bad alignment format.\n");
      break;
    case 'p':
      offset = get_arg_int(optarg);
      break;
    case 'n':
      offset = -1 * get_arg_int(optarg);
      break;
    case 'h':
      print_usage();
      exit(1);
    case '?':
      print_usage();
      exit(1);
    }
  }

  if (msa_fname == NULL || feat_fname == NULL) {
    print_usage();
    exit(1);
  }

  set_seed(-1);

  F = phast_fopen(feat_fname, "r");
  if ((gff = gff_read_set(F)) == NULL) { 
    die("ERROR: error reading %s.\n", feat_fname);
  }
  phast_fclose(F);

  /* handle case of local alignment specially -- avoid representing
     the alignment explicitly */
  F = phast_fopen(msa_fname, "r");
  if (format == UNKNOWN_FORMAT)
    format = msa_format_for_content(F, 1);
  if (format == LAV) {
    LocalPwAlignment *lpwa = NULL;
/*     int i; */

    fprintf(stderr, "WARNING: in local alignment mode, coordinates may only be mapped from query (reference) sequence to target (aligned) sequence.\n"); 

    lpwa = la_read_lav(F, 0);
    la_gff_transform(lpwa, gff);
/*     for (i = 0; i < lst_size(gff->features); i++) { */
/*       GFF_Feature *feat = lst_get_ptr(gff->features, i); */
/*       feat->start = la_get_target_coord(lpwa, feat->start); */
/*       feat->end = la_get_target_coord(lpwa, feat->end); */
/*     } */
  }

  else {                        /* normal alignment */
    msa = msa_new_from_file_define_format(F, format, NULL);
    phast_fclose(F);

    msa_map_gff_coords(msa, gff, src_ref, dest_ref, offset);
    msa_free(msa);
  }

  gff_print_set(stdout, gff);

  gff_free_set(gff);

  return 0;
}
