/* $Id: convert_coords.c,v 1.1.1.1 2004-06-03 22:43:12 acs Exp $
   Written by Adam Siepel, 2002
   Copyright 2002, Adam Siepel, University of California */

#include <stdlib.h>
#include <stdio.h>
#include <msa.h>
#include <gff.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include <local_alignment.h>

void print_usage() {
  fprintf(stderr, "USAGE: convert_coords -m <msa_fname> -f <feature_fname> [-s <src_frame>] [-d <dest_frame>] [-p] [-n] [-i PHYLIP|FASTA|PSU]\n\
\n\
Converts coordinates of features in a GFF file according to a multiple\n\
alignment.  Will map from the coordinate system of any sequence to any\n\
other sequence; can also map to or from the coordinate system of the\n\
entire alignment.  In addition, supports translation of coordinates by\n\
specified offset.\n\
\n\
Options:\n\
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
    -i PHYLIP|FASTA|PSU|LAV\n\
                    (default PHYLIP) Alignment format.  PSU is the\n\
                    format used by many of the\n\
                    tools developed at Penn. State University.  LAV is\n\
                    the format used by BLASTZ to represent local\n\
                    pairwise alignments.  If it is selected, such an\n\
                    alignment will be treated like a global alignment,\n\
                    with unaligned portions of the target sequence\n\
                    replaced by gaps (but the alignment will never be\n\
                    explicitly represented).\n\n");  
} 

int main(int argc, char* argv[]) {
  FILE* F;
  MSA *msa;
  msa_format_type format = PHYLIP;
  int src_ref = -1, dest_ref = 0, offset = 0;
  char *msa_fname = NULL, *feat_fname = NULL;
  GFF_Set *gff;
  char c;

  while ((c = getopt(argc, argv, "m:f:s:d:i:p:n:")) != -1) {
    switch(c) {
    case 'm':
      msa_fname = optarg;
      break;
    case 'f':
      feat_fname = optarg;
      break;
    case 's':
      src_ref = atoi(optarg);
      break;
    case 'd':
      dest_ref = atoi(optarg);
      break;
    case 'i':
      if (!strcmp(optarg, "PSU")) format = PSU;
      else if (!strcmp(optarg, "FASTA")) format = FASTA; 
      else if (!strcmp(optarg, "LAV")) format = LAV;
      else if (strcmp(optarg, "PHYLIP") != 0) { 
        fprintf(stderr, "ERROR: format must be \"PHYLIP,\" \"FASTA,\" or \"PSU.\"\n");
        print_usage(); 
        exit(1); 
      }
      break;
    case 'p':
      offset = atoi(optarg);
      break;
    case 'n':
      offset = -1 * atoi(optarg);
      break;
    case '?':
      print_usage();
      exit(1);
    }
  }

  if (msa_fname == NULL || feat_fname == NULL) {
    print_usage();
    exit(1);
  }

  if ((F = fopen(feat_fname, "r")) == NULL) {
    fprintf(stderr, "ERROR: cannot open %s.\n", feat_fname);
    exit(1);
  }
  if ((gff = gff_read_set(F)) == NULL) { 
    fprintf(stderr, "ERROR: error reading %s.\n", feat_fname);
    exit(1);
  }
  fclose(F);

  /* handle case of local alignment specially -- avoid representing
     the alignment explicitly */
  if (format == LAV) {
    LocalPwAlignment *lpwa = NULL;
/*     int i; */

    fprintf(stderr, "WARNING: in local alignment mode, coordinates may only be mapped from query (reference) sequence to target (aligned) sequence.\n"); 

    if ((F = fopen(msa_fname, "r")) == NULL) {
      fprintf(stderr, "ERROR: cannot open %s.\n", msa_fname);
      exit(1);
    }
    lpwa = la_read_lav(F, 0);
    la_gff_transform(lpwa, gff);
/*     for (i = 0; i < lst_size(gff->features); i++) { */
/*       GFF_Feature *feat = lst_get_ptr(gff->features, i); */
/*       feat->start = la_get_target_coord(lpwa, feat->start); */
/*       feat->end = la_get_target_coord(lpwa, feat->end); */
/*     } */
  }

  else {                        /* normal alignment */
    if ((F = fopen(msa_fname, "r")) == NULL) {
      fprintf(stderr, "ERROR: cannot open %s.\n", msa_fname);
      exit(1);
    }
    msa = msa_new_from_file(F, format, NULL);
    fclose(F);

    msa_map_gff_coords(msa, gff, src_ref, dest_ref, offset, NULL);
    msa_free(msa);
  }

  gff_print_set(stdout, gff);

  gff_free_set(gff);

  return 0;
}
