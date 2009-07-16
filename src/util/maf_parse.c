/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2009 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: maf_parse.c,v 1.42 2009-01-09 22:01:00 mt269 Exp $ */

#include <stdlib.h>
#include <stdio.h>
#include <msa.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include <string.h>
#include <ctype.h>
#include <sufficient_stats.h>
#include <local_alignment.h>
#include <maf.h>
#include <maf_block.h>

void print_usage() {
    printf("\n\
USAGE: maf_parse [OPTIONS] <infile>\n\
\n\
DESCRIPTION:\n\
\n\
    Reads a MAF file and perform various operations on it.\n\
    Performs parsing operations block-by-block whenever possible,\n\
    rather than storing entire alignment in memory.  \n\
    Can extract a sub-alignment from an alignment (by row\n \
    or by column).  Can extract features given GFF, BED, or \n\
    genepred file.  Can also extract sub-features such as CDS1,2,3\n\
    or 4d sites.  Can perform various functions such as gap\n\
    stripping or re-ordering of sequences.  Capable of reading and\n\
    writing in a few common formats, but will not load input or output\n\
    alignments into memory if output format is MAF.\n\
\n\
OPTIONS:\n\
\n\
 (Obtaining sub-alignments and re-ordering rows)\n\
    --start, -s <start_col>\n\
        Starting column of sub-alignment (indexing starts with 1).\n\
        Default is 1.\n\
\n\
    --end, -e <end_col>\n\
        Ending column of sub-alignment.  Default is length of\n\
        alignment.\n\
\n\
    --seqs, -l <seq_list>\n\
        Comma-separated list of sequences to include (default)\n\
        exclude (if --exclude).  Indicate by sequence number or name\n\
        (numbering starts with 1 and is evaluated *after* --order is\n\
        applied).\n\
\n\
    --exclude, -x\n\
        Exclude rather than include specified sequences.\n\
\n\
    --order, -O <name_list>\n\
        Change order of rows in alignment to match sequence names\n\
        specified in name_list.  The first name in the alignment becomes\n\
        the reference sequence.\n\
\n\
    --no-refseq, -n\n\
        Do not assume first sequence in MAF is refseq.  Instead, use coordinates\n\
        given by absolute position in alignment (starting from 1).\n\
\n\
(Splitting into multiple MAFs by length)\n\
    --split, -S length \n\
        Split MAF into pieces by length, and puts output in outRootX.maf, where\n\
        X=1,2,...,numPieces.  outRoot can be modified with --out-root, and the\n\
        minimum number of digits in X can be modified with --out-root-digits.\n\
        Splits between blocks, so that each output file does not exceed specified\n\
        length.  By default, length is counted by distance spanned in alignment\n\
        by refseq, unless --no-refseq is specified.\n\
\n\
   --out-root, -r <name>\n\
        Filename root for output files produced by --split (default \"maf_parse\").\n\
\n\
   --out-root-digits, -d <numdigits>\n\
        (for use with --split).  The minimum number of digits used to index each\n\
        output file produced by split.\n\
\n\
(Extracting features from MAF)\n\
    --features, -g <fname>\n\
        Annotations file.  May be GFF, BED, or genepred format.  Coordinates assumed\n\
        to be in frame of first sequence of alignment (reference sequence).  By\n\
        default, outputs subset of MAF which are labeled in annotations file.  But\n\
        can be used with --by-category, --by-group, and/or --do-cats to split MAF\n\
        by annotation type.  Implies --strip-i-lines, --strip-e-lines\n\
    --by-category, -L\n\
        (Requires --features).  Split by category, as defined by annotations file\n\
        and (optionally) category map (see --catmap). (NOT IMPLEMENTED YET)\n\
    --do-cats, -C <cat_list>\n\
        (For use with --by-category) Output sub-alignments for only the specified\n\
        categories. (NOT IMPLEMENTED YET)\n\
    --catmap, -c <fname>|<string>\n\
        (Optionally use with --by-category) Mapping of feature types to category\n\
        numbers.  Can either give a filename or an \"inline\" description of a\n\
        simple category map, e.g.,\n --catmap \"NCATS = 3 ; CDS 1-3\" or\n\
        --catmap \"NCATS = 1; UTR 1\". (NOT IMPLEMENTED YET)\n\
    --by-group, -P <tag>\n\
        (Requires --features).  Split by groups in annotation file, as defined\n\
        by specified tag. (NOT IMPLEMENTED YET)\n\
\n\
(Masking by quality score)\n\
    --mask-bases, -b <qscore>\n\
        Mask all bases with quality score <= n.  Note that n is in the same units\n\
        as displayed in the MAF (ranging from 0-9), and represents\n\
        min(9, floor(PHRED_score/5)).  Bases without any quality score will not be\n\
        masked.\n\
\n\
 (Other)\n\
    --strip-i-lines, -I\n\
        Remove lines in MAF starting with i.\n\
    --strip-e-lines, -E\n\
        Remove lines in MAF starting with e.\n\
    --help, -h\n\
        Print this help message.\n\n");
}

FILE *get_outfile(List *outfileList, Hashtable *outfileHash, String *name, char *out_root) {
  int idx;
  FILE *outfile;
  char *fname = malloc((strlen(out_root)+name->length+2)*sizeof(char));
  sprintf(fname, "%s.%s", out_root, name->chars);
  idx = ptr_to_int(hsh_get(outfileHash, fname));
  if (idx == -1) {
    hsh_put(outfileHash, fname, int_to_ptr(lst_size(outfileList)));
    outfile = mafBlock_open_file(fname);
    lst_push_ptr(outfileList, (void*)outfile);
    return outfile;
  }
  return (FILE*)lst_get_ptr(outfileList, idx);
}

void close_outfiles(List *outfileList, Hashtable *outfileHash) {
  int i;
  FILE *f;
  for (i=0; i<lst_size(outfileList); i++) {
    f = (FILE*)lst_get_ptr(outfileList, i);
    mafBlock_close_file(f);
  }
  lst_free(outfileList);
  hsh_free(outfileHash);
}


int main(int argc, char* argv[]) {
  char *maf_fname = NULL, *out_root_fname = "maf_parse";
  String *refseq = NULL, *currRefseq;
  int opt_idx, startcol = 1, endcol = -1, include = 1, splitInterval = -1;
  char c, outfilename[1000], splitFormat[100]="%s%.1i.maf", *group_tag = NULL;
  List *order_list = NULL, *seqlist_str = NULL, *cats_to_do_str=NULL, *cats_to_do=NULL;
  MafBlock *block;
  FILE *mfile, *outfile=NULL;
  int useRefseq=TRUE, currLen=-1, blockIdx=0, currSize, sortWarned=0;
  int lastIdx = 0, currStart=0, by_category = FALSE, i;
  GFF_Set *gff = NULL, *gffSub;
  CategoryMap *cm = NULL;
  int base_mask_cutoff = -1, stripILines=FALSE, stripELines=FALSE;
  List *outfileList=NULL;
  Hashtable *outfileHash=NULL;

  struct option long_opts[] = {
    {"start", 1, 0, 's'},
    {"end", 1, 0, 'e'},
    {"seqs", 1, 0, 'l'},
    {"exclude", 0, 0, 'x'},
    {"order", 1, 0, 'O'},
    {"split", 1, 0, 'S'},
    {"out-root", 1, 0, 'r'},
    {"out-root-digits", 1, 0, 'd'},
    {"no-refseq", 0, 0, 'n'},
    {"features", 1, 0, 'g'},
    {"by-category", 0, 0, 'L'},
    {"do-cats", 1, 0, 'C'},
    {"catmap", 1, 0, 'c'},
    {"by-group", 1, 0, 'P'},
    {"mask-bases", 1, 0, 'b'},
    {"strip-i-lines", 0, 0, 'I'},
    {"strip-e-lines", 0, 0, 'E'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "s:e:l:O:r:S:d:g:c:P:b:LnxEIh", long_opts, &opt_idx)) != -1) {
    switch(c) {
    case 's':
      startcol = get_arg_int(optarg);
      break;
    case 'e':
      endcol = get_arg_int(optarg);
      break;
    case 'l':
      seqlist_str = get_arg_list(optarg);
      break;
    case 'O':
      order_list = get_arg_list(optarg);
      break;
    case 'x':
      include = FALSE;
      break;
    case 'S':
      splitInterval = atoi(optarg);
      break;
    case 'r':
      out_root_fname = optarg;
      break;
    case 'd':
      sprintf(splitFormat, "%%s%%.%si.maf", optarg);
      break;
    case 'n':
      useRefseq = FALSE;
      break;
    case 'g':
      gff = gff_read_set(fopen_fname(optarg, "r"));
      stripILines=TRUE;
      stripELines=TRUE;
      break;
    case 'c':
      cm = cm_new_string_or_file(optarg);
      break;
    case 'C':
      cats_to_do_str = get_arg_list(optarg);
      break;
    case 'L':
      by_category = TRUE;
      break;
    case 'P':
      group_tag = optarg;
      break;
    case 'b':
      base_mask_cutoff = atoi(optarg);
      break;
    case 'E':
      stripELines=TRUE;
      break;
    case 'I':
      stripILines=TRUE;
      break;
    case 'h':
      print_usage();
      exit(0);
    case '?':
      fprintf(stderr, "Bad argument.  Try 'maf_parse -h' for help.\n");
      exit(1); 
    }
  }

  if (optind >= argc) 
    die("Missing alignment filename.  Try 'maf_parse -h' for help.\n");
  else if (optind == argc - 1) 
    maf_fname = argv[optind];
  else 
    die("ERROR: Too many arguments.  Try 'maf_parse -h' for help.\n");

  if (startcol < 1 || (endcol != -1 && endcol < startcol))
    die("ERROR: must have 1 <= start <= end <= [msa_length]\n");

  if ((group_tag != NULL || by_category) && gff == NULL)
    die("ERROR: --by-category and --by-group require --features.  Try \"maf_parse -h\""
	" for help.\n");

  if (group_tag != NULL && by_category) 
    die("ERROR: --by-category and --by-group cannot be used together.  Try \"maf_parse -h\""
	" for help.\n");
  
  if (splitInterval != -1 && gff != NULL)
    die("ERROR: can't use --split and --features together.  Try \"maf_parse -h\""
	"for help\n");

  if (group_tag != NULL || by_category) {
    outfileList = lst_new_ptr(10);
    outfileHash = hsh_new(100);
  }

  if (gff != NULL && cm == NULL) 
    cm = cm_new_from_features(gff);

  if (cats_to_do_str != NULL) {
    cats_to_do = cm_get_category_list(cm, cats_to_do_str, FALSE);
    if (gff != NULL) 
      gff_filter_by_type(gff, cats_to_do, 0, NULL);
  }

  /* Check to see if --do-cats names a feature which is length 1. 
     If so, set output_format to SS ? or FASTA ? */
  
  mfile = fopen_fname(maf_fname, "r");
  block = mafBlock_read_next(mfile, NULL, NULL);

  if (splitInterval == -1 && gff==NULL) {
    //TODO: do we want to copy header from original MAF in this case?
    mafBlock_open_file(NULL);
  }

  while (block != NULL) {
    if (order_list != NULL)
      mafBlock_reorder(block, order_list);
    if (seqlist_str != NULL)
      mafBlock_subSpec(block, seqlist_str, include);
    if (mafBlock_numSpec(block)==0 || mafBlock_all_gaps(block)) 
      goto get_next_block;
    if (stripILines)
      mafBlock_strip_iLines(block);
    if (stripELines)
      mafBlock_strip_eLines(block);
    if (base_mask_cutoff != -1)
      mafBlock_mask_bases(block, base_mask_cutoff);

    if (useRefseq) {  //get refseq and check that it is consistent in MAF file
      currRefseq = mafBlock_get_refSpec(block);
      if (refseq == NULL) 
	refseq = str_new_charstr(currRefseq->chars);
      else if (str_compare(refseq, currRefseq)!=0)
	die("Error: refseq not consistent in MAF (got %s, %s)\n",
	    refseq->chars, currRefseq->chars);
    }
    
    if (startcol != 1 || endcol != -1) 
      if (0 == mafBlock_trim(block, startcol, endcol, refseq, useRefseq ? 0 : lastIdx))
	goto get_next_block;

    currSize = mafBlock_get_size(block, refseq);
    if (useRefseq) {
      currStart = mafBlock_get_start(block, refseq);
      if (currStart < lastIdx && sortWarned == 0) {
	fprintf(stderr, "Warning: input MAF not sorted with respect to refseq.  Output files may not represent contiguous alignments. (%i, %i)\n", lastIdx, currStart);
	sortWarned = 1;
      }
    }
    else currStart = lastIdx;
    
    lastIdx = currStart + currSize;

    //split by length
    if (splitInterval != -1) {
      if (currLen == -1 || currLen+currSize > splitInterval) {
	sprintf(outfilename, splitFormat, out_root_fname, ++blockIdx);
	if (outfile != NULL) {
	  mafBlock_close_file(outfile);
	}
	outfile = mafBlock_open_file(outfilename);
	currLen = 0;
      }
      currLen += currSize;
    }
    else outfile = stdout;

    if (gff != NULL) {
      gffSub = gff_subset_range_overlap(gff, currStart, lastIdx);
      if (lst_size(gffSub->features) != 0) {
	if (by_category) gff_group_by_feature(gffSub);
	else if (group_tag != NULL) gff_group(gffSub, group_tag);
	gff_sort(gffSub);
	gff_flatten_within_groups(gffSub);
	for (i=0; i<lst_size(gffSub->features); i++) {
	  GFF_Feature *feat = (GFF_Feature*)lst_get_ptr(gffSub->features, i);
	  MafBlock *subBlock = mafBlock_copy(block);
	  mafBlock_trim(subBlock, feat->start, feat->end, refseq, 0);
	  if (by_category) 
	    outfile = get_outfile(outfileList, outfileHash, feat->feature, out_root_fname);
	  else if (group_tag != NULL) 
	    outfile = get_outfile(outfileList, outfileHash, 
				  gff_group_name(gffSub, feat), out_root_fname);
	  else outfile = stdout;
	  mafBlock_print(outfile, subBlock);
	  mafBlock_free(subBlock);
	}
      }
      gff_free_set(gffSub);
    } 
    else mafBlock_print(outfile, block);

  get_next_block:
    mafBlock_free(block);
    block = mafBlock_read_next(mfile, NULL, NULL);
  }
  if (by_category || group_tag != NULL)
    close_outfiles(outfileList, outfileHash);
  else mafBlock_close_file(outfile);

  fclose(mfile);
  return 0;
}
