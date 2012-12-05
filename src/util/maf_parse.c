/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: maf_parse.c,v 1.42 2009-01-09 22:01:00 mt269 Exp $ */

#include <stdlib.h>
#include <stdio.h>
#include <msa.h>
#include <getopt.h>
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
 (Output format)\n\
    --out-format, -o MAF|PHYLIP|FASTA|MPM|SS\n\
        (Default MAF).  Output file format.  SS format is only\n\
        available un-ordered.  Note that some options, which involve\n\
        reversing alignments based on strand, or stripping gaps, \n\
        cannot be output in MAF format and use FASTA by default.\n\
        Also note that when output format is not MAF, the entire\n\
        output must be loaded into memory.\n\
\n\
    --pretty, -p\n\
        Pretty-print alignment (use '.' when character matches\n\
        corresponding character in first sequence).  Ignored if\n\
        --out-format SS is selected.\n\
\n\
 (Obtaining sub-alignments and re-ordering rows)\n\
    --start, -s <start_col>\n\
        Start index of sub-alignment (indexing starts with 1).\n\
        Coordinates are in terms of the reference sequence unless\n\
        the --no-refseq option is used, in which case they are in\n\
        terms of alignment columns.  Default is 1.\n\
\n\
    --end, -e <end_col>\n\
        End index of sub-alignment.  Default is length of alignment.\n\
        Coordinates defined as in --start option, above.\n\
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
        Do not assume first sequence in MAF is refseq.  Instead, use\n\
        coordinates  given by absolute position in alignment (starting\n\
        from 1).\n\
\n\
(Splitting into multiple MAFs by length)\n\
    --split, -S length \n\
        Split MAF into pieces by length, and puts output in \n\
        outRootX.maf, where X=1,2,...,numPieces.  outRoot can be\n\
        modified with --out-root, and the minimum number of digits in X\n\
        can be modified with --out-root-digits.\n\
        Splits between blocks, so that each output file does not exceed\n\
        specified length.  By default, length is counted by distance\n\
        spanned in alignment by refseq, unless --no-refseq is specified.\n\
\n\
   --out-root, -r <name>\n\
        Filename root for output files produced by --split (default\n\
        \"maf_parse\").\n\
\n\
   --out-root-digits, -d <numdigits>\n\
        (for use with --split).  The minimum number of digits used to \n\
        index each output file produced by split.\n\
\n\
(Extracting features from MAF)\n\
    --features, -g <fname>\n\
        Annotations file.  May be GFF, BED, or genepred format.\n\
        Coordinates assumed to be in frame of first sequence of\n\
        alignment (reference sequence).  By default, outputs subset of \n\
        MAF which are labeled in annotations file.  But can be used with\n\
        --by-category, --by-group, and/or --do-cats to split MAF by\n\
        annotation type.  Or if used with --mask-features, is only used\n\
        to determine regions to mask.  Implies --strip-i-lines, \n\
        --strip-e-lines\n\
\n\
    --by-category, -L\n\
        (Requires --features).  Split by category, as defined by\n\
        annotations file and (optionally) category map (see --catmap).\n\
\n\
    --do-cats, -C <cat_list>\n\
        (For use with --by-category) Output sub-alignments for only the\n\
        specified categories.\n\
\n\
    --catmap, -c <fname>|<string>\n\
        (Optionally use with --by-category) Mapping of feature types to\n\
        category numbers.  Can either give a filename or an \"inline\"\n\
        description of a simple category map, e.g.,\n\
         --catmap \"NCATS = 3 ; CDS 1-3\" or\n\
         --catmap \"NCATS = 1; UTR 1\".\n\
\n\
    --by-group, -P <tag>\n\
        (Requires --features).  Split by groups in annotation file, as \n\
        defined by specified tag.\n\
\n\
(Masking by quality score)\n\
    --mask-bases, -b <qscore>\n\
        Mask all bases with quality score <= n.  Note that n is in the\n\
        same units as displayed in the MAF (ranging from 0-9), and\n\
        represents min(9, floor(PHRED_score/5)).  Bases without any\n\
        quality score will not be masked.\n\
\n\
    --masked-file, -m <filename>\n\
       (For use with --mask-bases).  Write a file containing all the\n\
       regions masked for low quality.  The file will be in 0-based\n\
       coordinates relative to the refseq, with an additional column\n\
       giving the name of the species masked.  Note that low-quality bases\n\
       masked at alignment columns with a gap in the reference sequence\n\
       may not be represented in the output file.\n\
\n\
    --mask-features -M <spec>\n\
      (Requires --features).  Mask all bases annotated in features in the\n\
      given species (can be a comma-delimited list of species).  Note that\n\
      coordinates are always in terms of refseq, even if a different species\n\
      is being masked.\n\
\n\
 (Other)\n\
    --strip-i-lines, -I\n\
        Remove lines in MAF starting with i.\n\
    --strip-e-lines, -E\n\
        Remove lines in MAF starting with e.\n\
    --help, -h\n\
        Print this help message.\n\n");
}



/*
  open a file with name out_root.name.maf, or returns it if already open.
  This is a bit messy because in some cases (splitting by feature) there may
  be more output files than the OS can handle.  But it would be computationally
  expensive to check and see which files are finished, assuming that the MAF is
  sorted.  

  So, if it tries to open a file and fails, it the goes through the list of
  filehandles, finds an open one, closes it, and tries to open the new one 
  again.  Repeat until successful.

  Then, if a filehandle needs to be re-opened, it is opened with append.  Again,
  if this is not successful, it looks for another file to close.  If it can't
  find one the program reports an error and dies.

  Finally, close_outfiles below checks and makes sure that all files
  are closed with mafBlock_close_file in the end, so that they get the #eof
  closer.
 */
FILE *get_outfile(List *outfileList, Hashtable *outfileHash, String *name, char *out_root,
		  int argc, char *argv[]) {
  int idx, i;
  FILE *outfile;
  char *fname = smalloc((strlen(out_root)+name->length+7)*sizeof(char));
  sprintf(fname, "%s.%s.maf", out_root, name->chars);
  idx = ptr_to_int(hsh_get(outfileHash, fname));
  if (idx == -1) {
    hsh_put(outfileHash, fname, int_to_ptr(lst_size(outfileList)));
    outfile = mafBlock_open_outfile(fname, argc, argv);
    while (outfile==NULL) {  //too many files are open, close one first
      for (i=0; i<lst_size(outfileList); i++) {
	outfile = (FILE*)lst_get_ptr(outfileList, i);
	if (outfile != NULL) break;
      }
      if (i == lst_size(outfileList)) {
	die("ERROR: too many files open in maf_parse\n");
      } else {
	phast_fclose(outfile);
	lst_set_ptr(outfileList, i, NULL);
      }
      outfile = mafBlock_open_outfile(fname, argc, argv);
    }
    lst_push_ptr(outfileList, (void*)outfile);
    sfree(fname);
    return outfile;
  }
  outfile = (FILE*)lst_get_ptr(outfileList, idx);
  if (outfile == NULL) { //has already been opened but then closed.
    outfile = phast_fopen_no_exit(fname, "a");
    while (outfile == NULL) {
      for (i=0; i<lst_size(outfileList); i++) {
	outfile = (FILE*)lst_get_ptr(outfileList, i);
	if (outfile != NULL) break;
      }
      if (i == lst_size(outfileList)) {
	die("ERROR: too many files open in maf_parse\n");
      } else {
	phast_fclose(outfile);
	lst_set_ptr(outfileList, i, NULL);
      }
      outfile = phast_fopen_no_exit(fname, "a");
    }
    lst_set_ptr(outfileList, idx, (void*)outfile);
  }
  sfree(fname);
  return outfile;
}


/* Closes all outfiles.  If already closed, reopen with append, add #eof 
   closer, and close again.  see comment above at get_outfile */
void close_outfiles(List *outfileList, Hashtable *outfileHash) {
  List *keys = hsh_keys(outfileHash);
  int *done, idx, i;
  char *fname;
  FILE *outfile;
  done = smalloc(lst_size(keys)*sizeof(int));
  for (i=0; i<lst_size(keys); i++) {
    done[i]=0;
    fname = (char*)lst_get_ptr(keys, i);
    idx = hsh_get_int(outfileHash, fname);
    outfile = (FILE*)lst_get_ptr(outfileList, idx);
    if (outfile != NULL) {
      mafBlock_close_outfile(outfile);
      done[i]=1;
    }
  }
  for (i=0; i<lst_size(keys); i++) {
    if (done[i]) continue;
    fname = (char*)lst_get_ptr(keys, i);
    outfile = phast_fopen(fname, "a");
    mafBlock_close_outfile(outfile);
  }
  sfree(done);
  lst_free(keys);
  lst_free(outfileList);
  hsh_free(outfileHash);
}


int main(int argc, char* argv[]) {
  char *maf_fname = NULL, *out_root_fname = "maf_parse", *masked_fn = NULL;
  String *refseq = NULL, *currRefseq;
  int opt_idx, startcol = 1, endcol = -1, include = 1, splitInterval = -1;
  char c, outfilename[1000], splitFormat[100]="%s%.1i.maf", *group_tag = NULL;
  List *order_list = NULL, *seqlist_str = NULL, *cats_to_do_str=NULL, *cats_to_do=NULL;
  MafBlock *block;
  FILE *mfile, *outfile=NULL, *masked_file=NULL;
  int useRefseq=TRUE, currLen=-1, blockIdx=0, currSize, sortWarned=0;
  int lastIdx = 0, currStart=0, by_category = FALSE, i, pretty_print = FALSE;
  int lastStart = -1, gffSearchIdx=0;
  GFF_Set *gff = NULL, *gffSub;
  GFF_Feature *feat;
  CategoryMap *cm = NULL;
  int base_mask_cutoff = -1, stripILines=FALSE, stripELines=FALSE;//, numspec=0;
  List *outfileList=NULL;
  Hashtable *outfileHash=NULL;//, *specNameHash=NULL;
  msa_format_type output_format = MAF;
  MSA *msa = NULL;//, **catMsa;
  char *mask_features_spec_arg=NULL;
  List *mask_features_spec=NULL;
  

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
    {"masked-file", 1, 0, 'm'},
    {"strip-i-lines", 0, 0, 'I'},
    {"strip-e-lines", 0, 0, 'E'},
    {"mask-features", 1, 0, 'M'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };


  while ((c = getopt_long(argc, argv, "s:e:l:O:r:S:d:g:c:P:b:o:m:M:pLnxEIh", long_opts, &opt_idx)) != -1) {
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
      sprintf(splitFormat, "%%s%%.%si.%%s", optarg);
      break;
    case 'n':
      useRefseq = FALSE;
      break;
    case 'g':
      gff = gff_read_set(phast_fopen(optarg, "r"));
      gff_sort(gff);
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
    case 'm':
      masked_fn = optarg;
      break;
    case 'M':
      mask_features_spec_arg = optarg;
      break;
    case 'E':
      stripELines=TRUE;
      break;
    case 'I':
      stripILines=TRUE;
      break;
    case 'o':
      output_format = msa_str_to_format(optarg);
      if (output_format == UNKNOWN_FORMAT) 
	die("ERROR: bad output format.  Try \"maf_parse -h\" for help.\n");
      if (output_format != MAF)
	die("Sorry, only MAF format output has been implemented right now.\n");
      break;
    case 'p':
      pretty_print = TRUE;
      break;
    case 'h':
      print_usage();
      exit(0);
    case '?':
      die("Bad argument.  Try 'maf_parse -h' for help.\n");
    }
  }

  if (optind >= argc) 
    die("Missing alignment filename.  Try 'maf_parse -h' for help.\n");
  else if (optind == argc - 1) 
    maf_fname = argv[optind];
  else 
    die("ERROR: Too many arguments.  Try 'maf_parse -h' for help.\n");
  
  set_seed(-1);

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
    cats_to_do = cm_get_category_str_list(cm, cats_to_do_str, FALSE);
    if (gff != NULL) 
      gff_filter_by_type(gff, cats_to_do, 0, NULL);
  }

  if (masked_fn != NULL) {
    if (base_mask_cutoff == -1)
      die("ERROR: need to use --mask-bases with --masked-file");
    masked_file = phast_fopen(masked_fn, "w");
  }

  if (mask_features_spec_arg != NULL) {
    if (gff==NULL)
      die("ERROR: need --features with --mask-features");
    mask_features_spec = lst_new_ptr(10);
    str_split(str_new_charstr(mask_features_spec_arg), ",", mask_features_spec);
    for (i=0; i < lst_size(mask_features_spec); i++) {
      fprintf(stderr, "masking species %s within features\n", 
	      ((String*)lst_get_ptr(mask_features_spec, i))->chars);
    }
  }

  /* Check to see if --do-cats names a feature which is length 1. 
     If so, set output_format to SS ? or FASTA ? */
  
  mfile = phast_fopen(maf_fname, "r");
  block = mafBlock_read_next(mfile, NULL, NULL);

  if (splitInterval == -1 && gff==NULL) {
    //TODO: do we want to copy header from original MAF in this case?
    mafBlock_open_outfile(NULL, argc, argv);
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
      mafBlock_mask_bases(block, base_mask_cutoff, masked_file);
    //TODO: still need to implement (either here or elsewhere)
    //    if (indel_mask_cutoff != -1) 
    //      mafBlock_mask_indels(block, indel_mask_cutoff, mfile);

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

    if (currStart < lastStart) gffSearchIdx = 0;
    lastStart = currStart;
    
    lastIdx = currStart + currSize;

    //split by length
    if (splitInterval != -1) {
      if (currLen == -1 || currLen+currSize > splitInterval) {
	sprintf(outfilename, splitFormat, out_root_fname, ++blockIdx,
		msa_suffix_for_format(output_format));
	if (output_format == MAF) {
	  if (outfile != NULL) mafBlock_close_outfile(outfile);
	  outfile = mafBlock_open_outfile(outfilename, argc, argv);
	}
	else if (output_format != MAF && msa != NULL) {
	  //	  msa_print_to_filename(msa, outfilename, output_format, pretty_print);
	  msa_free(msa);
	  msa = NULL;
	}
	currLen = 0;
      }
      currLen += currSize;
    }
    else outfile = stdout;
    if (gff != NULL && mask_features_spec != NULL) {
      gffSub = gff_subset_range_overlap_sorted(gff, currStart+1, lastIdx,
					       &gffSearchIdx);
      if (gffSub != NULL) {
	mafBlock_mask_region(block, gffSub, mask_features_spec);
	gff_free_set(gffSub);
      }
      mafBlock_print(outfile, block, pretty_print);


    } else if (gff != NULL) {
      gffSub = gff_subset_range_overlap_sorted(gff, currStart+1, lastIdx, 
					       &gffSearchIdx);
      if (gffSub != NULL) {
	if (by_category) gff_group_by_feature(gffSub);
	else if (group_tag != NULL) gff_group(gffSub, group_tag);
	gff_sort(gffSub);
	gff_flatten_within_groups(gffSub);
	for (i=0; i<lst_size(gffSub->features); i++) {
	  feat = (GFF_Feature*)lst_get_ptr(gffSub->features, i);
	  MafBlock *subBlock = mafBlock_copy(block);
	  mafBlock_trim(subBlock, feat->start, feat->end, refseq, 0);
	  if (by_category) 
	    outfile = get_outfile(outfileList, outfileHash, feat->feature, out_root_fname,
				  argc, argv);
	  else if (group_tag != NULL) 
	    outfile = get_outfile(outfileList, outfileHash, 
				  gff_group_name(gffSub, feat), out_root_fname,
				  argc, argv);
	  else outfile = stdout;
	  if (output_format == MAF)
	    mafBlock_print(outfile, subBlock, pretty_print);
	  //	  else msa_add_mafBlock(msa);
	  mafBlock_free(subBlock);
	}
	gff_free_set(gffSub);
      }
    }
    else {
      if (output_format == MAF) 
	mafBlock_print(outfile, block, pretty_print);
      //      else msa = msa_add_mafBlock(mafBlock, msa, );
    }
    
  get_next_block:
    mafBlock_free(block);
    block = mafBlock_read_next(mfile, NULL, NULL);
  }

  if (masked_file != NULL) fclose(masked_file);

  if (output_format == MAF) {
    if (by_category || group_tag != NULL)
      close_outfiles(outfileList, outfileHash);
    else if (outfile!=NULL) mafBlock_close_outfile(outfile);
  } else {
    msa_print(stdout, msa, output_format, pretty_print);
    msa_free(msa);
  }
  if (gff != NULL) gff_free_set(gff);
  phast_fclose(mfile);
  return 0;
}
