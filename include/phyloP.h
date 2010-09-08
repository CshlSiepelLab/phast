/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: phyloP.h,v 1.10 2008-11-12 02:07:59 acs Exp $ */

/* Functions that output data computed by phyloP */

#ifndef PHYLOP_H
#define PHYLOP_H

#include <list_of_lists.h>

void print_prior_only(FILE *outfile, int nsites, char *mod_fname, 
		      Vector *prior_distrib, ListOfLists *result);
void print_post_only(FILE *outfile, char *mod_fname, char *msa_fname, 
		     Vector *post_distrib, double ci, double scale,
		     ListOfLists *result);
void print_p(FILE *outfile, char *mod_fname, char *msa_fname, 
	     Vector *prior_distrib, double post_mean, 
	     double post_var, double ci, double scale,
	     ListOfLists *result);
void print_prior_only_joint(FILE *outfile, char *node_name, 
			    int nsites, char *mod_fname, 
                            Matrix *prior_distrib,
			    ListOfLists *result);
void print_post_only_joint(FILE *outfile, char *node_name, char *mod_fname, 
                           char *msa_fname, Matrix *post_distrib, 
                           double ci, double scale, double sub_scale,
			   ListOfLists *result);
void print_p_joint(FILE *outfile, char *node_name, char *mod_fname, 
		   char *msa_fname, 
                   double ci, Matrix *prior_joint, 
                   double post_mean, double post_var, 
                   double post_mean_sup, double post_var_sup, 
                   double post_mean_sub, double post_var_sub,
                   double scale, double sub_scale,
		   ListOfLists *result);
void print_feats_sph(FILE *outfile, p_value_stats *stats, GFF_Set *gff, 
                     mode_type mode, double epsilon, int output_gff,
		     ListOfLists *result);
void print_feats_sph_subtree(FILE *outfile, p_value_joint_stats *stats, 
			     GFF_Set *gff, mode_type mode, double epsilon, 
			     int output_gff, ListOfLists *result);
void print_quantiles(FILE *outfile, Vector *distrib, ListOfLists *result);
void print_wig(FILE *outfile, MSA *msa, double *tuple_pvals, char *chrom, 
	       int refidx, int log_trans, ListOfLists *result);
void print_base_by_base(FILE *outfile, char *header, char *chrom, MSA *msa, 
                        char **formatstr, int refidx, ListOfLists *result,
			int log_trans_outfile, int log_trans_results, int ncols, ...);
void print_feats_generic(FILE *outfile, char *header, GFF_Set *gff, 
			 char **formatstr, ListOfLists *result, 
			 int log_trans_outfile, int log_trans_results, int ncols, ...);
void print_gff_scores(FILE *outfile, GFF_Set *gff, double *pvals, 
		      int log_trans);

#endif
