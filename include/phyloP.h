/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2009 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: phyloP.h,v 1.10 2008-11-12 02:07:59 acs Exp $ */

/* Functions that output data computed by phyloP */

#ifndef PHYLOP_H
#define PHYLOP_H

void print_prior_only(int nsites, char *mod_fname, Vector *prior_distrib);
void print_post_only(char *mod_fname, char *msa_fname, Vector *post_distrib,
                     double ci, double scale);
void print_p(char *mod_fname, char *msa_fname, Vector *prior_distrib,
             double post_mean, double post_var, double ci, double scale);
void print_prior_only_joint(char *node_name, int nsites, char *mod_fname, 
                            Matrix *prior_distrib);
void print_post_only_joint(char *node_name, char *mod_fname, 
                           char *msa_fname, Matrix *post_distrib, 
                           double ci, double scale, double sub_scale);
void print_p_joint(char *node_name, char *mod_fname, char *msa_fname, 
                   double ci, Matrix *prior_joint, 
                   double post_mean, double post_var, 
                   double post_mean_sup, double post_var_sup, 
                   double post_mean_sub, double post_var_sub,
                   double scale, double sub_scale);
void print_feats_sph(p_value_stats *stats, GFF_Set *gff, 
                     mode_type mode, double epsilon, int output_gff);
void print_feats_sph_subtree(p_value_joint_stats *stats, GFF_Set *gff, 
                             mode_type mode, double epsilon, int output_gff);
void print_quantiles(Vector *distrib);
void print_wig(MSA *msa, double *tuple_pvals, char *chrom, int refidx,
               int log_trans);
void print_base_by_base(char *header, char *chrom, MSA *msa, 
                        char **formatstr, int refidx, int ncols, ...);
void print_feats_generic(char *header, GFF_Set *gff, char **formatstr, 
                         int ncols, ...);
void print_gff_scores(GFF_Set *gff, double *pvals, int log_trans);

#endif
