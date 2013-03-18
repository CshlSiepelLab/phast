/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** \file msa.h
   Reading and writing of sequence data to Multiple Sequence Alignment (MSA) file types (fasta, phylip, etc) and manipulating data:  Reverse complement, handle missing/gap data, identify informative sites
   \ingroup msa
*/

#ifndef MSA_H
#define MSA_H

#include <stdio.h>
#include <gff.h>
#include <category_map.h>
#include <markov_matrix.h>
#include <vector.h>
#include <hashtable.h>

/** Maximum "order" allowed when indexing columns (see
   msa_index_cols) */
#define MSA_MAX_ORDER 5         

/** Using an incomplete type in lieu of
  including "sufficient_stats.h"
  avoids problems caused by
  reciprocal dependencies in header
  files */
struct msa_ss_struct;

/** Multiple sequence alignment object */
typedef struct {
  int nseqs;                    /**< Number of sequences */
  unsigned int length;          /**< Number of columns */
  char *alphabet;               /**< Alphabet (see #DEFAULT_ALPHABET) */
  int inv_alphabet[NCHARS];     /**< Inverse of 'alphabet' (maps characters to their index in 'alphabet') */
  char *missing;                /**< Recognized missing data characters */
  int is_missing[NCHARS];       /**< Fast lookup of whether character
                                   is missing data char */
  char **names;			/**< Names of each sequence (2D char array) */
  char **seqs;			/**< Sequence data (2D char array) */
  int *categories;		/**< Categories for each coordinate */
  struct msa_ss_struct *ss;	/**< Sufficient Statistics */
  int ncats;			/**< Number of categories */
  int alloc_len;		/**< Length of memory allocated for sequence */
  int idx_offset;		/**< Index offset */
  int *is_informative;          /**< If non-NULL, indicates which
                                   sequences are to be considered
                                   "informative", e.g., for
                                   phylogenetic analysis */
} MSA;

/** Size of lookup tables */
#define NCHARS 256
/** Maximum line length to read in from a file at a time */
#define MAX_LINE_LEN 10000
/** Length of line to write out at a time */
#define OUTPUT_LINE_LEN 70
/** Default alphabet, assumed throughout PHAST */
#define DEFAULT_ALPHABET "ACGT"
/** Gap character, assumed throughout PHAST */
#define GAP_CHAR '-'
/** Default missing data character */
#define DEFAULT_MDATA_CHARS "*N"

/** Format types */
typedef enum {PHYLIP,           /**< PHYLIP format  */
              MPM,              /**< Format used by MultiPipMaker and
                                   some of Webb Miller's older
                                   tools  */
              FASTA,            /**< Standard FASTA format */
              SS,               /**< "Sufficient statistics" format,
                                   in which each unique alignment
                                   column (or tuple of columns) is
                                   represented only once, and a count
                                   is maintained of how many times it
                                   occurs */
              LAV,              /**< lav format, used by BLASTZ */
              MAF,              /**< Multiple Alignment Format (MAF)
				    used by MULTIZ and TBA  */
	      UNKNOWN_FORMAT    /**< Format unknown */
} msa_format_type; 

/** Strip all gaps (used with msa_strip_gaps) */
#define STRIP_ALL_GAPS -2 
/** Strip any gaps (used with msa_strip_gaps) */
#define STRIP_ANY_GAPS -1
/** Don't strip gaps (used with msa_strip_gaps) */
#define NO_STRIP 0


/** Coordinate map, defined by a sequence/alignment pair.  Allows fast
    conversion between the coordinate frame of a multiple alignment
    and the coordinate frame of one of the sequences in the
    alignment.  */
typedef struct {
  List *seq_list;               /**< list of indexes in sequence
                                   immediately following gaps,
                                   expressed in the coordinate frame
                                   of the sequence (starting with
                                   position 1) */
  List *msa_list;               /**< list of corresponding indices in
                                   the frame of the MSA */
  int seq_len;                  /**< length of sequence */
  int msa_len;                  /**< length of alignment */
} msa_coord_map;


/** \name MSA allocation functions
 \{ */
/** Creates a new MSA object.  
   @param seqs Two-dimensional character array to hold sequences
   @param names Two-dimensional character array to hold names
   @param nseqs Number of sequences
   @param length Length of sequences
   @param alphabet (Optional) Possible chars in sequences, if NULL then default alphabet is used
   @result Newly allocated MSA with provided sequences
   @note No new memory is allocated for seqs or names 
 */
MSA *msa_new(char **seqs, char **names, int nseqs, int length, 
             char *alphabet);

/** Re-allocate if sequence length increases
If the number of sequences increases, call msa_add_seq.
    @param msa MSA to reallocate
    @param new_len New length of sequences
    @param new_alloc_len unused
    @param do_cats reallocate to accommodate more categories
    @param store_order reallocate to accommodate more tuple orders
 */
void msa_realloc(MSA *msa, int new_length, int new_alloclen, int do_cats, int store_order);

/** Create a copy of an MSA.  
    @param msa MSA to copy data from
    @param suff_stats_only If set sequences aren't copied
    @result Newly allocated MSA with contents of provided file
   aren't copied 
*/
MSA *msa_create_copy(MSA *msa, int suff_stats_only);


/** \} \name MSA read/write file access functions
 \{ */

/** Creates a new alignment from the contents of the specified file and file type.
   @param F File descriptor of file containing alignment data
   @param format File format type
   @param alphabet (Optional) Chars allowed in sequences, if NULL default alphabet for DNA will be used.
   @result Newly allocated MSA with contents of provided file
   @todo { Does not handle interleaved PHYLIP files } 
*/
MSA *msa_new_from_file_define_format(FILE *F, msa_format_type format, char *alphabet);

/** Creates a new alignment from the contents of the specified file.
   @param F File descriptor of file containing alignment data
   @param alphabet (Optional) Chars allowed in sequences, if NULL default alphabet for DNA will be used.
   @result Newly allocated MSA with contents of provided file
   @note Does not handle MAF files or interleaved PHYLIP files.
   @see maf_read for MAF file reading
   @todo { Does not handle interleaved PHYLIP files } 
*/
MSA *msa_new_from_file(FILE *F, char *alphabet);

/** Create a single MSA from multiple files.
    @param msa_fname_list List of filenames containing MSA data
    @param seqnames List of sequence names to define order of sequences in combined MSA
    @param alphabet (Optional) Alphabet used in MSA files, if NULL Default Alphabet is used
    @result Newly created MSA populated with MSA data from multiple files    
 */
MSA *msa_concat_from_files(List *msa_fname_list, 
                           List *seqnames, char *alphabet);

/** Create a new alignment from a FASTA file.
   @param F File descriptor of FASTA file
   @param alphabet (Optional) Chars allowed in sequences, if NULL default alphabet for DNA will be used.
   @result Newly allocated MSA with contents of provided FASTA file
*/
MSA *msa_read_fasta(FILE *F, char *alphabet);

/** Read and return a single sequence from a FASTA file.
    @param F FASTA File to read sequence from
    @result Single sequence
 */
String *msa_read_seq_fasta(FILE *F);

/** Save an MSA to a file descriptor.
    @param F File descriptor to save MSA to
    @param msa MSA to save to file F
    @param format File format to output as i.e. FASTA
    @param pretty_print If set prints periods ('.') in place of characters identical to corresponding characters in the first sequence
*/
void msa_print(FILE *F, MSA *msa, msa_format_type format, int pretty_print);

/** Save an MSA to a file path
    @param filename Full path of where to save file
    @param msa MSA to save to file F
    @param format File format to output as i.e. FASTA
    @param pretty_print If set prints periods ('.') in place of characters identical to corresponding characters in the first sequence
*/
void msa_print_to_file(const char *filename, MSA *msa, msa_format_type format, 
		       int pretty_print);

/** Write a single line of summary statistics for alignment.
    Base freqs | Total number of columns | Number of columns containing any gaps | number of columns containing all gaps
    @param msa MSA containing sufficient statistics
    @param F File to write sufficient statistics to
    @param header Prints header line instead of summary statistics
    @param start (Optional) If not -1, Statistics only printed after at this index
    @param end (Optional) If not -1, Statistics only printed before this index
    @note Start and end are half-open, 0-based
    @note All following MSAs must have the same alphabet
*/
void msa_print_stats(MSA *msa, FILE *F, char *label, int header, int start,
                     int end);

/** \} \name MSA cleanup functions 
   \{ */

/** Free all categories in MSA but keep everything else
    @param msa MSA to free categories from
*/
void msa_free_categories(MSA *msa);

/** Free all sequences in MSA but keep everything else
    @param msa MSA to free sequences from
*/
void msa_free_seqs(MSA *msa);

/** Free MSA and all of its components
    @param msa MSA to free
    @warning Names and Sequences are also freed even if they have been allocated externally
*/
void msa_free(MSA *msa);

/** \} */

/** Reduce alignment to unordered Sufficient Statistics (SS) with tuple size 3
   containing only 4d sites.  
   @pre Assumes that msa_label_categories has already been called using a GFF where the CDS regions on the + strand have been given feature type CDSplus and CDS regions on - have type CDSminus.
   @param msa MSA to reduce
   @param cm Category Map
*/
void reduce_to_4d(MSA *msa, CategoryMap *cm);

/** Update msa-> length.
    Used after sites are removed to keep length statistic up to date.
    @param msa MSA to update
*/
void msa_update_length(MSA *msa);


/** "Project" alignment on specified sequence, by eliminating all
   columns in which that sequence has a gap.  
   @param msa MSA to project
   @param refseq Reference Sequence to project onto
   @note Indexing of sequences starts with 1 
*/
void msa_project(MSA *msa, int refseq);

/** Creates a sub-alignment consisting of the specified sequences
   within the specified range of columns.  
   @param seqlist (Optional) list of sequences to use in sub-alignment, 
          if NULL all sequences used
   @param include Listed sequence will be included if "include" == TRUE 
          and excluded otherwise (include is ignored if seqlist == NULL).  
   @param start_col Starting column of sub alignment [start_col, end_col)
   @param end_col Ending column of sub alignment [start_col, end_col)
   @note All memory is copied
   @note First character at index 0
*/
MSA* msa_sub_alignment(MSA *msa, List *seqlist, int include, int start_col, 
                       int end_col);
/** Creates a "coordinate map" object with respect to the designated sequence.  
   @param msa MSA containing reference sequence to create coordinate map from
   @param refseq Index of sequence within MSA to treat as refseq
   @result Newly created Coordinate Map from MSA
   @note Indexing begins with 1. 
*/
msa_coord_map* msa_build_coord_map(MSA *msa, int refseq);

/** Saves a Coordinate Map to a file.
    @param F File descriptor to save Coordinate Map to
    @param map Coordinate Map to save to file
*/ 
void msa_coord_map_print(FILE *F, msa_coord_map *map);

/** \name Convert coordinates functions
  \{ */

/**  Converts a sequence coordinate to an MSA coordinate using a coordinate map object.  
   @param map Coordinate Map
   @param seq_pos Sequence position on coordinate map
   @result Returns -1 if sequence coordinate is out of bounds, otherwise msa coordinate
   @note Indexing begins with 1. 
*/
int msa_map_seq_to_msa(msa_coord_map *map, int seq_pos);
/** Converts an MSA coordinate to a sequence coordinate using a coordinate map object.
    @param map Coordinate Map
    @param pos MSA sequence position
    @result Returns -1 if MSA coordinate is out of bounds, otherwise sequence coordinate
*/
int msa_map_msa_to_seq(msa_coord_map *map, int pos);

/**  Converts coordinates of all features in a GFF_Set from one frame of
   reference to another. 
   @param msa MSA 
   @param set Feature Set to extract features from 
   @param from_seq Starting part of an index between 1 and nseqs, or 0 (for the frame of the entire alignment), or -1 (to infer feature by feature by sequence name)
   @param to_seq Ending part of an index between 1 and nseqs, or 0 (for the frame of the entire alignment)
   @param offset Added to start and end of each feature coords
   @note Features whose start and end coords are out of range will be dropped; if only the start or the end is out of range, they will be truncated.
*/
void msa_map_gff_coords(MSA *msa, GFF_Set *set, int from_seq, int to_seq, 
                        int offset);


/** Returns an array of msa objects, one for each feature.
    @param msa MSA object
    @param gff Features object.  Will be modified by this function!  (See note)
    @return An array of newly allocated MSA objects, one for each feature in the gff which overlaps with the msa.  NULL if no features overlap.  The length of the array is equal to the number of features in the gff when the function returns.
    @note The gff is heavily modified by this function; elements that do not overlap with the MSA are removed, and others have their coordinates changed to the MSA frame of alignment.
 */
MSA **msa_split_by_gff(MSA *msa, GFF_Set *gff);


/** Converts coordinate from one map to coordinate on another.
    For convenience when going from one sequence to another.  
    @param from_map (Optional) use NULL to indicate frame of entire alignment
    @param to_map (Optional) use NULL to indicate frame of entire alignment
    @param coord Coordinate to map from from_map to to_map
    @result coordinate on to_map OR returns -1 if out of range.
*/
int msa_map_seq_to_seq(msa_coord_map *from_map, msa_coord_map *to_map, 
                       int coord);

/** \} */

/** Free a Coordinate Map object.
    @param map MSA Coordinate Map to free
*/
void msa_map_free(msa_coord_map *map);

/** 
   Add labels to MSA categories.
   @pre Coordinates of GFF_Set must be in frame of ref of entire alignment
   @param msa MSA to set categories for (msa->categories may be null)
   @param gff Feature Set to map categories to
   @param cm Category Map from features to categories
   @todo Document "MSA" convention; Provide option to specify source sequence
   explicitly; What to do with overlapping categories? for now, just rely on external code to order them appropriately ...FIXME
*/
void msa_label_categories(MSA *msa, GFF_Set *gff, CategoryMap *cm);

/** Return sequence index of given sequence name or -1 if not found.
    @param msa MSA containing desired sequence
    @param name Sequence name to get index of
    @result Sequence Index or -1 if not found
*/
int msa_get_seq_idx(MSA *msa, const char *name);

/** Adds a sequence name to the msa, and allocates space for the sequence.
   @pre Sequence must not already be present in MSA. 
   @param msa MSA to add another sequence to
   @param name Name of sequence to add
   @result new sequence index 
   @note Fills in new sequence with missing data, except for columns where all other
  sequences contain a gap, then it fills in a gap instead. 
*/
int msa_add_seq(MSA *msa, char *name);
/** Allocate space in col_tuples for more sequences.  
   @param msa MSA to add sequence(s) to
   @param new_nseq New number of sequences that MSA holds (new_nseq should be > msa->nseq)
   @note Fills in new sequence with missing data, except for columns where all other
         sequences contain a gap, then it fills in a gap instead. */
void msa_add_seq_ss(MSA *msa, int new_nseq);

/** \name MSA Reverse / Reverse Complement
  \{ */

/** Returns complement of a single DNA base.  
    @param c Single DNA base to complement
    @return Complement of single DNA base
    @note Leaves all bases but A,C,G,T unchanged. 
*/
char msa_compl_char(char c);

/** Reverse complements a DNA sequence.
    @param seq Sequence data to reverse complement
    @param length Length of sequence data in chars 
*/
void msa_reverse_compl_seq(char *seq, int length);

/** Reverse complements a segment of a sequence 
    @param seq Sequence data containing segment to reverse complement
    @param start Starting index of sequence to reverse complement
    @param end Ending index of sequence to reverse complement
*/
void msa_reverse_compl_seq_segment(char *seq, int start, int end);

/** Reverse complements an entire MSA. 
    @param msa MSA to reverse complement
*/
void msa_reverse_compl(MSA *msa);

/** Reverse complements a segment of an MSA.
   @param msa MSA containing segments to reverse complement
   @param start Starting index of sequences in MSA to reverse complement
   @param end Ending index of sequences in MSA to reverse complement
 */
void msa_reverse_compl_segment(MSA *msa, int start, int end);

/** Reverse complements all segments of an MSA corresponding to "groups"
   in a GFF that appear to be completely on the reverse strand.  
   @param msa MSA containing segments to reverse complement
   @param gff Feature Set to create groups from
   @param aux_data Auxiliary array of site-specific integers (e.g., gap patterns) to be kept in sync with the alignment and/or GFF_Set
   @note The GFF is partitioned into "groups" using gff_partition_by_group,
   and groups are tested using gff_reverse_strand_only (see docs for these
   functions). */
void msa_reverse_compl_gff(MSA *msa, GFF_Set *gff, int *aux_data);

/** Reverse complement segments of an MSA corresponding to groups of
   features on the reverse strand. 
   This function can be used to ensure that
   sites in strand-specific categories (e.g., 1st codon position,
   intron) are oriented consistently.  It can be useful in phylogenetic
   analysis and in training the transition probabilities of a phylo-HMM.
   Strandedness is tested using
   gff_reverse_strand_only. 
   @pre Features have been grouped as desired
   @pre Groups are non-overlapping (see gff_group and gff_remove_overlaps)
   @pre GFF_Set uses coordinate frame of the alignment
   @param msa (Optional) Alignment object.  If NULL, only the GFF_Set (and optionally aux_data) will be altered 
   @param feats Set of features
   @param aux_data Auxiliary array of site-specific integers (e.g., gap patterns) to be kept in sync with the alignment and/or GFF_Set
   @note Adjusts the coordinates in the GFF_Set accordingly.
*/
void msa_reverse_compl_feats(MSA *msa, GFF_Set *feats, int *aux_data);

/** \} */

/** Split a single MSA into multiple MSAs by category.
   @param msa MSA containing sequence data
   @param submsas List of MSAs, one entry per category
   @param cats_to_do (Optional) Categories to use when splitting (default is all categories)
   @param tuple_size Size of tuples in sub MSAs, if > 1 then tuple_size-1 columns of "missing data" characters (Ns) will be inserted between columns that were not adjacent in the original alignment.  */
void msa_partition_by_category(MSA *msa, List *submsas, List *cats_to_do, 
                               int tuple_size);


/** Get counts for each base.
   @param start (Optional) If not -1, frequencies calculated starting from this interval
   @param end (Optional) If not -1, frequencies calculated stopping at this interval
   @note start and end are half-open, 0-based
   @result Newly allocated vector containing base counts (size of strlen(alphabet) listed in order of the alphabet)
 */
Vector *msa_get_base_counts(MSA *msa, int start, int end);


/** Get frequencies for each base.
   @param start (Optional) If not -1, frequencies calculated starting from this interval
   @param end (Optional) If not -1, frequencies calculated stopping at this interval
   @note start and end are half-open, 0-based
   @result Newly allocated vector containing base frequencies (size of strlen(alphabet) listed in order of the alphabet)
 */
Vector *msa_get_base_freqs(MSA *msa, int start, int end);

/** Get frequencies of k-tuples of bases, rather than of individual bases.  
   By convention, it is the *last* character in each
   tuple whose category is considered (this convention makes sense for
   models that consider the *predecessors* of each base). 
   @param msa MSA containing Sequence Data or Sufficient Statistics (with tuple_size equal k)
   @param cat (Optional) If not -1, Only bases of this category are considered
   @note A gap anywhere in a k-tuple causes it to be ignored.  
*/
void msa_get_base_freqs_tuples(MSA *msa, Vector *freqs, int k, int cat);

/** Get background frequencies for codon data based on 3x4 model of codon
    frequencies.
    @param backgd (output) The background frequencies.  On input should be initialized to a vector of size 64.
    @param msa (input) An alignment object, assumed to represent in-frame codons.
 */
void msa_get_backgd_3x4(Vector *backgd, MSA *msa);

/** Get length of a particular sequence (alignment length minus gaps).
    @param msa MSA containing sequence to measure
    @param seqidx Index of sequence to measure
    @result length of sequence (alignment length minus gaps)
*/
int msa_seqlen(MSA *msa, int seqidx);


/** \name MSA Informative sites functions
  \{  */

/** Get number of informative sites for a specified category.
    @param msa MSA
    @param cat If not -1, Only consider columns of this category
    @result Number of informative sites for category cat 
    @note Informative sites contain at least two non-gaps (or non-missing chars).
*/
unsigned int msa_ninformative_sites(MSA *msa, int cat);

/** Get coordinates of informative sites in a GFF_Set in refseq reference frame.
    @param msa MSA 
    @param min_informative Must have at least this many non-missing (and non-gap if gaps_are_informative==0) sites to be Informative   
    @param speclist (Optional) Integer list of indices limiting which species to consider in determining informativeness OR all if NULL
    @result Set of informative features
*/
GFF_Set *msa_get_informative_feats(MSA *msa, int min_informative,
				   List *speclist, int refseq, 
				   int gaps_are_informative);


/** Set up array indicating which sequences are to be considered
    "informative".
   This function can be useful for phylogenetic analysis. 
   @param msa Alignment
   @param not_informative List of sequence names *NOT* to be considered informative 
*/
void msa_set_informative(MSA *msa, List *not_informative);

/** \} */
/** Index the columns.
    @param msa MSA
    @param If the columns should be in order 
*/
void msa_index_cols(MSA *msa, int order);

/** \name MSA Clean data functions
  \{ */
/**
Clean an alignment in preparation for codon-based analysis.
   First remove all gaps in refseq (unless refseq is NULL).  
   Returns an error if resulting sequence does not have multiple-of-three length.  If strand == '-',
   get the reverse complement of entire sequence.  Then go through
   each sequence and if any stop codons are encountered, mask out
   the stop codon and the rest of that sequence to the end of the MSA.
   @param MSA containing codons to clean
   @param refseq Reference sequence 
   @param strand what strand type the reference sequence is
   @result Always 0
 */ 
int msa_codon_clean(MSA *msa, const char *refseq, char strand);
/** Clean an alignment of coding sequences (CDS exons from genomic DNA
   or mRNAs).  
   Remove sites with gaps and short blocks of ungapped
   sites, also look for frame shifts.  
   @param msa MSA containing coding sequences
   @param refseq Sequence known to be a CDS, start with start codon, end with stop codon.  
   @param min_ncodons Minimum number of ungapped codons allowed.  
   @param errstr Error string populated in case an error occurs
   @param keep_stop_codons If 1 keep stop codons, otherwise remove
   @result  0 on success and 1 if the alignment is rejected
 */
int msa_coding_clean(MSA *msa, int refseq, int min_ncodons, 
                     String *errstr, int keep_stop_codons);

/** Clean an alignment of indel artifacts.  
   Replace the chars
   adjacent to each indel and gap-less subsequences of insufficient
   length with missing data characters.  Also eliminate any sites at
   which bases are present from too few species.  If a model is to be
   used that considers tuples of sites of size greater than 1, then an
   appropriate number of columns of missing data will be left between
   sites that were not adjacent in the original data set.  This
   routine ignores issues of frame (cf. msa_coding_clean, link below). 
   @param msa MSA to clean
   @param indel_border Number of chars adjacent to each indel to discard
   @param min_nbases Minimum number of consecutive gap-less bases (per sequence)
   @param min_nseqs Minimum number of sequences
   @param tuple_size Size of tuples to be considered; tuple_size-1 columns of missing data will be maintained between non-adjacent columns
   @param mdata_char Missing data character to use
   @see msa_coding_clean
*/
void msa_indel_clean(MSA *msa, int indel_border, int min_nbases, 
                     int min_nseqs, int tuple_size, char mdata_char);
/** \} */


/** Concatenate one MSA onto another.
    @param aggregate_msa MSA to add source_msa onto
    @param source_msa MSA to concatenate onto another  
*/
void msa_concatenate(MSA *aggregate_msa, MSA *source_msa);


/** Randomly permute the columns of a multiple alignment.  
   @param msa MSA to randomly permute columns of
*/
void msa_permute(MSA *msa);
/** Reorder rows of MSA so that names match specified target order.
   @param msa MSA containing rows to re-order
   @param target_order List of names in the desired order
   @warning All names in the msa must be present in target_order.   Rows of missing data will be added for names that are in target_order but
   not in msa 
*/
void msa_reorder_rows(MSA *msa, List *target_order);

/** Get a single character for specified sequence and position.
    Provides a layer of indirection to handle cases where sufficient statistics are and
   are not used.
    @param msa MSA 
    @param seq Sequence index within MSA to get data from
    @param pos Position within sequence to get data from
    @result Char at position 'pos' of sequence 'seq' in MSA 'msa'
 */
char msa_get_char(MSA *msa, int seq, int pos);

/** \name MSA File Format functions 
   \{ */

/** Translate format type into char*.
    @param format An msa format
    @result A char* describing the format (either "SS", "MAF", "FASTA", "PHYLIP", "MPM", or "UNKNOWN")
 */
char *msa_format_to_str(msa_format_type format);

/** Translate string to msa_format_type.
    @param str String to translate into msa_format_type
    @result Of type msa_format_type that corresponds to str
 */
msa_format_type msa_str_to_format(const char *str);

/** Get file format based on filename suffix.
    @param fname Filename to determine file format for
    @result File format determined by filename suffix
    @seealso msa_format_for_content
 */
msa_format_type msa_format_for_suffix(const char *fname);

/** Get file format based on file contents 
    @param F File descriptor to file (or stdin) containing file data
    @param die_if_unknown If TRUE, then exit with an error message if the format cannot be detected
    @result File format determined by file contents
*/
msa_format_type msa_format_for_content(FILE *F, int die_if_unknown);

/** Get the appropriate filename suffix depending on msa_format_type.
    @param t Format type 
    @result filename suffix
 */
char *msa_suffix_for_format(msa_format_type t);

/** \} \name MSA alphabet functions
  \{ */

/** Remove 'N' from alphabet.
    Sometimes useful when fitting tree models.
    @param msa MSA containing alphabet to modify
 */
void msa_remove_N_from_alph(MSA *msa);


/** Test if alphabet has lower case letters
    @param msa MSA
    @return TRUE if alphabet has lower case letters, FALSE otherwise
 */
int msa_alph_has_lowercase(MSA *msa);

/** Replace all lower case characters in alignment and alphabet with upper case characters.
    @param msa MSA containing alignment and sequence data
 */
void msa_toupper(MSA *msa);

/** Reset alphabet of MSA
    @param msa Alignment to reset alphabet for
    @param newalph New alphabet to use for MSA
 */
void msa_reset_alphabet(MSA *msa, char *newalph);

 /** \} \name MSA Missing Data / Gaps functions
  \{ */


/** Convert all missing data characters to gaps, except for refseq.
   Ns in the reference sequence (if refseq > 0) will be replaced by
   randomly chosen bases 
   @param msa MSA
   @param refseq Index of reference sequence
   @warning only use if Ns in reference sequence are rare 
*/
void msa_missing_to_gaps(MSA *msa, int refseq);

/** Mask out all alignment gaps of length greater than k by changing
    gap characters to missing data characters.  
    This function is useful when modeling micro-indels.  
    @param msa MSA containing alignment gaps to mask
    @param k Largest length of gaps allowed without being masked
    @param refseq Index of reference sequence
    @note If refseq is > 0, the designated sequence will not be altered. 
    @warning If MSA is stored only in
    terms of sufficient statistics, an explicit alignment will be
    created */
void msa_mask_macro_indels(MSA *msa, int k, int refseq);


/** Identify sites which do not appear to have any real alignment information.
   @pre noaln must be preallocated to size msa->length.
   @param msa MSA
   @param refseqidx Index of the reference sequence within the list of sequences in the MSA
   @param min_block_size Sites belonging to blocks at least this size with only reference sequence and gaps or missing data in all other sequences have no real alignment information
   @param noaln Array filled with 1s (indicating no alignment information) or 0s.  
   @note This routine now can typically be
   replaced by the simpler one (msa_missing_col), because of better handling of
   missing data 
   @see msa_mising_col
*/
void msa_find_noaln(MSA *msa, int refseqidx, int min_block_size, int *noaln);

/** Test for missing data at specified column (all sequences except refseq).
   @param msa MSA containing sequences
   @param ref Index of reference sequence in MSA
   @param pos Column to test for missing data
   @result TRUE if all sequences (except refseq) in column are missing data, false otherwise
*/
int msa_missing_col(MSA *msa, int ref, int pos);

/** Strip ANY or ALL gaps or perform projection.
   @param msa MSA to strip gaps from or projection
   @param gap_strip_mode How to strip gaps; if STRIP_ALL_GAPS or STRIP_ANY_GAPS, removes all columns with ALL or ANY gaps, respectively.  Otherwise, assumes a *projection* is desired onto the sequence whose
   index is gap_strip_mode (indexing starts with 1).
   @note Gaps are expected to be represented by GAP_CHAR.  
   @warning If msa->categories is non-NULL, will be adjusted accordingly. 
*/
void msa_strip_gaps(MSA *msa, int gap_strip_mode);

/** Get number of gapped columns.
   @param msa MSA containing sequence to measure gapped columns
   @param gap_strip_mode How to handle gaps when counting columns
   @param start Site at which to start counting
   @param end Site at which to stop counting
   @result Number of gapped columns
   @note If mode == STRIP_ANY_GAPS, a gapped column is one containing at least one gap
   @note if mode ==
   STRIP_ALL_GAPS, a gapped column is one containing only gaps
*/
int msa_num_gapped_cols(MSA *msa, int gap_strip_mode, int start, int end);


/**\} */
/** Translate a list of sequence names into MSA indexes and/or Translate a list of 1-based indices into 0-based indices.  
    Useful in converting command-line arguments
    @param msa MSA containing sequence names
    @param seqnames List of sequence names and/or 1-based indices
    @result Indices of sequences in MSA and/or 0-based indices
    @note Warn if a name has no match.  
    @code
      msa { 'panTro', 'felis', 'canis'}	
      msa_seq_indices(MSA, 'panTro','canis')  -> '0','2'
      msa_seq_indices(MSA, '1','2','3')       -> '0','1','2'
    @endcode
*/
List *msa_seq_indices(MSA *msa, List *seqnames);



/**  Delete specified columns from an alignment.
     @param msa Alignment
     @param delete_cols Columns to delete are specified
   by an array of FALSEs and TRUEs of size msa->length (TRUE means
   delete) 
*/
void msa_delete_cols(MSA *msa, int *delete_cols);


/** Calculate pairwise differences between two sequences
    @param msa MSA containing both sequences 
    @param idx1 Index of the first sequence 
    @param idx2 Index of the second sequence
    @param ignore_missing If 1 skip sites with missing data during computation
    @param ignore_gaps If 1 skip sites with gaps during computation
    @result Pairwise differences result
*/
double msa_fraction_pairwise_diff(MSA *msa, int idx1, int idx2, 
				  int ignore_missing, int ignore_gaps);
				  
/** Translate alignment.
   @param msa Alignment to translate
   @param oneframe If TRUE,
   then every third base is a new codon, regardless of gaps.  Otherwise
   translate each species separately taking gaps into consideration.  
   @param frame Indicates where to start the translation, is an offset from the first alignment column.  
   @result Pointer to translated alignment sequences data
   @note If oneframe is FALSE, frame should have length msa->nseqs.  Otherwise only the first value is accessed and applies to all species.
*/
char **msa_translate(MSA *msa, int oneframe, int *frame);

#endif
