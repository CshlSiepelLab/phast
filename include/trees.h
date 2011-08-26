/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file trees.h
    Functions and structures to hold, manipulate, and test trees, branches, and nodes.  
    @ingroup phylo
*/

#ifndef TREES_H
#define TREES_H

/** Maximum length of tree described as a string */
#define MAX_TREESTR_LEN 10000
/** Maximum line length that can be read from file at a time */
#define MAX_LINE_LEN 10000

#include <stdio.h>
#include <lists.h>
#include <stringsplus.h>

typedef struct tree_node TreeNode;

/* TODO: starting to accumulate a lot of data about the entire tree at
   each node; perhaps time to adopt a less purely recursive structure,
   with a higher level struct describing the tree as a whole. */
/** Structure describing a single node of a tree. 

  @warning The lists
   should be accessed only via the functions tr_preorder,
   tr_inorder, tr_postorder.  
  @warning Most of these auxiliary
   attributes (nnodes, height, nodes, preorder, inorder) are not
   guaranteed to remain correct if the structure of a tree is
   altered after initialization
*/
struct tree_node {
  TreeNode *parent;		/**< Node that is a parent to this node */
  TreeNode *lchild, 		/**< Node that is a child of this node, drawn to the left in a diagram */
   *rchild;   			/**< Node that is a child of this node, drawn to the right in a diagram */
  double dparent;		/**< Distance to parent node */
  char name [STR_MED_LEN];	/**< Name of this node i.e. 'Drosophila 23' usually one of the sequence names */
  void *data;                   /**< Allows generic data to be 
                                   associated with tree node */ 
  int id;			/**< Uniquely identifying id number for this node */
  int nnodes;                   /**< Number of nodes in subtree defined
                                   by this node */
  int height;                   /**< Height of this node in the tree
                                   (maximum distance to a leaf, in
                                   terms of number of edges) */
  char *label;                  /**< paml-style node label indicated by # 
                                   sign and a character string.  Used 
				   to indicate which lineage-specific 
				   model to use */
  int hold_constant;            /**< If 1, do not optimize this branch length. 
                                   Indicated by an exclamation point after the
                                   branch length */
  List *nodes;                  /**< List of nodes: ith element is a
                                   pointer to the node with id = i (Only guaranteed to be defined for
     				   the TreeNode at the root of a tree, defined on initialization) */
  List *preorder;               /**< List of nodes in the order of a
                                   preorder traversal from the root (Only guaranteed to be defined for
     				   the TreeNode at the root of a tree, defined on demand)*/
  List *inorder;                /**< List of nodes in the order of an
                                   inorder traversal from the root (Only guaranteed to be defined for
     				   the TreeNode at the root of a tree)*/
  List *postorder;              /**< List of nodes in the order of a
                                   postorder traversal from the
                                   root (Only guaranteed to be defined for
     				   the TreeNode at the root of a tree)*/
};

/** \name Tree allocation functions 
\{ */

/** Parse a tree from a file in Newick (New Hampshire) format 
   @param f File descriptor containing data to make new TreeNode from
   @result Newly allocated tree node with data from file
*/
TreeNode *tr_new_from_file(FILE *f);

/** Parse a Newick-formatted tree from a character string 
   @param s String containing data to make new TreeNode from
   @result Newly allocated tree node with data from string
*/
TreeNode *tr_new_from_string(const char *s);

/** Create and initialize a new tree node
   @result Newly allocated tree node
*/
TreeNode *tr_new_node();

/** \} */

/** Add specified child to specified parent, creating all requisite
   links.  If the parent already has two children, add a new node
   (to simulate an nary tree with a binary tree)
   @param parent Parent of the child node should be added below
   @param child Node to be added as a child below parent
 */
void tr_add_child(TreeNode *parent, TreeNode *child);

/** Returns a newly allocated char * with newick representation of tree.
  @param root Tree Node that is to be represented in a newick format string along with all its descendants
  @param show_branch_lengths Whether to include branch lengths in the newick format
  @result Newick tree string representing all nodes below supplied root node
 */
char *tr_to_string(TreeNode *root, int show_branch_lengths);
void tr_to_string_recur(char *str, TreeNode *node, int show_branch_lengths);

/** \name Tree print functions
\{ */

/** Print tree in New Hampshire / Newick format.
   @param f File descriptor of where to print/save tree in newick format
   @param root Tree node that will be printed to a file, along with all its descendants
   @param show_branch_lengths Whether to show branch lengths in the file
 */
void tr_print(FILE* f, TreeNode *root, int show_branch_lengths);
void tr_print_recur(FILE* f, TreeNode *n, int show_branch_lengths);

/** Print tree in New Hampshire / Newick format (Ordered).
   @param f File descriptor of where to print/save tree in newick format
   @param root Tree node that will be printed to a file, along with all its descendants
   @param show_branch_lengths Whether to show branch lengths in the file
 */
void tr_print_ordered(FILE* f, TreeNode *root, int show_branch_lengths);
void tr_print_ordered_recur(FILE* f, TreeNode *n, int *left_right,
                            int show_branch_lengths);


/**  Print a (very basic!) postscript rendering of a tree.
    @param F destination file
    @param tree Tree root
    @param Whether to print branch lengths by edges
    @param if TRUE, branches will be right-angled, otherwise will be diagonal
    @param use_branch_lens If TRUE, time will be laid out such that edges are proportional to branch lengths (dparent attributes)
    @param horizontal_layout If TRUE, tree will be laid out with root on left and leaves on right; otherwise, root will be at top and leaves on bottom
 */
void tr_print_ps(FILE *F, TreeNode *tree, int show_branch_lens,
                 int square_branches, int use_branch_lens, 
                 int horizontal_layout);

/** Print verbose description of each node
   @param F File descriptor to write to
   @param tree Tree to traverse in preorder and write out verbose description of each node
 */
void tr_print_nodes(FILE *F, TreeNode *tree);


/** \} */

/** Free a tree node and all of its descendants
   @param tree TreeNode to start the freeing at
*/
void tr_free(TreeNode *n);

/** Traverse tree to set nnodes at each node (postorder); also set
   height at each node, and create "nodes" list.
   @param tree Tree Node to start traversal at
 */
void tr_set_nnodes(TreeNode *tree);

/** Set the ID counter to zero */
void tr_reset_id();

/** \name Tree copy functions 
\{  */

/** Copy a tree onto an already existing destination tree node.
   @param src Source root tree node
   @param dest Destination root tree node
*/
void tr_cpy(TreeNode *dest, TreeNode *src);

/** Copy a tree creating newly allocated TreeNode.
   @param src Source root tree node to traverse to copy
   @result Newly allocated Tree Node with copied data
*/
TreeNode *tr_create_copy(TreeNode *src);

/** Copy contents of a single tree node (ignore pointers) 
    @param dest Destination Tree node to copy contents to
    @param src Source Tree node to copy contents from
*/
void tr_node_cpy(TreeNode *dest, TreeNode *src);
/** \} \name Tree order of traversal list functions
\{ */

/** Obtain a list representing a preorder traversal of the tree.
    Allows for simplified loops (no stacks!) and improved
    efficiency.
    @param tr First node in the pre-order list
    @result list representing a preorder traversal of the tree
 */
List *tr_preorder(TreeNode *tr);

/** Obtain a list representing an in-order traversal of the tree.
    Allows for simplified loops (no stacks!) and improved
    efficiency.
    @param tr First node in the pre-order list
    @result List representing an in-order traversal of the tree
 */
List *tr_inorder(TreeNode *tr);

/** Obtain a list representing a postorder traversal of the tree.
    Allows for simplified loops (no stacks!) and improved
    efficiency.
    @param tr First node in the post-order list   
    @result List representing post-order traversal of the tree
 */
List *tr_postorder(TreeNode *tr);

/** \} */

/** Provide x-y coordinates for layout. 
    @param x0 Upper left x bound
    @param y1 Upper left y bound
    @param x1 Lower right x bound
    @param y1 Lower right y band
    @param x On return, will contain x-coordinates for nodes in order of tree->nodes. Must be preallocated
    @param y On return, will contain y-coordinates for nodes in order of tree->nodes. Must be preallocated
    @param use_branch_lens IfTRUE, tree will be laid out such that edges are proportional to branch lengths (dparent attribute)
    @param horizontal If TRUE, tree will be laid out with root on left and leaves on right; otherwise, root will be at top and leaves at bottom
*/
void tr_layout_xy(TreeNode *tree, int x0, int y0, int x1, int y1, 
                  int *x, int *y, int use_branch_lens, int horizontal);

/** \name Tree branch length functions 
\{ */
/** Compute and return sum of lengths at all edges in a single node 
    @param t Tree Node to test lengths of edges with
    @result Total length
*/
double tr_total_len(TreeNode *t);

/** Compute and return sum of lengths of edges in subtree below
    given node
     @param sub_root Starting node of subtree for counting edges
     @result Number of edges found in subtree
*/
double tr_total_len_subtree(TreeNode *sub_root);

/** Compute and return maximum branch length in subtree.
   @param sub_root Starting node of subtree for determining largest branch length
   @result Largest branch length in tree/subtree specified
 */
double tr_max_branchlen(TreeNode *sub_root);

/** Compute and return distance from given node to root of tree.
   @param node Node to start measuring distance from.
   @result Distance from node 'node' to root
*/
double tr_distance_to_root(TreeNode *node);

/** \} */

/** Reurn a node from a tree given its name.
  @param t Tree containing node to return
  @param name Name of node to return
  @result NULL if not found, otherwise tree node specified by name
*/
TreeNode *tr_get_node(TreeNode *t, const char *name);

/** \name Tree scale functions 
\{*/
/** Scale all branch lengths by constant factor.
   @param t Root node of tree to modify
   @param scale_const Amount to scale all branch lengths by
*/
void tr_scale(TreeNode *t, double scale_const);

/** Scale all branch lengths by constant factor in subtree beneath given node 
   @param t Root node of tree containing subtree
   @param sub First node in subtree to modify
   @param scale_const Amount to scale all branch lengths by
   @param include_leading Whether to also scale the leading node
*/
void tr_scale_subtree(TreeNode *t, TreeNode *sub, double scale_const,
		      int include_leading);


/** Scale a tree so that the total branch length of some subtree is as
   defined by a second tree.  
   @param tree  Tree to scale
   @param sub Subtree of 'tree' with leaf
   names a subset of those in the
   first.
   @result Scale factor
   @note This function is similar to tr_extrapolate, but the branch
   length proportions of the larger tree are used without change and
   the smaller tree need not be a proper subtree of the larger
   tree. 
*/
double tr_scale_by_subtree(TreeNode *tree, TreeNode *sub);


/** \} \name Tree prune functions 
\{ */

/**  Prune away all leaves whose names are in (or not in) the specified
    list.  Nodes will be removed and branches combined (branch lengths
    added) to restore as a proper binary tree. 
    @param[in,out] t Tree to prune (may be altered because root can change)
    @param[in,out] names List of names.  On return, will contain list of names of leaves
           that were pruned away.
    @param[in] all_but if FALSE, prune leaves *in* 'names'; if TRUE, prune leaves *not in* 'names'
    @param[out] id_map if not NULL, should be allocated to the number of nodes 
    in original tree. On return, will be filled in with the new id for each 
    node
*/
void tr_prune(TreeNode **t, List *names, int all_but, int *id_map);

/** Prune away all nodes not in the specified subtree.
    @param t Root node of tree to prune
    @param n Root node of subtree to keep
*/
void tr_prune_supertree(TreeNode **t, TreeNode *n);

/** Prune away all nodes in the specified subtree.
    @param t Root node of tree to prune
    @param n Root node of subtree to prune away
 */
void tr_prune_subtree(TreeNode **t, TreeNode *n);

/** \} */

/**  Return the Least Common Ancestor (LCA) of the given species. 
     @pre Ids must be numbered in preorder (a node's parent always has a smaller id than it does and
    left descendants have smaller ids than right descendants)
     @param tree Tree containing species specified in list 'names'
     @param names Names of sequences to compute LCA for
     @result Least common ancestor of nodes specified in 'names'
  */
TreeNode *tr_lca(TreeNode *tree, List *names);

/** Given two trees, one of which is a (proper) subtree of the
    other, create a hybrid tree composed of the smaller tree and a
    scaled version of the larger tree.  
    @param sub Subtree of tree 'super'
    @param super Tree that contains the subtree 'sub'
    @result Hybrid tree composed of the smaller tree and a scaled version of the larger tree
    @note This
    function can be used to extrapolate from a small phylogeny for
    which accurate branch length estimation is possible (e.g., of
    eutherian mammals) to a larger phylogeny for which approximate
    branch length proportions are available, but absolute branch
    length estimates are not (e.g., of more distant vertebrates).
 */
TreeNode *tr_hybrid(TreeNode *sub, TreeNode *super);

/** Partition leaves of tree at (branch above) given node.  All
    descendant leaves of 'sub' will be added to 'inside' list and all
    non-descendants will be added to 'outside' list.  Lists must be
    pre-allocated. 
    @pre Lists 'inside' and 'outside' must be pre-allocated
    @param tree Tree containing subtree
    @param sub Node defining subtree
    @param inside List containing descendant leaves 
    @param outside List containing non-descendant leaves
*/
void tr_partition_leaves(TreeNode *tree, TreeNode *sub, List *inside, 
                         List *outside);

/** Partition leaves of tree at all nodes.  All
    descendant leaves of 'sub' will be added to 'inside' list and all
    non-descendants will be added to 'outside' list.  Lists must be
    pre-allocated. 
    @pre Lists 'inside' and 'outside' must be pre-allocated
    @param tree Tree containing subtree
    @param sub Node defining subtree
    @param inside (Optional) List containing descendant leaves 
    @param outside (Optional) List containing non-descendant leaves
 */
void tr_partition_nodes(TreeNode *tree, TreeNode *sub, List *inside, 
			List *outside);

/** Return a list of the leaf names in a given tree 
    @param tree Tree containing leaf names
    @result List of leaf names
*/
List *tr_leaf_names(TreeNode *tree);

/** Name all un-named ancestral nodes.
   If a node is unnamed, give
   it a name that is a concatenation of the name of a leaf from its
   left subtree and the name of a leaf from its right subtree.
   Leftmost descendants are selected, for lack of any better
   criterion.  
   @param tree Tree containing nodes that might not all have names
*/
void tr_name_ancestors(TreeNode *tree);

/** Re-root tree.

    Subtree originally beneath selected node will become
    right subtree of root, and remainder of tree will be left
    subtree.  
    @warning ids will not be altered, so they will no longer be 
    consistent with a preorder traversal of the tree 
    @warning will not maintain memory of which branches are to be held
    constant
  @param tree Tree to re-root
  @param newroot The new root for specified tree
  @param include_branch If include_branch == FALSE, the selected node will become
    the new root, and a zero-length branch to its right will connect
    it to its original subtree.  If instead include_branch == TRUE,
    then the branch above the selected node will also be included in the
    right subtree.  In this case, the selected node will become the
    right child of the new root and the branch in question will become
    the right branch beneath the new root.  The left branch beneath
    the new root will have length zero and will connect to the former
    parent of the selected node. 
 */
void tr_reroot(TreeNode *tree, TreeNode *newroot, int include_branch);

/** Return an array indicating whether each node is in the designated subtree
    @param t Tree containing subtree 'sub'
    @param sub Subtree of tree 't'
    @result Array of int indicating whether each node is (1) in he designated subtree or not (0). Indexed by id
 */
int* tr_in_subtree(TreeNode *t, TreeNode *sub);

/** \name Tree label functions 
\{ */

/** Set the label for a tree node
    @param t Tree Node to label
    @param label Label to apply to specified tree node
 */
void tr_label(TreeNode *t, const char *label);

/** Set label for a tree node (found by node name).
   @param tree Tree containing node to label
   @param nodename Name of the node to set label for
   @param label New label to apply to tree node
 */
void tr_label_node(TreeNode *t, const char *nodename, const char *label);

/** Set label for an entire subtree.
   @param tree Tree containing subtree to label
   @param subtreeNode Root node of the subtree to label
   @param include_leading_branch Whether to include node subtreeNode in the subtree or treat it as a parent
   @param label Label to apply to set  
*/
void tr_label_subtree(TreeNode *tree, const char *subtreeNode, 
		      int include_leading_branch,
		      const char *label);

/** Create a list of nodes with a given label.
    @param[in] tree Tree containing labeled nodes
    @param[in] label Label that identifies nodes to be added to list
    @param[out] rv List of nodes with given label
 */
void tr_get_labelled_nodes(TreeNode *tree, const char *label, List *rv);
/** \} */

#endif 
