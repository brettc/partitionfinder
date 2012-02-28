/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include <config.h>

#ifndef FREE_H
#define FREE_H

#include "utilities.h"


void Free_All_Nodes_Light(t_tree *tree);
void Free_All_Edges_Light(t_tree *tree);
void Free_Mat(matrix *mat);
void Free_Partial_Lk(phydbl *p_lk, int len, int n_catg);
void Free_Tree(t_tree *tree);
void Free_Edge(t_edge *b);
void Free_Node(t_node *n);
void Free_Cseq(calign *cdata);
void Free_Seq(align **d, int n_otu);
void Free_All(align **d, calign *cdata, t_tree *tree);
void Free_SubTree(t_edge *b_fcus, t_node *a, t_node *d, t_tree *tree);
void Free_Tree_Ins_Tar(t_tree *tree);
void Free_Tree_Lk(t_tree *tree);
void Free_NNI(t_tree *tree);
void Free_Edge_P_Lk_Struct(t_edge *b, t_tree *tree);
void Free_Node_Lk(t_node *n);
void Free_Edge_Lk(t_tree *tree, t_edge *b);
void Free_Model(model *mod);
void Free_Model_Complete(model *mod);
void Free_Model_Basic(model *mod);
void Free_Custom_Model(model *mod);
void Free(void *p);
void Free_Input(option *input);
void Free_Reachable(t_tree *tree);
void Free_St(supert_tree *st);
void Free_Eigen(eigen *eigen_struct);
void Free_Triplet(triplet *t);
void Free_Tree_Pars(t_tree *tree);
void Free_Edge_Pars(t_edge *b, t_tree *tree);
void Free_One_Spr(spr *this_spr);
void Free_Spr_List(t_tree *tree);
void Free_Actual_CSeq(calign *data);
void Free_Prefix_Tree(pnode *n, int size);
void Free_Pnode(pnode *n);
void Free_Edge_Labels(t_edge *b);
void Free_Nexus_Com(nexcom **com);
void Free_Nexus_Parm(nexparm *parm);
void Free_Nexus(option *io);
void Free_Optimiz(optimiz *s_opt);
void Free_Tree_List(t_treelist *list);
void Free_Bip(t_tree *tree);

#endif
