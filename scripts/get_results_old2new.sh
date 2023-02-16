#!/usr/bin/env bash

CP="$1"
PTISSUE="$2"
SPATH="$3"
cat ${CP}_colors.txt| sed "s/^/${CP} color /"
best="R"
echo "${CP} model ${best} NA"

grep "\->"  ${CP}_split/T-${PTISSUE}-${best}-binarized.dot| cut -d' ' -f1,3|sed 's/^\t//' >  ${CP}_split/T-${PTISSUE}-${best}-binarized_parent_child.txt
grep label  ${CP}_split/T-${PTISSUE}-${best}-binarized.dot| tr ',' ' '| cut -d' ' -f1,4,5| sed 's/^\t//'| sed 's/color=//'| sed 's/label="//'| sed 's/"//g'| sed 's/]//'| sed 's/\\n/__/' >  ${CP}_split/T-${PTISSUE}-${best}-binarized_nodes.txt
grep label  ${CP}_split/G-${PTISSUE}-${best}-binarized.dot| tr ',' ' '| cut -d' ' -f1,6| sed 's/\t//g'| sed 's/color=//'| sed 's/label=//'| sed 's/"//g'| sed 's/]//' | sed 's/^\t//'>  ${CP}_split/G-${PTISSUE}-${best}-binarized_labels.txt
grep "\->"  ${CP}_split/G-${PTISSUE}-${best}-binarized.dot| cut -d' ' -f1,3|sed 's/^\t//' >  ${CP}_split/G-${PTISSUE}-${best}-binarized_parent_child.txt

# migration
python ${SPATH}/print_migration_tree_from_dot.py ${CP}_split/G-${PTISSUE}-${best}-binarized_labels.txt  ${CP}_split/G-${PTISSUE}-${best}-binarized_parent_child.txt | sed "s/^/${CP} migration /" | tr '\t' ' '

python ${SPATH}/print_label_tree_from_dot.py ${CP}_split/T-${PTISSUE}-${best}-binarized_nodes.txt ${CP}_split/T-${PTISSUE}-${best}-binarized_parent_child.txt | sed "s/^/${CP} /"  | tr '\t' ' ' 

