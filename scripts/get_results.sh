#!/usr/bin/env bash
CP="$1"
PTISSUE="$2"
cat ${CP}_colors.txt| sed "s/^/${CP} color /"
best=`grep "PS" "${CP}_split/results.txt" | tail -1|cut -f6| sed 's/[a-z]//g'`
best_full=`grep "PS" "${CP}_split/results.txt" | tail -1|cut -f6`
# output data is the same even if file is labeled "R", "S", etc. as this only reflects the -m option, not the result
# fix best to "R" to ensure compatibility with '-m 3' option or "S" for '-m 1'
best="R"

echo "${CP} model ${best_full} NA"
cat ${CP}_split/${PTISSUE}-G-${PTISSUE}-${best}.tree | sed "s/^/${CP} migration /" | tr '\t' ' '
grep label ${CP}_split/${PTISSUE}-T-${PTISSUE}-${best}.dot| sed  's/.*color=//'| sed 's/,label="/\t/'| sed 's/".*//'| sed 's/\\n.*//'| sed "s/^/${CP} label /"  | tr '\t' ' '
cat ${CP}_split/${PTISSUE}-T-${PTISSUE}-${best}.tree | sed "s/^/${CP} tree /"  | tr '\t' ' '

