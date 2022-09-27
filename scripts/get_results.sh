CP="$1"
cat ${CP}_colors.txt| sed "s/^/${CP} color /"
best=`grep "PS" "${CP}_split/results.txt" | tail -1|cut -f6| sed 's/[a-z]//g'`
best_full=`grep "PS" "${CP}_split/results.txt" | tail -1|cut -f6`

echo "${CP} model ${best_full} NA"
cat ${CP}_split/PRL-G-PRL-${best}.tree | sed "s/^/${CP} migration /"
grep label ${CP}_split/PRL-T-PRL-${best}.dot| sed  's/.*color=//'| sed 's/,label="/\t/'| sed 's/".*//'| sed 's/\\n.*//'| sed "s/^/${CP} label /"
cat ${CP}_split/PRL-T-PRL-${best}.tree | sed "s/^/${CP} tree /"

