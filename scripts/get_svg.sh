CP="$1"
cat ${CP}_colors.txt| sed "s/^/${CP} color /"
best=`grep "PS" "${CP}_split/${CP}_results.txt" | tail -1|cut -f6| sed 's/[a-z]//g'`
echo "${CP} model ${best} NA"
cat ${CP}_split/${CP}_PRL-G-PRL-${best}.tree | sed "s/^/${CP} migration /"

dot -Tsvg ${CP}_split/${CP}_PRL-T-PRL-${best}.dot -o ${CP}_PRL-T-PRL-${best}.svg
