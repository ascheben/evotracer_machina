CP="$1"
sed -i '1s/^/0 ASVXXX\n/' ${CP}_tree_split.txt
sed -i -e '$aPRL 3' ${CP}_colors.txt
sed -i '1s/^/ASVXXX PRL\n/' ${CP}_labels_split.txt
