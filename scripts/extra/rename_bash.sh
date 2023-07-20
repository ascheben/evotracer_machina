ls -d *_split| while read l; do ls $l| while read m; do prefix=`echo $l| sed 's/_split//'`; echo mv ${l}/$m
${l}/${prefix}_${m};done;done > rename.sh
