grep label $1| sed  's/.*color=//'| sed 's/,label="/\t/'| sed 's/".*//'| sed 's/\\n.*//'
