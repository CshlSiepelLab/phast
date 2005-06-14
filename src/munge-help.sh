#!/bin/sh 

# preprocess help files for inclusion in executables

# let's just assume sed is in path
sed 's/\\$/\\\\/ ; s/$/\\n\\/ ; s/"/\\"/g ; s/%/%%/g ; 1s/^/char* HELP = "\\n/ ; $s/$/\n";/' $1


