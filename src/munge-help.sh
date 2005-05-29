#!/bin/sh 

# preprocess help files for inclusion in executables

# let's just assume sed is in path
sed 's/$/\\n\\/ ; s/"/\\"/g ; s/%/%%/g ; 1s/^/char* HELP = "\\n\\n/ ; $s/$/\n";/' $1


