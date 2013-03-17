#!/bin/sh 

# preprocess help files for inclusion in executables

# Note: any line in .help_src file beginning with a pound sign is discarded

# let's just assume sed is in path
sed '/^#.*$/d ; s/\\$/\\\\/ ; s/$/\\n\\/ ; s/"/\\"/g ; s/%/%%/g ; 1s/^/char* HELP = "\\n/ ; $s/$/\n";/' $1


