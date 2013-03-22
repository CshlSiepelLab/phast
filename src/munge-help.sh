#!/bin/bash 

# preprocess help files for inclusion in executables

# Note: any line in .help_src file beginning with a pound sign is discarded

# assume wc, cut, and sed is in the path

function mungehelp {
  file=$1
  if [ "$2" ]; then
     numchar=$2
  else 
     numchar=10000
  fi
  sed '/^#.*$/d ; s/\\$/\\\\/ ; s/$/\\n\\/ ; s/"/\\"/g ; s/%/%%/g ; 1s/^/char HELP['$numchar'] = "\\n/ ; $s/$/\n";/' $file
}

numchar=`mungehelp $1 | wc -c`
mungehelp $1 $numchar
