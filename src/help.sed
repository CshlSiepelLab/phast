#!/usr/bin/sed -f

# preprocess help files for inclusion in executables

s/$/\\n\\/
s/"/\\"/g
s/%/%%/g
1s/^/char* HELP = "\\n\\n/
$s/$/\n";/
