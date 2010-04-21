#!/usr/bin/perl -w
use Getopt::Long;

# script to test phast programs.  Rather than keeping around 
# copies of "good" results, run commands on version in two 
# different directories, and compare results.
# An old version can be checked out with:
# svn co -r 2000 $SVNROOT/phast/trunk phast-oldButGood/
# (make sure to change the phast directory in the old version's
# make-include.mk file, then compile).

# for each program it reads in a list of commands from the file 
# test_program.sh

# Commands beginning with @ are run twice (after stripping the !), 
# once from each bin directory, and the output from each is stored
# in a temporary file (temp_$filename_1.txt and temp_$filename_2.txt).  
# These files are compared, and an error is reported if there is any 
# difference.
#
# For programs (such as phyloFit) which do not write their output to
# stdout, file names can be provided to compare in place of stdout.
# A list of file names to compare is given before the command, each 
# preceded by a !.  If stdout should be included, add !stdout to the
# list.
# For example,
# !stdout !features.bed @phastCons --most-conserved features.bed align.ss model.mod 
# will run phastCons from the bin1 dir, copy stdout to temp_stdout_1 and
# features.bed to temp_features.bed_1, then run phastCons from the bin2
# dir, copy stdout to temp_stdout_2 and features.bed to temp_features.bed_2,
# then compare temp_stdout_1 to temp_stdout_2, and temp_features.bed_1 to
# temp_features.bed_2.  If there are any differences in either comparison
# it will report an error.
#
# Commands which do not begin with @ are treated as simple shell commands.
# Output is not captured or compared to anything.  (These can be used
# to produce temporary files which are used as input for later tests).
#
# WARNING: this program counts "errors" whenever the two binaries produce
# different output.  So if there is a problem with the command and they both
# produce the same error message, it will still be a "passed" test.  But
# in this case there is generally a message written to stderr, which should
# be fairly clear by watching the screen as the tests are run.

# TODO: add numerical comparison of files which may have slightly different
# numbers?

my $tempPrefix="temp";
my $help=0;
my $startLine=-1;
my $endLine=-1;
my $noQuitOnError=0;
my $doLineStr="";
my $scriptName="test_phast.sh";

Getopt::Long::Configure('no_ignore_case');
GetOptions('h|help' => \$help,
           't|temp=s' => \$tempPrefix,
	   'l|line=s' => \$doLineStr,
	   's|start-line=s' => \$startLine,
	   'i|script=s' => \$scriptName,
           'n|no-quit' => \$noQuitOnError) 
    or die "Invalid option.  Try ./testPhast.pl --help";

if ($help) {
    print "usage: perl testPhast.pl binDir1 binDir2 [progName1 progName2 ...]\
\
binDir1 and binDir2 are two different directories containing\
executables to be compared.  By default the commands in test_phast.sh\
are run.  Each \"testing\" command in this file should be preceded by\
a \"@\".  Then stdout/stderr will be captured and compared after\
calling the command from each binDir.  Additional results to be\
compared from the commands can come before the command as filenames\
preceded by a \"!\".\  Commands which do not begin with a \"@\" will\
be run just once, using the default PATH.  (This is useful for\
producing input files for the testing commands.)  If you want to avoid\
comparing stdout/stderr, you can precede the command with -stdout or\
-stderr.  Programs are labelled in the script file by a line like\
\"*** program ***\".  If program names are provided only commands\
within the blocks given are run.  Otherwise all commands are run.
\
example testing command:\
!tree.cons.mod !tree.noncons.mod \@phastCons --estimate-trees tree align.ss model.mod\
will run phastCons from each binDir and report an error if either\
command produces a different stdout, stderr, tree.cons.mod, or\
tree.noncons.mod file\
\
options:\
--temp,-t <tempPrefix>\
  prefix for temporary files used.  Default is \"temp\"\. Temporary
  files are cleaned if all tests are passed.\
--line,-l <line(s)>\
  only perform tests on specified line of script.  Lines are\
  1-based.  Can give several lines, ie, 1,4-6,8.  Commands which\
  do not begin with @ are always performed.
--start-line,-s <line>\
  only perform tests on lines of the script starting from this line.\
  Commands which do not begin with @ are always performed.\
--script-name,-i <scriptName>
  name of script with testing commands. Default is test_phast.sh.
--no-quit,-n\
  do not quit if error detected\
--help,-h\
  print this help\n";
    exit(0);
}

if (scalar(@ARGV) < 2) {
    die "usage: perl testPhast.pl bin1 bin2 [program1 program2 ...].  Try perl testPhast.pl --help";
}

my @programs=();
if (@ARGV >= 3) {
    @programs = @ARGV[2..(scalar(@ARGV)-1)];
}

my $bin1=$ARGV[0];
my $bin2=$ARGV[1];
my $numerror = 0;
my $numgood=0;


my @doLines=();
if ($doLineStr) {
    my @lineStrs=split(',', $doLineStr);
    foreach $lineStr (@lineStrs) {
	my @temp=split('-', $lineStr);
	if (scalar(@temp)==1) {
	    push(@doLines, $lineStr);
	} elsif (scalar(@temp)==2) {
	    die "format error $doLineStr $temp[0] > $temp[1]" if ($temp[0] > $temp[1]);
	    for (my $i=$temp[0]; $i<=$temp[1]; $i++) {
		push(@doLines, $i);
	    }
	} else {
	    die "format error $doLineStr";
	}
    }
#    print "doing lines " . join(',',@doLines) . "\n";
}
my $errorFlag=0;

sub compare_files {
    my $file1=$_[0];
    my $file2=$_[1];
    my $diffout = `diff --brief $file1 $file2`;
    if ($diffout) {
	print STDERR "ERROR: $file1 and $file2 differ\n";
	exit(-1) if (!$noQuitOnError);
	$errorFlag=1;
    } 
    system("rm -f $file1 $file2");
}


open(INFILE, $scriptName) or die "error opening $scriptName";

my $line=0;
my $currProgram="";
my $doThisProgram = !@programs;

while (<INFILE>) {
    $line++;
    my $cmd=$_;
    chomp($cmd);
    next if (!$cmd || $cmd =~ /^#/);
    if ($cmd =~ m/^[\s]*[\*]+[\s]*([A-Za-z][^\s]*)[\s]*[\*]*[\s]*$/) {
	$currProgram=$1;
	if (@programs) {
	    $doThisProgram=0;
	    foreach my $tempProgramName (@programs) {
		if ($tempProgramName eq $currProgram) {
		    $doThisProgram=1;
		}
	    }
	}
	next;
    }
    next if (!$doThisProgram);
    my $doTest=1;
    if (@doLines) {
	my $i;
	for ($i=0; $i < scalar(@doLines); $i++) {
	    last if ($doLines[$i] ==  $line);
	}
	$doTest=0 if ($i == scalar(@doLines));
    }
    $doTest=0 if ($line < $startLine);
    my @compareFiles=();
    my $compareStderr=1;
    my $compareStdout=1;
    while ($cmd =~ /^!/ || $cmd =~ /^-/) {
	my @fields=split(' ', $cmd);
	my $tempfile = substr($fields[0], 1);
	$compareStderr=0 if ($fields[0] eq "-stderr");
	$compareStdout=0 if ($fields[0] eq "-stdout");
	push(@compareFiles, $tempfile) if ($cmd =~ /^!/);
	$cmd = join(' ', @fields[1..(scalar(@fields)-1)]);
    }
    die if (!$cmd);
    if ($cmd =~ /^@/) {
	next if (!$doTest);
	$errorFlag=0;
	$cmd = substr($cmd, 1);
	print "$scriptName:" if (scalar(@ARGV) > 3);
	print "$line: $cmd\n";
	system("$bin1/$cmd 1>$tempPrefix.1.stdout 2>$tempPrefix.1.stderr");
	print `grep -iE "error|abort|fail|assertion" $tempPrefix.1.stderr`;
	foreach $file (@compareFiles) {
	    system("mv $file $tempPrefix.1.$file") if (! ($file eq "stdout"));
	}
	system("$bin2/$cmd >$tempPrefix.2.stdout 2>$tempPrefix.2.stderr");
	print `grep -iE "error|abort|fail|assertion" $tempPrefix.1.stderr`;
	foreach $file (@compareFiles) {
	    system("mv $file $tempPrefix.2.$file") if (! ($file eq "stdout"));
	}
	compare_files("$tempPrefix.1.stdout", "$tempPrefix.2.stdout") 
	    if ($compareStdout);
	compare_files("$tempPrefix.1.stderr", "$tempPrefix.2.stderr") 
	    if ($compareStderr);
	foreach $file (@compareFiles) {
	    compare_files("$tempPrefix.1.$file", "$tempPrefix.2.$file");
	}
	if ($errorFlag) {
	    $numerror++;
	}
	else {
	    $numgood++;
	}
    } else {
	die "can't compare files unless command preceded by @ at $cmd"
	    if (@compareFiles);
	system($cmd);
    }
}
close(INFILE);

print "passed $numgood tests\n";
if ($numerror > 0) {
    print "failed $numerror tests\n";
}
