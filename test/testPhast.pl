#!/usr/bin/perl -w

# script to test phast programs.  Rather than keeping around 
# copies of "good" results, run commands on version in two 
# different directories, and compare results.
# An old version can be checked out with:
# svn co -r 2000 $SVNROOT/phast/trunk phast-oldButGood/
# (make sure to change the phast directory in the old version's
# make-include.mk file, then compile).

# for each program it reads in a list of commands from the file 
# test_program.sh

# Commands beginning with @ are run only once, and the output is not 
# captured. (these may be used, for example, to produce temporary 
# files which are used as input for later tests).
#
# Other commands are run twice, once from each bin directory, and
# the output from each is stored in a temporary file (by default,
# "temp1.txt" and "temp2.txt".  These files are compared, and an 
# error is reported if there is any difference.
#
# For programs (such as phyloFit) which do not write their output to
# stdout, an additional type of line is provided.  If the command is
# preceded by a ! and a file name, then the command is run twice, but
# after each run, the file given is moved to the temporary file, and
# these files are compared.  stdout is NOT compared in this case.
# For example,
# !phyloFit.mod phyloFit align.ss  --tree "(spec1, (spec2, spec3))" -i SS
# will run the phyloFit command from the bin1 dir, then copy phyloFit.mod
# to temp1.txt, then run the command again from bin2, copy phyloFit.mod to
# temp2.txt, and then do a diff on temp1.txt and temp2.txt
#
# WARNING: this program counts "errors" whenever the two binaries produce
# different output.  So if there is a problem with the command and they both
# produce the same error message, it will still be a "passed" test.  But
# in this case there is generally a message written to stderr, which should
# be fairly clear by watching the screen as the tests are run.


if (scalar(@ARGV) != 2 && 
    scalar(@ARGV) != 3) {
    die "usage: perl testPhast.pl bin1 bin2 [progname], where bin1 and bin2 are directories containing phast executables to be compared";
}

my $bin1=$ARGV[0];
my $bin2=$ARGV[1];
my @progList=qw(phyloP);  # for now 
my $doProg="";
if (scalar(@ARGV) != 3) {
    $doProg = $ARGV[2];
}
my $tempfile1="temp1.txt";
my $tempfile2="temp2.txt";

for my $prog (@progList) {
    next if ($doProg && !($doProg eq $prog));
    open(INFILE, "test_$prog.sh") or die "error opening $prog.sh";
    my $numerror=0;
    my $numpassed=0;
    while (<INFILE>) {
	my $cmd=$_;
	chomp($cmd);
	next if (!$cmd || $cmd =~ /^#/);
	if ($cmd =~ /^@/) {
	    $cmd = substr($cmd, 1);
	    system($cmd);
	    next;
	} elsif ($cmd =~ /^!/) {
	    my @strings = split(' ', substr($cmd, 1));
	    my $outfile=$strings[0];
	    $cmd = join(' ', @strings[1..(scalar(@strings)-1)]);
	    print "testing $outfile from \"$cmd\"\n";
	    system("$bin1/$cmd > $tempfile1");
	    system("mv $outfile $tempfile1");
	    system("$bin2/$cmd > $tempfile2");
	    system("mv $outfile $tempfile2");
	} else {
	    print "testing \"$cmd\"\n";
	    system("$bin1/$cmd > $tempfile1");
	    system("$bin2/$cmd > $tempfile2");
	}
	my $diffout = `diff --brief $tempfile1 $tempfile2`;
	if ($diffout) {
	    $numerror++;
	    print STDERR "ERROR: \"$cmd\" results differed\n";
	    exit(-1);
	} else {
	    $numgood++;
	}
    }
    close(INFILE);
    print "passed $numgood $prog tests\n";
    if ($numerror > 0) {
	print "failed $numerror $prog tests\n";
    }
}
system("rm -f $tempfile1 $tempfile2");
