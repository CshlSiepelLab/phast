#!/usr/bin/perl

sub usage()
{
  print "\nUSAGE: perl $0 [binary] [man output directory]\n";
  print "     [binary] refers to an executable ELF binary that accepts -h as a flag to retrieve the help information\n";
  print "     [man output directory] refers to the directory where the generated man page files will be saved\n \n";
  exit 1
}

$ProgramToRun = $ARGV[0];
if ($ProgramToRun eq "")
{
  usage();
}
#print "Generating manual page for $ProgramToRun\n";

$dir = $ARGV[1];
if (! (-e $ProgramToRun))
{
  print "File does not exist\n";
  exit 1;
}

sub basename($) {
 my $file = shift;
 @parts = split(/\//, $file);
 return $parts[-1];
}

  $ProgramName = basename($ProgramToRun);

 @months = qw(January Febuary March April May June July August September October November Dececember);
 ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
 $year = 1900 + $yearOffset;
 $theTime = "$months[$month] $year";
 $outputname =  $dir . basename($ARGV[0]) . '.1';
 unlink($dir . basename($ARGV[0]) . '.1');
 open (OUTPUT, '>'. $dir . basename($ARGV[0]) . '.1');
sub trim($)
{
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}

$header = <<END;
.de  CW
.sp
.nf
.ft CW
..
.\" Process this file with
.\" groff -man -Tascii foo.1
.\"
.\" "verbatim" environment (from strace.1)
.de  CE
.ft
.fi
.sp
..
.\"
.TH $ProgramName 1 "$theTime" "Phast Package" "Phast Program"
END
print OUTPUT $header;

open FILE, "$ProgramToRun -h 2>&1 |";
while(<FILE>) {
  $LineFromFile = $_;
  if ($LineFromFile =~ m/(^[A-Z]+):/)
  {  
    $label = trim($1);
    print OUTPUT ".SH " . $label . "\n";
    $labelLength = length($label);
    for($i=0;$i<$labelLength;$i++)
    {  print OUTPUT " ";  }
    if ($LineFromFile =~ m/^[A-Z]+:(.*)/)
    { print OUTPUT basename($1) . "\n" };
  }
  else
  {
     #$origional = $LineFromFile;
     $LineFromFile =~ s/\\/\n/;
     $LineFromFile =~ s/\t/      /g;
     #if ($origional ne $LineFromFile)
       chomp($LineFromFile);
       chomp($LineFromFile); 

     print OUTPUT $LineFromFile . "\n";
  }
}

`gzip -f $outputname`
