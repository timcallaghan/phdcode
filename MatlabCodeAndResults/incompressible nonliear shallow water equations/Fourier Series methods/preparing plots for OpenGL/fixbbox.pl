# Script to replace all occurances of BoundingBox xmin ymim xmax ymax
# with new xmin ymin xmax ymax values in an EPS file. The script 
# accepts an EPS file and the four (space) separated new bounding
# box values and adjust the EPS file accordingly.

die "Usage: $0 input-filename\n"
  unless @ARGV == 1;

# Read in the filename from the command line
my ($infilename) = @ARGV;
$infilename = \*STDIN if $infilename eq '-';
	
# Declare the search and replace strings
$search = "BoundingBox:   139   227   495   583";
$replace = "BoundingBox:   140   228   490   578";
# Declare a string to modify the output file name from the input filename
$mod = "new";
# Declare the output file name to be the same
$outfilename = $mod.$infilename;
	
open(IN,$infilename) ||
    die "cannot open $infilename for reading: $!";
## optional test for overwrite...
#die "will not overwrite $outfilename" if -e $outfilename;
open(OUT,">$outfilename") ||
    die "cannot create $outfilename: $!";
while (<IN>) {    # read a line from file IN into $_
    s/$search/$replace/g; # change the lines
    print OUT $_; # print that line to file OUT
}
close(IN);
close(OUT);