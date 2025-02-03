#!/usr/bin/perl -w

# This script computes statistics 
# from a DOE experiment.

# Parse program arguments and set up output filenames.
use Cwd;
use File::Basename;

#system("source set_env_vars.sh");
if ( scalar(@ARGV) < 2 ) {
  print "Usage: ComputeDOEStatistics.pl <DOEFileList> <OutputFile>\n";
  exit;
}

# Initialize variables.
$outputdatapathPrefix = "./";


# Read command line.
$doeFileList = shift @ARGV;
$outputFile = shift @ARGV;
$csvfileSuffix = ".csv";

# Read inputs. 
print "\nReading $doeFileList\n";
open IN, "<$doeFileList"; # Open output file.
open OUT, ">$outputdatapathPrefix$outputFile$csvfileSuffix";

# Read the subject list file.
$firstlineString = "ExperimentID,";
$firstlineString .= "Ave. Volume_groundtruth(mm^3),Ave. Volume_detection(mm^3),".
    "Ave. Volume_overlap(mm^3),Ave. Volume_union(mm^3),Ave. Overlap_score,".
    "Ave. DICE_score,Ave. Sensitivity(%),Ave. Specificity(%),".
    "Ave. Classification_accuracy(%),Ave. Volume_difference(%),";
$firstlineString .= "Stdev. Volume_groundtruth(mm^3),Stdev. Volume_detection(mm^3),".
    "Stdev. Volume_overlap(mm^3),Stdev. Volume_union(mm^3),Stdev. Overlap_score,".
    "Stdev. DICE_score,Stdev. Sensitivity(%),Stdev. Specificity(%),".
    "Stdev. Classification_accuracy(%),Stdev. Volume_difference(%)";

$firstlineString .= "Parameters\n";

print OUT $firstlineString;

# $experimentCount = 1;
while(<IN>) {
    
    # For each experiment in the list,
    $line = $_;
    #
    chomp($line);
    print "\n$line\n";

    # Extract the directory name.
    $dirname = dirname($line);
    printf OUT ("%s,",$dirname);

    # Read the experiment list.
    open IN2, "<$line"; # Open output file.
    @scorecard = <IN2>;
    close IN2;

    # Separate header from measurements.
    $header =  shift @scorecard;
    chomp($header);
    (@headerTitles) = 
	     split(',', $header);
    chomp($line);

    # print ("$header\n");
    # print ("@scorecard\n");

    $experimentParameters=$headerTitles[$#headerTitles];
    # print $experimentParameters;

    # Compute averages and append parameters.
    @mymeanArray = &average(@scorecard);
    print ("@mymeanArray\n");

    # Compute stdevs and append parameters.
    @mystdevArray = &stdev(@scorecard);
    print ("@mystdevArray\n");

    # Write results to file.
    foreach(@mymeanArray) {
	printf OUT ("%.3f,", $_);
    }
    foreach(@mystdevArray) {
	printf OUT ("%.3f,", $_);
    }
    printf OUT ("%s",$experimentParameters);
    print OUT ("\n");
}

close OUT;
close IN;



# Compute average.
sub average{
    @input = @_;
    @meanArray = ();

    # First pass: compute averages.
    foreach(@input) {
	$line2 = $_;
	chomp($line2);
	(@parameters) = 
	     split(',', $line2);
	shift @parameters;
	$index = 0;
	foreach $value (@parameters) {
	    $meanArray[$index] += $value;
	    $index++;
	}
    }
    
    foreach(@meanArray) {
	$_ = $_ / ($#input + 1);
    }

    return @meanArray;
}


# Compute stdev.
sub stdev{
    @input = @_;
    @stdArray = ();

    # Compute average.
    @meanArray = &average(@input);
 
    # Squared error.
    foreach(@input){
	$line2 = $_;
	chomp($line2);
	(@parameters) = 
	     split(',', $line2);
	shift @parameters;
	$index = 0;
	foreach $value (@parameters) {
	    $tempValue = $value - $meanArray[$index];
	    $tempValue *= $tempValue;
	    $stdArray[$index] += $tempValue;
	    $index++;
	}
    }

    # Mean squared error.
    # Root mean squared error.
    foreach(@stdArray){
	$_ = $_ / ($#input);
	$_ = sqrt( $_ );
    }
    #print ("\n@stdArray\n");

    return @stdArray;
}
