#!/usr/bin/perl -w

# This script recursively executes our volumetric image concatenatenation program
# using a participant file list as input.


# Parse program arguments and set up output filenames.
use File::Basename;

#system("source set_env_vars.sh");

if ( scalar(@ARGV) < 2 ) {
  print "Usage: ThymusDixonConcatenation.pl <InputDataPath> <SubjectFile> \n";
  exit;
}

# Initialize variables.

$datapathPrefix = shift @ARGV;
$subjectlistFileName = shift @ARGV;
$outputdatapathPrefix = "./";

# $dicomconversionProgram = "dcm2nii";
# $parameterString = "-c N -d Y -i Y -n Y -p Y -f N -o " . $outputdatapathPrefix;
$dixonconcatenationProgram = "ConcatenateVolumes";
$formatExtension = ".nii";
$commentString = "#";

# Read inputs. 

print "Reading $subjectlistFileName\n";

open IN, "<$subjectlistFileName"; # Open output file.

# Read the subject list file.

while(<IN>) {
    
    # For each subject in the list,
    $line = $_;
    chomp($line);
    print "\n$line\n";
	
    # Parse the command line arguments.
    ($subjectID,$seriesID,@arguments) = split(',', $line);
	# $subjectID = pop(@arguments);
	# $seriesID = pop(@arguments);
	
    if( $subjectID ne $commentString)
    {
	# Put together the filenames and directories.
	$input = "";
	foreach $argument (@arguments) {
	$input = "$input ${datapathPrefix}/$subjectID/$argument";
	}
	$output = "$outputdatapathPrefix/$subjectID.$seriesID$formatExtension";
    print "\n$subjectID\n"; 
	print "\n$dixonconcatenationProgram --in $input --out $output\n";
	system("$dixonconcatenationProgram --in $input --out $output");

    }
}


close IN;

