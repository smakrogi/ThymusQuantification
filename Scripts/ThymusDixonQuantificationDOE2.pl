#!/usr/bin/perl -w

# This script implements a DOE experiment  
# over a subject list using Dixon fat and water images as input.
# Main steps:
# Read csv file with the test parameters and their ranges 
# and generate output text file with many commandlines.
# Read the output text file line-by-line:
#    make directory and cd to it
#    append the command line suffix to the prefix
#    run the experiment
#    cd up.


# Parse program arguments and set up output filenames.
use File::Basename;
use Cwd;

#system("source set_env_vars.sh");
if ( scalar(@ARGV) < 3 ) {
  print "Usage: ThymusDixonQuantificationDOE2.pl <DOE file> <SubjectFile> <Input Data Path> <Extra Arguments>\n";
  print "Thymus quantification using original Dixon fat and water images for segmentation.\n";
  exit;
}

# Initialize variables.
$outputdatapathPrefix = "./";
$parameterParser = "DOEParametersParser";
$thymusQuantifier = "ThymusDixonQuantificationCL2.pl";
$experimentPrefix = "Experiment_";

# Read command line.
$doeFile = shift @ARGV;
$subjectlistFileName = shift @ARGV;
$datapathPrefix = shift @ARGV;    #"/home/makrogianniss/Data/MRIConvertOutput";
$additionalArguments = " ";
foreach(@ARGV) {
 $additionalArguments =
 $additionalArguments . " " . $_;
}

# Run java code to generate command line suffices.
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$timestamp = sprintf("_%4d%02d%02d_%02d%02d%02d",
		     $year+1900,$mon+1,$mday,$hour,$min,$sec);
$textfileSuffix = ".txt";
$commandlineFile = "CommandLineParameters$timestamp". $textfileSuffix;

print "$parameterParser $doeFile $commandlineFile";
system("$parameterParser $doeFile $commandlineFile");

# Read inputs. 
print "Reading $commandlineFile\n";
open IN, "<$commandlineFile"; # Open output file.

# Read the subject list file.
$experimentCount = 1;
while(<IN>) {
    
    # For each subject in the list,
    $line = $_;
    #
    chomp($line);
    print "\n$line\n";

    # Make and cd to new directory.
    ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    $timeStamp = sprintf("_%4d%02d%02d_%02d%02d%02d",
			 $year+1900,$mon+1,$mday,$hour,$min,$sec);
    $outputDir = "$outputdatapathPrefix$experimentPrefix$timeStamp";
    print"\nmkdir $outputDir\n";
    system("mkdir $outputDir");
    print "\ncd $outputDir\n";
    # system("cd $outputDir");
    chdir("$outputDir") or die "$!";

    # Make info file with parameter values.
    open INFO, ">$outputDir.info.txt"; #or die $!;
    print INFO "Parameters: $line\n";
    close INFO;

    # Execute the CL script for all subjects with the defined parameters.
    $experimentCommand = "$thymusQuantifier $subjectlistFileName";
    $experimentCommand = $experimentCommand . " $datapathPrefix $line";
    $experimentCommand = $experimentCommand . $additionalArguments;
    print "\n$experimentCommand\n";
    system("$experimentCommand");

    # Go to original directory.
    chdir("..") or die "$!";;

    $experimentCount++;
}
