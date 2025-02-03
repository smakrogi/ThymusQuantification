#!/usr/bin/perl -w

# This script implements a pipeline for  
# registering all members of a subject list
# to a common template.
# The steps are:
# 1. Read subject list
# 2. Read template volume (original, mask)
# 3. For each subject in list:
#    4. read water, fat and mask volumes of the subject
#    5. crop volumes around seeds
#    6. linearly align subject to template (flirt, itk, etc)
#    7. store transformation matrix to file
#    8. apply transformation to intensity and mask volumes
#    9. add current registered volume to average intensity volume.
#   10. add current registered mask to average mask volume.
# 11. divide average intensity and mask volumes by the subject list length.


# Parse program arguments and set up output filenames.
use File::Basename;
use Cwd;

#system("source set_env_vars.sh");
if ( scalar(@ARGV) < 3 ) {
  print "Usage: ThymusDixonQuantificationCL.pl <SubjectFile> <Input Data Path> <Command line parameters>\n";
  print "Thymus quantification pipeline using fat and water ratio images as input to segmentation.\n";
   exit;
}

# Initialize variables.
$outputdatapathPrefix = "./";
# $headerinfoProgram = "fslinfo";
$cropProgram = "CropVolume";
$croppedSuffix = "_cropped";
$halfROI_Length = "65.0";  # in mm,
# $imagearithmeticProgram = "/home/makrogianniss/Software/ITK_Apps_x86_64_dyn_rel/ImageCalculator/ImageCalculator --add";
$fatwatercalculatorProgram = "FatWaterCalculator";
$thymusquantificationProgram = "ThymusQuantification";
$segmentationValidationProgram = "MaskOverlapAnalysis";
$formatExtension = ".nii";
$commentString = "#";

# Read command line.
$subjectlistFileName = shift @ARGV;
$datapathPrefix = shift @ARGV;    #"/home/makrogianniss/Data/MRIConvertOutput";

$thymusquantificationParameters = " ";
foreach(@ARGV) {
 $thymusquantificationParameters =
 $thymusquantificationParameters . " " . $_;
}
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
    # ($subjectID,$WSimageFile,$FSimageFile,$coordR,$coordA,$coordI ) = split(/ /, $line, 7);
    ($subjectID,$WSimageFile,$FSimageFile,$refMaskImageFile,$coordR,$coordA,$coordI) = 
	split(',', $line);
    #$subjectID = $line;

    if( $subjectID ne $commentString)
    {
	# Put together the filenames and directories.
	$originalWSFilename = "${datapathPrefix}/${subjectID}/${WSimageFile}";
	$originalWSBasename  = basename($WSimageFile, $formatExtension);
	$originalFSFilename = "${datapathPrefix}/${subjectID}/${FSimageFile}";
	$originalFSBasename  = basename($FSimageFile, $formatExtension);
	$refMaskFilename = "${datapathPrefix}/${subjectID}/${refMaskImageFile}";
	print "\n$subjectID\n"; 

	# Make directory named after the subject (if it doesn't exist).
	$outputPath = "${outputdatapathPrefix}${subjectID}";
	system("mkdir $outputPath");

	# Run (in the subject directory):

	# ## Display header info.
	# print "$headerinfoProgram $originalWSFilename\n";
	# system("$headerinfoProgram $originalWSFilename");
	# print "$headerinfoProgram $originalFSFilename\n";
	# system("$headerinfoProgram $originalFSFilename");
    
	# Crop FS and WS volumes.
	$croppedWSFilename = "$outputPath" . "/" . 
	    "${originalWSBasename}$croppedSuffix$formatExtension";
	$WSCroppingCommand = "$cropProgram  --outputdir $outputPath --output $croppedWSFilename --subjectid ${subjectID} --halfROIEdgeLength $halfROI_Length --seed $coordR,$coordA,$coordI $originalWSFilename";
	print "\n$WSCroppingCommand\n";
	system($WSCroppingCommand);
	$croppedFSFilename = "$outputPath" . "/" . 
	    "${originalFSBasename}$croppedSuffix$formatExtension";
	$FSCroppingCommand = "$cropProgram --output $croppedFSFilename --outputdir $outputPath --subjectid ${subjectID} --halfROIEdgeLength $halfROI_Length --seed $coordR,$coordA,$coordI $originalFSFilename";
	print "\n$FSCroppingCommand\n";
	system($FSCroppingCommand);

	# Compute Fat & water ratio images.
	$fatratioImage = "$outputPath" . "/" .
	    "${subjectID}$croppedSuffix"."FatRatio"."$formatExtension";
	$waterratioImage = "$outputPath" . "/" .
	    "${subjectID}$croppedSuffix"."WaterRatio"."$formatExtension";
	$fatandwaterImage = "$outputPath" . "/" .
	    "${subjectID}$croppedSuffix"."FatandWater"."$formatExtension";
	
	$fatwatercalculatorCommand = "$fatwatercalculatorProgram -r --outputdir $outputPath --subjectid ${subjectID} --outfatwaterimage $fatandwaterImage --outfatimage $fatratioImage --outwaterimage $waterratioImage $croppedWSFilename $croppedFSFilename";
	print "\n$fatwatercalculatorCommand\n";
	system("\n$fatwatercalculatorCommand\n");
 

    }

}

close IN;
