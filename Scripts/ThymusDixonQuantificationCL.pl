#!/usr/bin/perl -w

# This script implements an image processing 
# pipeline on a subject list of thymus Dixon MRI scans.
# The preprocessing stages are:
# 1. Display F and W volume info.
# 2. Compute fat and water ratio images.
# 3. Crop all volumes around seed voxel.
# (4. F+W)
# 5. Thymus segmentation, volumetrics and fat-water ratios.
# 6. Validation against manually segmented reference mask.
# S. Makrogiannis, Mar. 2011.

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


# $thymusquantificationParameters = " ";
# if($segmentationMethod eq "GEOD_ACT_CONT_3D")
# {
#     if($gradientMethod eq "PARZEN")
#    {
#	$thymusquantificationParameters = 
#	    $thymusquantificationParameters .
#	    " --parzband 80 "; # 0.1:fr&wr, 150:ws&fs, 0.075:fr, 80:ws
#    }
#    if($gradientMethod eq "CONVENTIONAL")
#    {
#	$thymusquantificationParameters = 
#	    $thymusquantificationParameters .
#	    "--lssigmoidbeta 6000 --lssigmoidalpha -80 " . #:ws&fs
#	    "--lsmaxrmserror 0.015";
#	    "--lssigmoidbeta 0.001 --lssigmoidalpha -0.001 " . #:fr&wr
#	    "--lsmaxrmserror 0.015";
#	    "--lssigmoidbeta 0.02 --lssigmoidalpha -0.008 " . #:fr
#	    "--lsmaxrmserror 0.015";
#	    "--lssigmoidbeta 40 --lssigmoidalpha -8 " . #:ws
#	    "--lsmaxrmserror 0.015";
#    }
# }

# Generate overall score files for all participants.
$scorefilePrefix = "AllScores";
$scorefilePrefix2 = "AllScoresValidation";
$scorefileSuffix = ".csv";
# ($temp1,$temp2,$temp3) = fileparse($subjectlistFileName, "\.txt");
fileparse_set_fstype("VMS");
$scorefilenameBase = basename($subjectlistFileName, ".txt");
($sec,$min,$hour,$mday,$mon,$year) = localtime(time);
$timestamp = sprintf("_%4d%02d%02d_%02d%02d%02d",
		     $year+1900,$mon+1,$mday,$hour,$min,$sec);
$scorefileName = "$outputdatapathPrefix$scorefilePrefix$scorefilenameBase$timestamp$scorefileSuffix";
# print "$scorefileName\n";
open OUTPUT, ">$scorefileName";
$firstlineString = "SubjectId,ThymusSize(px),ThymusVol(mm^3),WaterRatioMean,WaterRatioRMSE,FatRatioMean,FatRatioRMSE,$thymusquantificationParameters\n";
print OUTPUT $firstlineString;
$scorefileName2 = "$outputdatapathPrefix$scorefilePrefix2$scorefilenameBase$timestamp$scorefileSuffix";
# print "$scorefileName\n";
open OUTPUT2, ">$scorefileName2";
$firstlineString2 = "Subject_Id,Volume_groundtruth(mm^3),Volume_detection(mm^3),Volume_overlap(mm^3),Volume_union(mm^3),Overlap_score,DICE_score,Sensitivity(%),Specificity(%),Classification_accuracy(%),Volume_Differences(%),$thymusquantificationParameters\n";
print OUTPUT2 $firstlineString2;

$index = 1;


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
    
	# Segment the Thymus.
	# Compute volumes and water-fat ratios.
	$autosegmentationMask = "$outputPath" . "/" . "${subjectID}AutoMask$formatExtension";
	$thymusQuantificationCommand = "$thymusquantificationProgram  --subjectid ${subjectID} --seed $coordR,$coordA,$coordI --outputdir $outputPath --outmask $autosegmentationMask $thymusquantificationParameters $fatratioImage $waterratioImage";
	print "\n$thymusQuantificationCommand\n";
	system($thymusQuantificationCommand );
    
	# Compute overlap scores for segmentation validation.
	$validationCSVFile = "$outputPath" . "/" .
	    "${subjectID}OverlapScore.csv";
	$validationCommand = "$segmentationValidationProgram --subjectid ${subjectID} --csvOutputFile $validationCSVFile $autosegmentationMask $refMaskFilename ";
	print "\n$validationCommand\n";
	system($validationCommand);

	# Read volumetrics from each case and 
	# write them to a general file.
	$scorecardFile = "$outputPath"."/"."${subjectID}".".stats.csv";
	open SCORE, "<$scorecardFile"; #or die $!;
	@scorecard = <SCORE>;
	print "$scorecard[$index]\n";
	$localScore = $scorecard[$index];
	chomp($localScore);
	$outputString = "$localScore" . "\n";
	print OUTPUT "$outputString";
	print "\n";
	close SCORE;

	# Read overlap scores from each case and 
	# write them to a general file.
	$scorecardFile2 = $validationCSVFile;
	open SCORE2, "<$scorecardFile2"; #or die $!;
	@scorecard2 = <SCORE2>;
	print "$scorecard2[$index]\n";
	$localScore2 = $scorecard2[$index];
	chomp($localScore2);
	$outputString = "$localScore2" . "\n";
	print OUTPUT2 "$outputString";
	print "\n";
	close SCORE2;
	
    }

}

# print OUTPUT "$thymusquantificationParameters\n";
# print OUTPUT2 "$thymusquantificationParameters\n";

close OUTPUT;
close OUTPUT2;

close IN;

