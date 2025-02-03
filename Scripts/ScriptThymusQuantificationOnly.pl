#!/usr/bin/perl -w

#/*===========================================================================
#
# Program:   Quantification of thymus tissues from MRI.
# Module:    $RCSfile: ThymusQuantification.cxx,v $
# Language:  Perl
# Date:      $Date: 2010/04/23 10:42:32 $
# Version:   $Revision: 0.5 $
# Author:    S. K. Makrogiannis
# 3T MRI Facility National Institute on Aging/National Institutes of Health.
#
# =============================================================================*/


# Thymus volumetrics and fat-water ratios.
# S. Makrogiannis, Apr. 2010.

# Parse program arguments and set up output filenames.
use File::Basename;
use Cwd;

#system("source set_env_vars.sh");
if ( scalar(@ARGV) < 4 ) {
  print "Usage: ScriptThymusQuantificationOnly.pl <SubjectFile> 
<Input Data Path> <Gradient Algo> <Segmentation Algo>\n";
  exit;
}

# Initialize variables.
$outputdatapathPrefix = "./";
# $headerinfoProgram = "fslinfo";
# $cropProgram = "/home/makrogianniss/Software/CropVolume_Rel_Dyn_Bin/bin/CropVolume";
# $n3shadingcorrectionProgram = "N3BiasCorrector";
$N3correctionSuffix = "_n3_corrected";
# $n3scaleFactor = 2; # previous: 2
# $n3defaultparameters = "mask.img 50 4 ";
# $biasfieldSuffix = "_bias_field";
# $registrationProgram = "flirt";
$croppedSuffix = "_cropped";
$registeredSuffix = "_NS_aligned";
# $registrationParameters = " -cost normmi -dof 3 -nosearch"; # -cost mutualinfo,normmi -dof 6 -bins 128 -nosearch
# $halfROI_Length = "40.0";  # in mm, 20.0
# $imagearithmeticProgram = "/home/makrogianniss/Software/ITK_Apps_x86_64_dyn_rel/ImageCalculator/ImageCalculator --add";
$thymusquantificationProgram = "/home/makrogianniss/Software/ThymusQuantification_Rel_Dyn_Bin/bin/ThymusQuantification";
$thymusquantificationParameters = " ";

# Read command line.
$subjectlistFileName = shift @ARGV;
$datapathPrefix = shift @ARGV;    #"/home/makrogianniss/Data/MRIConvertOutput";
$gradientMethod = shift @ARGV;
$segmentationMethod = shift @ARGV;

# Suggested settings for level sets.
# 3D Parzen.
$thymusquantificationParameters = " ";
if($segmentationMethod eq "LEVEL_SETS_3D")
{
    if($gradientMethod eq "PARZEN")
    {
	$thymusquantificationParameters = 
	    "--roiradius 1.5 --lscurvscale 0.15 --lspropscale 1.2 ".
	    " --parzrad 2 --parzband 75 ".
	    "--lssigmoidbeta 20 --lssigmoidalpha -4 ";
    }
    if($gradientMethod eq "CONVENTIONAL")
    {
	$thymusquantificationParameters = 
	    "--roiradius 1.5 --lscurvscale 0.15 --lspropscale 1.2 ".
	    "--lssigmoidbeta 600 --lssigmoidalpha -80 ";
    }
}
if($segmentationMethod eq "WS_LEVEL_SETS")
{
    
    # 3D Gradient.
    # WS Parzen or Gradient.
    $thymusquantificationParameters = 
	"--roiradius 1.5 --lscurvscale 0.15 --lspropscale 1.2 ".
	" --parzrad 2 --parzband 75 ".
	"--lssigmoidbeta 15 --lssigmoidalpha -3 ";
}

$formatExtension = ".nii";
$commentString = "#";


# Generate overall score file for all participants.
$scorefilePrefix = "AllScores";
$scorefileSuffix = ".csv";
# ($temp1,$temp2,$temp3) = fileparse($subjectlistFileName, "\.txt");
fileparse_set_fstype("VMS");
$scorefilenameBase = basename($subjectlistFileName, ".txt");
$scorefileName = "$outputdatapathPrefix$scorefilePrefix$scorefilenameBase".
    ".$gradientMethod.$segmentationMethod$scorefileSuffix";
# print "$scorefileName\n";
open OUTPUT, ">$scorefileName";
$firstlineString = "Gradient Method: $gradientMethod, Segmentation Method: $segmentationMethod\n";
$secondlineString = "SubjectId,ThymusSize(px),ThymusVol(mm^3),WaterRatioMean,".
    "WaterRatioRMS,FatRatioMean,FatRatioRMS\n";
print OUTPUT $firstlineString;
print OUTPUT $secondlineString;
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
    # ($subjectID, $NSimagePath, $NSimageFile, $FSimagePath, $FSimageFile, $WSimagePath, $WSimageFile) = split(/ /, $line, 7);
	($subjectID,$NSimageFile,$WSimageFile,$FSimageFile,$coordR,$coordA,$coordI) = split(',', $line);
    #$subjectID = $line;

    if( $subjectID ne $commentString)
    {

	# Put together the filenames and directories.
	# $dataPath = "${datapathPrefix}/${subjectID}";
	$originalNSFilename = "${datapathPrefix}/${subjectID}/${NSimageFile}";
	$originalWSFilename = "${datapathPrefix}/${subjectID}/${WSimageFile}";
	$originalFSFilename = "${datapathPrefix}/${subjectID}/${FSimageFile}";
	print "\n$subjectID\n"; 

	# Make directory named after the subject (if it doesn't exist).
	$outputPath = "${outputdatapathPrefix}${subjectID}";
	system("mkdir $outputPath");

	# Run (in the subject directory):

	# ## Display header info.
	# print "$headerinfoProgram $originalNSFilename\n";
	# system("$headerinfoProgram $originalNSFilename");
	# print "$headerinfoProgram $originalWSFilename\n";
	# system("$headerinfoProgram $originalWSFilename");
	# print "$headerinfoProgram $originalFSFilename\n";
	# system("$headerinfoProgram $originalFSFilename");

	# Segment the Thymus.
	# Compute volumes and water-fat ratios.
	# $thymusMaskVolume = "$outputPath"."/"."${subjectID}"."_ThymusMask"."${formatExtension}";
	# $thymusMaskVolumeNoPath = "${subjectID}"."_ThymusMask"."${formatExtension}";
	# $N3correctedNSFilenameNoPath = "${originalNSFilename}$croppedSuffix${N3correctionSuffix}$formatExtension";
	# $registeredWSvolumeNoPath = "${originalWSFilename}$croppedSuffix${N3correctionSuffix}$registeredSuffix$formatExtension";
	# $registeredFSvolumeNoPath = "${originalFSFilename}$croppedSuffix${N3correctionSuffix}$registeredSuffix$formatExtension";

	$N3correctedNSFilename = "${originalNSFilename}$croppedSuffix".
	    "${N3correctionSuffix}$formatExtension";
	$registeredWSvolume = "${originalWSFilename}".
	    "$croppedSuffix${N3correctionSuffix}$registeredSuffix$formatExtension";
	$registeredFSvolume = "${originalFSFilename}$croppedSuffix".
	    "${N3correctionSuffix}$registeredSuffix$formatExtension";

	$thymusquantificationCommand = 
	    "$thymusquantificationProgram --subjectid ${subjectID} " .
	    "--nonsuppressed $N3correctedNSFilename ".
	    "--watersuppressed $registeredWSvolume " .
	    "--fatsuppressed $registeredFSvolume --seedpoint $coordA $coordR ".
	    "$coordI --outputpath $outputPath --segmentalgo $segmentationMethod ".
	    "--gradalgo $gradientMethod $thymusquantificationParameters";
	print "$thymusquantificationCommand\n";
	system($thymusquantificationCommand);

	# Read classification scores from each case and 
	# write them to the general file.

	$scorecardFile = "$outputPath"."/"."${subjectID}"."_".
	    "$gradientMethod"."_"."$segmentationMethod".".stats.csv";
	open SCORE, "<$scorecardFile"; #or die $!;
	@scorecard = <SCORE>;
	print "$scorecard[$index]\n";
	$localScore = $scorecard[$index];
	chomp($localScore);

	# $outputString = "${outputString}" . "\n";
	$outputString = "$localScore" . "\n";
	# print "\n";
	# print "$firstlineString";
	# print "$outputString";
	print OUTPUT "$outputString";
	print "\n";
	close SCORE;
	#system ("rm $scorecardFile");
	}

}

close OUTPUT;

close IN;

