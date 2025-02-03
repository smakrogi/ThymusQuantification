#!/usr/bin/perl -w

# This script implements an image pre-processing 
# pipeline on a subject list.
# The preprocessing stages are:
# 1. Display NS,FS and WS volume info.
# 2. Crop all volumes around seed voxel.
# 3. Co-register FS and WS to NS.
# 4. WS+FS
# 5. Thymus volumetrics and fat-water ratios.
# S. Makrogiannis, Apr. 2010.

# Parse program arguments and set up output filenames.
use File::Basename;
use Cwd;

#system("source set_env_vars.sh");
if ( scalar(@ARGV) < 4 ) {
  print "Usage: ScriptThymusQuantification.pl <SubjectFile> <Input Data Path> <Gradient Algo> <Segmentation Algo>\n";
  exit;
}

# Initialize variables.
$outputdatapathPrefix = "./";
# $headerinfoProgram = "fslinfo";
$cropProgram = "/home/makrogianniss/Software/CropVolume_Rel_Dyn_Bin/bin/CropVolume";
$registrationProgram = "flirt";
$croppedSuffix = "_cropped";
$registeredSuffix = "_NS_aligned";
$registrationParameters = " -cost normmi -dof 3 -nosearch"; # -cost mutualinfo,normmi -dof 6 -bins 128 -nosearch
$halfROI_Length = "40.0";  # in mm, 20.0
$imagearithmeticProgram = "/home/makrogianniss/Software/ITK_Apps_x86_64_dyn_rel/ImageCalculator/ImageCalculator --add";
$thymusquantificationProgram = "/home/makrogianniss/Software/ThymusQuantification_Rel_Dyn_Bin/bin/ThymusQuantification";
# $thymusquantificationParameters = " --roiradius 1.0 ";
# $thymusquantificationParameters = " --gradalgo CONVENTIONAL 
# --segmentalgo LEVEL_SETS_3D  --roiradius 1.5 --lssigmoidbeta 300 
# --lscurvscale 0.01 --lsmaxit 250 --lspropscale 1.5 ";
# $thymusquantificationParameters = " --gradalgo PARZEN 
# --segmentalgo WS_LEVEL_SETS --roiradius 2.0 --lssigmoidbeta 50 
# --lssigmoidalpha -15 --lscurvscale 0.1 --lsmaxit 250 --lspropscale 1.5 ";
$formatExtension = ".nii";
$commentString = "#";

# Read command line.
$subjectlistFileName = shift @ARGV;
$datapathPrefix = shift @ARGV;    #"/home/makrogianniss/Data/MRIConvertOutput";
$gradientMethod = shift @ARGV;
$segmentationMethod = shift @ARGV;

# Generate overall score file for all participants.
$scorefilePrefix = "AllScores";
$scorefileSuffix = ".csv";
# ($temp1,$temp2,$temp3) = fileparse($subjectlistFileName, "\.txt");
fileparse_set_fstype("VMS");
$scorefilenameBase = basename($subjectlistFileName, ".txt");
$scorefileName = "$outputdatapathPrefix$scorefilePrefix$scorefilenameBase".".$gradientMethod.$segmentationMethod$scorefileSuffix";
# print "$scorefileName\n";
open OUTPUT, ">$scorefileName";
$firstlineString = "SubjectId,ThymusSize(px),ThymusVol(mm^3),WaterRatioMean,WaterRatioRMS,FatRatioMean,FatRatioRMS\n";
print OUTPUT $firstlineString;
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
    ($subjectID,$NSimageFile,$WSimageFile,$FSimageFile,$coordR,$coordA,$coordI) = 
	split(',', $line);
    #$subjectID = $line;

    if( $subjectID ne $commentString)
    {

	# Put together the filenames and directories.
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
    
    # Crop NS, FS and WS volumes.
    # CropVolume ~/Data/DCM2NII_Output/20091215_141023MRI-08002/20091215_141023WIPeTHRIVEFSSENSEMRI-08002s1901a1019.nii -0.369535029 -73.5488663 32.0585556 20.0 20091215_141023WIPeTHRIVEFSSENSEMRI-08002s1901a1019_cropped.nii
    $croppedNSFilename = "$outputPath" . "/" . 
	"${NSimageFile}$croppedSuffix$formatExtension";
    print "$cropProgram $originalNSFilename $coordA $coordR $coordI 
$halfROI_Length $croppedNSFilename";
    system("$cropProgram $originalNSFilename $coordA $coordR $coordI 
$halfROI_Length $croppedNSFilename");
    
    $croppedWSFilename = "$outputPath" . "/" . 
	"${WSimageFile}$croppedSuffix$formatExtension";
    system("$cropProgram $originalWSFilename $coordA $coordR $coordI 
$halfROI_Length $croppedWSFilename");
    
    $croppedFSFilename = "$outputPath" . "/" . 
	"${FSimageFile}$croppedSuffix$formatExtension";
    system("$cropProgram $originalFSFilename $coordA $coordR $coordI 
$halfROI_Length $croppedFSFilename");
    
    ## Register WS and FS images to the NS space.
    
    ## WS -> NS
    # flirt -cost mutualinfo -dof 3 -nosearch -v 
# -in 20091215_141023WIPeTHRIVEWSSENSEMRI-08002s1801a1018_cropped.nii 
# -ref 20091215_141023WIPeTHRIVENOSUPPSENSEMRI-08002s2001a1020_cropped.nii 
# -out 20091215_141023WIPeTHRIVEWSSENSEMRI-08002s1801a1018_cropped_NS_aligned.nii
    $registeredWSvolume = "$outputPath" . "/" . 
	"${WSimageFile}$croppedSuffix$registeredSuffix";
    $registrationWSMatrix = "$registeredWSvolume.matr";
    print "$registrationProgram $registrationParameters -in $croppedWSFilename 
-ref $croppedNSFilename -out $registeredWSvolume -omat $registrationWSMatrix\n";
    system("$registrationProgram $registrationParameters -in $croppedWSFilename 
-ref $croppedNSFilename -out $registeredWSvolume -omat $registrationWSMatrix");
    print "gunzip -f $registeredWSvolume.nii.gz\n";
    system("gunzip -f $registeredWSvolume.nii.gz");
    
    ## FS -> NS
    # flirt -cost mutualinfo -dof 3 -nosearch -v 
# -in 20091215_141023WIPeTHRIVEFSSENSEMRI-08002s1901a1019_cropped 
# -ref 20091215_141023WIPeTHRIVENOSUPPSENSEMRI-08002s2001a1020_cropped 
# -out 20091215_141023WIPeTHRIVEFSSENSEMRI-08002s1901a1019_cropped_NS_aligned.nii
    $registeredFSvolume = "$outputPath" . "/" . 
	"${FSimageFile}$croppedSuffix$registeredSuffix";
    $registrationFSMatrix = "$registeredFSvolume.matr";
    print "$registrationProgram $registrationParameters -in $croppedFSFilename 
-ref $croppedNSFilename -out $registeredFSvolume -omat $registrationFSMatrix\n";
    system("$registrationProgram $registrationParameters -in $croppedFSFilename 
-ref $croppedNSFilename -out $registeredFSvolume -omat $registrationFSMatrix");
    print "gunzip -f $registeredFSvolume.nii.gz\n";
    system("gunzip -f $registeredFSvolume.nii.gz");
	
    # Add registered WS and FS images.
    #  ~/Software/ITK_Apps_x86_64_dyn_rel/ImageCalculator/ImageCalculator 
# --add --in 20091215_141023WIPeTHRIVEWSSENSEMRI-08002s1801a1018_cropped_NS_aligned.nii 
# 20091215_141023WIPeTHRIVEFSSENSEMRI-08002s1901a1019_cropped_NS_aligned.nii 
# --out sum.nii
    $summedVolume = "$outputPath"."/"."${subjectID}"."_WS_FS_sum".
	"${formatExtension}";
    print "$imagearithmeticProgram --in $registeredWSvolume $registeredFSvolume 
--out $summedVolume";
    system("$imagearithmeticProgram --in $registeredWSvolume $registeredFSvolume 
--out $summedVolume");
    
    # Segment the Thymus.
    # Compute volumes and water-fat ratios.
    # $thymusMaskVolume = "$outputPath"."/"."${subjectID}"."_ThymusMask"."${formatExtension}";
    # $thymusMaskVolumeNoPath = "${subjectID}"."_ThymusMask"."${formatExtension}";
    # $croppedNSFilenameNoPath = "${NSimageFile}$croppedSuffix$formatExtension";
    # $registeredWSvolumeNoPath = "${WSimageFile}$croppedSuffix$registeredSuffix";
    # $registeredFSvolumeNoPath = "${FSimageFile}$croppedSuffix$registeredSuffix";

    print "$thymusquantificationProgram --subjectid ${subjectID} 
--nonsuppressed $croppedNSFilename --watersuppressed $registeredWSvolume 
--fatsuppressed $registeredFSvolume --seedpoint $coordA $coordR $coordI 
--outputpath $outputPath --segmentalgo $segmentationMethod 
--gradalgo $gradientMethod $thymusquantificationParameters ";
    system("$thymusquantificationProgram --subjectid ${subjectID} 
--nonsuppressed $croppedNSFilename --watersuppressed $registeredWSvolume 
--fatsuppressed $registeredFSvolume --seedpoint $coordA $coordR $coordI 
--outputpath $outputPath --segmentalgo $segmentationMethod 
--gradalgo $gradientMethod $thymusquantificationParameters");
    
    # Read classification scores from each case and 
    # write them to the general file.

    # $scorecardFile = "$outputPath"."/"."${subjectID}".".stats.csv";
    $scorecardFile = "$outputPath"."/"."${subjectID}"."_".
	"$gradientMethod"."_"."$segmentationMethod".".stats.csv";
    open SCORE, "<$scorecardFile"; #or die $!;
    @scorecard = <SCORE>;
    print "$scorecard[$index]\n";
    $localScore = $scorecard[$index];
    chomp($localScore);
    # $outputString = sprintf '%s,%s,%s,%s,%s,%s', $subjectID, 
    # $lesionType, $coordX, $coordY, $coordZ, $localScore;
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

