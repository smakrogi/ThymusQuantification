#!/usr/bin/perl -w

# This script implements an image pre-processing 
# pipeline on a subject list.
# The preprocessing stages are:
# 1. Display NS,FS and WS volume info.
# 2. Crop around all volumes seed voxel.
# 3. Shading correction.
# 4. Co-register FS and WS to NS.
# S. Makrogiannis, Apr. 2010.

# Parse program arguments and set up output filenames.
use File::Basename;

#system("source set_env_vars.sh");
if ( scalar(@ARGV) < 2 ) {
  print "Usage: Preprocessing.pl <SubjectFile> <Input Data Path>\n";
  exit;
}

# Initialize variables.
$outputdatapathPrefix = "./";
$headerinfoProgram = "fslinfo";
$cropProgram = "/home/makrogianniss/Software/CropVolume_Rel_Dyn_Bin/bin/CropVolume";
$n3shadingcorrectionProgram = "/home/makrogianniss/Software/N4ITK_x86_64_dyn_rel/itkN3MRIBiasFieldCorrectionImageFilterTest";
$N3correctionSuffix = "_n3_corrected";
$n3scaleFactor = 2; # previous: 2
$n3defaultparameters = "mask.img 50 4 ";
$biasfieldSuffix = "_bias_field";
$registrationProgram = "flirt";
$croppedSuffix = "_cropped";
$registeredSuffix = "_NS_aligned";
$registrationParameters = " -cost normmi -dof 3 -nosearch "; # -cost mutualinfo,normmi -dof 6 -bins 128 -nosearch
$halfROI_Length = "40.0";  # in mm, 20.0
$imagearithmeticProgram = "/home/makrogianniss/Software/ITK_Apps_x86_64_dyn_rel/ImageCalculator/ImageCalculator --add";
$formatExtension = ".nii";
$commentString = "#";

# Read command line.
$subjectlistFileName = shift @ARGV;
$datapathPrefix = shift @ARGV;    #"/home/makrogianniss/Data/MRIConvertOutput";

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
	($subjectID,$NSimageFile,$FSimageFile,$WSimageFile,$coordR,$coordA,$coordI) = split(',', $line);
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

	## Display header info.
	print "$headerinfoProgram $originalNSFilename\n";
	system("$headerinfoProgram $originalNSFilename");
	print "$headerinfoProgram $originalWSFilename\n";
	system("$headerinfoProgram $originalWSFilename");
	print "$headerinfoProgram $originalFSFilename\n";
	system("$headerinfoProgram $originalFSFilename");

	# Crop NS, FS and WS volumes.
	# CropVolume ~/Data/DCM2NII_Output/20091215_141023MRI-08002/20091215_141023WIPeTHRIVEFSSENSEMRI-08002s1901a1019.nii -0.369535029 -73.5488663 32.0585556 20.0 20091215_141023WIPeTHRIVEFSSENSEMRI-08002s1901a1019_cropped.nii
	$croppedNSFilename = "$outputPath" . "/" . "${NSimageFile}$croppedSuffix$formatExtension";
	print "$cropProgram $originalNSFilename $coordA $coordR $coordI $halfROI_Length $croppedNSFilename";
	system("$cropProgram $originalNSFilename $coordA $coordR $coordI $halfROI_Length $croppedNSFilename");
	
	$croppedWSFilename = "$outputPath" . "/" . "${WSimageFile}$croppedSuffix$formatExtension";
	system("$cropProgram $originalWSFilename $coordA $coordR $coordI $halfROI_Length $croppedWSFilename");

	$croppedFSFilename = "$outputPath" . "/" . "${FSimageFile}$croppedSuffix$formatExtension";
	system("$cropProgram $originalFSFilename $coordA $coordR $coordI $halfROI_Length $croppedFSFilename");
	
    # Try to correct shading.
    # N3BiasCorrector 3 BLSA_5898_06_312283400_301_T1_3D_80_SLICE_20091125_cropped.nii.gz n3_corrected_cropped.nii
	$N3correctedNSFilename = "$outputPath" . "/" . "${NSimageFile}${croppedSuffix}${N3correctionSuffix}" . $formatExtension;
	$N3biasfieldNSFilename = "$outputPath" . "/" . "${NSimageFile}${croppedSuffix}${N3correctionSuffix}${biasfieldSuffix}" . $formatExtension;
    print "$n3shadingcorrectionProgram 3 $croppedNSFilename $N3correctedNSFilename $n3scaleFactor $n3defaultparameters $N3biasfieldNSFilename\n";
    system("$n3shadingcorrectionProgram 3 $croppedNSFilename $N3correctedNSFilename $n3scaleFactor $n3defaultparameters $N3biasfieldNSFilename");

    $N3correctedWSFilename = "$outputPath" . "/" . "${WSimageFile}${croppedSuffix}${N3correctionSuffix}" . ".nii";
	$N3biasfieldWSFilename = "$outputPath" . "/" . "${WSimageFile}${croppedSuffix}${N3correctionSuffix}${biasfieldSuffix}" . ".nii";
    print "$n3shadingcorrectionProgram 3 $croppedWSFilename $N3correctedWSFilename $n3scaleFactor $n3defaultparameters $N3biasfieldWSFilename\n";
    system("$n3shadingcorrectionProgram 3 $croppedWSFilename $N3correctedWSFilename $n3scaleFactor $n3defaultparameters $N3biasfieldWSFilename");

    $N3correctedFSFilename = "$outputPath" . "/" . "${FSimageFile}${croppedSuffix}${N3correctionSuffix}" . ".nii";
	$N3biasfieldFSFilename = "$outputPath" . "/" . "${FSimageFile}${croppedSuffix}${N3correctionSuffix}${biasfieldSuffix}" . ".nii";
    print "$n3shadingcorrectionProgram 3 $croppedFSFilename $N3correctedFSFilename $n3scaleFactor $n3defaultparameters $N3biasfieldFSFilename\n";
    system("$n3shadingcorrectionProgram 3 $croppedFSFilename $N3correctedFSFilename $n3scaleFactor $n3defaultparameters $N3biasfieldFSFilename");

	## Register WS and FS images to the NS space.

	## WS -> NS
	# flirt -cost mutualinfo -dof 3 -nosearch -v -in 20091215_141023WIPeTHRIVEWSSENSEMRI-08002s1801a1018_cropped.nii -ref 20091215_141023WIPeTHRIVENOSUPPSENSEMRI-08002s2001a1020_cropped.nii -out 20091215_141023WIPeTHRIVEWSSENSEMRI-08002s1801a1018_cropped_NS_aligned.nii
	$registeredWSvolume = "$outputPath" . "/" . "${WSimageFile}$croppedSuffix${N3correctionSuffix}$registeredSuffix";
	$registrationWSMatrix = "$registeredWSvolume.matr";
	print "$registrationProgram $registrationParameters -in $N3correctedWSFilename -ref $N3correctedNSFilename -out $registeredWSvolume -omat $registrationWSMatrix\n";
    system("$registrationProgram $registrationParameters -in $N3correctedWSFilename -ref $croppedNSFilename -out $registeredWSvolume -omat $registrationWSMatrix");
	print "gunzip -f $registeredWSvolume.nii.gz\n";
	system("gunzip -f $registeredWSvolume.nii.gz");

	## FS -> NS
	# flirt -cost mutualinfo -dof 3 -nosearch -v -in 20091215_141023WIPeTHRIVEFSSENSEMRI-08002s1901a1019_cropped -ref 20091215_141023WIPeTHRIVENOSUPPSENSEMRI-08002s2001a1020_cropped -out 20091215_141023WIPeTHRIVEFSSENSEMRI-08002s1901a1019_cropped_NS_aligned.nii
	$registeredFSvolume = "$outputPath" . "/" . "${FSimageFile}$croppedSuffix${N3correctionSuffix}$registeredSuffix";
	$registrationFSMatrix = "$registeredFSvolume.matr";
    print "$registrationProgram $registrationParameters -in $N3correctedFSFilename -ref $N3correctedNSFilename -out $registeredFSvolume -omat $registrationFSMatrix\n";
    system("$registrationProgram $registrationParameters -in $N3correctedFSFilename -ref $N3correctedNSFilename -out $registeredFSvolume -omat $registrationFSMatrix");
	print "gunzip -f $registeredFSvolume.nii.gz\n";
	system("gunzip -f $registeredFSvolume.nii.gz");
	
	# Add registered WS and FS images.
	#  ~/Software/ITK_Apps_x86_64_dyn_rel/ImageCalculator/ImageCalculator --add --in 20091215_141023WIPeTHRIVEWSSENSEMRI-08002s1801a1018_cropped_NS_aligned.nii 20091215_141023WIPeTHRIVEFSSENSEMRI-08002s1901a1019_cropped_NS_aligned.nii --out sum.nii
	$summedVolume = "$outputPath"."/"."${subjectID}"."_WS_FS_sum"."${formatExtension}";
	print "$imagearithmeticProgram --in $registeredWSvolume $registeredFSvolume --out $summedVolume";
	system("$imagearithmeticProgram --in $registeredWSvolume $registeredFSvolume --out $summedVolume");	
	
    }
}


close IN;

