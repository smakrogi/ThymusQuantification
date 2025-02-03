#!/bin/bash


# convert from DICOM to nifti.
dcm2nii -c N -d N -i Y -n Y -p Y -f N -o ./ /home/makrogianniss/Data/Imagebox/1.3.46.670589.11.24071.5.0.5536.2009120911084340070

# merge (fat-suppressed) slabs into one big volume.
fslmerge -z  WIPDIXONSENSEMRI-09201.FS.nii WIPDIXONSENSEMRI-09201s901a1009_1.nii WIPDIXONSENSEMRI-09201s1001a1010_1.nii WIPDIXONSENSEMRI-09201s1101a1011_1.nii WIPDIXONSENSEMRI-09201s1201a1012_1.nii WIPDIXONSENSEMRI-09201s1301a1013_1.nii
