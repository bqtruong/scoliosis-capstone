% Directory containing all patient sensitive DICOM images
dicomDirName = 'C:\Users\Brian\Documents\School\4th Year\BME 4063\adult 1 output';

% Folder created under dicomDirName with anonymized DICOM images
outputFolder = 'output';

% Make output folder if it doesn't exist
mkdir(dicomDirName, outputFolder);

% Get all files within specified directory
dicomDir = dir(dicomDirName);

% Only get the names of the files
dicomFiles = extractfield(dicomDir, 'name');

% Only files with a '.dcm' extension
dicomFiles = dicomFiles(find(contains(dicomFiles, '.dcm')));

parfor x = 1:length(dicomFiles)
    % For every file, anonymize the DICOM, preserving specified fields if
    % at all, and output to subfolder
    old_file = fullfile(dicomDirName, dicomFiles{x});
    new_file = fullfile(dicomDirName, outputFolder, dicomFiles{x});
    original_dicom = dicomread(old_file);
    metadata = dicominfo(old_file);
    
    metadata.AccessionNumber = '';
    metadata.ReferringPhysicianName = '';
    metadata.ReferringPhysicianAddress = '';
    metadata.ReferringPhysicianTelephoneNumbers = '';
    metadata.ReferringPhysicianIdentificationSequence = '';
    metadata.PhysiciansOfRecord = '';
    metadata.PatientID = '';
    metadata.PatientName = '';
    metadata.PatientBirthDate = '';
    metadata.OtherPatientIDsSequence = '';
    metadata.PatientBirthName = '';
    metadata.PatientAddress = '';
    metadata.PatientMotherBirthName = '';
    metadata.PatientTelephoneNumbers = '';
    
    dicomwrite(original_dicom, new_file, metadata);
end