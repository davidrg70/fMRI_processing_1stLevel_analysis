% normalise T1 images of all the subjects indicated...
% normalisation to MNI using SPM warping fields

% based on https://neuroimaging-core-docs.readthedocs.io/en/latest/pages/spm_anat_normalization-rev.html

clear all;
addpath('/work/users/d/g/dga/tools/spm12_7487/');

analysis_dir = '/work/users/d/g/dga/SeLECTS_BIDS/';
scripts_dir = '/work/users/d/g/dga/code/spm/';
cd(analysis_dir);

subjects_list = {'pseudonyms'};

%% loop

for i = 1:length(subjects_list)
    clear matlabbatch; clear t1_file;
    spm('defaults','fmri');
    spm_jobman('initcfg');

    t1_file = fullfile(analysis_dir,subjects_list{i},'anat',[subjects_list{i},'_T1w.nii,1']); % NOTE THE ,1

    matlabbatch{1}.spm.spatial.preproc.channel.vols = {t1_file}; % unprocessed T1
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'/users/d/g/dga/tools/spm12_7487/tpm/TPM.nii,1'};
    matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'/users/d/g/dga/tools/spm12_7487/tpm/TPM.nii,2'};
    matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'/users/d/g/dga/tools/spm12_7487/tpm/TPM.nii,3'};
    matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'/users/d/g/dga/tools/spm12_7487/tpm/TPM.nii,4'};
    matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'/users/d/g/dga/tools/spm12_7487/tpm/TPM.nii,5'};
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'/users/d/g/dga/tools/spm12_7487/tpm/TPM.nii,6'};
    matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write = [0 1];
    matlabbatch{2}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
    matlabbatch{2}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Segment: Bias Corrected (1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
    matlabbatch{2}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85]; % bounding box MNI dimensions, following https://www.fil.ion.ucl.ac.uk/spm/docs/tutorials/fmri_preprocessing/scripting/
    matlabbatch{2}.spm.spatial.normalise.write.woptions.vox = [1 1 1]; % By choosing [1 1 1] instead of [2 2 2], we are retaining the higher resolution and tissue contrast of the original image. This assumes the original image is higher resolution.
    matlabbatch{2}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{2}.spm.spatial.normalise.write.woptions.prefix = 'w';

    spm_jobman('run', matlabbatch);
end