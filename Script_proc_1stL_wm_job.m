% MULTIPLE-SUBJECTS SCRIPT WORKING MEMORY (2021) DG

spm('defaults','fmri');

analysis_dir = '/home/uni10/nmri/projects/dgarnica/MRI_EEG_PROSPECTIVE/ROLANDIC'; % It is necessary to write this as char (to concatenate it with the rest of data as char)
cd(analysis_dir);
subjects_list = {'pseudonyms'};

%% Processing

for i = 1:length(subjects_list)
    clear matlabbatch;
    %spm_jobman('initcfg');
    
    if ~exist([analysis_dir '/' subjects_list{i} '/func/' 'Preproc_wm'])
        mkdir([analysis_dir '/' subjects_list{i} '/func/' 'Preproc_wm']);
    end    
    
    % check the ART preprocessing in existing files:
    load(fullfile(analysis_dir, subjects_list{i}, ['/conf/', 'EPI_fMRI_WM_PA_mcf_topped_nii', '/subjectinfo_EPI_fMRI_WM_PA_mcf_topped_nii.mat']));
    
%     switch subject.ART_threshold
%         case 'liberal'
%             ARTthreshold = 'liberal';
%         case 'intermediate'
%             ARTthreshold = 'intermediate';
%         case 'conservative'
%             ARTthreshold = 'conservative';
%         otherwise
%             ARTthreshold = 'individual';
%             warning('The ART threshold is individual, and not a default threshold from ART-Conn')
%     end
    
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.name = 'T1';
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {{[analysis_dir '/' subjects_list{i} '/anat/' subjects_list{i} '.hrT1.nii']}};
    matlabbatch{2}.cfg_basicio.file_dir.file_ops.cfg_named_file.name = 'Session1';
    matlabbatch{2}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {{[analysis_dir '/' subjects_list{i} '/func/' subjects_list{i} '.EPI_fMRI.WM.PA.mcf.topped' '.nii']}};
    % matlabbatch{2}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {{[analysis_dir '/' subjects_list{i} '/func/' 'art_' subjects_list{i} '.EPI_fMRI.WM.PA.mcf.topped_' , ARTthreshold, '.nii']}};
    matlabbatch{3}.cfg_basicio.file_dir.dir_ops.cfg_named_dir.name = 'Basedir';
    matlabbatch{3}.cfg_basicio.file_dir.dir_ops.cfg_named_dir.dirs = {{[analysis_dir '/' subjects_list{i} '/func/' 'Preproc_wm']}};
    matlabbatch{4}.cfg_basicio.file_dir.dir_ops.cfg_cd.dir(1) = cfg_dep('Named Directory Selector: Basedir(1)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','dirs', '{}',{1}));
    matlabbatch{5}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent(1) = cfg_dep('Named Directory Selector: Basedir(1)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','dirs', '{}',{1}));
    matlabbatch{5}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'MeanImages';
    matlabbatch{6}.spm.util.exp_frames.files(1) = cfg_dep('Named File Selector: Session1(1) - Files', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
    matlabbatch{6}.spm.util.exp_frames.frames = Inf;
    matlabbatch{7}.spm.util.imcalc.input(1) = cfg_dep('Expand image frames: Expanded filename list.', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{7}.spm.util.imcalc.output = 'Session1_Mean';
    matlabbatch{7}.spm.util.imcalc.outdir(1) = cfg_dep('Make Directory: Make Directory ''MeanImages''', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','dir'));
    matlabbatch{7}.spm.util.imcalc.expression = 'mean(X)';
    matlabbatch{7}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{7}.spm.util.imcalc.options.dmtx = 1;
    matlabbatch{7}.spm.util.imcalc.options.mask = 0;
    matlabbatch{7}.spm.util.imcalc.options.interp = 1;
    matlabbatch{7}.spm.util.imcalc.options.dtype = 4;
    matlabbatch{8}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Named File Selector: T1(1) - Files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
    matlabbatch{8}.spm.spatial.coreg.estimate.source(1) = cfg_dep('Image Calculator: ImCalc Computed Image: Session1_Mean', substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{8}.spm.spatial.coreg.estimate.other(1) = cfg_dep('Expand image frames: Expanded filename list.', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{8}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{8}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{8}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{8}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    matlabbatch{9}.spm.spatial.preproc.channel.vols(1) = cfg_dep('Named File Selector: T1(1) - Files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
    matlabbatch{9}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{9}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{9}.spm.spatial.preproc.channel.write = [0 1];
    matlabbatch{9}.spm.spatial.preproc.tissue(1).tpm = {'/home/uni10/nmri/projects/dgarnica/CerebroMatic/MainStudy_Segmentation_parameters/mw_com_prior_Age_0120.nii,1'};
    matlabbatch{9}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{9}.spm.spatial.preproc.tissue(1).native = [1 1];
    matlabbatch{9}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{9}.spm.spatial.preproc.tissue(2).tpm = {'/home/uni10/nmri/projects/dgarnica/CerebroMatic/MainStudy_Segmentation_parameters/mw_com_prior_Age_0120.nii,2'};
    matlabbatch{9}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{9}.spm.spatial.preproc.tissue(2).native = [1 1];
    matlabbatch{9}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{9}.spm.spatial.preproc.tissue(3).tpm = {'/home/uni10/nmri/projects/dgarnica/CerebroMatic/MainStudy_Segmentation_parameters/mw_com_prior_Age_0120.nii,3'};
    matlabbatch{9}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{9}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{9}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{9}.spm.spatial.preproc.tissue(4).tpm = {'/home/uni10/nmri/projects/dgarnica/CerebroMatic/MainStudy_Segmentation_parameters/mw_com_prior_Age_0120.nii,4'};
    matlabbatch{9}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{9}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{9}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{9}.spm.spatial.preproc.tissue(5).tpm = {'/home/uni10/nmri/projects/dgarnica/CerebroMatic/MainStudy_Segmentation_parameters/mw_com_prior_Age_0120.nii,5'};
    matlabbatch{9}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{9}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{9}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{9}.spm.spatial.preproc.tissue(6).tpm = {'/home/uni10/nmri/projects/dgarnica/CerebroMatic/MainStudy_Segmentation_parameters/mw_com_prior_Age_0120.nii,6'};
    matlabbatch{9}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{9}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{9}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{9}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{9}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{9}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{9}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{9}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{9}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{9}.spm.spatial.preproc.warp.write = [0 0];
    matlabbatch{10}.spm.tools.dartel.warp1.images{1}(1) = cfg_dep('Segment: rc1 Images', substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','rc', '()',{':'}));
    matlabbatch{10}.spm.tools.dartel.warp1.images{2}(1) = cfg_dep('Segment: rc2 Images', substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','rc', '()',{':'}));
    matlabbatch{10}.spm.tools.dartel.warp1.settings.rform = 0;
    matlabbatch{10}.spm.tools.dartel.warp1.settings.param(1).its = 3;
    matlabbatch{10}.spm.tools.dartel.warp1.settings.param(1).rparam = [4 2 1e-06];
    matlabbatch{10}.spm.tools.dartel.warp1.settings.param(1).K = 0;
    matlabbatch{10}.spm.tools.dartel.warp1.settings.param(1).template = {'/home/uni10/nmri/projects/dgarnica/CerebroMatic/MainStudy_DARTEL_templates/mw_com_Template_1_Age_0120.nii'};
    matlabbatch{10}.spm.tools.dartel.warp1.settings.param(2).its = 3;
    matlabbatch{10}.spm.tools.dartel.warp1.settings.param(2).rparam = [2 1 1e-06];
    matlabbatch{10}.spm.tools.dartel.warp1.settings.param(2).K = 0;
    matlabbatch{10}.spm.tools.dartel.warp1.settings.param(2).template = {'/home/uni10/nmri/projects/dgarnica/CerebroMatic/MainStudy_DARTEL_templates/mw_com_Template_2_Age_0120.nii'};
    matlabbatch{10}.spm.tools.dartel.warp1.settings.param(3).its = 3;
    matlabbatch{10}.spm.tools.dartel.warp1.settings.param(3).rparam = [1 0.5 1e-06];
    matlabbatch{10}.spm.tools.dartel.warp1.settings.param(3).K = 1;
    matlabbatch{10}.spm.tools.dartel.warp1.settings.param(3).template = {'/home/uni10/nmri/projects/dgarnica/CerebroMatic/MainStudy_DARTEL_templates/mw_com_Template_3_Age_0120.nii'};
    matlabbatch{10}.spm.tools.dartel.warp1.settings.param(4).its = 3;
    matlabbatch{10}.spm.tools.dartel.warp1.settings.param(4).rparam = [0.5 0.25 1e-06];
    matlabbatch{10}.spm.tools.dartel.warp1.settings.param(4).K = 2;
    matlabbatch{10}.spm.tools.dartel.warp1.settings.param(4).template = {'/home/uni10/nmri/projects/dgarnica/CerebroMatic/MainStudy_DARTEL_templates/mw_com_Template_4_Age_0120.nii'};
    matlabbatch{10}.spm.tools.dartel.warp1.settings.param(5).its = 3;
    matlabbatch{10}.spm.tools.dartel.warp1.settings.param(5).rparam = [0.25 0.125 1e-06];
    matlabbatch{10}.spm.tools.dartel.warp1.settings.param(5).K = 4;
    matlabbatch{10}.spm.tools.dartel.warp1.settings.param(5).template = {'/home/uni10/nmri/projects/dgarnica/CerebroMatic/MainStudy_DARTEL_templates/mw_com_Template_5_Age_0120.nii'};
    matlabbatch{10}.spm.tools.dartel.warp1.settings.param(6).its = 3;
    matlabbatch{10}.spm.tools.dartel.warp1.settings.param(6).rparam = [0.25 0.125 1e-06];
    matlabbatch{10}.spm.tools.dartel.warp1.settings.param(6).K = 6;
    matlabbatch{10}.spm.tools.dartel.warp1.settings.param(6).template = {'/home/uni10/nmri/projects/dgarnica/CerebroMatic/MainStudy_DARTEL_templates/mw_com_Template_6_Age_0120.nii'};
    matlabbatch{10}.spm.tools.dartel.warp1.settings.optim.lmreg = 0.01;
    matlabbatch{10}.spm.tools.dartel.warp1.settings.optim.cyc = 3;
    matlabbatch{10}.spm.tools.dartel.warp1.settings.optim.its = 3;
    matlabbatch{11}.spm.tools.dartel.crt_warped.flowfields(1) = cfg_dep('Run Dartel (existing Templates): Flow Fields', substruct('.','val', '{}',{10}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '()',{':'}));
    matlabbatch{11}.spm.tools.dartel.crt_warped.images{1}(1) = cfg_dep('Segment: Bias Corrected (1)', substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
    matlabbatch{11}.spm.tools.dartel.crt_warped.jactransf = 0;
    matlabbatch{11}.spm.tools.dartel.crt_warped.K = 6;
    matlabbatch{11}.spm.tools.dartel.crt_warped.interp = 1;
    matlabbatch{12}.spm.tools.dartel.crt_warped.flowfields(1) = cfg_dep('Run Dartel (existing Templates): Flow Fields', substruct('.','val', '{}',{10}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '()',{':'}));
    matlabbatch{12}.spm.tools.dartel.crt_warped.images{1}(1) = cfg_dep('Named File Selector: Session1(1) - Files', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
    matlabbatch{12}.spm.tools.dartel.crt_warped.jactransf = 0;
    matlabbatch{12}.spm.tools.dartel.crt_warped.K = 6;
    matlabbatch{12}.spm.tools.dartel.crt_warped.interp = 1;
    matlabbatch{13}.spm.spatial.smooth.data(1) = cfg_dep('Create Warped: Warped Images (Image 1)', substruct('.','val', '{}',{12}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '()',{':', 1}));
    matlabbatch{13}.spm.spatial.smooth.fwhm = [9 9 9];
    matlabbatch{13}.spm.spatial.smooth.dtype = 0;
    matlabbatch{13}.spm.spatial.smooth.im = 0;
    matlabbatch{13}.spm.spatial.smooth.prefix = 's';
    
    spm_jobman('run', matlabbatch);
end

%% Create Brainmask (Binarization 'i1 > 0.1', for GM,WM,CSF)

for i = 1:length(subjects_list) % create warp (in batch)
    clear matlabbatch;
    
    % warp c1 GM
    matlabbatch{1}.spm.tools.dartel.crt_warped.flowfields = {[analysis_dir '/' subjects_list{i} '/anat/' 'u_rc1' subjects_list{i} '.hrT1.nii']};
    matlabbatch{1}.spm.tools.dartel.crt_warped.images = {{[analysis_dir '/' subjects_list{i} '/anat/' 'c1' subjects_list{i} '.hrT1.nii']}};
    matlabbatch{1}.spm.tools.dartel.crt_warped.jactransf = 0;
    matlabbatch{1}.spm.tools.dartel.crt_warped.K = 6;
    matlabbatch{1}.spm.tools.dartel.crt_warped.interp = 1;
    
    % warp c2 WM
    matlabbatch{2}.spm.tools.dartel.crt_warped.flowfields = {[analysis_dir '/' subjects_list{i} '/anat/' 'u_rc1' subjects_list{i} '.hrT1.nii']};
    matlabbatch{2}.spm.tools.dartel.crt_warped.images = {{[analysis_dir '/' subjects_list{i} '/anat/' 'c2' subjects_list{i} '.hrT1.nii']}};
    matlabbatch{2}.spm.tools.dartel.crt_warped.jactransf = 0;
    matlabbatch{2}.spm.tools.dartel.crt_warped.K = 6;
    matlabbatch{2}.spm.tools.dartel.crt_warped.interp = 1;
    
    % warp c3 CSF
    matlabbatch{3}.spm.tools.dartel.crt_warped.flowfields = {[analysis_dir '/' subjects_list{i} '/anat/' 'u_rc1' subjects_list{i} '.hrT1.nii']};
    matlabbatch{3}.spm.tools.dartel.crt_warped.images = {{[analysis_dir '/' subjects_list{i} '/anat/' 'c3' subjects_list{i} '.hrT1.nii']}};
    matlabbatch{3}.spm.tools.dartel.crt_warped.jactransf = 0;
    matlabbatch{3}.spm.tools.dartel.crt_warped.K = 6;
    matlabbatch{3}.spm.tools.dartel.crt_warped.interp = 1;
    
    spm_jobman('run', matlabbatch);
end

for i = 1:length(subjects_list) % image calculator (in batch)
    clear matlabbatch;
    
    matlabbatch{1}.spm.util.imcalc.input = {
        [analysis_dir '/' subjects_list{i} '/anat/' 'wc1' subjects_list{i} '.hrT1.nii' ',1'] % GM (ALREADY WARPED)
        [analysis_dir '/' subjects_list{i} '/anat/' 'wc2' subjects_list{i} '.hrT1.nii' ',1'] % WM (ALREADY WARPED)
        [analysis_dir '/' subjects_list{i} '/anat/' 'wc3' subjects_list{i} '.hrT1.nii' ',1'] % CSF within skull (ALREADY WARPED)
        };
    matlabbatch{1}.spm.util.imcalc.output = 'brainmask_binarized';
    matlabbatch{1}.spm.util.imcalc.outdir = {[analysis_dir '/' subjects_list{i} '/anat/']};
    matlabbatch{1}.spm.util.imcalc.expression = 'i1 > 0.1';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = -7;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    
    spm_jobman('run', matlabbatch);
end

%% Model specification and 1st-level Analysis
cd(analysis_dir);
folder_scripts = '/home/uni10/nmri/projects/dgarnica/MRI_EEG_PROSPECTIVE/ROLANDIC/scripts';
addpath(genpath(folder_scripts));
fprintf('Current ROLANDIC scripts path added \n'); % IT MUST EXCLUDE scripts/utilities/afni_matlab

spm('defaults','fmri');

for i = 1:length(subjects_list)
    clear matlabbatch;
    
    % Calls function to extract the task and scan data (including movement params) - IT NEEDS THE SUBJECT STRUCTURE LOADED!
    structure = [analysis_dir '/' subjects_list{i} '/conf/' 'EPI_fMRI_WM_PA_mcf_topped_nii/' 'subjectinfo_EPI_fMRI_WM_PA_mcf_topped_nii.mat'];
    load(structure);
    [subject] = dg_child_taskdata_trim(subject); % BRINGS TASK DATA AND TRIMS DATASET (TRIMMING ALSO THE 6 MOVEMENT PARAMETERS)
    
    if ~exist([analysis_dir '/' subjects_list{i} '/func/' '1stLevel_wm'])
        mkdir([analysis_dir '/' subjects_list{i} '/func/' '1stLevel_wm']);
    end
    
    cd([analysis_dir '/' subjects_list{i} '/func/' '1stLevel_wm']);    
    matlabbatch{1}.spm.stats.fmri_spec.dir = {[analysis_dir '/' subjects_list{i} '/func/' '1stLevel_wm']};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.4;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    
    dataset = [analysis_dir '/' subjects_list{i} '/anat/' 'sw' subjects_list{i} '.EPI_fMRI.WM.PA.mcf.topped.nii'];
    volumes = spm_select('Expand',dataset);
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(volumes);
    
    if ~isempty(subject.scan_trimmed_length)
        warning('Total length of dataset is greater than the last scanned task timepoint! \n');
        fprintf(['Setting a trimmed length of the scan dataset for ', subjects_list{i},'! \n']);
        trimmed_dataset = volumes(subject.scan_trimmed_length,:); % it must index the whole width of chars with colon (:)
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(trimmed_dataset);
    elseif isempty(subject.scan_trimmed_length)
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(volumes);
    end
    % matlabbatch{1}.spm.stats.fmri_spec.sess.scans(1) =  cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{13}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    
    % BRINGING RESULTS_Encoding
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name = 'Encoding';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset = [subject.WMtask.Onsets_encoding];
    
    % BRINGING RESULTS_Duration_Encoding
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = 2;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).orth = 1;
    
    % BRINGING RESULTS_Retrieval
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name = 'Retrieval';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset = [subject.WMtask.Onsets_retrieval];
    
    % BRINGING RESULTS_Duration_Retrieval
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = 3;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).orth = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
    
    % BRINGING REGRESSORS
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    % matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {[analysis_dir '/' subjects_list{i} '/func/' 'mov_params_wm.mat']};
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    
    % BRINGING BRAINMASK
    matlabbatch{1}.spm.stats.fmri_spec.mask = {[analysis_dir '/' subjects_list{i} '/anat/' 'brainmask_binarized' '.nii']};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    % MODEL REVIEW
    matlabbatch{2}.spm.stats.review.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.review.display.matrix = 1;
    matlabbatch{2}.spm.stats.review.print = 'jpg';
    
    % MODEL ESTIMATION
    matlabbatch{3}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.fmri_est.write_residuals = 1;
    matlabbatch{3}.spm.stats.fmri_est.method.Classical = 1;
    
    % OPTIMIZED CENSORING TOOLBOX (M. WILKE) it does scrubbing and motion censoring, and brings the corresponding regressors to the model
    matlabbatch{4}.spm.tools.optcens_cfg.ospms(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{4}.spm.tools.optcens_cfg.approach = 1;
    matlabbatch{4}.spm.tools.optcens_cfg.docalc = 3;
    matlabbatch{4}.spm.tools.optcens_cfg.o_minrem = 0;
    matlabbatch{4}.spm.tools.optcens_cfg.o_maxrem = 20;
    matlabbatch{4}.spm.tools.optcens_cfg.zeropad = 0;
    matlabbatch{4}.spm.tools.optcens_cfg.adapt = 0;
    
    % CONTRAST MANAGER
    matlabbatch{5}.spm.stats.con.spmmat(1) = cfg_dep('Optimized Censoring Toolbox: Optimally censored & interpolated SPM.mat(s)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','both'));
    matlabbatch{5}.spm.stats.con.consess{1}.tcon.name = 'Encoding > Retrieval';
    matlabbatch{5}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
    matlabbatch{5}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{5}.spm.stats.con.consess{2}.tcon.name = 'Encoding < Retrieval';
    matlabbatch{5}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
    matlabbatch{5}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{5}.spm.stats.con.delete = 0;
    
    % RESULTS REPORT
    matlabbatch{6}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{6}.spm.stats.results.conspec(1).titlestr = '';
    matlabbatch{6}.spm.stats.results.conspec(1).contrasts = Inf;
    matlabbatch{6}.spm.stats.results.conspec(1).threshdesc = 'FWE';
    matlabbatch{6}.spm.stats.results.conspec(1).thresh = 0.05;
    matlabbatch{6}.spm.stats.results.conspec(1).extent = 0;
    matlabbatch{6}.spm.stats.results.conspec(1).conjunction = 1;
    matlabbatch{6}.spm.stats.results.conspec(1).mask.none = 1;
    matlabbatch{6}.spm.stats.results.conspec(2).titlestr = '';
    matlabbatch{6}.spm.stats.results.conspec(2).contrasts = Inf;
    matlabbatch{6}.spm.stats.results.conspec(2).threshdesc = 'none';
    matlabbatch{6}.spm.stats.results.conspec(2).thresh = 0.001;
    matlabbatch{6}.spm.stats.results.conspec(2).extent = 0;
    matlabbatch{6}.spm.stats.results.conspec(2).conjunction = 1;
    matlabbatch{6}.spm.stats.results.conspec(2).mask.none = 1;
    matlabbatch{6}.spm.stats.results.units = 1;
    matlabbatch{6}.spm.stats.results.export{1}.jpg = true;
    matlabbatch{6}.spm.stats.results.export{2}.ps = true;
    matlabbatch{6}.spm.stats.results.export{3}.csv = true;
    
    spm_jobman('run', matlabbatch);
end
