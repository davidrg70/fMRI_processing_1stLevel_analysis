% Written by David Garnica, dga@email.unc.edu
% 2025, UNC Psychiatry 2025
% ALL-SUBJECTS SCRIPT TO PROCESS AND RUN 1ST-LEVEL ANALYSIS OF --REST--
clear all;
close all;
addpath('/work/users/d/g/dga/tools/spm12_7487/');

analysis_dir = '/work/users/d/g/dga/SeLECTS_BIDS/';
cd(analysis_dir);
subjects_list = {'subs-ids'};
subjects_ids = {'pseudonyms'};

TR = 1.4;
FWHM = 5; % mm - smoothing kernel
Trim = 1; % 1 for trimming (when dataset is longer than task timepoints)
          % 0 for no trimming
          % I KEEP 1 AS DEFAULT (REASON: trim whenever it's necessary to avoid movement outliers to be included in postprocessing but out of the task range!)
FD_threshold = 0.2; % 0.3mm used ORIGINALLY by mw_optcens.m script, lines 100 and 101!
    % (0.2 mm is stringent after Fair et al. 2013 https://doi.org/10.3389/fnsys.2012.00080)
r = 50; % mm % following the fMRI-Quality-Checker 'calculateQCmeasures.m' script
% it is a scaling factor that converts rotational displacements from radians to millimeters
% (typically, 𝑟 = 50 mm, assuming an average head radius of 50 mm).
% following https://wiki.cam.ac.uk/bmuwiki/FMRI, citing Power et al 2012 https://doi.org/10.1016/j.neuroimage.2011.10.018 and Patel et al 2014 https://doi.org/10.1016/j.neuroimage.2014.03.012)
scrubbing_threshold = 0.4; % 40%

% spm; % OPEN SPM GUI

%% Preproc after FSL TopUp and McFLIRT

for i = 1:length(subjects_list)
    clear matlabbatch;
    spm('defaults','fmri');
    spm_jobman('initcfg');

    t1_file = fullfile(analysis_dir,subjects_list{i},'anat',[subjects_list{i},'_T1w.nii,1']); % NOTE THE ,1
    epi_file = fullfile(analysis_dir,subjects_list{i},'func', [subjects_ids{i},'.EPI_fMRI.rs.PA.mcf.topped.nii']);
    volumes = spm_select('Expand',epi_file);
    volumes_cell = cellstr(volumes);
    V = spm_vol(epi_file);
    pixdim = sqrt(sum(V(1).mat(1:3,1:3).^2)); % [pixdim1 pixdim2 pixdim3] GET THE ORIGINAL VOXELS DIMENSION
                                              % e.g. [3.0 3.0 3.3] NOT ISOTROPIC
    iso_dim = min(pixdim);                    % use the smallest dimension to the define an isotropic voxel dimension
    isovox  = [iso_dim iso_dim iso_dim];      % e.g. [3.0 3.0 3.0]

    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {t1_file}; % unprocessed T1
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = volumes_cell(1);       % CELL: first volume of EPI-BOLD
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = volumes_cell(2:end);    % CELL: remainder of volumes of EPI-BOLD
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
    matlabbatch{2}.spm.spatial.preproc.channel.vols = {t1_file}; % unprocessed T1
    matlabbatch{2}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{2}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{2}.spm.spatial.preproc.channel.write = [0 1];
    matlabbatch{2}.spm.spatial.preproc.tissue(1).tpm = {'/work/users/d/g/dga/tools/spm12_7487/tpm/TPM.nii,1'};
    matlabbatch{2}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{2}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{2}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{2}.spm.spatial.preproc.tissue(2).tpm = {'/work/users/d/g/dga/tools/spm12_7487/tpm/TPM.nii,2'};
    matlabbatch{2}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{2}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{2}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{2}.spm.spatial.preproc.tissue(3).tpm = {'/work/users/d/g/dga/tools/spm12_7487/tpm/TPM.nii,3'};
    matlabbatch{2}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{2}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{2}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{2}.spm.spatial.preproc.tissue(4).tpm = {'/work/users/d/g/dga/tools/spm12_7487/tpm/TPM.nii,4'};
    matlabbatch{2}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{2}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{2}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{2}.spm.spatial.preproc.tissue(5).tpm = {'/work/users/d/g/dga/tools/spm12_7487/tpm/TPM.nii,5'};
    matlabbatch{2}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{2}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{2}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{2}.spm.spatial.preproc.tissue(6).tpm = {'/work/users/d/g/dga/tools/spm12_7487/tpm/TPM.nii,6'};
    matlabbatch{2}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{2}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{2}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{2}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{2}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{2}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{2}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{2}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{2}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{2}.spm.spatial.preproc.warp.write = [0 1];
    matlabbatch{3}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
    matlabbatch{3}.spm.spatial.normalise.write.subj.resample = volumes_cell;
    matlabbatch{3}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85]; % bounding box MNI dimensions, following https://www.fil.ion.ucl.ac.uk/spm/docs/tutorials/fmri_preprocessing/scripting/
    matlabbatch{3}.spm.spatial.normalise.write.woptions.vox = isovox;
    matlabbatch{3}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{3}.spm.spatial.normalise.write.woptions.prefix = 'w';
    matlabbatch{4}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    matlabbatch{4}.spm.spatial.smooth.fwhm = [FWHM FWHM FWHM];
    matlabbatch{4}.spm.spatial.smooth.dtype = 0;
    matlabbatch{4}.spm.spatial.smooth.im = 0;
    matlabbatch{4}.spm.spatial.smooth.prefix = 's';

    spm_jobman('run', matlabbatch);
end

%% Binarized Masks
% Process:
% 1) Normalise.write takes the native tissue maps (c1…c3) and deforms them into
% the same MNI grid I used for the EPI (so they will align with the wP.*.nii).
% It creates out wc1…, wc2…, wc3… files
% 2) Imcalc reads those three, sums them, then thresholds at 0.1 to create a
% single binary "brainmask.nii" in MNI space.

clear matlabbatch;
for i = 1:length(subjects_list)
    % 1) Normalise the GM/WM/CSF maps into MNI (to match the w*EPI)
    clear matlabbatch;
    spm('defaults','fmri');
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {fullfile(analysis_dir,subjects_list{i},'anat',['y_' subjects_list{i} '_T1w.nii'])};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {
        fullfile(analysis_dir,subjects_list{i},'anat',['c1' subjects_list{i} '_T1w.nii'])
        fullfile(analysis_dir,subjects_list{i},'anat',['c2' subjects_list{i} '_T1w.nii'])
        fullfile(analysis_dir,subjects_list{i},'anat',['c3' subjects_list{i} '_T1w.nii'])};
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [3 3 3];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w'; % produces wc1*, wc2*, wc3*

    spm_jobman('run',matlabbatch);

    % 2) Imcalc to merge & threshold into a single binarized mask
    clear matlabbatch;
    spm('defaults','fmri');
    matlabbatch{1}.spm.util.imcalc.input = {
        fullfile(analysis_dir,subjects_list{i},'anat',['wc1' subjects_list{i} '_T1w.nii'])
        fullfile(analysis_dir,subjects_list{i},'anat',['wc2' subjects_list{i} '_T1w.nii'])
        fullfile(analysis_dir,subjects_list{i},'anat',['wc3' subjects_list{i} '_T1w.nii'])};
    matlabbatch{1}.spm.util.imcalc.output = 'brainmask';
    matlabbatch{1}.spm.util.imcalc.outdir = {fullfile(analysis_dir, subjects_list{i},'anat')};
    % combine GM+WM+CSF and threshold at 0.1
    % matlabbatch{1}.spm.util.imcalc.expression = 'i1 + i2 + i3 > 0.1'; whole-brain tissue mask (GM,WM,CSF combined>0.1)
    matlabbatch{1}.spm.util.imcalc.expression = 'i1 > 0.1'; % GM-only mask (that is,only voxels with GM probability>0.1)
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = -7;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;  % float

    spm_jobman('run',matlabbatch);
end

%% 1st-Level of analysis
behavioral_dir = '/work/users/d/g/dga/SeLECTS_behavioral/';
par_files_dir = '/work/users/d/g/dga/SeLECTS_par_files/';
scripts_dir = '/work/users/d/g/dga/code/spm/';
cd(scripts_dir);

for i = 1:length(subjects_list)
    clear matlabbatch;
    clear subject;
    spm('defaults','fmri');

    % dataset = [analysis_dir subjects_list{i} '/func/' 'sw' subjects_ids{i} '.EPI_fMRI.rs.PA.mcf.topped.nii']; % load only preprocessed data
    % ---NOTE: warped but unsmoothed images for future graph theory analysis...
    dataset = [analysis_dir subjects_list{i} '/func/' 'w' subjects_ids{i} '.EPI_fMRI.rs.PA.mcf.topped.nii']; % load only preprocessed data
    volumes = spm_select('Expand',dataset);

    % Extract the task and scan data (including 6 movement params) - LOADS SUBJECT STRUCTURE!
    structure = fullfile(behavioral_dir,subjects_list{i},'subjectinfo_EPI_fMRI_rs_PA_mcf_topped_nii.mat');
    load(structure);
    new_par_filepath = fullfile(par_files_dir,subjects_list{i}, [subjects_ids{i},'.EPI_fMRI.rs.PA.mcf.par']);
    subject.new_par_filepath = new_par_filepath;

    % copy, paste, and rename the .PAR FILE, so that OPTCENS can read it and calculate FD/STS
        % REMOVE THIS FILE LATER, TO AVOID FD DATA BE CROSSED WITH WL AND WM CONDITIONS
    rp_filename = [analysis_dir subjects_list{i} '/func/rp_motion.txt'];
    [status,message,messageId] = copyfile(new_par_filepath, rp_filename);

    % sanity check, compare subject struct motion data to current .par file data
    rowsSame = []; colsSame = [];
    motion_params = dlmread(new_par_filepath);
    TF = isfield(subject,'R_movement_params'); % check if the R movementparams field exists!
    if TF == 1 % if exists, check that movement parameters have same size
        try % checking number of volumes!
            rowsSame = all(subject.R_movement_params==motion_params,2);
        catch
            difference = abs(length(subject.R_movement_params) - length(motion_params));
            % determine which is greater, and make an array accordingly for later logical checking
            if length(motion_params) > length(subject.R_movement_params)
                rowsSame = cat(1, ones(length(subject.R_movement_params),1), zeros(difference,1));
            elseif length(subject.R_movement_params) > length(motion_params)
                rowsSame = cat(1, ones(length(motion_params),1), zeros(difference,1));
            end
            warning('Number of volumes between previous saved movement parameters are different from those in pasted .par file');
            warning('Taking only copied data, from .par file');
            subject.R_movement_params = motion_params;
        end
    elseif TF == 0 % if it doesn't exist, it was not created before...then create it now
        subject.R_movement_params = motion_params;
    end

    % determine 24P, following Friston et al (1996)
    cd(scripts_dir); % go to that dir to allow the scripts running 'in-place'
    [regs24P] = dg_calculate_24P(subject);
    subject.TwentyFourP = regs24P;

    % get aCompCor components, after running my script "run_acompcor.sh"
    aCompCor_file = fullfile(par_files_dir,subjects_list{i},'rs_acompcor_all.tsv');
    if ~exist(aCompCor_file)
        error('aCompCor components not computed yet for REST task of %s', subject.id);
    end
    aCompCor_all = readtable(aCompCor_file,'FileType','text','Delimiter','\t');
    aCompCor_six = aCompCor_all(:,1:6);
    subject.Six_aCompCor = table2array(aCompCor_six);

    % calculate FD to determine volumes exceeding the FD threshold
    [out_dir,~,~] = fileparts(subject.new_par_filepath);
    [~,filename_func,~] = fileparts(dataset);
    [mean_FD, FD] = dg_calculate_FD(motion_params, filename_func, out_dir, r, true);
    std_FD = std(FD);
    num_vols_to_scrub = length(FD(FD > FD_threshold));

    % print mean FD and SD of FD, informative
    if mean_FD >= FD_threshold
        warning('%s Mean FD = %.2f, greater than/equal to FD threshold (%.1f)', subjects_ids{i}, mean_FD, FD_threshold);
        fprintf('%s SD of FD = %.2f \n', subjects_ids{i}, std_FD);
    else 
        fprintf('%s Mean FD = %.2f \n', subjects_ids{i}, mean_FD);
        fprintf('%s SD of FD = %.2f \n', subjects_ids{i}, std_FD);
    end

    % determine a scrubbing-spike regression matrix -> GLM!
    idx = FD > FD_threshold;
    spikes_to_regress = find(idx); % identifies "bad" volumes
    spike_regression_matrix = eye(length(volumes));
    spike_regression_matrix = spike_regression_matrix(:, spikes_to_regress);
    subject.spike_regression_matrix = spike_regression_matrix;

    % CONFOUNDS_REGRESSORS = horzcat(subject.TwentyFourP,...
    %                                subject.Six_aCompCor,...
    %                                subject.spike_regression_matrix); % ALL REGRESSORS FOR SPM GLM!

    % DO NOT INCLUDE SPIKE REGRESSION MATRIX NOW, ALLOW OPTCENS TO CALCULATE THIS FOR INTERPOLATION...
    CONFOUNDS_REGRESSORS = horzcat(subject.TwentyFourP,...
                                   subject.Six_aCompCor); % ALL REGRESSORS FOR SPM GLM!    

    % check whether scrubbing vector > scrubbing threshold, and print informative
    scrubbing_proportion = length(spikes_to_regress) / length(volumes);
    scrubbing_percent = scrubbing_proportion*100;
    if scrubbing_proportion >= scrubbing_threshold
        % warning('%.2f%% of volumes will be scrubbed, greater than/equal to scrubbing threshold (%.0f%%)', scrubbing_percent, scrubbing_threshold_percent);
    else
        % fprintf('%.2f%% of volumes will be scrubbed, lower to scrubbing threshold (%.0f%%) \n', scrubbing_percent, scrubbing_threshold_percent);
    end

    % save new brainmask filepath in subject struct
    if ~exist(fullfile(analysis_dir,subjects_list{1},'anat','brainmask.nii'))
        error('Brainmask creation for %s was unsuccessful or not performed. Please, do check.', subjects_ids{i});
    elseif exist(fullfile(analysis_dir,subjects_list{1},'anat','brainmask.nii'))
        subject.brainmask = fullfile(analysis_dir,subjects_list{1},'anat','brainmask.nii');
    end

    % update and sabe subject struct, including mean FD, SD of FD, and % of vols scrubbed
    subject.raw_dataset = [analysis_dir subjects_list{i} '/func/' subjects_ids{i} '.EPI_fMRI.rs.PA.mcf.topped.nii'];
    subject.analysis_dir = analysis_dir;
    subject.mri_T1 = fullfile(analysis_dir,subjects_list{i},'anat',[subjects_list{i},'_T1w.nii,1']);
    subject.QCdir = fullfile(par_files_dir,subjects_list{i});
    subject.fmri_par = new_par_filepath;
    subject.FD_mean = mean_FD;
    subject.FD_SD = std_FD;
    subject.scrubbing_percent = scrubbing_percent;
    toRemove = {'montage','valid_channels','Scans_rejected','fname_art_reg',...
                'ART_threshold','scrubbing_ARTarray'}; % list of fields to remove
    present = toRemove(isfield(subject,toRemove)); % find which of those are actually in the struct
    subject   = rmfield(subject,present); % remove only the present ones
    save(structure,'subject'); % finish editing and save the subject struct

    % create a dir to write 1stLevel results
    if Trim == 1
        firstLevel_dir = [analysis_dir subjects_list{i} '/func/' '1stLevel_rs_02FD'];
        if ~exist(firstLevel_dir)
            mkdir(firstLevel_dir);
        end
    elseif Trim == 0
        firstLevel_dir = [analysis_dir subjects_list{i} '/func/' '1stLevel_rs_untrimmed_optcens'];
        if ~exist(firstLevel_dir)
            mkdir(firstLevel_dir);
        end
    end
    cd(firstLevel_dir);

    % START SPM MODEL SPECIFICATION
    matlabbatch{1}.spm.stats.fmri_spec.dir = {firstLevel_dir};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

    % Trim dataset: Total scan length - task length (after data saved in the Subject struct!!!!)
    % Also trim the regressors matrix accordingly!!!!
    if Trim == 1
        if ~isempty(subject.scan_trimmed_length)
            warning('Total length of dataset is greater than the last scanned task timepoint! \n');
            warning('This functional dataset was trimmed before! \n');
            fprintf(['Setting a trimmed length of the scan dataset for ', subjects_list{i},'! \n']);
            trimmed_dataset = volumes(subject.scan_trimmed_length,:); % it must index the whole width of chars with colon (:)
            CONFOUNDS_REGRESSORS = CONFOUNDS_REGRESSORS(subject.scan_trimmed_length,:); % trims regressors matrix
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(trimmed_dataset);
        elseif isempty(subject.scan_trimmed_length)
            fprintf(['RS scan was not trimmed for ', subjects_ids{i},'! \n']);
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(volumes);
        end
    elseif Trim == 0
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(volumes);
    end

    % save regressors matrix as a tsv file, because SPM needs it as fname!
    % regressors_matrix_file = fullfile(out_dir, 'rs_all_regressors.tsv');
    regressors_matrix_file = fullfile(firstLevel_dir, 'rs_all_regressors.mat');
    R = CONFOUNDS_REGRESSORS;
    % dlmwrite(regressors_matrix_file,CONFOUNDS_REGRESSORS,'delimiter','\t','precision',  '%.6f');
    save(regressors_matrix_file, "R");

    % ONSET AND DURATION OF REST
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name = 'Rest';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = length(volumes) * TR;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).orth = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};

    % BRINGING REGRESSORS
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    % IMPORTANT NOTE: CONFOUNDS_REGRESSORS HAS TO BE CALLED 'R' SO THAT SPM READS ALL REGRESSORS HERE
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {[firstLevel_dir '/rs_all_regressors.mat']};
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;

    % BRINGING BRAINMASK
    matlabbatch{1}.spm.stats.fmri_spec.mask = {[analysis_dir '/' subjects_list{i} '/anat/' 'brainmask' '.nii']};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

    % MODEL REVIEW
    matlabbatch{2}.spm.stats.review.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.review.display.matrix = 1;
    matlabbatch{2}.spm.stats.review.print = 'jpg';
    
    % MODEL ESTIMATION
    matlabbatch{3}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.fmri_est.write_residuals = 1;
    matlabbatch{3}.spm.stats.fmri_est.method.Classical = 1;

    % OPTIMIZED CENSORING TOOLBOX (M. WILKE) it does motion censoring and interpolation, and brings the corresponding regressors to the model
    %  -TAKE INTERPOLATION RESULTS ONLY-
    matlabbatch{4}.spm.tools.optcens_cfg.ospms(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{4}.spm.tools.optcens_cfg.approach = 2; % approach 2 is for resting-state
    matlabbatch{4}.spm.tools.optcens_cfg.docalc = 3;
    matlabbatch{4}.spm.tools.optcens_cfg.o_minrem = 0;
    % adapt the max % of volumes to censor/interpolate (case-by-case decision: % vols with FD > FD threshold)
    matlabbatch{4}.spm.tools.optcens_cfg.o_maxrem = scrubbing_percent;
    matlabbatch{4}.spm.tools.optcens_cfg.zeropad = 0;
    matlabbatch{4}.spm.tools.optcens_cfg.adapt = 0;

    % CONTRAST MANAGER
    matlabbatch{5}.spm.stats.con.spmmat(1) = cfg_dep('Optimized Censoring Toolbox: Optimally censored & interpolated SPM.mat(s)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','int'));
    matlabbatch{5}.spm.stats.con.consess{1}.tcon.name = 'Rest > Baseline';
    matlabbatch{5}.spm.stats.con.consess{1}.tcon.weights = 1; % only one condition
    matlabbatch{5}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{5}.spm.stats.con.delete = 0;  % keep previous contrasts if they exist

    % RESULTS REPORT
    % FWE
    matlabbatch{6}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{6}.spm.stats.results.conspec(1).titlestr = 'Rest > Baseline';
    matlabbatch{6}.spm.stats.results.conspec(1).contrasts = 1;
    matlabbatch{6}.spm.stats.results.conspec(1).threshdesc = 'FWE';
    matlabbatch{6}.spm.stats.results.conspec(1).thresh = 0.05;
    matlabbatch{6}.spm.stats.results.conspec(1).extent = 0;
    matlabbatch{6}.spm.stats.results.conspec(1).conjunction = 1;
    matlabbatch{6}.spm.stats.results.conspec(1).mask.none = 1;
    % Unc
    matlabbatch{6}.spm.stats.results.conspec(2).titlestr = 'Rest (p<0.001 uncorrected)';
    matlabbatch{6}.spm.stats.results.conspec(2).contrasts = 1;
    matlabbatch{6}.spm.stats.results.conspec(2).threshdesc = 'none';
    matlabbatch{6}.spm.stats.results.conspec(2).thresh = 0.001;
    matlabbatch{6}.spm.stats.results.conspec(2).extent = 0;
    matlabbatch{6}.spm.stats.results.conspec(2).conjunction = 1;
    matlabbatch{6}.spm.stats.results.conspec(2).mask.none = 1;
    matlabbatch{6}.spm.stats.results.units = 1; % 1 = slices
    matlabbatch{6}.spm.stats.results.export{1}.jpg = true;
    matlabbatch{6}.spm.stats.results.export{2}.ps = true;
    matlabbatch{6}.spm.stats.results.export{3}.csv = true;

    spm_jobman('run', matlabbatch);
end

% lastly, remove all "rp_motion.txt" files copied! SEE WHY ABOVE!
for i = 1:length(subjects_list)
    rp_filename = [analysis_dir subjects_list{i} '/func/rp_motion.txt'];
    if exist(rp_filename)
        delete(rp_filename);
        fprintf('rp_motion.txt file removed for RS and %s \n', subjects_list{i});
    elseif ~exist(rp_filename)
        warning('No rp_motion copied before for RS and %s. Do check why...', subjects_list{i});
    end
end