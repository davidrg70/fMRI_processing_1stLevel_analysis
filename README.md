# fMRI_processing_1stLevel_analysis
Scripts for fMRI processing and 1st-Level analysis (SPM)

Script_proc_1stL_wl_job / NEW.m
Preprocesses fMRI data for the Ebner et al (2011) task (block-design), which is a phonological task, self-paced. Later it runs a 1st-Level analysis (2-sample T-test to compare the two conditions: Phonological >/< Fractals). Anat and Func normalization to MNI space.

Script_proc_1stL_wm_job / NEW.m
Preprocesses fMRI data for the Siffredi et al (2011) task (block-design), which is a working memory task, fixed-length blocks. Later it runs a 1st-Level analysis (2-sample T-test to compare the two conditions: Encoding >/< Retrieval). Anat and Func normalization to MNI space.

Script_proc_1stL_rs_job / NEW.m
Does the same but for a resting-state run, no task events, 1st-level of analysis compares Rest to baseline (constant term).

The 3 scripts perform interpolation based on DVARS and FD threshold = 0.2 mm, no scrubbing results taken, using OPTCENS (Wilke & Baldeweg 2019, https://doi.org/10.1016/j.jneumeth.2019.02.008).
