% 2025-05-20 AndyP
% batch script for a conn toolbox gPPI for Explore (testing..)

NSUBJECTS=1;
cwd=pwd;
FUNCTIONAL_FILE=cellstr(conn_dir('nfas-sub-202200_task-clockRev_run-*_space-MNI152NLin2009cAsym_desc-preproc_bold.nii'));
STRUCTURAL_FILE=cellstr(conn_dir('sub-202200_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii'));
if rem(length(FUNCTIONAL_FILE),NSUBJECTS),error('mismatch number of functional files %n', length(FUNCTIONAL_FILE));end
if rem(length(STRUCTURAL_FILE),NSUBJECTS),error('mismatch number of anatomical files %n', length(FUNCTIONAL_FILE));end
nsessions=2;
FUNCTIONAL_FILE=reshape(FUNCTIONAL_FILE,[NSUBJECTS,nsessions]);
STRUCTURAL_FILE={STRUCTURAL_FILE{1:NSUBJECTS}};
disp([num2str(size(FUNCTIONAL_FILE,1)),' subjects']);
disp([num2str(size(FUNCTIONAL_FILE,2)),' sessions']);
TR=0.6; % Repetition time

feedback_regressor = cellstr(conn_dir('sub-*-run-*feedback-regressor.mat'));
clock_regressor = cellstr(conn_dir('sub-*-run-*clock-onset-regressor.mat'));
entropy_PM_regressor = cellstr(conn_dir('sub-*-run-*clock-v_entropy_wi-regressor.mat'));

%% CONN-SPECIFIC SECTION: RUNS PREPROCESSING/SETUP/DENOISING/ANALYSIS STEPS
%% Prepares batch structure
clear batch;
batch.filename=fullfile(cwd,'Explore_gPPI.mat');            % New conn_*.mat experiment name

%% SETUP & PREPROCESSING step (using default values for most parameters, see help conn_batch to define non-default values)
% CONN Setup                                            % Default options (uses all ROIs in conn/rois/ directory); see conn_batch for additional options
% CONN Setup.preprocessing                               (realignment/coregistration/segmentation/normalization/smoothing)
batch.Setup.isnew=1;
batch.Setup.nsubjects=NSUBJECTS;
batch.Setup.RT=TR;                                        % TR (seconds)
batch.Setup.functionals=repmat({{}},[NSUBJECTS,1]);       % Point to functional volumes for each subject/session
for nsub=1:NSUBJECTS,for nses=1:nsessions,batch.Setup.functionals{nsub}{nses}{1}=FUNCTIONAL_FILE{nsub,nses}; end; end %note: each subject's data is defined by 2 sessions and one single (4d) file per session
batch.Setup.structurals=STRUCTURAL_FILE;                  % Point to anatomical volumes for each subject
nconditions=3;                                

batch.Setup.conditions.names=[{'clock_onset'}, {'feedback_onset'},{'clock_onset_v_entropy_wi_PM'}];
for ncond=1:nconditions
    for nsub=1:NSUBJECTS
        for nses=1:nsessions
            load(feedback_regressor{nses})
            load(clock_regressor{nses})
            load(entropy_PM_regressor{nses})
            task_regressors0{1+ncond}{nsub}{nses}=clock_onset_rdmU;
            task_regressors0{1+ncond}{nsub}{nses}=feedback_onset_rdmU; 
            task_regressors0{1+ncond}{nsub}{nses}= cPM;
        end
    end
end


batch.Setup.preprocessing.steps='default_mni';
batch.Setup.preprocessing.sliceorder='interleaved (Siemens)';
batch.Setup.done=1;
batch.Setup.overwrite='Yes';

batch.Denoising.done = 1; % skip denoising

batch.Analysis = 2;

batch.Analyses.done = 0; %
batch.Analyses.overwrite = 1;
batch.Analyses.name = 'gPPI';
batch.Analyses.measure = 3;
batch.Analyses.weight = 2;
batch.Analyses.modulation = 2;
batch.Analyses.conditions = [];
batch.Analyses.type = 3;
batch.Analyses.modulation = 2;
batch.Analyses.conditions = [];
batch.Analyses.sources = [{'entropy'}    {'atlas.Hippocampus r'}    {'atlas.Hippocampus l'}];

conn_batch(batch);
