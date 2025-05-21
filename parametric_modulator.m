% 2025-05-15 AndyP
% Convolve a dmUBLOCK(1) output from AFNI with a parametric modulator for
% entropy
% save feedback, clock onset, and entropy PM regressors for import into
% CONN toolbox SPM wrapper.
homedir = '/Users/dnplserv/gPPI/Explore_HC_only';
cd(homedir)
D = dir('sub-*');
TR = 0.6;

for iD = 1:length(D)
    cd(strcat(D(iD).name,'/func'));
    
    disp(D(iD).name);
    
    clock_files = dir('*clock_onset_regressor_dmU_*.1D');
    feedback_files = dir('*feedback_regressor_dmU_*.1D');
    df_file = dir('*-entropy-PM.csv');
    
    clock_file0 = nan;
    feedback_file0 = nan;
    df_file0 = nan;
    for iR = 1:2
        
        if  regexp(clock_files(iR).name,strcat('_',mat2str(iR))) > 0
            clock_file0 = clock_files(iR).name;
        end
        if  regexp(feedback_files(iR).name,strcat('_',mat2str(iR))) > 0
            feedback_file0 = feedback_files(iR).name;
        end
        if  regexp(df_file(iR).name,strcat('run-',mat2str(iR))) > 0
            df_file0 = df_file(iR).name;
        end
        
        if ~all(isnan(clock_file0)) & ~all(isnan(feedback_file0)) & ~all(isnan(df_file0))
        
        clock_onset_rdmU = table2array(readtable(clock_file0,'FileType','text'));
        feedback_onset_rdmU = table2array(readtable(feedback_file0,'FileType','text'));
        df = table2array(readtable(df_file0));
        
        nTR = length(clock_onset_rdmU(:,1));
        onsets = interp1(1:nTR,1:nTR,df(:,3)./TR,'nearest'); % convert timing of clock onset into nearest TR, which is nearest index (time/TR)
        parameter = nanzscore(df(:,4));
        clock_onset_rdmUc = clock_onset_rdmU(:,1) - nanmean(clock_onset_rdmU(:,1));
        PM = zeros(nTR,1);
        PM(onsets) = parameter / max(parameter);
        cPM = conv(PM,clock_onset_rdmUc(:,1));
        cPM = cPM(nTR:end);
        
        save(strcat(D(iD).name,'-run-',mat2str(iR),'feedback-regressor.mat'),'feedback_onset_rdmU');
        save(strcat(D(iD).name,'-run-',mat2str(iR),'clock-onset-regressor.mat'),'clock_onset_rdmU');
        save(strcat(D(iD).name,'-run-',mat2str(iR),'clock-v_entropy_wi-regressor.mat'),'cPM');
        
        else
            warning(strcat('subject ', D(iD).name, 'run- ',iR, 'missing'));
        end
        
    end
    
    cd(homedir);
end
