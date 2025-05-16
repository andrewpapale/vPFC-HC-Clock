% 2025-05-15 AndyP
% Convolve a dmUBLOCK(1) output from AFNI with a parametric modulator for
% entropy

% first import test1.csv (entropy for 202200, run 1)
TR = 0.6;
nTR = length(rdmU1(:,3)); % 3 is dmUBLOCK(1), for now
onsets = interp1(1:nTR,1:nTR,test(:,4)./TR,'nearest');
parameter = test1;
PM = zeros(nTR,1); 
PM(onsets) = parameter / max(parameter);
cPM = conv(PM,rdmU1(:,3));
cPM = cPM(nTR:end);

