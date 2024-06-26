% 2024-06-26 AndyP
% RT-RT correlation matrix

uS = unique(X.id);
nS = length(uS);
Ratt = [];
Rnon = [];
for iS=1:nS
    
    k = zeros(size(X,1),1);
    k(X.id==uS(iS) & (X.group_bin=='HL' | X.group_bin=='LL')) = 1; % current subject index to keep
    RT = X.rt_csv(k==1);
    B = X.block(k==1);
    uB = unique(B);
    nB = length(uB);
    for iB=1:nB
        RT0 = RT(B==uB(iB));
        if (length(RT0)==40)
            %T0 = squareform(pdist(RT0,'euclidean'));
            Ratt = cat(2,Ratt,RT0);
        end
    end
    
    k(X.id==uS(iS) & X.group_bin=='NON') = 1; % current subject index to keep
    RT = X.rt_csv(k==1);
    B = X.block(k==1);
    run_trial0 = X.run_trial(k==1);
    uB = unique(B);
    nB = length(uB);
    for iB=1:nB
        RT0 = RT(B==uB(iB));
        if (length(RT0)==40)
            %T0 = squareform(pdist(RT0,'euclidean'));
            Rnon = cat(2,Rnon,RT0);
        end
    end
    
    fprintf('%d/%d\n',iS,nS)
    
end