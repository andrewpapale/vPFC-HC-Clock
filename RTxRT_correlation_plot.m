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
    
    k(X.id==uS(iS) & (X.group_bin=='NON')) = 1; % current subject index to keep
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


% rng(0,'twister');
% Tdiff = [];
% Pdiff = [];
% for i=1:1000
%     r = randi([1,size(Rnon,2)],1,size(Ratt,2));
%     [Tdiff0,Pdiff0] = corrcoef(Ratt' - Rnon(:,r)');
%     Tdiff = cat(3,Tdiff,Tdiff0);
%     Pdiff = cat(3,Pdiff,Pdiff0);
%     fprintf('%d/%d\n',i,1000)
% end
[Tatt,Patt] = corrcoef(Ratt');
[Tnon,Pnon] = corrcoef(Rnon');
x = 1:1:40;
y = 1:1:40;
pcolor(smooth2a(Tatt - Tnon,1,1)); shading flat; axis square;
hold on
plot(x,y,'color','r','linewidth',2);
grid on
shading flat
colorbar
set(gca,'fontsize',24)
caxis([-0.1 0.1])
title('ATT - NON Pearson Correlation Coefficients','fontsize',24)

% Tdiff1 = Tdiff;
% Tdiff1(Pdiff > 1/(40*40)) = NaN;
% pcolor((smooth2a(squeeze(nanmean(Tdiff1,3)),2,2))); shading flat;
