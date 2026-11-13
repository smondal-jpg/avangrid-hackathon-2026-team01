% Segmentation: k-means load-shape clusters and rule-based business segments -> seg.mat

wd = 'E:\Avangrid Hackathon 2026\matlab';
load(fullfile(wd,'features.mat')); load(fullfile(wd,'class.mat'));
fd = fullfile(wd,'figs'); if ~isfolder(fd), mkdir(fd); end
nC = height(F); rng(1);

P = F{:,compose('p%02d',0:23)}; Mo = F{:,compose('m%02d',1:12)};
Z = normalize([P, Mo]); Z(~isfinite(Z)) = 0;
ok = ~C.Vacant & F.hasData;

sub = find(ok); sub = sub(randperm(numel(sub),min(6000,numel(sub))));
ev = evalclusters(Z(sub,:),'kmeans','silhouette','KList',4:10);
fprintf('silhouette by k: %s (best k = %d); using k = 6 for business-usable segments\n',mat2str(round(ev.CriterionValues,3)),ev.OptimalK);
k = 6;
[idx,~] = kmeans(Z(ok,:),k,'Replicates',5,'MaxIter',500);
cl = zeros(nC,1); cl(ok) = idx;
Pc = splitapply(@(x) mean(x,1),P(ok,:),idx); Mc = splitapply(@(x) mean(x,1),Mo(ok,:),idx);
cen = [Pc Mc];
met = [mean(Mc(:,[12 1 2]),2)./mean(Mc(:,[6 7 8]),2), Mc(:,7), mean(Pc(:,17:22),2), ...
       mean(Pc(:,11:16),2), mean(Pc(:,1:5),2), mean(Pc(:,7:10),2)];
lbl = ["Winter-peaking (electric heat)","Summer-peaking (cooling)","Evening-peaking (commuter)", ...
       "Daytime-active (home all day)","Night-heavy / daytime-low","Morning-peaking"];
Dv = normalize(met); nm = repmat("Typical / balanced",k,1); used = false(1,6); done = false(k,1);
for it = 1:k
    Dm = Dv; Dm(done,:) = -Inf; Dm(:,used) = -Inf;
    [mx,ii] = max(Dm(:)); if mx < 0.5, break, end
    [j,m] = ind2sub(size(Dm),ii); nm(j) = lbl(m); used(m) = true; done(j) = true;
end
sz = splitapply(@median,F.ann_kwh(ok),idx);
nameK = compose("C%d %s (median %.1f MWh/yr)",(1:k)',nm,sz/1000);
G = table(F.id,cl,'VariableNames',{'id','cluster'});
G.cluster_name = repmat("Vacant / very low use",nC,1); G.cluster_name(ok) = nameK(idx); G.cluster_name(~F.hasData) = "No data";

cpTop = F.cp_kw >= prctile(F.cp_kw(ok),90) & ok;
G.HighElectrification = (C.EV & C.EHeat) | (C.EV & F.ann_kwh>15000) | (C.HP & C.WaterHeat);
G.PeakContributor     = cpTop;
G.Seasonal            = C.Seasonal & ok;
G.FlexibleLoad        = C.EV | C.Pool | C.WaterHeat | C.HotTub;
G.EmergingCandidate   = ~C.EV & ~C.EHeat & ~C.Vacant & (C.Solar | (C.AC & F.ann_kwh>prctile(F.ann_kwh(ok),75)));
G.LowUse              = C.Vacant;

segs = {'HighElectrification','PeakContributor','Seasonal','FlexibleLoad','EmergingCandidate','LowUse'};
T2 = table(segs',zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,1),'VariableNames', ...
    {'segment','customers','pct','median_kwh','median_cp_kw','share_of_system_peak_pct'});
for i = 1:6
    s = G.(segs{i});
    T2{i,2:6} = [sum(s), 100*sum(s)/sum(F.hasData), median(F.ann_kwh(s)), median(F.cp_kw(s)), 100*sum(F.cp_kw(s))/sum(F.cp_kw(F.hasData))];
end
disp(T2)
K2 = groupsummary(table(G.cluster_name(F.hasData),F.ann_kwh(F.hasData),F.cp_kw(F.hasData),double(C.EV(F.hasData)),double(C.EHeat(F.hasData)),double(C.Solar(F.hasData)),'VariableNames', ...
    {'cluster','kwh','cp_kw','EV','EHeat','Solar'}),'cluster',{'median','mean'},{'kwh','cp_kw','EV','EHeat','Solar'});
disp(K2(:,{'cluster','GroupCount','median_kwh','median_cp_kw','mean_EV','mean_EHeat','mean_Solar'}))

f = figure('Position',[100 100 1300 450]);
subplot(1,2,1); plot(0:23,cen(:,1:24)','LineWidth',2); xlabel('Hour of day'); ylabel('Load / daily mean');
title('Cluster daily shapes'); grid on; xlim([0 23]); legend(nameK,'Location','northwest','FontSize',7);
subplot(1,2,2); plot(1:12,cen(:,25:36)','LineWidth',2); xlabel('Month'); ylabel('Load / annual mean');
title('Cluster seasonal shapes'); grid on; xlim([1 12]);
exportgraphics(f,fullfile(fd,'fig_clusters.png'),'Resolution',150);

save(fullfile(wd,'seg.mat'),'G','T2','K2','cen','nameK');
