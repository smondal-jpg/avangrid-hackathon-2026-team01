% Figures: appliance signatures, technology counts, segment share of system peak

wd = 'E:\Avangrid Hackathon 2026\matlab';
load(fullfile(wd,'ami.mat'),'kwh');
load(fullfile(wd,'features.mat')); load(fullfile(wd,'class.mat')); load(fullfile(wd,'seg.mat'));
fd = fullfile(wd,'figs');

pick = @(mask,score) find(mask & score==max(score(mask)),1);
ids = [pick(C.EV & F.isEV & F.n_runh>100, F.plat_kw.*(F.plat_kw<9)), ...
       pick(C.ResistHeat | C.HP, C.pEHeat.*F.r2), ...
       pick(C.Pool & ~C.EHeat & F.winsum<0.9, F.pool_kw.*(F.pool_kw<4).*F.r2), ...
       pick(F.isPV, F.solar_idx)];
ttl = ["EV charging (L2 blocks, evenings)","Electric heating (winter, all hours)", ...
       "Pool pump (summer daytime block)","Solar PV (midday import near zero)"];
f = figure('Position',[50 50 1400 700]);
for i = 1:4
    subplot(2,2,i);
    x = reshape(double(kwh(ids(i),:)),24,[]);
    imagesc(1:size(x,2),0:23,x); axis xy; colormap(gca,turbo); cb = colorbar; cb.Label.String = 'kWh';
    caxis([0 prctile(x(:),99.5)]);
    set(gca,'XTick',[1 92 182 274 365],'XTickLabel',{'Sep','Dec','Mar','Jun','Aug'});
    ylabel('Hour of day'); title(sprintf('%s - customer %d',ttl(i),ids(i)));
end
exportgraphics(f,fullfile(fd,'fig_signatures.png'),'Resolution',150);

S2 = S(~ismember(S.technology,{'EV_L1','Vacant'}),:);
[~,o] = sort(S2.customers,'descend'); S2 = S2(o,:);
f = figure('Position',[50 50 1000 420]);
n = height(S2); v = flip(S2.customers); pc = flip(S2.pct);
barh(v); grid on; set(gca,'YTick',1:n,'YTickLabel',flip(S2.technology),'TickLabelInterpreter','none');
xlabel('Customers'); title(sprintf('Detected technologies (%d customers with data)',sum(F.hasData)));
text(v,1:n,compose('  %d (%.1f%%)',v,pc),'FontSize',8);
exportgraphics(f,fullfile(fd,'fig_tech_counts.png'),'Resolution',150);

ok = F.hasData & G.cluster>0;
[g,nmK] = findgroups(G.cluster_name(ok));
shC = splitapply(@numel,F.cp_kw(ok),g)/sum(ok)*100;
shP = splitapply(@sum,F.cp_kw(ok),g)/sum(F.cp_kw(ok))*100;
f = figure('Position',[50 50 1100 420]);
bar(categorical(nmK),[shC shP]); ylabel('%'); grid on
legend('Share of customers','Share of system peak load','Location','northwest');
title('Load-shape segments: who drives the system peak?');
exportgraphics(f,fullfile(fd,'fig_segment_peak_share.png'),'Resolution',150);
disp(table(nmK,shC,shP,'VariableNames',{'segment','pct_customers','pct_system_peak'}))
