% Classification: logistic regression for EV and heat pump, threshold rules for the rest -> class.mat

wd = 'E:\Avangrid Hackathon 2026\matlab';
load(fullfile(wd,'features.mat'));
nC = height(F); rng(1);

feat = setdiff(F.Properties.VariableNames,{'id','isEV','isHP','isPV','hasData','missing','peak_month'},'stable');
hd = find(F.hasData); nd = numel(hd);
Xf = F{:,feat};
sig = @(x,x0,s) 1./(1+exp(-(x-x0)/s));                   % 0.5 at threshold x0
cLab = 0.25;                                             % labelled share of EV / HP customers

C = table(F.id,'VariableNames',{'id'}); tic
Xs = normalize(sign(Xf(hd,:)).*log1p(abs(Xf(hd,:)))); Xs(~isfinite(Xs)) = 0;
lam = 0.004;
for tech = ["EV","HP"]
    y = F.("is"+tech)(hd);
    cvm = fitclinear(Xs,y,'Learner','logistic','Regularization','lasso','Lambda',lam,'Prior','uniform','KFold',5);
    [~,sc] = kfoldPredict(cvm); g = sc(:,2);                % out-of-fold scores
    mdl = fitclinear(Xs,y,'Learner','logistic','Regularization','lasso','Lambda',lam,'Prior','uniform');
    [~,o] = sort(g,'descend'); ys = double(y(o)); k = (1:nd)';
    Lk = cumsum(ys); prec = min(1,Lk./(cLab*k)); rec = Lk/sum(ys);
    N = round(sum(y)/cLab);                                 % top N = labels / 0.25
    flag = false(nC,1); flag(hd(o(1:N))) = true;
    p = zeros(nC,1); p(hd(o)) = 1 - (k-1)/nd;               % score percentile
    [~,~,~,auc] = perfcurve(y,g,true);
    C.("p"+tech) = p; C.(tech+"_flag") = flag; PR.(tech) = [rec prec];
    fprintf('%s (%.0f s): CV AUC %.3f | top N = %d (labels/0.25): est. precision %.2f, recall %.2f\n', ...
        tech,toc,auc,N,prec(N),rec(N));
    [~,kb] = sort(abs(mdl.Beta),'descend'); kb = kb(mdl.Beta(kb)~=0);
    M.(tech) = mdl; M.(tech+"_coef") = table(feat(kb)',mdl.Beta(kb),'VariableNames',{'feature','coef'});
    fprintf('   %d features kept | top: %s\n',numel(kb),strjoin(compose('%s(%+.2f)',string(feat(kb(1:min(8,end))))',mdl.Beta(kb(1:min(8,end)))).',', '));
end
fprintf('Label check: %.0f%% of labelled EVs show <20 h/yr of high-power charging blocks (others: %.0f%%)\n', ...
    100*mean(F.n_runh(F.isEV)<20),100*mean(F.n_runh(F.hasData & ~F.isEV)<20));

C.EV       = C.EV_flag | F.isEV;
C.EV_L2    = C.EV & F.ev_kw>=3;
C.EV_L1    = C.EV & ~C.EV_L2 & F.l1_nights>=0.1;
C.pSolar   = sig(F.solar_idx,0.2,0.05);
C.Solar    = C.pSolar>=0.5 | F.isPV;
C.SolarUnreported = C.pSolar>=0.5 & ~F.isPV;
C.pSolarBatt = C.Solar.*sig(F.zero_eve,0.05,0.02);
C.SolarBatt  = C.pSolarBatt>=0.5;
C.pAC      = sig(F.cool_slope,0.5,0.15).*sig(F.corrT_sum,0.3,0.08);
C.AC       = C.pAC>=0.5;
C.pEHeat   = sig(F.heat_share,0.25,0.05).*sig(F.heat_slope,0.6,0.15);
C.EHeat    = C.pEHeat>=0.5 | C.HP_flag | F.isHP;
fprintf('HP labels inside electric-heat rule: %.0f%% (all customers %.0f%%)\n',100*mean(C.pEHeat(F.isHP)>=0.5),100*mean(C.pEHeat(hd)>=0.5));
C.HP       = C.HP_flag | F.isHP;
C.ResistHeat = C.EHeat & ~C.HP;
C.pWaterHeat = sig(F.wh_kw,0.45,0.12).*sig(F.b0,12,4);
C.WaterHeat  = C.pWaterHeat>=0.5;
C.pPool    = sig(F.pool_kw,1.0,0.2);
C.Pool     = C.pPool>=0.5;
C.pHotTub  = sig(F.ht_sum,0.3,0.07).*sig(F.ht_win,0.3,0.07);
C.HotTub   = C.pHotTub>=0.5;
C.Seasonal = F.low_months>=3;
C.Vacant   = F.ann_kwh<1000;
C.AlwaysOn   = F.base5>=1.0 & ~C.Vacant;
C.ShortCycle = F.spike1>=150 & ~C.Vacant;
C.SuppHeat   = ~C.EHeat & F.heat_slope>=0.2 & F.heat_share>=0.08 & F.heat_share<0.25 & ~C.Vacant;
C.TimerLoad  = F.n_runh>=50 & F.start_conc>=0.6 & ~C.Vacant;
lg = C.Properties.VariableNames(varfun(@islogical,C,'OutputFormat','uniform'));
for i = 1:numel(lg), C.(lg{i}) = C.(lg{i}) & F.hasData; end

techs = {'EV','EV_L2','EV_L1','AC','EHeat','HP','ResistHeat','WaterHeat','Pool','HotTub','Solar','SolarUnreported','SolarBatt','SuppHeat','ShortCycle','AlwaysOn','TimerLoad','Seasonal','Vacant'};
S = table(techs',zeros(numel(techs),1),zeros(numel(techs),1),zeros(numel(techs),1),zeros(numel(techs),1), ...
    'VariableNames',{'technology','customers','pct','median_kwh','median_peak_kw'});
for i = 1:numel(techs)
    k = C.(techs{i});
    S.customers(i) = sum(k); S.pct(i) = 100*sum(k)/nd;
    S.median_kwh(i) = median(F.ann_kwh(k)); S.median_peak_kw(i) = median(F.max_kw(k));
end
disp(S)
fprintf('Hot tub homes that also have an EV: %d\n',sum(C.HotTub & C.EV));
main = [C.AC C.EHeat C.WaterHeat C.Pool C.HotTub C.EV C.Solar]; nt = sum(main(hd,:),2);
fprintf('Main technologies per home (AC, electric heat, water heater, pool, hot tub, EV, solar): none %d | one %d | two %d | three+ %d\n', ...
    sum(nt==0),sum(nt==1),sum(nt==2),sum(nt>=3));
fprintf('Solar rule vs complete PV list: recall %.2f, precision %.2f\n', ...
    mean(C.pSolar(F.isPV)>=0.5), sum(C.pSolar>=0.5 & F.isPV)/sum(C.pSolar>=0.5 & F.hasData));

f = figure('Position',[100 100 900 400]);
for i = 1:2
    tk = ["EV","HP"]; subplot(1,2,i); plot(PR.(tk(i))(:,1),PR.(tk(i))(:,2),'LineWidth',2); grid on
    xlabel('Recall (of labelled)'); ylabel('Estimated precision'); ylim([0 1]); xline(PR.(tk(i))(round(sum(F.("is"+tk(i)))/cLab),1),'--','top N'); title(tk(i)+": PU-estimated precision-recall");
end
fd = fullfile(wd,'figs'); if ~isfolder(fd), mkdir(fd); end
exportgraphics(f,fullfile(fd,'fig_pu_pr.png'),'Resolution',150);
save(fullfile(wd,'class.mat'),'C','S','M','PR','feat');
