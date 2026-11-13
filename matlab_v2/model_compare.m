% Model types, feature families and weather tests for the EV and heat-pump models

wd = 'E:\Avangrid Hackathon 2026\matlab';
tic
load(fullfile(wd,'features.mat')); load(fullfile(wd,'class.mat'),'C');
has = F.hasData==1; hd = find(has);
feat = setdiff(F.Properties.VariableNames,{'id','isEV','isHP','isPV','hasData','missing','peak_month'},'stable');
X = F{hd,feat}; X(~isfinite(X)) = 0;
fam = {'size of the home',      {'ann_kwh','mean_kw','max_kw','p99','p95','base5','load_factor','cp_kw','cp_ratio'}; ...
       'weather response',      {'b0','heat_slope','cool_slope','r2','heat_share','cool_share','cold_night','mild_night','cold_ratio','vcold_kw','c30_kw','vc_ratio','corrT_win','corrT_sum','pool_kw','wh_kw'}; ...
       'blocks and spikes',     {'jumps5','evday3','evday3_n','evday5','spike1','ev_kw','l1_nights','plat_kw','n_runh','start_conc','ht_sum','ht_win'}; ...
       'zero readings',         {'zero_all','zero_mid','zero_night','solar_idx','zero_eve','mid_night'}; ...
       'shape of the day',      [compose('p%02d',0:23) {'eve_share'}]; ...
       'shape of the year',     [compose('m%02d',1:12) {'winsum','low_months'}]};
fprintf('%d measurements: %s\n',numel(feat),strjoin(cellfun(@(a,b) sprintf('%s %d',a,sum(ismember(feat,b))),fam(:,1),fam(:,2),'UniformOutput',false)',' | '));

for tech = ["EV","HP"]
    y = F.("is"+tech)(hd)==1; N = round(sum(y)/0.25);
    fprintf('\n=== %s: %d listed homes, top N = %d ===\n',tech,sum(y),N);
    fprintf('1) Model types (same 5 folds)\n');
    show('logistic regression with lasso (our model)',cvlin(X,y,3),y,N);
    show('random forest (150 trees)',cvmodel('forest',X,y),y,N);
    show('boosted trees (150 rounds)',cvmodel('boost',X,y),y,N);
    if exist('fitcnet','file')>0, show('neural network (one hidden layer, 16 units)',cvmodel('net',X,y),y,N); end
    fprintf('2) Families of measurements (logistic regression)\n');
    for i = 1:size(fam,1)
        m = ismember(feat,fam{i,2});
        show(sprintf('without %s',fam{i,1}),cvlin(X(:,~m),y,1),y,N);
        show(sprintf('only %s',fam{i,1}),cvlin(X(:,m),y,1),y,N);
    end
    a = zeros(1,numel(feat)); for j = 1:numel(feat), [~,~,~,a(j)] = perfcurve(y,X(:,j),true); end
    [~,o] = sort(abs(a-0.5),'descend');
    fprintf('   best single measurements (AUC alone; below 0.5 = lower in listed homes): %s\n',strjoin(compose('%s %.2f',string(feat(o(1:10)))',a(o(1:10))')',', '));
    fprintf('   (%.0f s)\n',toc);
end

fprintf('\n=== Weather ===\n');
y = F.isPV(hd)==1; m = ismember(feat,fam{2,2});
fprintf('Solar list (%d homes): AUC with all measurements, without weather, weather only\n',sum(y));
show('all measurements',cvlin(X,y,1),y,sum(y)); show('without weather response',cvlin(X(:,~m),y,1),y,sum(y)); show('only weather response',cvlin(X(:,m),y,1),y,sum(y));
yH = F.isHP==1 & has; A = C.pEHeat>=0.5 & has; k = sum(A);
ws = F.winsum; ws(~has | ~isfinite(ws)) = -Inf; [~,o] = sort(ws,'descend'); B = false(size(A)); B(o(1:k)) = true;
fprintf('Electric heating without weather data: the %d homes with the highest winter-to-summer ratio (calendar only)\n',k);
fprintf('   weather rule (bH and heating share): %d of the %d listed heat-pump homes are in it | calendar rule: %d | homes in both rules: %d of %d\n',sum(yH & A),sum(yH),sum(yH & B),sum(A & B),k);
kA = sum(C.AC & has);
sr = F.m07./max((F.m05+F.m10)/2,0.05); sr(~has | ~isfinite(sr)) = -Inf; [~,o2] = sort(sr,'descend'); B2 = false(size(A)); B2(o2(1:kA)) = true;
fprintf('Air conditioning without weather data: the %d homes with the highest July-to-spring/fall ratio agree with our AC flag in %d homes (%.0f%%)\n',kA,sum(B2 & C.AC & has),100*sum(B2 & C.AC & has)/kA);

load(fullfile(wd,'ami.mat'),'kwh'); nD = 365;
rs = hd(round(linspace(1,numel(hd),6000))); Xk = double(kwh(rs,:)); Xk(isnan(Xk)) = 0;
Dk = squeeze(sum(reshape(Xk,numel(rs),24,nD),2));
dm = @(x) mean(reshape(double(x(1:8760)),24,nD),1)';
Td = dm(temp); Hu = dm(humid); Wi = dm(wind); hddv = max(60-Td,0); cddv = max(Td-65,0);
D1 = [ones(nD,1) hddv cddv]; D2 = [D1 cddv.*(Hu-mean(Hu))/std(Hu) hddv.*(Wi-mean(Wi))/std(Wi)];
r2 = @(D) 1 - sum((Dk - (Dk*pinv(D)')*D').^2,2)./sum((Dk-mean(Dk,2)).^2,2);
ra = r2(D1); rb = r2(D2); ac = C.AC(rs); eh = C.pEHeat(rs)>=0.5;
fprintf('Humidity and wind (6,000 homes): share of day-to-day change explained by temperature alone, then with humidity (on hot days) and wind (on cold days) added\n');
fprintf('   AC homes: %.1f%% -> %.1f%% | electric-heating homes: %.1f%% -> %.1f%% | all homes: %.1f%% -> %.1f%%\n', ...
    100*median(ra(ac),'omitnan'),100*median(rb(ac),'omitnan'),100*median(ra(eh),'omitnan'),100*median(rb(eh),'omitnan'),100*median(ra,'omitnan'),100*median(rb,'omitnan'));
fprintf('\nDone in %.0f s\n',toc);

function g = cvlin(X,y,nrep)
% out-of-fold score, logistic regression
Xs = normalize(sign(X).*log1p(abs(X))); Xs(~isfinite(Xs)) = 0; g = zeros(size(y));
for s = 1:nrep
    rng(s); cvm = fitclinear(Xs,y,'Learner','logistic','Regularization','lasso','Lambda',0.004,'Prior','uniform','KFold',5);
    [~,sc] = kfoldPredict(cvm); g = g + normalize(sc(:,2));
end
end

function g = cvmodel(kind,X,y)
% out-of-fold score, other model types
rng(1); cv = cvpartition(y,'KFold',5); g = zeros(size(y)); Xs = sign(X).*log1p(abs(X));
for f = 1:5
    tr = find(training(cv,f)); te = test(cv,f);
    pos = tr(y(tr)); neg = tr(~y(tr)); neg = neg(randperm(numel(neg),min(numel(neg),6000))); ii = [pos; neg];
    switch kind
        case 'forest', mdl = fitcensemble(Xs(ii,:),y(ii),'Method','Bag','NumLearningCycles',150,'Learners',templateTree('MinLeafSize',10));
        case 'boost',  mdl = fitcensemble(Xs(ii,:),y(ii),'Method','LogitBoost','NumLearningCycles',150,'LearnRate',0.1,'Learners',templateTree('MaxNumSplits',6));
        case 'net',    mdl = fitcnet(Xs(ii,:),y(ii),'LayerSizes',16,'Standardize',true,'Lambda',1e-3);
    end
    [~,sc] = predict(mdl,Xs(te,:)); g(te) = sc(:,2);
end
end

function show(name,g,y,N)
[~,o] = sort(g,'descend'); L = cumsum(double(y(o))); [~,~,~,auc] = perfcurve(y,g,true); k = [min(200,N) N];
fprintf('   %-46s AUC %.3f | top %d: %.0f%% (+/-%.0f) | top %d: %.0f%% (+/-%.0f)\n',name,auc,k(1),100*L(k(1))/(0.25*k(1)),400*sqrt(L(k(1)))/k(1),k(2),100*L(k(2))/(0.25*k(2)),400*sqrt(L(k(2)))/k(2));
end
