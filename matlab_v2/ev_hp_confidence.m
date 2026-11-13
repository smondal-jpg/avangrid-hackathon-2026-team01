% EV charging sessions and heat-pump checks -> evhp.mat

wd = 'E:\Avangrid Hackathon 2026\matlab'; nd = 'E:\Avangrid Hackathon 2026\matlab_v2';
tic
load(fullfile(wd,'ami.mat'),'kwh'); load(fullfile(wd,'features.mat'),'F','t','Tday'); load(fullfile(wd,'class.mat'),'C','PR');
n = size(kwh,1); H = 8760; has = F.hasData==1; hd = find(has);
hodv = hour(t(1:H)); hodv = hodv(:); monv = month(t(1:H)); monv = monv(:); Td = Tday(:)';
vc = Td<20; if nnz(vc)<5, vc = Td<=prctile(Td,5); end
sh = Td>=45 & Td<55;
S = zeros(n,8); Hf = zeros(n,4);                          % S: nS pm cons nNight nGood en ln unif
for c = 1:2000:n
    r = c:min(c+1999,n); X = double(kwh(r,1:H)); med = median(X,2,'omitnan'); m = isnan(X); X(m) = 0; X = X + m.*med; X(isnan(X)) = 0;
    S(r,:) = sessions(X,3,hodv,monv);
    Dk = squeeze(sum(reshape(X,numel(r),24,[]),2)); b0 = F.b0(r); bH = F.heat_slope(r);
    Hf(r,1) = (mean(Dk(:,vc),2)-b0)./max(bH*mean(60-Td(vc)),1);
    Hf(r,2) = (mean(Dk(:,sh),2)-b0)./max(bH*mean(60-Td(sh)),1);
    Hf(r,3) = F.cool_slope(r)./(abs(bH)+0.05);
end
if exist(fullfile(nd,'tod.mat'),'file')
    L = load(fullfile(nd,'tod.mat'),'P'); Hf(:,4) = mean(L.P(:,[1 2 8],2),2,'omitnan')./(abs(mean(L.P(:,3:7,2),2,'omitnan'))+1e-3);
end
Hf(~isfinite(Hf)) = 0; S(~isfinite(S)) = 0;
nn = S(:,4); pm = S(:,2);
band = (pm>=2.8 & pm<=3.9) | (pm>=5.5 & pm<=7.8) | (pm>=8.8 & pm<=12);   % charger bands: 3.3-3.8, 5.8-7.7, 9.6-11.5 kW
save(fullfile(nd,'evhp.mat'),'S','band','Hf');
fprintf('Charging sessions and heat-pump physics for %d homes in %.0f s\n',sum(has),toc);
feat = setdiff(F.Properties.VariableNames,{'id','isEV','isHP','isPV','hasData','missing','peak_month'},'stable');
Xold = F{hd,feat};

fprintf('\n=== EV ===\n');
yE = F.isEV==1 & has; oth = has & ~yE;
fprintf('Night charging sessions per year (2+ hours, 3 kW or more above the home''s rolling level):\n');
fprintf('   homes on Avangrid''s EV list: %.0f%% have fewer than 10, %.0f%% fewer than 40 | all other homes: %.0f%% and %.0f%%\n', ...
    100*mean(nn(yE)<10),100*mean(nn(yE)<40),100*mean(nn(oth)<10),100*mean(nn(oth)<40));
ks = [100 200 400 772];
fprintf('Current model (saved in class.mat), estimated precision: %s\n',strjoin(compose('top %d: %.0f%%',[ks' 100*PR.EV(ks,2)])',' | '));
fprintf('Lists made from physical rules:\n');
rep('40+ night sessions a year',has & nn>=40,yE);
rep('40+ night sessions at a standard charger power',has & nn>=40 & band,yE);
rep('80+ night sessions at a standard charger power',has & nn>=80 & band,yE);
rep('120+ night sessions at a standard charger power',has & nn>=120 & band,yE);
rep('Current EV flags (model top 772)',C.EV_flag,yE);
rep('Current flags AND 20+ sessions at charger power',C.EV_flag & nn>=20 & band,yE);
rep('Current flags AND 40+ sessions at charger power',C.EV_flag & nn>=40 & band,yE);
rep('Current flags, NOT 20+ sessions at charger power',C.EV_flag & ~(nn>=20 & band),yE);
fprintf('Retrained models (same method as step 3, 5-fold, 3 runs):\n');
y = yE(hd); Xnew = [Xold S(hd,:) double(band(hd)) nn(hd).*band(hd)];
pk('current 84 measurements',cvscore(Xold,y),y,ks);
g1 = cvscore(Xnew,y); pk('84 + charging-session measurements',g1,y,ks);
pk('same, trained only on listed homes that charge',cvscore(Xnew,y & nn(hd)>=20),y,ks);
pk('physical score alone (night sessions at charger power)',nn(hd).*band(hd),y,ks);

fprintf('\n=== HEAT PUMP ===\n');
yH = F.isHP==1 & has; A = C.pEHeat>=0.5; B = C.SuppHeat; D = C.AC & ~A & ~B; N0 = has & ~A & ~B & ~D;
fprintf('Where the %d homes on Avangrid''s heat-pump list fall (all homes in brackets):\n',sum(yH));
fprintf('   pass electric-heating test %.0f%% (%.0f%%) | supplemental heat %.0f%% (%.0f%%) | AC only %.0f%% (%.0f%%) | none of these %.0f%% (%.0f%%)\n', ...
    100*[mean(A(yH)) mean(A(has)) mean(B(yH)) mean(B(has)) mean(D(yH)) mean(D(has)) mean(N0(yH)) mean(N0(has))]);
g = {'pass electric-heating test',A & yH; 'supplemental heat',B & yH; 'AC only',D & yH; 'none of these',N0 & yH};
for i = 1:4
    fprintf('   listed homes, %-28s %3d homes | median b0 %.1f, bH %.2f, bC %.2f, yearly kWh %.0f\n',g{i,1},sum(g{i,2}), ...
        median(F.b0(g{i,2})),median(F.heat_slope(g{i,2})),median(F.cool_slope(g{i,2})),median(F.ann_kwh(g{i,2})));
end
ks = [50 100 200 372];
fprintf('Current model (saved in class.mat), estimated precision: %s\n',strjoin(compose('top %d: %.0f%%',[ks' 100*PR.HP(ks,2)])',' | '));
fprintf('Lists made from rules:\n');
rep('Electric-heating test',A & has,yH);
rep('Electric-heating test AND passes the AC test',A & C.AC & has,yH);
rep('Supplemental heat AND passes the AC test',B & C.AC & has,yH);
rep('Current heat-pump flags (model top 372)',C.HP_flag,yH);
fprintf('Retrained models (same method as step 3, 5-fold, 3 runs):\n');
y = yH(hd); Xh = [Xold Hf(hd,:)];
pk('current 84 measurements',cvscore(Xold,y),y,ks);
pk('84 + heat-pump physics',cvscore(Xh,y),y,ks);
pool = A(hd) | B(hd); gp = -99*ones(numel(hd),1); gp(pool) = cvscore(Xh(pool,:),y(pool));
pk(sprintf('two-stage: only the %d homes with a heating signal',sum(pool)),gp,y,ks);
pool = A(hd); gp = -99*ones(numel(hd),1); gp(pool) = cvscore(Xh(pool,:),y(pool));
pk(sprintf('two-stage: only the %d electric-heating homes',sum(pool)),gp,y,ks);
fprintf('\nDone in %.0f s. Saved %s\n',toc,fullfile(nd,'evhp.mat'));

function S = sessions(X,thr,hodv,monv)
% runs of 2+ hours at thr kW above the rolling 25-h median
[n,H] = size(X); e = X - movmedian(X,[12 12],2); on = e>=thr;
d = diff(int8([false(n,1) on false(n,1)]),1,2);
[cs,rs] = find((d(:,1:H)==1)'); [ce,~] = find((d(:,2:H+1)==-1)');   % start hour, home, end hour (by home, then time)
ln = ce-cs+1; k = ln>=2; rs = rs(k); cs = cs(k); ce = ce(k); ln = ln(k);
E = cumsum(e,2); i1 = sub2ind([n H],rs,ce); tot = E(i1); m = cs>1; tot(m) = tot(m) - E(sub2ind([n H],rs(m),cs(m)-1));
ef = e(sub2ind([n H],rs,cs)); el = e(i1);
p = max(ef,el); m3 = ln>=3; p(m3) = (tot(m3)-ef(m3)-el(m3))./(ln(m3)-2);
h = hodv(cs); mo = monv(cs); night = h>=18 | h<=4;
nS = accumarray(rs,1,[n 1]); pm = accumarray(rs,p,[n 1],@median);
cons = abs(p-pm(rs))<=0.15*pm(rs);
R = [accumarray(rs,double(ismember(mo,[12 1 2])),[n 1])/90, accumarray(rs,double(ismember(mo,[6 7 8])),[n 1])/92, ...
     accumarray(rs,double(ismember(mo,[3 4 5 9 10 11])),[n 1])/183];
S = [nS, pm, accumarray(rs,double(cons),[n 1])./max(nS,1), accumarray(rs,double(night),[n 1]), accumarray(rs,double(night & cons),[n 1]), ...
     accumarray(rs,tot,[n 1])./max(nS,1), accumarray(rs,ln,[n 1])./max(nS,1), min(R,[],2)./(max(R,[],2)+1e-9)];
end

function g = cvscore(X,y)
% out-of-fold score, 3 random splits
Xs = normalize(sign(X).*log1p(abs(X))); Xs(~isfinite(Xs)) = 0; g = zeros(size(y));
for s = 1:3
    rng(s); cvm = fitclinear(Xs,y,'Learner','logistic','Regularization','lasso','Lambda',0.004,'Prior','uniform','KFold',5);
    [~,sc] = kfoldPredict(cvm); g = g + normalize(sc(:,2));
end
end

function pk(name,g,y,ks)
% AUC and estimated precision of the top k
ks = ks(:); [~,o] = sort(g,'descend'); L = cumsum(double(y(o))); [~,~,~,auc] = perfcurve(y,g,true);
s = compose('top %d: %.0f%% (+/-%.0f)',[ks, 100*L(ks)./(0.25*ks), 400*sqrt(L(ks))./ks]);
fprintf('   %-54s AUC %.3f | %s\n',name,auc,strjoin(s',' | '));
end

function rep(name,m,y)
k = sum(m); L = sum(y(m));
fprintf('   %-52s %5d homes | %3d on Avangrid''s list | est. precision %3.0f%% (+/-%.0f) | finds %2.0f%% of the list\n', ...
    name,k,L,100*L/(0.25*max(k,1)),400*sqrt(L)/max(k,1),100*L/sum(y));
end
