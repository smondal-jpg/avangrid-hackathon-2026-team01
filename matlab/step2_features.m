% Load and weather features per customer -> features.mat

wd = 'E:\Avangrid Hackathon 2026\matlab';
load(fullfile(wd,'ami.mat'));
[nC,nH] = size(kwh); nD = nH/24;

hod  = hour(t)';  mon = month(t)';
Tday = mean(reshape(temp,24,nD),1)';
Td   = repelem(Tday',24);
monD = mon(1:24:end);
hdd  = max(60-Tday,0);  cdd = max(Tday-65,0);
A    = [ones(nD,1) hdd cdd];

sysLoad = sum(kwh,1,'omitnan');                          % system load
[~,ord] = sort(sysLoad,'descend'); topH = ord(1:100);     % 100 system peak hours

win = ismember(mon,[12 1 2]); smr = ismember(mon,[6 7 8]); s2 = ismember(mon,4:9);
nightH = ismember(hod,0:5); mid = ismember(hod,10:14);
mildD  = ismember(mon,[4 5 9 10]) & Td>50 & Td<68;
poolS  = smr & Td<70 & hod>=8 & hod<=18;
poolR  = ismember(mon,[4 5 10]) & Td>=55 & Td<70 & hod>=8 & hod<=18;
nhD    = ismember(monD(2:end),4:10);
evPair = [1 2 20 21 22 23];                              % blocks starting 19:00-01:00

names = [{'missing','ann_kwh','mean_kw','max_kw','p99','p95','base5','load_factor','zero_all', ...
    'jumps5','evday3','evday3_n','evday5','spike1','ev_kw','l1_nights','plat_kw','n_runh','start_conc', ...
    'b0','heat_slope','cool_slope','r2','heat_share','cool_share','cold_night','mild_night','cold_ratio', ...
    'vcold_kw','c30_kw','vc_ratio','corrT_win','corrT_sum', ...
    'zero_mid','zero_night','solar_idx','zero_eve','mid_night','winsum', ...
    'pool_kw','wh_kw','ht_sum','ht_win','low_months','cp_kw','cp_ratio','eve_share','peak_month'}, ...
    compose('p%02d',0:23), compose('m%02d',1:12)];
F = zeros(nC,numel(names));

chunk = 2000; tic
for c0 = 1:chunk:nC
    r = c0:min(c0+chunk-1,nC); n = numel(r);
    X = double(kwh(r,:));
    nanM = isnan(X); miss = mean(nanM,2);
    med = median(X,2,'omitnan'); med(isnan(med)) = 0;
    X(nanM) = 0; X = X + nanM.*med;

    D  = reshape(X,n,24,nD);                              % customer x hour x day
    Ex = D - median(D,2);                                 % excess over the day's median

    hi3 = Ex>=3; two3 = hi3(:,1:23,:) & hi3(:,2:24,:); evd = any(two3,2);
    hi5 = Ex>=5; two5 = hi5(:,1:23,:) & hi5(:,2:24,:);
    evday3 = mean(evd,3); evday3n = mean(any(two3(:,evPair,:),2),3); evday5 = mean(any(two5,2),3);
    spike1 = sum(sum(hi3 & ~circshift(hi3,1,2) & ~circshift(hi3,-1,2),2),3);
    mx = max(Ex,[],2); evkw = sum(mx.*evd,3)./max(sum(evd,3),1);
    z = false(n,1,nD); inrun = cat(2,two3,z) | cat(2,z,two3);   % hours inside blocks
    pl = Ex; pl(~inrun) = NaN; platkw = median(reshape(pl,n,[]),2,'omitnan'); platkw(isnan(platkw)) = 0;
    nrun = sum(reshape(inrun,n,[]),2);
    [~,fh] = max(two3,[],2); sc = zeros(n,23);
    for h = 1:23, sc(:,h) = sum((fh==h) & evd,3); end
    sconc = max(sc,[],2)./max(sum(sc,2),1);

    nightSeq = cat(2,D(:,19:24,1:end-1),D(:,1:8,2:end));   % 18:00 -> 08:00
    l1 = nightSeq - prctile(reshape(nightSeq,n,[]),10,2);
    l1 = l1>=0.9 & l1<=2.3; run = zeros(n,1,nD-1); cnt = run;
    for h = 1:14, run = (run+1).*l1(:,h,:); cnt = max(cnt,run); end
    l1n = sum(squeeze(cnt>=6) & nhD,2)/sum(nhD);

    dk = squeeze(sum(D,2));                               % daily kWh
    cf = (A\dk')'; pr = cf*A';
    r2 = 1 - sum((dk-pr).^2,2)./sum((dk-mean(dk,2)).^2,2);

    cn = mean(X(:,nightH & temp'<25),2); mn = mean(X(:,nightH & temp'>50 & temp'<65),2);
    vc = mean(X(:,temp'<10),2); c30 = mean(X(:,temp'>25 & temp'<35),2);
    cw = corr(X(:,win)',temp(win)); cs = corr(X(:,smr)',temp(smr));

    zm = mean(X(:,mid&s2)<0.05,2); zn = mean(X(:,ismember(hod,0:4)&s2)<0.05,2);
    ze = mean(X(:,ismember(hod,18:22))<0.05,2);
    mnr = mean(X(:,mid&s2),2)./max(mean(X(:,ismember(hod,0:4)&s2),2),0.02);

    nt = D(:,1:6,:); ht = any((nt - median(nt,2))>=1.5 & ~inrun(:,1:6,:),2);   % heater cycles, EV blocks excluded
    mo = zeros(n,12); for m = 1:12, mo(:,m) = mean(X(:,mon==m),2); end
    [~,pm] = max(mo,[],2);
    prof = squeeze(mean(D,3)); prof = prof./max(mean(prof,2),1e-3);

    F(r,:) = [miss, sum(X,2), mean(X,2), max(X,[],2), prctile(X,99,2), prctile(X,95,2), prctile(X,5,2), ...
        mean(X,2)./max(max(X,[],2),1e-3), mean(X<0.02,2), sum(diff(X,1,2)>=5,2), ...
        evday3, evday3n, evday5, spike1, evkw, l1n, platkw, nrun, sconc, ...
        cf, r2, cf(:,2)*sum(hdd)./sum(dk,2), cf(:,3)*sum(cdd)./sum(dk,2), cn, mn, cn./max(mn,0.05), ...
        vc, c30, vc./max(c30,0.05), cw, cs, ...
        zm, zn, zm-zn, ze, mnr, mean(X(:,win),2)./max(mean(X(:,smr),2),0.05), ...
        mean(X(:,poolS),2)-mean(X(:,poolR),2), ...
        mean(X(:,mildD & ismember(hod,6:8)),2)-mean(X(:,mildD & ismember(hod,2:4)),2), ...
        mean(any(ht(:,:,ismember(monD,6:8)),2),3), mean(any(ht(:,:,ismember(monD,[12 1 2])),2),3), ...
        sum(mo<0.25*median(mo,2),2), mean(X(:,topH),2), mean(X(:,topH),2)./max(mean(X,2),1e-3), ...
        sum(X(:,ismember(hod,16:20)),2)./max(sum(X,2),1e-3), pm, prof, mo./max(mean(mo,2),1e-3)];
    fprintf('features %5d-%5d  %.0f s\n',r(1),r(end),toc);
end
F(~isfinite(F)) = 0;
F = array2table(F,'VariableNames',names);
F.id = (1:nC)';
F.hasData = F.missing < 0.5;
F.isEV = ismember(F.id,evID) & F.hasData; F.isHP = ismember(F.id,hpID) & F.hasData; F.isPV = ismember(F.id,pvID) & F.hasData;

save(fullfile(wd,'features.mat'),'F','sysLoad','topH','temp','humid','wind','t','Tday','hdd','cdd');
fprintf('Saved features.mat: %d customers with data x %d features | labels with data: EV %d HP %d PV %d\n',sum(F.hasData),numel(names),sum(F.isEV),sum(F.isHP),sum(F.isPV));
