% Demand flexibility and battery (2 x Enphase IQ 10C: 20.16 kWh, 7.68 kW, 89.4%, 10% reserve)
% rtDispatch: hour-by-hour controller on past data only; dispatch: day-ahead plan for comparison
% Output: bess.mat, Avangrid_results.xlsx, figures

wd = 'E:\Avangrid Hackathon 2026\matlab';
load(fullfile(wd,'ami.mat'),'kwh','t');
load(fullfile(wd,'features.mat')); load(fullfile(wd,'class.mat')); load(fullfile(wd,'seg.mat'));
fd = fullfile(wd,'figs'); [nC,nH] = size(kwh); nD = nH/24;
hod = hour(t)'; mon = month(t)'; monD = mon(1:24:end);

bp.Ecap = 20.16; bp.Pmax = 7.68; bp.eta = 0.894; bp.reserve = 0.10; pf = 0.95; bp.Pch = 3*pf;
bp.Edel = bp.Ecap*(1-bp.reserve)*sqrt(bp.eta);          % kWh delivered per day
Pon = 0.28; Poff = 0.10; Pexp = 0.08; capVal = 100;     % $/kWh on-peak, off-peak, export; $/kW-yr
capexKWh = 1000; upgCost = 8000; life = 15; rate = 0.07;   % $/kWh installed, service upgrade $, years, discount rate
annF = (1-(1+rate)^-life)/rate;                          % annuity factor

[~,o8] = sort(mean(reshape(sysLoad,24,nD),2),'descend');
winH = sort(o8(1:8))';                                   % event window: 8 highest hours of the average day
evDay = false(1,nD); evDay(unique(ceil(topH/24))) = true;
bp.dm = ones(1,24,nD); bp.dm(1,:,evDay) = 0; bp.dm(1,winH,evDay) = 1;   % discharge allowed
bp.cm = ones(1,24,nD); bp.cm(1,:,evDay) = 0; bp.cm(1,1:7,evDay) = 1;   % charging allowed
bp.evDay = evDay; bp.solarDay = reshape(ismember(monD,4:9),1,1,[]);
bp.s = 1.1; bp.c = 0.8; bp.emerg = 0.85; bp.HC = 10; bp.hold = 21; bp.inWin = ismember(1:24,winH); bp.winEnd = max(winH);
[~,o6] = sort(mean(reshape(sysLoad,24,nD),2),'descend'); touH = sort(o6(1:6))';
price = Poff*ones(1,24,nD); price(1,touH,:) = Pon;          % 6-h on-peak window
[~,pkDay] = max(max(reshape(sysLoad,24,nD),[],1)); pkHrs = (pkDay-1)*24+(1:24);
fprintf('Event window hours (ending): %s | %d event days | system peak day %s\n', ...
    mat2str(winH),sum(evDay),string(t(pkHrs(1)),'dd-MMM-yyyy'));
summerPeak = mean(ismember(mon(topH),6:8))>0.5;
mildRef = ismember(mon,[4 5 9 10]) & ismember(hod,unique(hod(topH)));

R = zeros(nC,17); dPk = zeros(nC,24); dPkDA = zeros(nC,24); lam = 0.8;
chunk = 2000; tic
for c0 = 1:chunk:nC
    r = c0:min(c0+chunk-1,nC); n = numel(r);
    X = double(kwh(r,:)); nanM = isnan(X);
    med = median(X,2,'omitnan'); med(isnan(med)) = 0; X(nanM) = 0; X = X + nanM.*med;
    D = reshape(X,n,24,nD);
    Fc = dayForecast(D,lam);
    [Dn,disT,chT,bl] = rtDispatch(D,C.Solar(r),bp,Fc,price,Pexp);   % real-time controller
    [Dp,disp_,chp,solCh] = dispatch(D,C.Solar(r),bp,Fc);    % day-ahead plan
    Xn = reshape(Dn,n,nH); Xp = reshape(Dp,n,nH); Edp = sum(disp_,2);
    pk = squeeze(max(D,[],2)); fcErr = 100*mean(abs(squeeze(max(Fc,[],2))-pk)./max(pk,0.1),2);

    L = D(:,winH,:); medPk = median(max(L,[],2),3);
    needE = prctile(squeeze(sum(max(L-medPk,0),2)),90,2);
    dmed = reshape(repmat(median(D,2),1,24,1),n,nH);

    R(r,:) = [max(X,[],2), max(Xn,[],2), mean(pk-squeeze(max(Dn,[],2)),2), ...
        mean(X(:,topH),2), mean(Xn(:,topH),2), disT, ...
        bl, ...
        needE, mean(max(X(:,topH)-dmed(:,topH),0),2), 0.3*max(mean(X(:,topH),2)-mean(X(:,mildRef),2),0), ...
        mean(pk >= 0.6*max(X,[],2),2), chT, ...
        max(Xp,[],2), mean(Xp(:,topH),2), mean(pk-squeeze(max(Dp,[],2)),2), ...
        sum(sum((disp_-chp).*price,2),3) - sum(Edp/bp.eta.*solCh,3)*Pexp, fcErr];
    dPk(r,:) = X(:,pkHrs)-Xn(:,pkHrs); dPkDA(r,:) = X(:,pkHrs)-Xp(:,pkHrs);
    fprintf('bess %5d-%5d  %.0f s\n',r(1),r(end),toc);
end

B = array2table(R,'VariableNames',{'peak_before_kw','peak_after_kw','daily_peak_red_kw','cp_before_kw','cp_after_kw', ...
    'shifted_kwh','bill_savings_usd','need_kwh_p90','flex_ev_kw','flex_thermo_kw','persistent_peak_frac','grid_charge_kwh', ...
    'peak_after_da_kw','cp_after_da_kw','daily_peak_red_da_kw','bill_savings_da_usd','fc_peak_err_pct'});
B.id = F.id;
B.cp_red_kw    = B.cp_before_kw - B.cp_after_kw;
B.cp_red_da_kw = B.cp_before_kw - B.cp_after_da_kw;
B.peak_red_da_pct = 100*(B.peak_before_kw-B.peak_after_da_kw)./max(B.peak_before_kw,0.1);
B.peak_red_pct = 100*(B.peak_before_kw-B.peak_after_kw)./max(B.peak_before_kw,0.1);
B.units_rec    = min(2,max(1,ceil(B.need_kwh_p90/(10.08*(1-bp.reserve)*sqrt(bp.eta)))));
B.capacity_value_usd = B.cp_red_kw*capVal;

flexPool = C.Pool.*max(F.pool_kw,0)*summerPeak;
flexWH   = C.WaterHeat.*0.5.*max(F.wh_kw,0);
B.dr_kw  = C.EV.*B.flex_ev_kw + (C.AC|C.EHeat).*B.flex_thermo_kw + flexPool + flexWH;

incrEV   = C.EV.*F.ev_kw;
incrHeat = C.EHeat.*max(F.vcold_kw-F.mild_night,0);
incrNew  = G.EmergingCandidate*7.2;                      % new Level-2 charger
B.incr_kw   = max([incrEV, incrHeat, incrNew],[],2);
B.over3kVA  = B.incr_kw/pf > 3;
B.batt_cover_kw = min(bp.Pmax,max(B.incr_kw-3*pf,0));
B.units_rec = max(B.units_rec,1+(B.batt_cover_kw>bp.Pmax/2));   % one unit = 3.84 kW
B.new_load  = G.EmergingCandidate & B.over3kVA;             % new load

pr = @(x) tiedrank(x)/numel(x);
elig = F.hasData & ~C.Vacant & ~G.Seasonal & F.max_kw>3*pf;
B.enable = B.batt_cover_kw/bp.Pmax;                      % share of battery power needed
B.score = 0.25*pr(B.cp_red_kw) + 0.20*pr(B.peak_red_pct) + 0.10*pr(B.persistent_peak_frac) + ...
          0.25*B.enable + 0.10*pr(B.bill_savings_usd) + 0.10*C.Solar;
B.capex_usd   = bp.Ecap*capexKWh*ones(nC,1);
B.avoided_upg = upgCost*(B.new_load & B.incr_kw-3*pf <= bp.Pmax);
B.npv_benefit = annF*(B.bill_savings_usd + B.capacity_value_usd) + B.avoided_upg;
B.bcr         = B.npv_benefit./B.capex_usd;
B.score(~elig) = 0;
q = prctile(B.score(elig),[95 85 70]);
B.tier = repmat("Not prioritised",nC,1);
B.tier(B.score>=q(3)) = "Tier 3"; B.tier(B.score>=q(2)) = "Tier 2"; B.tier(B.score>=q(1)) = "Tier 1";
B.tier(~elig) = "Not eligible";

tiers = ["Tier 1","Tier 2","Tier 3"];
P = array2table(zeros(3,14),'VariableNames',{'customers','cp_relief_MW','avg_cp_red_kw','avg_peak_red_pct','shifted_MWh_yr', ...
    'bill_savings_k_usd','capacity_value_k_usd','new_load_enabled','avg_units','capex_M_usd','npv_benefit_M_usd','bcr', ...
    'cp_relief_da_MW','avg_peak_red_da_pct'});
P.tier = tiers';
for i = 1:3
    s = B.tier==tiers(i);
    P{i,1:14} = [sum(s), sum(B.cp_red_kw(s))/1000, mean(B.cp_red_kw(s)), mean(B.peak_red_pct(s)), ...
        sum(B.shifted_kwh(s))/1000, sum(B.bill_savings_usd(s))/1000, sum(B.capacity_value_usd(s))/1000, sum(B.new_load(s)), ...
        mean(B.units_rec(s)), sum(B.capex_usd(s))/1e6, sum(B.npv_benefit(s))/1e6, sum(B.npv_benefit(s))/sum(B.capex_usd(s)), ...
        sum(B.cp_red_da_kw(s))/1000, mean(B.peak_red_da_pct(s))];
end
P = movevars(P,'tier','Before',1);
disp(P)
for i = 1:3
    s = B.tier==tiers(i); v = [sum(B.avoided_upg(s)) annF*sum(B.capacity_value_usd(s)) annF*sum(B.bill_savings_usd(s))]/1e6; bc = sum(v)*1e6/sum(B.capex_usd(s));
    fprintf('%s: benefit $%.1fM = avoided upgrade %.1f + peak relief %.1f + bill savings %.1f | cost $%.1fM | benefit/cost %.2f at $1,000/kWh, %.2f at $700/kWh | breaks even at $%.0f/kWh\n', ...
        tiers(i),sum(v),v,sum(B.capex_usd(s))/1e6,bc,bc/0.7,1000*bc);
    fprintf('        %d homes | new load kept within 3 kVA: %d | load above 3 kVA already there (no upgrade credit): %d | two units needed: %d | one unit enough: %d\n', ...
        sum(s),sum(B.new_load(s)),sum(B.over3kVA(s) & ~B.new_load(s)),sum(B.units_rec(s)==2),sum(B.units_rec(s)==1));
end
s1 = B.tier=="Tier 1"; s12 = ismember(B.tier,tiers(1:2));
pk0 = max(sysLoad); pk1 = max(sysLoad(pkHrs)-sum(dPk(s1,:),1)); pk2 = max(sysLoad(pkHrs)-sum(dPk(s12,:),1));
pk2p = max(sysLoad(pkHrs)-sum(dPkDA(s12,:),1));
fprintf('System peak day (real-time controller): %.1f MW -> %.1f MW (Tier 1) -> %.1f MW (Tier 1+2) | day-ahead plan Tier 1+2: %.1f MW\n',pk0/1000,pk1/1000,pk2/1000,pk2p/1000);
fprintf('Home''s own highest hour, Tier 1+2: real-time %.1f%% lower (rises in %.1f%% of homes) | day-ahead plan %.1f%% lower (rises in %.1f%% of homes)\n', ...
    mean(B.peak_red_pct(s12)),100*mean(B.peak_red_pct(s12)<-0.5),mean(B.peak_red_da_pct(s12)),100*mean(B.peak_red_da_pct(s12)<-0.5));
fprintf('Home''s own highest hour, all eligible homes: real-time %.1f%% lower (rises in %.1f%%) | day-ahead plan %.1f%% lower (rises in %.1f%%)\n', ...
    mean(B.peak_red_pct(elig)),100*mean(B.peak_red_pct(elig)<-0.5),mean(B.peak_red_da_pct(elig)),100*mean(B.peak_red_da_pct(elig)<-0.5));
fprintf('Day-ahead forecast error of daily peak (median customer): %.0f%%\n',median(B.fc_peak_err_pct(elig)));
fprintf('Demand flexibility (no battery) at system peak: %.1f MW (EV %.1f, thermostat %.1f, pool %.1f, water heat %.1f)\n', ...
    sum(B.dr_kw)/1000, sum(C.EV.*B.flex_ev_kw)/1000, sum((C.AC|C.EHeat).*B.flex_thermo_kw)/1000, sum(flexPool)/1000, sum(flexWH)/1000);

f = figure('Position',[100 100 1000 420]); hold on
plot(0:23,sysLoad(pkHrs)/1000,'k','LineWidth',2);
plot(0:23,(sysLoad(pkHrs)-sum(dPk(s1,:),1))/1000,'LineWidth',2);
plot(0:23,(sysLoad(pkHrs)-sum(dPk(s12,:),1))/1000,'LineWidth',2);
plot(0:23,(sysLoad(pkHrs)-sum(dPkDA(s12,:),1))/1000,'--','LineWidth',1.5);
xlabel('Hour of day'); ylabel('System load (MW)'); grid on; xlim([0 23]);
title(sprintf('System peak day %s: batteries at Tier 1 / Tier 1+2 customers',string(t(pkHrs(1)),'dd-MMM-yyyy')));
legend('Today','Tier 1 (real-time controller)','Tier 1+2 (real-time controller)','Tier 1+2 (day-ahead plan)','Location','northwest');
exportgraphics(f,fullfile(fd,'fig_peakday.png'),'Resolution',150);

t1 = find(B.tier=="Tier 1"); [~,j] = max(B.cp_red_kw(t1)); ex = t1(j);
x = double(kwh(ex,:)); x(isnan(x)) = median(x,'omitnan');
x3 = reshape(x,1,24,nD); xn = reshape(rtDispatch(x3,C.Solar(ex),bp,dayForecast(x3,lam),price,Pexp),1,nH);
h = max(1,(pkDay-4)*24+1):min(nH,(pkDay+3)*24);
f = figure('Position',[100 100 1000 380]);
area(t(h),x(h)','FaceColor',[0.8 0.8 0.8],'EdgeColor','none'); hold on
plot(t(h),xn(h)','LineWidth',1.6);
ev = t(h(evDay(ceil(h/24)) & ismember(mod(h-1,24)+1,winH)));
plot(ev,zeros(size(ev)),'k.','MarkerSize',6);
ylabel('kW'); grid on; legend('Without battery','With battery','System-event window','Location','northwest');
title(sprintf('Tier 1 customer %d: week around the system peak (%.1f kW less at system peak hours)',ex,B.cp_red_kw(ex)));
exportgraphics(f,fullfile(fd,'fig_bess_example.png'),'Resolution',150);

save(fullfile(wd,'bess.mat'),'B','P','winH','pkHrs','dPk','dPkDA','bp');

out = fullfile(wd,'Avangrid_results.xlsx');
T = [C(:,{'id','pEV','EV','EV_L2','EV_L1','pHP','HP','EHeat','ResistHeat','pAC','AC','pWaterHeat','WaterHeat', ...
    'pPool','Pool','pHotTub','HotTub','pSolar','Solar','SolarUnreported','SolarBatt','SuppHeat','ShortCycle','AlwaysOn','TimerLoad','Seasonal','Vacant'}), ...
    G(:,{'cluster_name','HighElectrification','PeakContributor','FlexibleLoad','EmergingCandidate'}), ...
    B(:,{'tier','score','units_rec','cp_red_kw','peak_red_pct','shifted_kwh','bill_savings_usd','dr_kw','incr_kw','over3kVA','batt_cover_kw','capex_usd','npv_benefit','bcr','cp_red_da_kw','peak_red_da_pct'}), ...
    F(:,{'ann_kwh','max_kw','cp_kw','ev_kw','heat_slope','cool_slope','pool_kw'})];
writetable(T,out,'Sheet','Customers');
writetable(S,out,'Sheet','Technology summary');
writetable(T2,out,'Sheet','Business segments');
writetable(P,out,'Sheet','BESS tiers');
fprintf('Saved %s\n',out);

function [Dn,disT,chT,bill] = rtDispatch(D,isPV,bp,Fc,price,Pexp)
[n,~,nD] = size(D); se = sqrt(bp.eta); U = bp.Ecap*(1-bp.reserve);
Efull = repmat(bp.Edel,n,1,nD);
[~,a] = topfill(D,Efull,bp.Pmax);        ths  = reshape(a,n,nD);   % level the battery could have held each day
[~,a] = topfill(Fc.*bp.dm,Efull,bp.Pmax); thev = reshape(a,n,nD);   % event-window level from the forecast
dpk = reshape(max(D,[],2),n,nD);
M = movmax(ths,[27 0],2); tha = bp.s*[ths(:,1) M(:,1:end-1)];   % target level, past 28 days
M = movmax(dpk,[27 0],2); P28 = [dpk(:,1) M(:,1:end-1)];    % recent peak, past 28 days
evnext = [bp.evDay(2:end) false];
soc = U*ones(n,1); Dn = D; disT = zeros(n,1); chT = zeros(n,1); bill = zeros(n,1);
for d = 1:nD
    th = tha(:,d); capE = max(th,bp.c*P28(:,d)); em = bp.emerg*P28(:,d); sol = isPV & bp.solarDay(d);
    for h = 1:24
        L = D(:,h,d); ch = zeros(n,1);
        if bp.evDay(d)
            if bp.inWin(h)
                p = min(min(bp.Pmax,max(L-thev(:,d),0)),soc*se);
            else
                p = min(min(bp.Pmax,max(L-em,0)),soc*se);
                if h<=bp.HC || h>bp.winEnd, ch = min(min(bp.Pch,max(capE-L,0)),(U-soc)/se).*(p==0); end
            end
        elseif evnext(d) && h>=bp.hold
            p = min(min(bp.Pmax,max(L-em,0)),soc*se);
            ch = min(min(bp.Pch,max(capE-L,0)),(U-soc)/se).*(p==0);
        else
            p = min(min(bp.Pmax,max(L-th,0)),soc*se);
            ch = min(min(bp.Pch,max(th-L,0)),(U-soc)/se).*(p==0);
        end
        soc = soc + ch*se - p/se;
        g = ch.*~sol;                                     % grid charging
        Dn(:,h,d) = L - p + g;
        disT = disT + p; chT = chT + g;
        bill = bill + (p-g)*price(1,h,d) - (ch-g)*Pexp;
    end
end
end

function [Dn,dis,ch,solCh] = dispatch(D,isPV,bp,Fc)
if nargin < 4, Fc = D; end                               % plan on the actual load
mu = mean(Fc,2);
E = min(bp.Edel,sum(max(Fc-mu,0),2));
E(:,:,bp.evDay) = bp.Edel;
E = min(E,bp.Pch*bp.eta*sum(bp.cm,2));
[~,th] = topfill(Fc.*bp.dm,E,bp.Pmax);                   % plan on the forecast
raw = min(max(D-th,0),bp.Pmax).*bp.dm;                   % applied to the actual load
cum = cumsum(raw,2); dis = max(min(raw,E-(cum-raw)),0);   % stop when the battery is empty
Ed = sum(dis,2);
solCh = isPV & bp.solarDay;
ch = valfill(Fc-dis+(1-bp.cm)*1e3,Ed/bp.eta,bp.Pch).*~solCh;   % charge in valley hours
Dn = D-dis+ch;
end

function Fc = dayForecast(D,lam)
nD = size(D,3); Fc = D; W = D;                           % fading average of previous same weekdays
for d = 2:nD
    if d <= 7, Fc(:,:,d) = D(:,:,d-1); continue, end
    if d > 14, W(:,:,d) = (1-lam)*D(:,:,d-7) + lam*W(:,:,d-7); else, W(:,:,d) = D(:,:,d-7); end
    Fc(:,:,d) = 0.5*D(:,:,d-1) + 0.5*W(:,:,d);
end
end

function [q,hi] = topfill(L,E,P)
lo = zeros(size(E)); hi = max(L,[],2);
for it = 1:35
    th = (lo+hi)/2; big = sum(min(max(L-th,0),P),2) > E;
    lo(big) = th(big); hi(~big) = th(~big);
end
q = min(max(L-hi,0),P);
end

function q = valfill(L,E,P)
lo = min(L,[],2); hi = max(L,[],2)+P;
for it = 1:45
    th = (lo+hi)/2; big = sum(min(max(th-L,0),P),2) > E;
    hi(big) = th(big); lo(~big) = th(~big);
end
q = min(max(lo-L,0),P);
end
