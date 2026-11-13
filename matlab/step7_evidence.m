% One example home per technology with the rule drawn on it -> figs\ev_*.png, figs\hm_*.png

wd = 'E:\Avangrid Hackathon 2026\matlab';
load(fullfile(wd,'ami.mat'),'kwh'); load(fullfile(wd,'features.mat')); load(fullfile(wd,'class.mat'),'C');
fd = fullfile(wd,'figs'); nD = 365;
hod = hour(t)'; mon = month(t)'; monD = mon(1:24:end); Td = repelem(Tday',24);
teal = [14 124 134]/255; org = [232 159 0]/255; gry = [92 107 115]/255; red = [200 65 75]/255;
set(groot,'defaultAxesFontName','Times New Roman','defaultTextFontName','Times New Roman','defaultAxesFontSize',10);
pick   = @(mask,score) find(mask & score==max(score(mask)),1);
getx   = @(id) fillmissing(double(kwh(id,:)),'constant',median(double(kwh(id,:)),'omitnan'));
prof   = @(x,m) arrayfun(@(h) mean(x(m & hod==h)),0:23);
newfig = @() figure('Position',[100 100 370 255],'Color','w');
done   = @(f,name) exportgraphics(f,fullfile(fd,name),'Resolution',300);
plain  = ~C.Solar & ~C.EHeat & ~C.Vacant;
Ts = linspace(min(Tday),max(Tday),100)';
fitT = @(id) F.b0(id) + F.heat_slope(id)*max(60-Ts,0) + F.cool_slope(id)*max(Ts-65,0);
dayk = @(id) sum(reshape(getx(id),24,nD),1)';

% EV Level 2
id = pick(C.EV_L2 & F.isEV & plain & F.plat_kw>=5 & F.plat_kw<=9, F.evday3);
x = getx(id); D = reshape(x,24,nD); thr = median(D,1)+3;
blk = D>=thr; blk = blk & ([blk(2:end,:); false(1,nD)] | [false(1,nD); blk(1:end-1,:)]);
s = movsum(any(blk,1).*ismember(monD,[4 5 9 10]),[0 6]); [~,d0] = max(s(1:nD-6));
h = (d0-1)*24+(1:168); b = blk(:,d0:d0+6); b = b(:)';
f = newfig(); hold on
h1 = bar(0:167,x(h).*b,1,'FaceColor',org,'EdgeColor','none');
h2 = plot(0:167,x(h),'Color',teal,'LineWidth',1.2);
h3 = plot(0:167,repelem(thr(d0:d0+6),24),'--','Color',gry);
set(gca,'XTick',0:24:168,'XTickLabel',[]); grid on; xlim([0 167]); ylim([0 1.4*max(x(h))]); ylabel('kW');
xlabel('7 days in a row (each grid line = midnight)');
legend([h1 h2 h3],{'charging block','home load','normal level + 3 kW'},'Location','northeast','Box','off','FontSize',9);
done(f,'ev_l2.png'); heatfig(getx(id),fd,'hm_l2.png');
fprintf('EV L2: customer %d | block power %.1f kW | %.0f days a year with a block\n',id,F.plat_kw(id),365*F.evday3(id));

% EV Level 1
ns = ismember(hod,[18:23 0:7]); ok = ismember(mon,4:10); best = -inf;
m = C.EV_L1 & ~C.Solar & F.zero_mid<0.1 & F.missing<0.02; if ~any(m), m = C.EV_L1; end
for k = find(m)'
    xk = getx(k); bk = prctile(xk(ns),10);
    sk = mean(xk(ns&ok)-bk>=0.9 & xk(ns&ok)-bk<=2.3) - 2*mean(xk(~ns&ok)-bk>=0.9);
    if sk>best, best = sk; id = k; end
end
x = getx(id); base = prctile(x(ns),10);
inb = ns & x-base>=0.9 & x-base<=2.3;
s = movsum((inb - 2*(~ns & x-base>=0.9)).*ok,[0 119]); [~,h0] = max(s(1:end-119));
h0 = min(floor((h0-1)/24)*24+13,numel(x)-119); h = h0:h0+119;
f = newfig(); hold on
h1 = patch([0 119 119 0],base+[0.9 0.9 2.3 2.3],org,'FaceAlpha',0.25,'EdgeColor','none');
h2 = plot(0:119,x(h),'Color',teal,'LineWidth',1.2);
h3 = yline(base,'--','Color',gry);
set(gca,'XTick',0:24:120,'XTickLabel',[]); grid on; xlim([0 119]); ylim([0 max(4.5,1.5*max(x(h)))]); ylabel('kW');
xlabel('5 days in a row (each grid line = noon)');
legend([h1 h2 h3],{'Level 1 band: 0.9-2.3 kW above base','home load','night base'},'Location','northeast','Box','off','FontSize',9);
done(f,'ev_l1.png'); heatfig(getx(id),fd,'hm_l1.png');
fprintf('EV L1: customer %d | %.0f%% of Apr-Oct nights have a 6+ hour low step\n',id,100*F.l1_nights(id));

% air conditioning
id = pick(C.AC & plain & ~C.EV & ~C.Pool & F.cool_slope>=1 & F.cool_slope<=2, F.r2);
dk = dayk(id);
f = newfig(); hold on
scatter(Tday,dk,8,gry,'filled','MarkerFaceAlpha',0.5); plot(Ts,fitT(id),'Color',org,'LineWidth',2);
xline(65,':'); xlabel(sprintf('Daily mean outdoor temperature (%cF)',176)); ylabel('kWh that day');
text(min(Tday)+2,0.93*max(dk),sprintf('cooling slope = %.1f kWh per %cF',F.cool_slope(id),176),'Color',org,'FontWeight','bold');
done(f,'ev_ac.png'); heatfig(getx(id),fd,'hm_ac.png');
fprintf('AC: customer %d | cooling slope %.2f | heating slope %.2f\n',id,F.cool_slope(id),F.heat_slope(id));

% electric heating
i1 = pick(C.ResistHeat & ~C.AC & ~C.EV & ~C.Solar & F.heat_slope>=1.5, F.r2);
m2 = C.HP & C.pEHeat>=0.5 & C.AC & ~C.Solar; if any(m2 & F.isHP), m2 = m2 & F.isHP; end
i2 = pick(m2,F.r2);
f = newfig(); hold on
scatter(Tday,dayk(i1),8,red,'filled','MarkerFaceAlpha',0.4); scatter(Tday,dayk(i2),8,teal,'filled','MarkerFaceAlpha',0.4);
h1 = plot(Ts,fitT(i1),'Color',red,'LineWidth',2); h2 = plot(Ts,fitT(i2),'Color',teal,'LineWidth',2);
xline(60,':'); xlabel(sprintf('Daily mean outdoor temperature (%cF)',176)); ylabel('kWh that day');
legend([h1 h2],{sprintf('resistance heat: %.1f kWh per %cF colder',F.heat_slope(i1),176), ...
    sprintf('heat pump: %.1f heating, %.1f cooling',F.heat_slope(i2),F.cool_slope(i2))},'Location','northeast','Box','off','FontSize',9);
done(f,'ev_heat.png'); heatfig(getx(i1),fd,'hm_heat.png');
fprintf('Resistance heat: customer %d | heating slope %.2f | heating share %.0f%%\n',i1,F.heat_slope(i1),100*F.heat_share(i1));
fprintf('Heat pump: customer %d | heating slope %.2f | cooling slope %.2f | heating share %.0f%% | on Avangrid list %d\n',i2,F.heat_slope(i2),F.cool_slope(i2),100*F.heat_share(i2),F.isHP(i2));

% water heater
mildD = ismember(mon,[4 5 9 10]) & Td>50 & Td<68;
cand = find(C.WaterHeat & plain & ~C.EV & ~C.Pool & ~C.TimerLoad & F.missing<0.02 & F.wh_kw>=0.6 & F.wh_kw<=1.8);
cand = cand(round(linspace(1,numel(cand),min(numel(cand),300)))); best = -inf;
for k = cand'
    xk = getx(k); mr = arrayfun(@(mm) mean(xk(mon==mm & ismember(hod,6:8))) - mean(xk(mon==mm & ismember(hod,2:4))),1:12);
    if min(mr)>best, best = min(mr); id = k; end
end
p = prof(getx(id),mildD);
rs = find(F.hasData & ~C.Vacant); rs = rs(round(linspace(1,numel(rs),1500)));
Xs = double(kwh(rs,:)); typ = arrayfun(@(h) median(mean(Xs(:,mildD & hod==h),2,'omitnan'),'omitnan'),0:23);
f = newfig(); hold on; yl = 1.15*max(p);
patch([1.5 4.5 4.5 1.5],[0 0 yl yl],gry,'FaceAlpha',0.15,'EdgeColor','none');
patch([5.5 8.5 8.5 5.5],[0 0 yl yl],org,'FaceAlpha',0.25,'EdgeColor','none');
h1 = plot(0:23,typ,'Color',gry,'LineWidth',1.5); h2 = plot(0:23,p,'Color',teal,'LineWidth',2);
xlim([0 23]); ylim([0 yl]); xlabel('Hour of day (mild spring and fall days)'); ylabel('kW');
text(9,0.93*yl,sprintf('6-9 AM minus 2-5 AM = %.2f kW',F.wh_kw(id)),'Color',org,'FontWeight','bold');
legend([h2 h1],{'this home','typical home'},'Location','east','Box','off','FontSize',9);
done(f,'ev_wh.png'); heatfig(getx(id),fd,'hm_wh.png');
fprintf('Water heater: customer %d | morning rise %.2f kW | mild-day use %.1f kWh\n',id,F.wh_kw(id),F.b0(id));

% pool
id = pick(C.Pool & plain & ~C.EV & F.pool_kw<4 & F.winsum<0.9, F.pool_kw.*F.r2);
x = getx(id); ps = prof(x,ismember(mon,[6 7 8]) & Td<70); pr = prof(x,ismember(mon,[4 5 10]) & Td>=55 & Td<70);
f = newfig(); hold on; hh = 8:18;
h1 = patch([hh fliplr(hh)],[ps(hh+1) fliplr(pr(hh+1))],org,'FaceAlpha',0.35,'EdgeColor','none');
h2 = plot(0:23,pr,'Color',gry,'LineWidth',1.5); h3 = plot(0:23,ps,'Color',teal,'LineWidth',2);
xlim([0 23]); ylim([0 1.5*max(ps)]); xlabel('Hour of day'); ylabel('kW');
legend([h3 h2 h1],{sprintf('cool summer days (below 70 %cF)',176),sprintf('spring / fall days (55-70 %cF)',176),sprintf('extra daytime load = %.1f kW',F.pool_kw(id))}, ...
    'Location','northwest','Box','off','FontSize',9);
done(f,'ev_pool.png'); heatfig(getx(id),fd,'hm_pool.png');
fprintf('Pool: customer %d | extra summer daytime load %.2f kW\n',id,F.pool_kw(id));

% hot tub
m = C.HotTub & ~C.EV & ~C.EHeat & F.missing<0.02; m2 = m & F.ht_sum>=0.4 & F.ht_sum<=0.7 & F.ht_win>=0.4 & F.ht_win<=0.7; if any(m2), m = m2; end
id = pick(m,-abs(F.ht_sum-0.55)-abs(F.ht_win-0.55));
D = reshape(getx(id),24,nD); nt = D(1:6,:); ex = max(nt-median(nt,1),[],1);
f = newfig(); hold on
h1 = bar(1:nD,ex.*(ex>=1.5),1,'FaceColor',org,'EdgeColor','none'); h2 = bar(1:nD,ex.*(ex<1.5),1,'FaceColor',gry,'EdgeColor','none');
h3 = yline(1.5,'--','Color',red,'LineWidth',1.2);
set(gca,'XTick',[1 92 182 274 365],'XTickLabel',{'Sep','Dec','Mar','Jun','Aug'}); xlim([0 366]); ylim([0 1.45*max(ex)]);
xlabel('One bar per night, for the whole year'); ylabel('kW above the rest of the night');
legend([h1 h2 h3],{'night with a heater cycle','other nights','1.5 kW line'},'Location','northwest','Box','off','FontSize',9);
done(f,'ev_hottub.png'); heatfig(getx(id),fd,'hm_hottub.png');
fprintf('Hot tub: customer %d | heater cycle on %.0f%% of summer nights and %.0f%% of winter nights\n',id,100*F.ht_sum(id),100*F.ht_win(id));

% solar and solar + battery
s2 = ismember(mon,4:9); zs = @(x) arrayfun(@(h) 100*mean(x(s2 & hod==h)<0.05),0:23);
i1 = pick(C.Solar & F.isPV & ~C.SolarBatt & ~C.Vacant, F.solar_idx);
m2 = C.SolarBatt & F.solar_idx>=0.3; if any(m2 & F.isPV), m2 = m2 & F.isPV; end
i2 = pick(m2,F.zero_eve);
f = newfig(); hold on
patch([9.5 14.5 14.5 9.5],[0 0 100 100],gry,'FaceAlpha',0.12,'EdgeColor','none');
patch([17.5 22.5 22.5 17.5],[0 0 100 100],org,'FaceAlpha',0.12,'EdgeColor','none');
h1 = plot(0:23,zs(getx(i1)),'Color',teal,'LineWidth',2); h2 = plot(0:23,zs(getx(i2)),'Color',org,'LineWidth',2);
xlim([0 23]); ylim([0 135]); set(gca,'YTick',0:20:100); xlabel('Hour of day (April to September)'); ylabel('% of hours meter reads near zero');
text(10.3,105,'midday','Color',gry,'FontSize',9); text(18.3,105,'evening','Color',org,'FontSize',9);
legend([h1 h2],{'solar home','solar + battery home'},'Location','north','Orientation','horizontal','Box','off','FontSize',9);
done(f,'ev_solar.png'); heatfig(getx(i1),fd,'hm_solar.png');
fprintf('Solar: customer %d | near zero at midday %.0f%%, at night %.0f%%\n',i1,100*F.zero_mid(i1),100*F.zero_night(i1));
fprintf('Solar + battery: customer %d | near zero in the evening %.0f%% of hours | on Avangrid solar list %d\n',i2,100*F.zero_eve(i2),F.isPV(i2));
close all
fprintf('Saved 8 rule plots and 8 heat maps in %s\n',fd);

function heatfig(x,fd,name)
% whole-year heat map of one home
f = figure('Position',[100 100 235 255],'Color','w');
x = reshape(x,24,[]); imagesc(1:size(x,2),0:23,x); axis xy; colormap(gca,turbo); caxis([0 prctile(x(:),99.5)]);
set(gca,'XTick',[1 92 182 274],'XTickLabel',{'Sep','Dec','Mar','Jun'},'YTick',0:6:18); ylabel('Hour of day');
exportgraphics(f,fullfile(fd,name),'Resolution',300);
end
