% Summary tables of the weather model: kWh(day) = b0 + bH*HDD60 + bC*CDD65

wd = 'E:\Avangrid Hackathon 2026\matlab';
load(fullfile(wd,'features.mat'),'F'); load(fullfile(wd,'class.mat'),'C'); load(fullfile(wd,'seg.mat'),'G');
ok = F.hasData==1 & ~C.Vacant; sig = @(x,x0,s) 1./(1+exp(-(x-x0)/s));
V = [F.b0 F.heat_slope F.cool_slope 100*F.heat_share 100*F.cool_share];
row = @(name,m) fprintf('%-34s %6d   %6.1f   %6.2f   %6.2f   %5.0f%%   %5.0f%%\n',name,sum(m),median(V(m,:),1,'omitnan'));
head = @() fprintf('%-34s %6s   %6s   %6s   %6s   %6s   %6s\n','','homes','b0','bH','bC','heat%','cool%');

q = [5 25 50 75 95];
fprintf('\n1) Range over the %d occupied homes (percentiles)\n',sum(ok));
disp(array2table(round(prctile(V(ok,1:3),q)',2),'VariableNames',compose('p%d',q), ...
    'RowNames',{'b0 (kWh on a mild day)','bH (kWh per F colder)','bC (kWh per F hotter)'}));

fprintf('2) Median values per technology (heat%% / cool%% = share of yearly use that follows cold / hot weather)\n'); head();
none = ok & ~(C.AC | C.EHeat | C.WaterHeat | C.Pool | C.HotTub | C.EV | C.Solar | C.SuppHeat);
T = {'All occupied homes',ok; 'No technology found',none; 'Air conditioning',C.AC & ok; 'Electric heating (any)',C.EHeat & ok; ...
     '  heat pump',C.HP & ok; '  resistance heat',C.ResistHeat & ok; 'Supplemental electric heat',C.SuppHeat & ok; ...
     'EV',C.EV & ok; 'Water heater',C.WaterHeat & ok; 'Pool',C.Pool & ok; 'Hot tub',C.HotTub & ok; 'Solar',C.Solar & ok};
for i = 1:size(T,1), row(T{i,1},T{i,2}); end

fprintf('\n3) Median values per group of homes\n'); head();
cn = unique(G.cluster_name(ok));
for i = 1:numel(cn), row(char(cn(i)),ok & G.cluster_name==cn(i)); end
S = {'High electrification',G.HighElectrification; 'Peak contributor',G.PeakContributor; 'Flexible load',G.FlexibleLoad; 'Likely next adopter',G.EmergingCandidate};
for i = 1:size(S,1), row(S{i,1},S{i,2} & ok); end
P = F{:,compose('p%02d',0:23)}; Mo = F{:,compose('m%02d',1:12)};
fprintf('\n   What defines each group (1.00 = the home''s own average)\n');
fprintf('%-58s %7s %7s %7s %7s %7s %7s\n','','0-5h','6-9h','10-15h','17-21h','Dec-Feb','Jun-Aug');
for i = 1:numel(cn)
    m = ok & G.cluster_name==cn(i);
    fprintf('%-58s %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n',char(cn(i)),mean(mean(P(m,1:5),2)),mean(mean(P(m,7:10),2)), ...
        mean(mean(P(m,11:16),2)),mean(mean(P(m,18:22),2)),mean(mean(Mo(m,[12 1 2]),2)),mean(mean(Mo(m,6:8),2)));
end

fprintf('\n4) What happens if a threshold moves\n');
fprintf('   Air conditioning, cooling slope bC at least ...\n');
for x = [0.3 0.4 0.5 0.6 0.7]
    fprintf('      %.1f : %5d homes (%.0f%% of occupied homes)\n',x,sum(ok & sig(F.cool_slope,x,0.15).*sig(F.corrT_sum,0.3,0.08)>=0.5), ...
        100*mean(sig(F.cool_slope(ok),x,0.15).*sig(F.corrT_sum(ok),0.3,0.08)>=0.5));
end
fprintf('   Electric heating, heating slope bH at least ... (heating share at least 25%%)\n');
hp = F.isHP==1 & ok;
for x = [0.4 0.5 0.6 0.7 0.8]
    p = sig(F.heat_share,0.25,0.05).*sig(F.heat_slope,x,0.15)>=0.5;
    fprintf('      %.1f : %5d homes (%.1f%% of occupied homes) | %.0f%% of Avangrid heat-pump list passes\n',x,sum(p & ok),100*mean(p(ok)),100*mean(p(hp)));
end
fprintf('   Solar, midday-minus-night zero gap at least ...\n');
pv = F.isPV==1 & ok;
for x = [0.1 0.2 0.3 0.4]
    p = F.solar_idx>=x & ok;
    fprintf('      %2.0f points : %4d homes | finds %.0f%% of Avangrid solar list | %.0f%% of our flags are on the list\n',100*x,sum(p),100*mean(p(pv)),100*mean(pv(p)));
end
fprintf('   Pool, extra summer daytime load at least ...\n');
for x = [0.6 0.8 1.0 1.2 1.5], fprintf('      %.1f kW : %5d homes\n',x,sum(ok & F.pool_kw>=x)); end
fprintf('   Water heater, morning rise at least ... (mild-day use at least 12 kWh)\n');
for x = [0.30 0.45 0.60], fprintf('      %.2f kW : %5d homes\n',x,sum(ok & sig(F.wh_kw,x,0.12).*sig(F.b0,12,4)>=0.5)); end
fprintf('   Hot tub, share of nights with a heater cycle at least ... (summer and winter)\n');
for x = [0.2 0.3 0.4], fprintf('      %.0f%% : %5d homes\n',100*x,sum(ok & sig(F.ht_sum,x,0.07).*sig(F.ht_win,x,0.07)>=0.5)); end
fprintf('   EV Level 2, block power at least ... (EV homes only)\n');
for x = [2 3 4], fprintf('      %d kW : %5d homes\n',x,sum(C.EV & F.ev_kw>=x)); end
