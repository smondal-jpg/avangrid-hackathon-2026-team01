% Weather model per 3-hour block of the day -> tod.mat, figs\fp_*.png

wd = 'E:\Avangrid Hackathon 2026\matlab'; nd = 'E:\Avangrid Hackathon 2026\matlab_v2';
tic
load(fullfile(wd,'ami.mat'),'kwh'); load(fullfile(wd,'features.mat'),'F','temp','t'); load(fullfile(wd,'class.mat'),'C');
if ~exist(fullfile(nd,'figs'),'dir'), mkdir(fullfile(nd,'figs')); end
n = size(kwh,1); nD = 365; has = F.hasData==1;
blk = {'0003','0306','0609','0912','1215','1518','1821','2124'};
mn = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
pname = cell(1,36); for k = 1:36, pname{k} = sprintf('%s%d',mn{ceil(k/3)},k-3*(ceil(k/3)-1)); end
Tb = squeeze(mean(reshape(double(temp(1:8760)),3,8,nD),1));   % temperature of each block
dt = t(1:24:8760); per = (month(dt)-1)*3 + min(ceil(day(dt)/10),3);
A = full(sparse((1:nD)',per(:),1,nD,36)); A = A./sum(A,1);   % 10-day period average

P = nan(n,8,3); R2 = nan(n,8); RES = nan(n,8,36,'single');
for c = 1:2000:n
    r = c:min(c+1999,n); X = double(kwh(r,:));
    med = median(X,2,'omitnan'); m = isnan(X); X(m) = 0; X = X + m.*med;
    Y = squeeze(mean(reshape(X,numel(r),3,8,nD),2));        % homes x block x day, mean kW
    for b = 1:8
        D = [ones(nD,1) max(60-Tb(b,:)',0) max(Tb(b,:)'-65,0)];
        y = squeeze(Y(:,b,:)); B = y*pinv(D)'; res = y - B*D';
        P(r,b,:) = reshape(B,[],1,3);
        R2(r,b) = 1 - sum(res.^2,2)./sum((y-mean(y,2)).^2,2);
        RES(r,b,:) = reshape(res*A,[],1,36);
    end
end
save(fullfile(nd,'tod.mat'),'P','R2','RES','blk','pname');
fprintf('Fingerprints for %d homes in %.0f s | median R2 per block: %s\n',sum(has),toc,mat2str(round(median(R2(has,:),1,'omitnan'),2)));

% median fingerprint per group
G = {'All homes',has; 'Avangrid EV list',has & F.isEV==1; 'Avangrid heat-pump list',has & F.isHP==1; ...
     'Avangrid solar list',has & F.isPV==1; 'Our pool flag',has & C.Pool; 'Our water-heater flag',has & C.WaterHeat};
nG = size(G,1); M = zeros(nG,8,3);
for g = 1:nG
    M(g,:,:) = median(P(G{g,2},:,:),1,'omitnan');
    fprintf('\n%s (%d homes)\n',G{g,1},sum(G{g,2}));
    disp(array2table(round(squeeze(M(g,:,:))',3),'VariableNames',strcat('h',blk), ...
        'RowNames',{'b0 (base kW)','bH (kW per F colder)','bC (kW per F hotter)'}));
end

% AUC with and without the fingerprint
feat = setdiff(F.Properties.VariableNames,{'id','isEV','isHP','isPV','hasData','missing','peak_month'},'stable');
hd = find(has); Xold = F{hd,feat}; Xfp = reshape(P(hd,:,:),numel(hd),24);
Rs = double(RES(hd,:,:)); Xres = [std(Rs,0,3) max(Rs,[],3)-min(Rs,[],3)];
fprintf('\nCross-validated AUC (mean of 3 runs)      %-14s %-14s %-14s %-14s\n','current 84','84 + 24 new','84 + 24 + res','24 new only');
for tech = ["EV","HP","PV"]
    y = F.("is"+tech)(hd)==1;
    a = [cvauc(Xold,y) cvauc([Xold Xfp],y) cvauc([Xold Xfp Xres],y) cvauc(Xfp,y)];
    fprintf('%-3s (%3d labels)                          %-14.3f %-14.3f %-14.3f %-14.3f\n',tech,sum(y),a);
end

% figure 1
set(groot,'defaultAxesFontName','Times New Roman','defaultTextFontName','Times New Roman','defaultAxesFontSize',10);
cm = [linspace(1,14/255,64)' linspace(1,124/255,64)' linspace(1,134/255,64)'];
ttl = {'b0: base load (kW)','bH: extra kW per F colder','bC: extra kW per F hotter'}; fmt = {'%.2f','%.3f','%.3f'};
f = figure('Position',[50 50 1250 330],'Color','w'); tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for p = 1:3
    ax = nexttile; V = M(:,:,p); imagesc(V); colormap(ax,cm); caxis([0 max(V(:))]);
    set(ax,'XTick',1:8,'XTickLabel',{'0-3','3-6','6-9','9-12','12-15','15-18','18-21','21-24'},'YTick',1:nG,'YTickLabel',[]);
    if p==1, set(ax,'YTickLabel',G(:,1)); end
    title(ttl{p}); xlabel('Hour of day');
    for g = 1:nG, for b = 1:8
        text(b,g,sprintf(fmt{p},V(g,b)),'HorizontalAlignment','center','FontSize',8,'Color',[1 1 1]*(V(g,b)>=0.6*max(V(:))));
    end, end
end
exportgraphics(f,fullfile(nd,'figs','fp_groups.png'),'Resolution',300);

% figure 2
cmd = [linspace(0.13,1,32)' linspace(0.40,1,32)' linspace(0.67,1,32)'; linspace(1,0.78,32)' linspace(1,0.25,32)' linspace(1,0.29,32)'];
f = figure('Position',[50 50 1250 420],'Color','w'); tiledlayout(1,nG-1,'TileSpacing','compact','Padding','compact');
Rall = squeeze(median(RES(has,:,:),1,'omitnan'));
for g = 2:nG
    ax = nexttile; V = double(squeeze(median(RES(G{g,2},:,:),1,'omitnan'))) - double(Rall);
    imagesc(1:36,1:8,V); axis xy; colormap(ax,cmd); caxis([-1 1]*max(abs(V(:))));
    set(ax,'XTick',2:6:36,'XTickLabel',mn(1:2:12),'YTick',1:8,'YTickLabel',[]);
    if g==2, set(ax,'YTickLabel',{'0-3','3-6','6-9','9-12','12-15','15-18','18-21','21-24'}); ylabel('Hour of day'); end
    title(G{g,1},'FontSize',9); colorbar('southoutside');
end
exportgraphics(f,fullfile(nd,'figs','fp_residual.png'),'Resolution',300);
close all
fprintf('\nSaved tod.mat and 2 figures in %s (%.0f s in total)\n',nd,toc);

function a = cvauc(X,y)
% 5-fold AUC, 3 random splits
Xs = normalize(sign(X).*log1p(abs(X))); Xs(~isfinite(Xs)) = 0; a = zeros(1,3);
for s = 1:3
    rng(s); cvm = fitclinear(Xs,y,'Learner','logistic','Regularization','lasso','Lambda',0.004,'Prior','uniform','KFold',5);
    [~,sc] = kfoldPredict(cvm); [~,~,~,a(s)] = perfcurve(y,sc(:,2),true);
end
a = mean(a);
end
