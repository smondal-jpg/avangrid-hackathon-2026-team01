% Load AMI csv files, weather and customer lists into ami.mat

root = 'E:\Avangrid Hackathon 2026\12 Nov';
out  = 'E:\Avangrid Hackathon 2026\matlab\ami.mat';

nC = 35725; nH = 8760; t0 = datenum(2025,9,1);
kwh = nan(nC,nH,'single');

files = dir(fullfile(root,'2026 Avangrid Hackathon Data_*','*.csv'));
fmt = '%f %s %f %*s %f %*s %*s %*s';      % customer, date, hour_ending, kWh
tic
for i = 1:numel(files)
    fid = fopen(fullfile(files(i).folder,files(i).name));
    c = textscan(fid,fmt,'Delimiter',',','HeaderLines',1);
    fclose(fid);
    d = char(c{2});
    dn = datenum((d(:,1:4)-'0')*[1000;100;10;1],(d(:,6:7)-'0')*[10;1],(d(:,9:10)-'0')*[10;1]);
    h = (dn-t0)*24 + c{3};
    ok = h>=1 & h<=nH & c{1}>=1 & c{1}<=nC;
    kwh(sub2ind([nC nH],c{1}(ok),h(ok))) = c{4}(ok);
    fprintf('%3d/%d  %s  rows %d  %.0f s\n',i,numel(files),files(i).name(1:10),numel(h),toc);
end

w = readtable(fullfile(root,'Weather Data_Avangrid_Hackathon .xlsx'));
[~,iu] = unique(string(w.date)+string(w.time),'stable');   % drop DST duplicate hour
w = w(iu(1:nH),:);
temp  = fillmissing(w.temp,'linear');                       % fill the 94-h gap in March
humid = fillmissing(w.humidity,'linear');
wind  = fillmissing(w.windspeed,'linear');

lf = fullfile(root,'EV_Heat Pump_Solar PV Customers.xlsx');
evID = readmatrix(lf,'Sheet','EV Customers');
hpID = readmatrix(lf,'Sheet','HP Customers');
pvID = readmatrix(lf,'Sheet','PV Customers');
evID = evID(~isnan(evID)); hpID = hpID(~isnan(hpID)); pvID = pvID(~isnan(pvID));

t = datetime(2025,9,1) + hours(0:nH-1)';
if ~isfolder(fileparts(out)), mkdir(fileparts(out)); end
save(out,'kwh','temp','humid','wind','t','evID','hpID','pvID','-v7.3');
fprintf('Saved %s | missing %.2f%% | EV %d HP %d PV %d\n',out,100*mean(isnan(kwh(:))),numel(evID),numel(hpID),numel(pvID));
