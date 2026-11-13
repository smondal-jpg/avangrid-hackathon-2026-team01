% Add the corrected csv files (10001-20000, parts 44-49) to ami.mat

wd  = 'E:\Avangrid Hackathon 2026\matlab';
src = 'E:\Avangrid Hackathon 2026\12 Nov';

load(fullfile(wd,'ami.mat'));
[nC,nH] = size(kwh); t0 = datenum(2025,9,1);
files = [dir(fullfile(src,'2026 Avangrid Hackathon Data_10001_20000_v2','*.csv')); ...
         dir(fullfile(src,'2026 Avangrid Hackathon Data_Parts 44 to 49','*.csv'))];
before = sum(mean(isnan(kwh),2)<0.5);
hdr = ['CUSTOMER_NUMBER,usage_date,hour_ending,day_type,hourly_kwh,interval_count,interval_type,is_complete_hour' newline];
tic
for i = 1:numel(files)
    s = fileread(fullfile(files(i).folder,files(i).name));
    s(s==char(13)) = [];
    s = strrep(s,hdr,'');
    s = strrep(s,'E-','E#'); s = strrep(s,'-',','); s = strrep(s,'E#','E-');
    s = strrep(s,',Weekday',''); s = strrep(s,',Weekend',''); s = strrep(s,',15,min','');
    s(s==newline) = ',';
    v = reshape(sscanf(s,'%f,'),8,[]);
    h = (datenum(v(2,:),v(3,:),v(4,:))-t0)*24 + v(5,:);
    ok = h>=1 & h<=nH & v(1,:)>=1 & v(1,:)<=nC;
    kwh(sub2ind([nC nH],v(1,ok),h(ok))) = v(6,ok);
    fprintf('%3d/%d  customers %d-%d  %.0f s\n',i,numel(files),min(v(1,:)),max(v(1,:)),toc);
end
save(fullfile(wd,'ami.mat'),'kwh','temp','humid','wind','t','evID','hpID','pvID','-v7.3','-nocompression');
fprintf('Updated ami.mat | customers with data: %d (before: %d)\n',sum(mean(isnan(kwh),2)<0.5),before);
