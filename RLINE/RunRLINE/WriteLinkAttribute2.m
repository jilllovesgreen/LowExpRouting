link = readtable('LinkAttribute.csv');
X_begin = link.X_begin;
Y_begin = link.Y_begin;
Z_begin = link.Z_begin;
Z_begin(Z_begin == -9999) = 0;
X_end = link.X_end;
Y_end = link.Y_end;
Z_end = link.Z_end;
Z_end(Z_end == -9999) = 0;
KPH =link.KPH;
METERS = sqrt((X_begin - X_end).^2 + (Y_begin-Y_end).^2 + (Z_begin - Z_end).^2);%zeros(len/2,1); 
NetworkID = link.NetworkID;
len = length(NetworkID);
hist_rel_10 = zeros(len,1);
hist_rel_22 = zeros(len,1);

spd_10 = readtable('..\..\Routing\HistSpeed\Monday_10_speed.csv'); % optional, the historical speed table
spd_22 = readtable('..\..\Routing\HistSpeed\Monday_22_speed.csv'); 

for i1=1:len
    id_i = NetworkID(i1); 
    if ~any(spd_10.network_id == id_i)
        hist_rel_10(i1) = 100; % if the link has no value in hist db, use 100% of KPH value
        hist_rel_22(i1) = 100;
        continue
    end
hist_rel_10(i1) = mean(spd_10.monday_10_speed_kph(spd_10.network_id == id_i));
hist_rel_22(i1) = mean(spd_22.monday_22_speed_kph(spd_22.network_id == id_i));
end
clear spd_10 spd_22

link = removevars(link, 'hist_speed_kph');
link = removevars(link, 't_sec' );
link.METERS = METERS;
link = addvars(link, hist_rel_10 );
link = addvars(link, hist_rel_22 );

%=== add freeflow speed metrics from hist speed data
if ~issorted(NetworkID)
   error('NetworkID is not sorted') 
end
FFS_table = readtable('..\..\Routing\HistSpeed\KPH_vs_FFS.csv');
FFS_table = sortrows(FFS_table);
FFS = FFS_table.SPFREEFLOW;
FFS(FFS==0) = link.KPH(FFS==0);  % replace 0 as KPH
FFS(isnan(FFS)) = link.KPH(isnan(FFS)); % replace nan as FFS
% for i1 = 1: len  % go through NetworkID in table link
%     id_i = NetworkID(i1); 
%     if ~any(FFS_NID == id_i) || isnan(mean(FFS_ffs(FFS_NID == id_i))) || mean(FFS_ffs(FFS_NID == id_i)) == 0
%         FFS(i1) = link.KPH(i1);
%         continue
%     end
%     FFS(i1) = mean(FFS_ffs(FFS_NID == id_i));
% end
link = addvars(link, FFS );
 writetable(link, ['LinkAttribute_' loc '.xlsx']);
