%Plot routes for each OD pair
% plot shortest routes based on truck routes, (shortest routes optional) 
% inhaled mass MY12 10am, inhaled mass MY07 10am, inhaled mass MY12 10 pm
ID = 1:1100;
%OD = ID(pm ==min(pm(t_change ==0)));  
OD = trips_check;
%OD  =  ; %[586, 934, 1059, 1079];% =id2(41:end);
linkAttrib = readtable('..\RLINE\RunRLINE_Script\LinkAttribute_LB.xlsx');
network_id = linkAttrib.NetworkID;
X1 = linkAttrib.X_begin;
Y1 = linkAttrib.Y_begin;
X2 = linkAttrib.X_end;
Y2 = linkAttrib.Y_end;

xmin = min([X1 ; X2]);
xmax = max([X1 ; X2]);
ymin = min([Y1 ; Y2]);
ymax = max([Y1 ; Y2]);

t = readtable('path_table_tier1_d_tier2_t_hr10.csv');  % fastest on truck route network, aka, tiered fastest
d = readtable('path_table_d.csv');  % shortest path
%d_TR = readtable('path_table_d_TR.csv'); % shortest on truck route network, aka, tiered shortest
im_my12_hr10 = readtable('path_table_tier1_d_tier2_inhale_MY12_hr10.csv'); % tiered LER 10 AM
im_my12_hr22 = readtable('path_table_tier1_d_tier2_inhale_MY12_hr22.csv');  % tiered LER 10 PM

for i1 = 1:length(OD)
    
   ODi = OD(i1);
   path4 = im_my12_hr10.(['R_' num2str(ODi) ]);
   path5 = im_my12_hr22.(['R_' num2str(ODi) ]); 
   
   if all(path4 == path5)
      continue 
   end
   
   
   figure, sgtitle(['Trip ' num2str(ODi)])
   
   pos = 1;% figure position
   path1 = t.(['R_' num2str(ODi) ]);
   path1 = path1(path1>0);
   subplot(2,2,pos);
   PlotPath(path1,network_id, X1, Y1, X2, Y2);
   title('Tiered Fastest Route')
   axis([xmin, xmax, ymin, ymax])
   axis equal
   
   pos = 2;   
   path2 = d.(['R_' num2str(ODi) ]);
   path2 = path2(path2>0);
   subplot(2,2,pos);
   PlotPath(path2, network_id, X1, Y1, X2, Y2);
   title('Shortest Route')
   axis([xmin, xmax, ymin, ymax])
   axis equal

%    pos = 2;
%    path3 = d_TR.(['R_' num2str(ODi) ]);
%    path3 = path3(path3>0); % there are trailing zeros 
%    subplot(2, 2, pos);
%    PlotPath(path3, network_id, X1, Y1, X2, Y2);
%    title('Shortest Route based on Truck Routes')
%    %axis([xmin, xmax, ymin, ymax])
%    axis equal
   
   pos = 3; 
   %path4 = im_my12_hr10.(['R_' num2str(ODi) ]);
   path4 = path4(path4>0);
   subplot(2, 2, pos);
   PlotPath(path4,network_id, X1, Y1, X2, Y2);
   axis([xmin, xmax, ymin, ymax])
   title('Tiered LER 10 AM')
   axis equal
   
   pos = 4 ;
   %path5 = im_my12_hr22.(['R_' num2str(ODi) ]);
   path5 = path5(path5>0);
   subplot(2,2, pos);
   PlotPath(path5,network_id, X1, Y1, X2, Y2);
   title('Tiered LER 10 PM')
   axis([xmin, xmax, ymin, ymax])
   axis equal
   
  
end
