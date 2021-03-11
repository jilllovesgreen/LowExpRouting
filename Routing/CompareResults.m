k = 1609.3; % meter to mile conversion factor
ta = 'hr22'; %22pm, 10am
MYa ='12';

tb =  ta;
MYb = MYa;

inputWeightA=['tier1_d_tier2_t_' ta];  %'';  % compare t, dist, IM, weight_t_10am
inputWeightB=['tier1_d_tier2_inhale_MY12_' tb];  % tier_routing_t_hr10,  tier_routing_inhale_MY12_hr10

cTableA = readtable(['cost_table_' inputWeightA '.csv']);
cTableB = readtable(['cost_table_' inputWeightB '.csv']);

%----
% pm_change= (cTableB.(['PM25_' MYb '_' tb])-cTableA.(['PM25_' MYa '_' ta]))./cTableA.(['PM25_' MYa '_' ta])*100;
% d_change= (cTableB.d - cTableA.d)./cTableA.d*100;
% %t_change=(cTableB.(['t_' tb])-cTableA.(['t_' ta]))./cTableA.(['t_' ta])*100;
% % pm = pm_change;
% % d = d_change;
% % t = t_change;  %
% % disp('[min(pm), mean(pm), max(pm), min(d), mean(d), max(d), min(t), mean(t), max(t)]')
% % [min(pm), mean(pm), max(pm), min(d), mean(d), max(d), min(t), mean(t), max(t)]
% 
% cTableB(pm_change > 0,:) = cTableA(pm_change > 0,:); %caution, use to
% %remove pm_change > 0 cases
% cTableB(d_change < 0,:) = cTableA(d_change < 0,:); % caution 

%----
pm_change= (cTableB.(['PM25_' MYb '_' tb])-cTableA.(['PM25_' MYa '_' ta]))./cTableA.(['PM25_' MYa '_' ta])*100;
d_change= (cTableB.d - cTableA.d)./cTableA.d*100;
t_change=(cTableB.(['t_' tb])-cTableA.(['t_' ta]))./cTableA.(['t_' ta])*100;
nox_change = (cTableB.(['NOx_' MYb '_' tb])-cTableA.(['NOx_' MYa '_' ta]))./cTableA.(['NOx_' MYa '_' ta])*100;
% original units: d meter, t second, PM2.5 ug, NOx ug, CO2 g
% updated units:   mile,    min,       ug,         ug,    kg
metricsA = [cTableA.d/k, cTableA.(['t_' ta])/60, ...
    cTableA.(['PM25_' MYa '_' ta]), cTableA.(['NOx_' MYa '_' ta]), cTableA.(['CO2_' MYa '_' ta])/1000];
    
metricsB = [cTableB.d/k, cTableB.(['t_' tb])/60, ...
    cTableB.(['PM25_' MYb '_' tb]), cTableB.(['NOx_' MYb '_' tb]), cTableB.(['CO2_' MYb '_' tb])/1000];

%a(1)=a(1)/1609.3; a(2) = a(2)/60; a(3:4) = a(3:4)/1000;a(5) = a(5)/1000;
metrics_change = (metricsB - metricsA)./metricsA*100;
totalMetrics=[metricsA metricsB metrics_change];

sumNumTrips=[];

x_axis = t_change; % x axis !!!! 

y_axis = pm_change(pm_change~=0);
x_axis = x_axis(x_axis ~=0);
%brk1=prctile(pm_change(pm_change~=0),0:10:100);   % y axis,  pm change
%brk2=prctile(metric_change(metric_change ~=0),0:10:100);  % x axis
brk1 = -100:10:0;  % y axis break points, usually pm/nox change
%brk2 = [min(metric_change) -30 0 2 5 7 10 20 40 50 max(metric_change)];
%brk2 = min(x_axis): (max(x_axis)-min(x_axis))/10 : max(x_axis);
%brk2 = prctile(x_axis, 0:10:100);

brk2 = [0 2	5	10	15	20	25	30	50	60	max(x_axis)];  % x axis break points, ususally time/dist increase

for i1 = 1:(length(brk1)-1)
    limit1 = y_axis>brk1(i1)& y_axis<=brk1(i1+1); % y axis
    for i2=1:(length(brk2)-1)
        limit2 = x_axis>brk2(i2)& x_axis<=brk2(i2+1); % x axis  
        sumNumTrips(i1,i2) = sum(limit1&limit2);
        if i1==1 || i2 ==1
        limit1 = y_axis>=brk1(i1)& y_axis<=brk1(i1+1);
        limit2 = x_axis >= brk2(i2)& x_axis <= brk2(i2+1);
        sumNumTrips(i1,i2)=sum(limit1&limit2);  % zero increase already removed
        continue
        end
    end
end
sumNumTripsMat=flip(sumNumTrips)/sum(sumNumTrips(:))*100;

test= metricsA;
x = t_change;

xA= [sum(x == 0) mean(test(x == 0, :), 1); ...
    sum(x >0 & x <=10) mean(test(x >0 & x <=10, :), 1);...
    sum(x >10 & x<=30)  mean(test(x >10 & x<=30, :), 1); ...
    sum(x > 30) mean(test(x > 30,:), 1); ...
    length(x)  mean(test, 1)];
v_mph = xA(:,2)./xA(:,3)*60;  % mile/min  to mile/hour
xA_with_speed = [xA(:,1:3), v_mph,  xA(:,4:end)];

% xA= [ sum(d<= -30) mean(test(d<= -30, :)); ...
%     sum(d> -30 & d <0) mean(test(d> -30 & d <0, :)); ...
%     sum(d >0 & d <=10) mean(test(d >0 & d <=10, :));...
%     sum(d >10 & d<=30)  mean(test(d >10 & d<=30, :)); ...
%     sum(d > 30) mean(test(d > 30,:))];
% 
% xA2 = [xA; sum(d~=0) mean(test(d~=0, :))];
% v_mph = xA2(:,2)./xA2(:,3)/1609.3*3600;  % meter/second/1069.3*3600
% xA_with_speed = [xA2(:,1:3), v_mph,  xA2(:,4:end)];

test = metricsB;
xB= [ sum(x == 0) mean(test(x == 0, :), 1); ...
    sum(x >0 & x <=10) mean(test(x >0 & x <=10, :), 1);...
    sum(x >10 & x<=30)  mean(test(x >10 & x<=30, :), 1); ...
    sum(x > 30) mean(test(x > 30,:), 1);...
    length(x)  mean(test, 1)]; % last row is the overall average

v_mph = xB(:,2)./xB(:,3)*60;  % mile/min*60
xB_with_speed = [xB(:,1:3), v_mph, xB(:,4:end)]; % insert the mph column to the table

compMatrix = [xA_with_speed  xB_with_speed(:,2:end)];
