% RUNROUTING Run the routing based on given OD pairs and required weight constrain
% The routing will be devided into two tiers 
% Tier 1: truck routes and connectors, find the route where the combined cost are lowest
% Tier 2: all other roads, find the fastest one
addpath('./weights')
load('storePOI_JID'); % the linkID of the stores, adjust as needed
load('borderPOI_JID'); % the linkID of the entry/exit points on the border, adjust as needed. The OD pair are actually link pairs, not node pairs
alpha = [0.25, 0.25, 0.5]; % the weight for NOx, PM2.5 and time/dist/truckroute, adjust as needed
inputWeight_tier1 ='d';% d, t_hr10, d_TR,  inhale_MY12_hr10,
inputWeight_tier2 = 't_hr22'; %'inhale_MY12_hr10'; %inhale_MY12_hr10

route_orig = zeros(length(borderPOI)*length(storePOI),1); % allocate orig vector
route_dest = route_orig;  % allocate destination vector
cnt = 0;                          
for j1 = 1:length(borderPOI)  % allocate origin to destination list
    for k1 = 1:length(storePOI)
       cnt = cnt+1;
       route_orig(cnt) = borderPOI(j1);
       route_dest(cnt) = storePOI(k1);
    end
end
% if there were 8 origins, 10 destinations, we will run all origins to
% destinations and all destinations to origins
% total 8x10x2 = 160 trips
route_start = [route_orig; route_dest]; 
route_end = [route_dest; route_orig];
m = length(route_start);
costTable_tier1 = zeros(m,17);

% all table headers of the following table must match the example tables
T = readtable('multinet_LB.xlsx');  % network table
Tw1 = readtable('./weights/length_duration_weights_unique_id.csv'); % attributes of length/duration/weighted truck route length(weighted_tr_len)
Tw2 = readtable('./weights/inhaled_mass_weights.csv'); % calculated inhaled mass ug/link
Tw3 = readtable('./weights/CO2_weights.csv'); % CO2 mass emission g/link

if ~ all(Tw1.network_id == Tw2.network_id) || ~ all(Tw2.network_id == Tw3.network_id)
   error('Network ID do not match') 
end

% prepare the graph based on input weights
% ==== Tier 1, find the section of route on non-truck roads
% ==== for first half of trips start from border points to store, break
% fastest route to start -> break1 (which we run in tier 2 routing); break2 -> end (which we keep);
% ==== for second half of trips start from store to border points, break
% fastest route to start -> break1 (which we keep); break2 -> end (run in tier 2 routing);

[G, node_map, node1,node2,edgeList, weight_for_routing] = PrepareGraph(inputWeight_tier1,T,Tw1,Tw2,Tw3,alpha);
recordPath_tier1 = cell(m,2);
break1 = zeros(m,1);
tr_percent = zeros(m,1);
%check_case_tier1 = zeros(m, 1);
for i1 =  1:m
    s_JID = route_start(i1); % start link
    e_JID = route_end(i1);   % end link
    half_label = i1 > m/2;  % half_label is 0 if it is the first half of OD pairs: entry points to store; 1 if it the second half: stores to entry points 
    [sum_cost_tier1, path, break1(i1), tr_percent(i1), path_link_keep] = TierRouting(s_JID, e_JID, G, node_map, T, Tw1, Tw2, Tw3,node1,node2,edgeList, weight_for_routing, half_label, i1); 
    costTable_tier1(i1,:) = [s_JID, e_JID, sum_cost_tier1];
    recordPath_tier1{i1, 1} = length(path_link_keep);
    recordPath_tier1{i1, 2} = path_link_keep;
end

%==== Tier 2 
% first half, start -> break1 (which we run in tier 2 routing); break2 -> end (which we keep from tier 1);
% second half,  start -> break1 (which we keep from tier1); break2 -> end (run in tier 2 routing);
% inputWeight_tier2 = 'inhale_MY12_hr22'; %'inhale_MY12_hr10'; %inhale_MY12_hr10
Tw1_tier2 = Tw1(Tw1.tr_label2 ==1, :);
Tw2_tier2 = Tw2(Tw1.tr_label2 ==1, :);
Tw3_tier2 = Tw3(Tw1.tr_label2 ==1, :);
costTable_tier2 = zeros(m, 17);
recordPath_tier2 = cell(m, 2);
check_case_tier2 = zeros(m, 1);
[G, node_map, node1,node2,edgeList, weight_for_routing] = PrepareGraph(inputWeight_tier2,T,Tw1_tier2,Tw2_tier2,Tw3_tier2,alpha);

for i1 = 1:m
    half_label = i1 > m/2;
    if half_label == 0 % first half, start -> break1 (which we run in tier 2 routing); break2 -> end (which we keep in last run);
        s_JID = route_start(i1); % start link
        e_JID = break1(i1);   % end link
    else  % second half,  start -> break1 (which we keep); break2 -> end (run in tier 2 routing);
        s_JID = break1(i1);
        e_JID = route_end(i1);
    end
    if s_JID == 0 || e_JID == 0
       error('Start or end link is wrong') 
    end
      % half_label is 0 if it is the first half of OD pairs: entry points to store; 1 if it the second half: stores to entry points 
    [sum_cost, path] = routing(s_JID, e_JID, G, node_map, T, Tw1_tier2, Tw2_tier2, Tw3_tier2, node1,node2,edgeList, weight_for_routing); 
    costTable_tier2(i1,:) = [s_JID, e_JID, sum_cost];
    recordPath_tier2{i1, 1} = length(path);
    recordPath_tier2{i1, 2} = path;
end

%====
s_JID=route_start; t_link=route_end;
costTable = costTable_tier1 + costTable_tier2;
d=costTable(:,3); t_hr10=costTable(:,4); t_hr22=costTable(:,5);
NOx_07_hr10=costTable(:,6); NOx_07_hr22=costTable(:,7); 
NOx_12_hr10=costTable(:,8); NOx_12_hr22=costTable(:,9);
PM25_07_hr10=costTable(:,10); PM25_07_hr22=costTable(:,11);
PM25_12_hr10=costTable(:,12); PM25_12_hr22=costTable(:,13);
CO2_07_hr10=costTable(:,14);  CO2_07_hr22=costTable(:,15);
CO2_12_hr10=costTable(:,16);  CO2_12_hr22=costTable(:,17);

results=table(s_JID, t_link, d, t_hr10, t_hr22,...   % weight table A. e.g., time, duration
        PM25_07_hr10, PM25_07_hr22,PM25_12_hr10, PM25_12_hr22,...    
        NOx_07_hr10, NOx_07_hr22, NOx_12_hr10, NOx_12_hr22,...
        CO2_07_hr10, CO2_07_hr22,CO2_12_hr10, CO2_12_hr22);  

writetable(results, ['cost_table_tier1_' inputWeight_tier1 '_tier2_' inputWeight_tier2 '.csv']);
%save(['cost_table_' inputWeight ], 'costTable') %'_tr_weightv2'

num_path = cell2mat(recordPath_tier1(:,1)) + cell2mat(recordPath_tier2(:,1));
R = zeros( max(num_path) , m); %path of the route in the form of every link network_id
for i2 = 1:m
   R(1: num_path(i2), i2)= [recordPath_tier1{i2,2} ; recordPath_tier2{i2,2}];
end

pathTable= table(R);
writetable(pathTable, ['path_table_tier1_' inputWeight_tier1 '_tier2_' inputWeight_tier2 '.csv' ]);
