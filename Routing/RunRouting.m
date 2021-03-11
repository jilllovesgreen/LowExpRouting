%RUNROUTING Run the routing based on given OD pairs and required weight constrain
%d: routing based on edge distance
%d_TR: routing based on truck route label and edge distance, a truck will try all truck routes first and then non-truck routes to get to destination   
%t_10am, routing based on Tuesday 10 am historical speed by TomTom Multinet  
%t_22pm, routing based on Tuesday 10 pm historical speed by TomTom Multinet  
%inhale_MY07_10am, inhaled mass from a Model Year 2007 truck at 10 am,
%inhale_MY07_22pm, inhaled mass from a Model Year 2007 truck at 10 pm
%inhale_MY12_10am, inhaled mass from a Model Year 2012 truck at 10 am
%inhale_MY12_22pm, inhaled mass from a Model Year 2012 truck at 10 pm
% for all inhalation weight, check PrepareGraph Line 99 to adjust the composite weight with
% time/distance/truck distance
addpath('./weights')
load('storePOI_JID'); % the linkID of the stores, adjust as needed
load('borderPOI_JID'); % the linkID of the entry/exit points on the border, adjust as needed. The OD pair are actually link pairs, not node pairs
inputWeight ='d_TR';% d, t_hr10, d_TR,  inhale_MY12_hr10,
alpha = [0.25, 0.25, 0.5]; % the weight for NOx, PM2.5 and time/dist/truckroute, adjust as needed

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
costTable = zeros(m,17);

% all table headers of the following table must match the example tables
T = readtable('multinet_LB.xlsx');  % network table
Tw1 = readtable('./weights/length_duration_weights_unique_id.csv'); % attributes of length/duration/weighted truck route length(weighted_tr_len)
Tw2 = readtable('./weights/inhaled_mass_weights.csv'); % calculated inhaled mass ug/link
Tw3 = readtable('./weights/CO2_weights.csv'); % CO2 mass emission g/link

% prepare the graph based on input weights
[G, node_map, node1,node2,edgeList, weight_for_routing] = PrepareGraph(inputWeight,T,Tw1,Tw2,Tw3,alpha);

recordPath = cell(m,2);
% break1 = zeros(m,1);
% break2 = zeros(m,1);
% tr_percent = zeros(m,1);
%tr_p_entry_stretch = zeros(m,1); % the longest stretch connected with entry/exit points
for i1 = 101%1:m
    s_link = route_start(i1); % start link
    e_link = route_end(i1);   % end link
    half_label = i1 > m/2;  % half_label is 0 if it is the first half of OD pairs: entry points to store; 1 if it the second half: stores to entry points 
    [sum_cost, path] = routing(s_link, e_link, G,node_map, T, Tw1,Tw2, Tw3,node1,node2,edgeList, weight_for_routing); 
    costTable(i1,:) = [ s_link, e_link, sum_cost];
    recordPath{i1, 1}= length(path);
    recordPath{i1, 2}= path;
end
% T_out = [path_weight_d, path_weight_t_10am, path_weight_t_22pm,...
%     path_weight_NOx_07_10am, path_weight_NOx_07_22pm, path_weight_NOx_12_10am,...
%     path_weight_NOx_12_22pm, path_weight_PM25_07_10am, path_weight_PM25_07_22pm, path_weight_PM25_12_10am,...
%     path_weight_PM25_12_22pm,...
%     path_weight_CO2_07_10am, path_weight_CO2_07_22pm, ...
%     path_weight_CO2_12_10am, path_weight_CO2_12_22pm]; 

s_link=costTable(:,1); t_link=costTable(:,2);
d=costTable(:,3); t_hr10=costTable(:,4); t_hr22=costTable(:,5);
NOx_07_hr10=costTable(:,6); NOx_07_hr22=costTable(:,7); 
NOx_12_hr10=costTable(:,8); NOx_12_hr22=costTable(:,9);
PM25_07_hr10=costTable(:,10); PM25_07_hr22=costTable(:,11);
PM25_12_hr10=costTable(:,12); PM25_12_hr22=costTable(:,13);
CO2_07_hr10=costTable(:,14);  CO2_07_hr22=costTable(:,15);
CO2_12_hr10=costTable(:,16);  CO2_12_hr22=costTable(:,17);

results=table(s_link, t_link, d, t_hr10, t_hr22,...   % weight table A. e.g., time, duration
        PM25_07_hr10, PM25_07_hr22,PM25_12_hr10, PM25_12_hr22,...    
        NOx_07_hr10, NOx_07_hr22, NOx_12_hr10, NOx_12_hr22,...
        CO2_07_hr10, CO2_07_hr22,CO2_12_hr10, CO2_12_hr22);  

writetable(results, ['cost_table_' inputWeight '.csv']);
%save(['cost_table_' inputWeight ], 'costTable') %'_tr_weightv2'

R = zeros( max(cell2mat(recordPath(:,1))), m); %path of the route in the form of every link network_id
for i2 = 1:m
   R(1:recordPath{i2,1}, i2)= recordPath{i2,2}';
end

pathTable= table(R);
writetable(pathTable, ['path_table_' inputWeight '.csv' ]);
