function [G, node_map,node1,node2,edge, weight_for_routing]  = PrepareGraph( inputWeight,T,Tw1,Tw2,Tw3, alpha)
%PREPAREGRAPH Calcaulate Graphs for multiple routing objectives
truck_links = Tw1.network_id;

kk = T.IDtxt == truck_links(1);

for ii = 2: length(truck_links)
    kk = kk | T.IDtxt == truck_links(ii);
end

T = T(kk,:);

node_map = containers.Map('KeyType','double','ValueType','int64');
node_map_reverse = containers.Map('KeyType','int64','ValueType','double');

%edge_map = containers.Map('KeyType','int64','ValueType','double');

kk = strcmp(T.ONEWAY, 'N');

T = T(~kk,:);
node_list = unique([T.F_JID; T.T_JID]);
n = 0;
for ii = 1:length(node_list)    
    n = n + 1;
    node_map(node_list(ii)) = n;
    node_map_reverse(n) = node_list(ii);
end

edge = T.IDtxt;
m0 = length(edge);

node1 = int64(zeros(m0,1));
node2 = int64(zeros(m0,1));
m = m0;
for ii = 1:m0
    if strcmp(T.ONEWAY(ii),'FT')
        node1(ii) = node_map(T.F_JID(ii));
        node2(ii) = node_map(T.T_JID(ii));
    else
        if strcmp(T.ONEWAY(ii), 'TF')
            node1(ii) = node_map(T.T_JID(ii));
            node2(ii) = node_map(T.F_JID(ii));
        
        else
            if strcmp(T.ONEWAY(ii), '')
                node1(ii) = node_map(T.T_JID(ii));
                node2(ii) = node_map(T.F_JID(ii));
            
                %add edge
                edge = [edge; edge(ii)];
                node1 = [node1; node2(ii)];
                node2 = [node2; node1(ii)];
                m = m + 1;
            end
        end
    end

end


%% genrate weights
%x1 = table(node1, node2, edge);
%writetable(x1, 'node_link.xlsx');
weight_t_hr10 = zeros(m,1); % duration 10:00 am
weight_t_hr22 = zeros(m,1); % duration 22:00 pm
weight_d = zeros(m,1);      % distance
weight_d_tr_w=zeros(m,1);       % distance weighted based on truck routes length
weight_CO2_07_hr10 = zeros(m,1);
weight_CO2_07_hr22 = zeros(m,1);
weight_CO2_12_hr10 = zeros(m,1);
weight_CO2_12_hr22 = zeros(m,1);

weight_inhale_07_hr10 = zeros(m,1);
weight_inhale_07_hr22 = zeros(m,1);
weight_inhale_12_hr10 = zeros(m,1);
weight_inhale_12_hr22 = zeros(m,1);

NOx_07_hr10_norm = my_norm(Tw2.NOx_07_hr10_ug);
NOx_07_hr22_norm = my_norm(Tw2.NOx_07_hr22_ug);
NOx_12_hr10_norm = my_norm(Tw2.NOx_12_hr10_ug);
NOx_12_hr22_norm = my_norm(Tw2.NOx_12_hr22_ug);
PM25_07_hr10_norm = my_norm(Tw2.PM25_07_hr10_ug);
PM25_07_hr22_norm = my_norm(Tw2.PM25_07_hr22_ug);
PM25_12_hr10_norm = my_norm(Tw2.PM25_12_hr10_ug);
PM25_12_hr22_norm = my_norm(Tw2.PM25_12_hr22_ug);

duration_hr10_norm = my_norm(Tw1.duration_hr10_s);
duration_hr22_norm = my_norm(Tw1.duration_hr22_s);

d_tr_w_norm = my_norm(Tw1.weighted_tr_len); %!!! caution!

for ii = 1:m
    kk = Tw1.network_id == edge(ii);
    if sum(kk) > 0
        weight_t_hr10(ii) = Tw1.duration_hr10_s(kk);
        weight_t_hr22(ii) = Tw1.duration_hr22_s(kk);
        weight_d(ii) = Tw1.length_m(kk); 
        weight_d_tr_w(ii)=Tw1.weighted_tr_len(kk);
        
        weight_inhale_07_hr10(ii) = nansum(alpha .*[NOx_07_hr10_norm(kk), PM25_07_hr10_norm(kk), duration_hr10_norm(kk)]); % duration_10am_norm(kk)
        weight_inhale_07_hr22(ii) = nansum(alpha .*[NOx_07_hr22_norm(kk), PM25_07_hr22_norm(kk), duration_hr22_norm(kk)]); % duration_22pm_norm(kk)
        weight_inhale_12_hr10(ii) = nansum(alpha .*[NOx_12_hr10_norm(kk), PM25_12_hr10_norm(kk), duration_hr10_norm(kk)]); % duration_10am_norm(kk)
        weight_inhale_12_hr22(ii) = nansum(alpha .*[NOx_12_hr22_norm(kk), PM25_12_hr22_norm(kk), duration_hr22_norm(kk)]); %duration_22pm_norm(kk)
        
        weight_CO2_07_hr10(ii) = Tw3.CO2_07_hr10_g(kk);
        weight_CO2_07_hr22(ii) = Tw3.CO2_07_hr22_g(kk);
        weight_CO2_12_hr10(ii) = Tw3.CO2_12_hr10_g(kk);
        weight_CO2_12_hr22(ii) = Tw3.CO2_12_hr22_g(kk);
    end
end

%% set start-edge and end-edge

%% generate directed graph with different weights
%G = digraph(node1, node2, inputWeight);

weight_for_routing = weight_d; % defalt
 
switch inputWeight
    case  'd'
        weight_for_routing = weight_d;
        %G = digraph(node1, node2, weight_d);  
    case  'd_TR' % truck route weight
        weight_for_routing=weight_d_tr_w;
    case  't_hr10'
        weight_for_routing = weight_t_hr10;
        %G = digraph(node1, node2, weight_t_10am);
    case  't_hr22'
        weight_for_routing = weight_t_hr22;
        %G = digraph(node1, node2, weight_t_22pm);
    case 'inhale_MY07_hr10' 
        weight_for_routing = weight_inhale_07_hr10;
        %G = digraph(node1, node2, weight_inhale_07_10am); %!!!
    case 'inhale_MY07_hr22'
        weight_for_routing = weight_inhale_07_hr22;
        %G = digraph(node1, node2, weight_inhale_07_22pm);
    case 'inhale_MY12_hr10'
        weight_for_routing = weight_inhale_12_hr10;
        %G = digraph(node1, node2, weight_inhale_12_10am);
    case 'inhale_MY12_hr22'
        weight_for_routing = weight_inhale_12_hr22;
        %G = digraph(node1, node2, weight_inhale_12_22pm);
end

G = digraph(node1, node2, weight_for_routing);

%% routing and output edge path
function y = my_norm(x)  
    y = x /max(x); % caution
end

end

