% 
% T = readtable('../streets_LB_TomTom_ID_v2.csv');
% Tw1 = readtable('../length_duration_weights.csv');
% Tw2 = readtable('../inhaled_mass_weights.csv');
% Tw3 = readtable('../CO2_weights.csv');

%% generate directed graph

% T_network = readtable('streets_LB_TomTom_UTM_gen20m_split_truck_roads.csv');

function [sumTable, path_link_id]=routing(s_link, t_link, G,node_map, T, Tw1,Tw2, Tw3,node1,node2,edge, weight_for_routing)

%disp(['Start link: ' num2str(s_link) ' Dest link: ' num2str(t_link)])
      
s_node = node_map(s_link); %node_map(T.F_JID(T.IDtxt==s_link));
t_node = node_map(t_link); %node_map(T.T_JID(T.IDtxt==t_link));
%plot(G);
%[nodepath, path_weight, edgepath] = shortestpath(G,s_node,t_node); %!!
[~, ~, edgepath] = shortestpath(G,s_node,t_node); %!!
e_path = G.Edges(edgepath,:);

if isempty(e_path)
    sumTable = zeros(1, 15);
    path_link_id = [];
    disp(['Empty path! Start and end link: ' num2str(s_link) ' ' num2str(t_link) ]) 
    return
end
% edgepath is the nominal edge path, not the correct network ID

path_link_id = zeros(size(e_path,1),1); % the true network ID of a path, use this to map trip on ArcMap

for i1 = 1:size(e_path,1)
    % match back the edge to the network_id
    kk = node1 == e_path.EndNodes(i1,1) & node2 == e_path.EndNodes(i1,2) & weight_for_routing == e_path.Weight(i1);
    link_kk = edge(kk);
    if sum(kk)> 1
        error('Link number error')
    end
    path_link_id(i1) = link_kk(1);
end

% if isempty(path_link_id) % if path is empty
%     disp(['Empty path! Start and end link: ' num2str(s_link) ' ' num2str(t_link) ]) 
%     path_link_id = [s_link; t_link];
%     check_case = 6;
% elseif s_link ~= path_link_id(1) && t_link ~= path_link_id(end)  % if neither path start nor end 
%     path_link_id = [s_link; path_link_id; t_link];
%     check_case = 1;
%     %disp('Case 1')
% elseif s_link == path_link_id(1) && t_link ~= path_link_id(end)  % if start equals but end not 
%     path_link_id = [path_link_id; t_link];
%     check_case = 2;
%     %disp('Case 2')
% %     disp(['Other cases: start link ' num2str(s_link)  ' and path link starts with '  num2str(path_link_id(1)) ])
% %     disp(['end link ' num2str(t_link) ' and path link ends with ' num2str(path_link_id(end))])
% elseif s_link ~= path_link_id(1) && t_link == path_link_id(end)  % if start does not equal but end equals
%     path_link_id = [s_link; path_link_id];
%     check_case = 3;
%     %disp('Case 3')
% else % path start and end both equal to start and end point 
%     check_case = 4;
% end
% 
% if s_link == t_link
%     disp(['Start and end links are the same: check case was ' num2str(check_case)])
%     disp(num2str(path_link_id))
%     disp(num2str(s_link))
%     path_link_id = s_link;
%     check_case = 5; 
% end

path_len = length(path_link_id);

path_weight_t_hr10 = zeros(path_len,1); % duration 10:00 am
path_weight_t_hr22 = zeros(path_len,1); % duration 22:00 pm
path_weight_d = zeros(path_len,1);      % distance

path_weight_NOx_07_hr10 = zeros(path_len,1); % NOx from MY2007 truck at 10am
path_weight_NOx_07_hr22 = zeros(path_len,1); 
path_weight_NOx_12_hr10 = zeros(path_len,1);
path_weight_NOx_12_hr22 = zeros(path_len,1);

path_weight_PM25_07_hr10 = zeros(path_len,1);
path_weight_PM25_07_hr22 = zeros(path_len,1);
path_weight_PM25_12_hr10 = zeros(path_len,1);
path_weight_PM25_12_hr22 = zeros(path_len,1);

path_weight_CO2_07_hr10 = zeros(path_len,1);
path_weight_CO2_07_hr22 = zeros(path_len,1);
path_weight_CO2_12_hr10 = zeros(path_len,1);
path_weight_CO2_12_hr22 = zeros(path_len,1);

truck_route_label = zeros(path_len,1);

for i1 = 1:size(path_link_id) 
    kk = Tw1.network_id == path_link_id(i1); % !!
    if sum(kk) > 0
        path_weight_t_hr10(i1) = Tw1.duration_hr10_s(kk);
        path_weight_t_hr22(i1) = Tw1.duration_hr22_s(kk);
        path_weight_d(i1) = Tw1.length_m(kk);        
        
        path_weight_NOx_07_hr10(i1) = Tw2.NOx_07_hr10_ug(kk);
        path_weight_NOx_07_hr22(i1) = Tw2.NOx_07_hr22_ug(kk);
        path_weight_NOx_12_hr10(i1) = Tw2.NOx_12_hr10_ug(kk);
        path_weight_NOx_12_hr22(i1) = Tw2.NOx_12_hr22_ug(kk);
        
        path_weight_PM25_07_hr10(i1) = Tw2.PM25_07_hr10_ug(kk);
        path_weight_PM25_07_hr22(i1) = Tw2.PM25_07_hr22_ug(kk);
        path_weight_PM25_12_hr10(i1) = Tw2.PM25_12_hr10_ug(kk);
        path_weight_PM25_12_hr22(i1) = Tw2.PM25_12_hr22_ug(kk);
             
        path_weight_CO2_07_hr10(i1) = Tw3.CO2_07_hr10_g(kk);
        path_weight_CO2_07_hr22(i1) = Tw3.CO2_07_hr22_g(kk);
        path_weight_CO2_12_hr10(i1) = Tw3.CO2_12_hr10_g(kk);
        path_weight_CO2_12_hr22(i1) = Tw3.CO2_12_hr22_g(kk);
        
        truck_route_label(i1) = Tw1.tr_label2(kk);
    end
end

T_out = [path_weight_d, path_weight_t_hr10, path_weight_t_hr22,...
    path_weight_NOx_07_hr10, path_weight_NOx_07_hr22, path_weight_NOx_12_hr10,...
    path_weight_NOx_12_hr22, path_weight_PM25_07_hr10, path_weight_PM25_07_hr22, path_weight_PM25_12_hr10,...
    path_weight_PM25_12_hr22,...
    path_weight_CO2_07_hr10, path_weight_CO2_07_hr22, ...
    path_weight_CO2_12_hr10, path_weight_CO2_12_hr22];

sumTable=sum(T_out,1); % must specify the sum dimention 1 here otherwise error when only one link is the solution

end
