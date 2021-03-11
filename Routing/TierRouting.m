% 
% T = readtable('../streets_LB_TomTom_ID_v2.csv');
% Tw1 = readtable('../length_duration_weights.csv');
% Tw2 = readtable('../inhaled_mass_weights.csv');
% Tw3 = readtable('../CO2_weights.csv');

%% generate directed graph

% T_network = readtable('streets_LB_TomTom_UTM_gen20m_split_truck_roads.csv');

function [sumTable, path_link_id, b1, tr_percent, path_link_keep]=TierRouting(s_link, t_link, G,node_map, T, Tw1,Tw2, Tw3,node1,node2,edge, weight_for_routing, half_label, routeID)

%disp(['Start link: ' num2str(s_link) ' Dest link: ' num2str(t_link)])
s_node = node_map(s_link);
t_node = node_map(t_link);
%plot(G);
%[nodepath, path_weight, edgepath] = shortestpath(G,s_node,t_node); %!!

[~, ~, edgepath] = shortestpath(G, s_node, t_node); %!!
e_path = G.Edges(edgepath,:);
% edgepath is the nominal edge path, not the correct network ID
if isempty(e_path)
    sumTable = zeros(1, 15);
    path_link_id = [];
    b1 = 0;
    tr_percent = 0;
    path_link_keep = [];
    disp(['Empty path! Start and end link: ' num2str(s_link) ' ' num2str(t_link) ]) 
    return
end

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
        truck_route_label(i1) = Tw1.tr_label2(kk);
    end
end
tr_percent = sum(truck_route_label)/length(path_link_id)*100;
% special case: entire path is on truck route already
if tr_percent == 100 && half_label == 0 
    b1 = t_link;
    path_link_keep = []; 
    %tr_p_entry_stretch = 100;
elseif tr_percent == 100 && half_label == 1 
    b1 = s_link;
    path_link_keep = []; 
elseif tr_percent ~= 100 && half_label == 0  % first half of OD pairs: entry points to store, entry points are already on truck routes 
    if any(truck_route_label)  % true if there is at least one truck route link in the path
        % last link of longest consecutive truck route, which will be
        % the break point of truck route and non-truck route
        temp_path = flip(path_link_id);  % flip the path id and find the last truck route link in the path
        temp_tr_label = flip(truck_route_label);
        ind = find(temp_tr_label==1, 1); % the last truck route link in the path
        link1 = temp_path(ind); % the break point of truck route, special case: b1 is starting point, b2 is second link
        %tr_p_entry_stretch = ind/length(path_link_id)*100;
        if ind > 1
            link2 =  temp_path(ind -1);   % the break point of nontruck route 
            path_link_keep = flip(temp_path(1: ind-1));
            
            JID1a = T.F_JID(T.IDtxt == link1);
            JID1b = T.T_JID(T.IDtxt == link1);
            JID2a = T.F_JID(T.IDtxt == link2);
            JID2b = T.T_JID(T.IDtxt == link2);
            b1 = intersect([ JID1a JID1b], [JID2a, JID2b]);
            
        else  % ind == 1 if truck route is actually at the end
            b1 = t_link;
            %disp(['Warning break point is at route ' num2str(routeID) '. Index of break point and length of path are: ' num2str(ind)  ' and '  num2str(length(path_link_id))])
            %num2str(ind) 
            %num2str(length(path_link_id))
            path_link_keep = [];
        end
        
        
    else  % if none of the path is truck route, then use the shortest path 
        %b1 = 0;
        b1 = 0;
        path_link_keep = path_link_id;
        disp('Warning no truck route link in the path' )
    end
elseif tr_percent ~= 100 && half_label == 1  % second half of OD pairs:  store back to entry points, stores are not always on truck routes, end point (entry/exit) is always on truck routes
    if any(truck_route_label)
        %ind = max(diff(find(diff([NaN temp_tr_label' NaN]))));
        ind = find(truck_route_label==1, 1) ; %max(diff(find(diff([NaN truck_route_label' NaN]))));
        %tr_p_entry_stretch = ind/length(path_link_id)*100;
        link2 = path_link_id(ind); % the break point of nontruck route, because switched start and end points, now break point two is truck route
        if ind >1  % spcial case: b2 is end point, b1 is second to last link
            link1 =  path_link_id(ind - 1);   % the break point of nontruck route 
            path_link_keep = path_link_id(1: ind-1);
            JID1a = T.F_JID(T.IDtxt == link1);
            JID1b = T.T_JID(T.IDtxt == link1);
            JID2a = T.F_JID(T.IDtxt == link2);
            JID2b = T.T_JID(T.IDtxt == link2);
            b1 = intersect([ JID1a JID1b], [JID2a, JID2b]);
        else %ind==1
            %disp(['Warning break point is at route ' num2str(routeID) '. Index of break point and length of path are: ' num2str(ind)  ' and '  num2str(length(path_link_id))])
            %num2str(ind) 
            %num2str(length(path_link_id))
            b1 = s_link;
            path_link_keep = [];
        end
        
    else
        b1 = 0;
        path_link_keep = path_link_id;
        disp('Warning no truck route link in the path' )
    end
    
else
    disp(['Truck route percent and half label are: ' num2str(tr_percent) , num2str(half_label)])
end

if isempty(path_link_keep)
    sumTable = zeros(1, 15);
else
    for i1 = 1:size(path_link_keep) 
        kk = Tw1.network_id == path_link_keep(i1); % !!
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
%writetable(sumTable, ['sum_weight_' inputWeight  '.csv']); %!!!!!!!!

end
