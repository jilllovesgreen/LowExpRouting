% some of Shorstest Dist based on City TR are much longer than Fastest Route 
% because it takes too much detour on truck routes
% this script find such routes and replace them with fastest routes
% therefore the SDCTR is more realistic
inputWeightA='t_hr10';  %'d_TR_tr_weightv2_syn';  % compare t, dist, IM, weight_t_10am
inputWeightB='tier_routing_d_TR';  %'inhale_MY12_10am_tr_weightv2';  % inhale_MY12_22pm, d_TR_tr_weightv2, inhale_MY12_10am_tr_weightv2
t = 'hr10';  % 10: am, 22: 10 pm
MYa ='12';
MYb ='12';
cTableA = readtable(['cost_table_' inputWeightA '.csv']);
cTableB = readtable(['cost_table_' inputWeightB '.csv']);

t_change=(cTableB.(['t_' t])-cTableA.(['t_' t]))./cTableA.(['t_' t])*100;
% --- replace d_TR routes whoes time is 30% more than fastest routes into
% fastest routes
table_syn = cTableB;  % synthesized table
table_syn(t_change > 30,:) = cTableA(t_change > 30,:);
writetable(table_syn, ['cost_table_d_TR_syn_' t '.csv'])