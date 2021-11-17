% compare results at system level

pr = 0.2; timethresh = 0.2; % all in decimal;
%t1 = zeros(length(pr)*length(timethresh)+1, 7); % caution
k = 1609.3; % meter to mile conversion factor
ta = 'hr10'; %22pm, 10am
%MYa ='12';

tb =  ta;
%MYb = MYa;

inputWeightA=['t_' ta];  %'';  % compare t, dist, IM, weight_t_10am
inputWeightB=['inhale_MY12_' tb];  % tier_routing_t_hr10,  tier_routing_inhale_MY12_hr10

cTableA = readtable(['cost_table_' inputWeightA '.csv']);  % fastest route
cTableB = readtable(['cost_table_' inputWeightB '.csv']);  % LER

%----random run1
[rs_yes_tech, baseline_trips] = RandomSample(pr); 
t_change=(cTableB.(['t_' tb])-cTableA.(['t_' ta]))./cTableA.(['t_' ta]);
take_LER_pot = t_change <= timethresh & t_change>0; % take LER potential
take_LER = rs_yes_tech.* take_LER_pot;
no_LER = baseline_trips - take_LER;
cTableA_syn = cTableA{:,:}.*baseline_trips;
cTableB_syn_1 = cTableA{:,:}.* no_LER + cTableB{:,:}.*take_LER;
%---- random run2
[rs_yes_tech, baseline_trips] = RandomSample(pr); 
t_change=(cTableB.(['t_' tb])-cTableA.(['t_' ta]))./cTableA.(['t_' ta]);
take_LER_pot = t_change <= timethresh & t_change>0; % take LER potential
take_LER = rs_yes_tech.* take_LER_pot;
no_LER = baseline_trips - take_LER;
cTableB_syn_2 = cTableA{:,:}.* no_LER + cTableB{:,:}.*take_LER;
%----random run3
[rs_yes_tech, baseline_trips] = RandomSample(pr); 
t_change=(cTableB.(['t_' tb])-cTableA.(['t_' ta]))./cTableA.(['t_' ta]);
take_LER_pot = t_change <= timethresh & t_change>0; % take LER potential
take_LER = rs_yes_tech.* take_LER_pot;
no_LER = baseline_trips - take_LER;
cTableB_syn_3 = cTableA{:,:}.* no_LER + cTableB{:,:}.*take_LER;

cTableB_syn = (cTableB_syn_1 + cTableB_syn_2 + cTableB_syn_3)/3;

sumA = sum(cTableA_syn(:, [3 4 8 12 18 19 16]), 1); %  time, dist, PM inhale, NOx inhale, PM emis, NOx emis
sumA(:,1) =sumA(:,1)/k; sumA(:,2) = sumA(:,2)/3600; sumA(:, 5:end)= sumA(:, 5:end)/1000;
sumB = sum(cTableB_syn(:, [3 4 8 12 18 19 16]), 1);
sumB(:,1) =sumB(:,1)/k; sumB(:,2) = sumB(:,2)/3600; sumB(:, 5:end)= sumB(:, 5:end)/1000;
