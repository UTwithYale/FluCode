% Para: PriorAge, PriorRisk
% Para: vaccin_PriorG and vaccin_Num
% Para: PertVacci Iy the probability of taking vaccine.


% Weekly allocation (proportional to state population) and 
% For first 125 million doses: 0.20 (0–4), 0.50 (5–17), 0.2 (18–49), 0.09 (50–64), 0.01 (65+)
% After 125 million doses (children should have achieved 80% coverage): 0.00 (0–4), 0.00 (5–17), 0.25 (18–49), 0.35 (50–64), 0.40 (65+)
% Coverage target 
% 80% for each age group 
% Once 80% of the age group is vaccinated, vaccination would end and remaining vaccine allocated to other age groups

Prop_assigns = [0.20, 0.50, 0.2,  0.09, 0.01; ...
    0, 0.00, 0.25, 0.35, 0.40];

if vac_num_used * Pop_US /Pop_217  <= 125*10^6
    Prop_assign = Prop_assigns(1,:);
else
    Prop_assign = Prop_assigns(2,:);
end

% vac_num_used
vaccin_NumLast = vaccin_Num;
if vaccin_Num>0
    for i = 1:LocationNum
        for a = 1:size(S,3)
            for r = 1:size(S,4)
                v = 1; % let v=1 denote the unvaccianted group
                tmp_Pop = S(t-1, i, a, r, v) + L(t-1, i, a, r, v)+ C(t-1, i, a, r, v) ...
                    + Ia(t-1, i, a, r, v)+Iy(t-1, i, a, r, v)+It(t-1, i, a, r, v) + R(t-1, i, a, r, v);
                % tmp_Pop = (tmp_Pop*PertVacci(a) );
                
                temp_vac = vaccin_Num * Pop_MetroAll(i) / Pop_217 * Prop_assign(a) /sum(Prop_assign(a:end));
                
                if t-1==91 && i==1 && a==1 && r==2 && v==1
                    1;
                end
                
                if isinf(temp_vac)
                    temp_vac = 0;
                end
                if temp_vac>0.8*tmp_Pop
                    V(t-1, i, a, r) = 0.8*tmp_Pop * PertVacci(a);
                    vaccin_Num = vaccin_Num-0.8*tmp_Pop;
                    vac_num_used = vac_num_used+0.8*tmp_Pop;
                else
                    V(t-1, i, a, r) = temp_vac*PertVacci(a);
                    vaccin_Num = vaccin_Num - temp_vac;
                    vac_num_used = vac_num_used+temp_vac;
                end
            end
        end
    end
   
end