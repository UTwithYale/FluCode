function res4 = plotByDL4_ageG(fitness, hourlyPerD, Pop_MetroAll, tmpR, shiftWeek, Pop_Metro, tmpAgeAttack38_52_propor)

DL = fitness{1}; % incidence of Symptomatic illnesses
% DH_shf = fitness{2}; % incidence of Hospitalizations
% D = fitness{3}; % incidence of Deaths
% YT = fitness{4}; % incidence of Treated

%% DL
temp_DL = DL;

temp = permute(sum(temp_DL, [5,4,2]) , [3, 1, 2]);
I_day = [];
for i = 1:size(temp, 1)
    temp1 = temp(i,:);
    I_day = [I_day aggHourlyToDay(hourlyPerD, temp1)];
end
I_dayS100 = I_day;
res4{1} = I_dayS100;

temp = reshape(sum(sum(sum(sum(DL,5),4),2),1),1,5);
res4{5} = temp./sum(sum(sum(Pop_Metro(:,:,:,:),4),3),1);

%% DH_shf
for ifit = 2:4
    temp_DL = fitness{ifit};
    temp = permute( sum(temp_DL, [5,4,2]) , [3, 1, 2]);
    I_day = [];
    for i = 1:size(temp, 1)
        temp1 = temp(i,:);
        I_day = [I_day aggHourlyToDay(hourlyPerD, temp1)];
    end
    I_dayS100 = I_day;
    res4{ifit} = I_dayS100;
end