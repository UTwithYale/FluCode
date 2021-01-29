
function fitness = model_SECIR_V20_sto_Jul(x, xi, Pop_Metro, vaccin,Humidity_mean_data, ageG2, ageG, Pop_MetroAll, Humidity_mean_day_data, R0City1,...
    LocationNumTest, tempCaseRateState, stateName, shiftWeek, Pop_Metro17, betasRatio, schoolClosure_Day, VX, SC, AV, flag_realVX, EPI, Tpeirod, sto1div5, humidityYesOrNot, HoustonFlag)
%% Parameter Calibration
temp = [betasRatio, 1]*x(1);
I0 = 0.001*x(2);
x3 = [1  temp];

%% reduction in local transmIysion rate for indivItuals vIyiting other nodes
Rho= x3(1);
%% beta = b1*exp(-180*temp_HumItity)+b2;
b1 = x3(2); b2 = x3(3);

if ~exist('humidityYesOrNot')
    humidityYesOrNot = 1;
end

if ~exist('HoustonFlag')
    HoustonFlag = 0;
end

% setTpeirod_hourlyPerD;
if ~exist('Tpeirod')
    Tpeirod = 7*53;
end

if ~exist('sto1div5')
    sto1div5 = 0;
end


% Tpeirod = 7*28;%7*14;
% Tpeirod = 7*40;%7*14;
% Tpeirod = 7*65;%7*14;
hourlyPerD = 5; % 5 parts each day

% Tpeirod = 7*30;
% hourlyPerD = 1; % 5 parts each day

% schoolImpactSettings;
% θ(t): impact of school calendar on transmIysion rates among children on day t; set to one on school days, reduced on weekends, holItays and school closures
theta = 0.451;
%%
% timeSetting;
% tNum_Begin = 20090406; % Week 17
tNum_Begin = 20090401; %tNum_End = 20191231; % tNum_End = 20100228;
% tNum_End = 20111015; 
tNum_End = 20100331;
% Week 38	September 14, 2009

if flag_realVX==0
    tNum_Begin_sim = 20090901; %+7*shiftWeek;
else
    tNum_Begin_sim = 20090914; %+7*shiftWeek;
end

if HoustonFlag==1
%     tNum_Begin_sim = 20090901;
    tNum_Begin_sim = 20091121;
%     tNum_Begin_sim = 20091201;
    
%     tNum_Begin_sim = 20091001;
%     tNum_Begin_sim = 20091015;
%     tNum_Begin_sim = 20091214;
end
DateNum_Begin = datenum( num2str(tNum_Begin_sim),'yyyymmdd')+7*shiftWeek;
tNum_Year = 20090101;
DateNum_BeginSinceYear = DateNum_Begin-datenum( num2str(tNum_Year),'yyyymmdd')+1;

%%
LocationNum = size(Pop_Metro,1);
ageNum = size(Pop_Metro,2); riskNum = size(Pop_Metro,3); vacNum = size(Pop_Metro,4);

DL = zeros(Tpeirod*hourlyPerD, LocationNum, ageNum, riskNum, vacNum); % (t,i,a,r,v)

VacationSetting;
Vac_2009;
% para_Mar_Jul_v2;% similar as parasettings_Sto_Preimmunity_xi_Mark
para_Mar_Jul_v22;

% save test.mat
% fitness = 100000*sum(sum(sum(sum(DL,5),4),3),2)';
% fitness = DL; 
% fitness = It;

% Risk of hospitalization among high-risk symptomatically infected individuals by age (overall) (%)
% Risk of hospitalization among non-high-risk symptomatically infected individuals by age (overall) (%):  
Risk_hosp = [4.2,1.6,2.9,2.9,16.6; ...
    4.2,1.6,2.9,2.9,16.6; ...
    0.4,0.1,0.2,0.2,2.4]*0.01;

if EPI==2
    Risk_hosp(:,1:4) = Risk_hosp(:,1:4)*2;
end

Day_hosp = [4.2,3.9,5.4,6.8,6.1; ...
    4.2,3.9,5.4,6.8,6.1; ...
    2.7,3.6,5.8,7.5,5.1];
% Day_hosp = zeros(3,5);
Hour_hosp = round(Day_hosp*hourlyPerD);

DH = DL(:,:,:,1:3,1);
for t=1:size(DH,1)
    for i=1:size(DH,2)
        for r=1:3
            DH(t,i,:,r) = permute(DH(t,i,:,r), [3,1,2]).*Risk_hosp(r,:)';
        end
    end
end




% save test.mat
DH_shf = zeros(Tpeirod*hourlyPerD, LocationNum, ageNum, 3);
for t=1:size(DH,1)
    for i=1:size(DH,2)
        for a=1:5
            for r=1:3
                temp = Hour_hosp(r,a);
                if t-temp>0
                    DH_shf(t,i, a, r) =  DH_shf(t, i,a, r)+DH(t-temp, i, a, r);
                end
            end
        end
    end
end


% Percentage of persons with influenza symptoms who seek healthcare and are diagnosed with flu (%):
% P_d = [67, 52, 37, 43, 56]*0.01;

% Reduction in all-cause hospital admission (adults)
% RR 0.37 (0.17–0.81)
% Dobson (attached above)

for t=1:size(DH,1)
    for i=1:size(DH,2)
        for a=1:5
            for r=1:3
                DH_shf(t,i, a, r) =  DH_shf(t, i,a, r)*P_d(a)*0.37 + DH_shf(t, i,a, r)*(1-P_d(a));
            end
        end
    end
end



tempD =  zeros(Tpeirod*hourlyPerD, LocationNum, ageNum, riskNum, vacNum); % (t,i,a,r,v)
tempD(2:end,:,:,:,:) = D(2:end,:,:,:,:) - D(1:(end-1),:,:,:,:);

fitness{1} = DL; % incidence of Symptomatic illnesses
fitness{2} = DH_shf; % incidence of Hospitalizations
fitness{3} = tempD; % incidence of Deaths
fitness{4} = YT; % incidence of Treated

end

% if plotFlag==100
%     Iasd_Single = 100000*sum( sum(sum(sum(DL,5),4),3),2)'./sum(Pop_MetroAll(LocationNumTest) );
%     Iasd_Single7 = aggDayToWeek(Iasd_Single, hourlyPerD);
%     fitness = (sum(tempCaseRateState( (shiftWeek+38):(shiftWeek+38-1+length(Iasd_Single7))) - Iasd_Single7).^2)/length(Iasd_Single7);
% end
