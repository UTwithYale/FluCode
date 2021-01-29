
function fitness = model_SECIR_V20_sto_Jul2(x, xi, Pop_Metro, vaccin,Humidity_mean_data, ageG2, ageG, Pop_MetroAll, Humidity_mean_day_data, R0City1,...
    LocationNumTest, tempCaseRateState, stateName, shiftWeek, Pop_Metro17, betasRatio, schoolClosure_Day, VX, SC, AV, flag_realVX, EPI, temps_startSeed)
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

if ~exist('sto1div5')
    sto1div5 = 0;
end

% setTpeirod_hourlyPerD;
% Tpeirod = 7*53;
Tpeirod = 7*15;%7*14;
% Tpeirod = 7*35;%7*14;
hourlyPerD = 5; % 5 parts each day

% schoolImpactSettings;
% θ(t): impact of school calendar on transmIysion rates among children on day t; set to one on school days, reduced on weekends, holItays and school closures
theta = 0.451;
tNum_Begin = 20090401;
tNum_End = 20100331;
tNum_Begin_sim = 20090914; %+7*shiftWeek;
DateNum_Begin = datenum( num2str(tNum_Begin_sim),'yyyymmdd')+7*shiftWeek;
tNum_Year = 20090101;
DateNum_BeginSinceYear = DateNum_Begin-datenum( num2str(tNum_Year),'yyyymmdd')+1;

%%
LocationNum = size(Pop_Metro,1);
ageNum = size(Pop_Metro,2); riskNum = size(Pop_Metro,3); vacNum = size(Pop_Metro,4);

DL = zeros(Tpeirod*hourlyPerD, LocationNum, ageNum, riskNum, vacNum); % (t,i,a,r,v)

VacationSetting;
Vac_2009;
para_Mar_Jul_v22;

fitness = DL; 

end

