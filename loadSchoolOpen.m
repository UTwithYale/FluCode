% tNum_Begin = 20090401;
% tNum_End = 20091231;
function VacRatioM = loadSchoolOpen(tNum_Begin, tNum_End, VacRatio)
%% Names of 500 cities
MyFileInfo = strcat('500_Cities_List.xlsx');
%500_Cities_List.xlsx: StateAbbr	StateDesc	CityName	UniqueID	PopulationCount	GeoLocation	CityFIPS
[NUM,TXT,RAW] = xlsread( strcat('DataSource/', MyFileInfo)); %
State500 = TXT(2:end,2);
CityName = TXT(2:end,3);
StateFullName= (TXT(2:end,2));
TXT_Mas6 = TXT(2:end,[6]);
Num_Mas6 = []; % for airport
for i=1:length(TXT_Mas6)
    temp = TXT_Mas6{i};
    temp = strsplit(temp,{'(', ')', ','});
    Num_Mas6 = [Num_Mas6; str2num(temp{2}) str2num(temp{3})];
end
FIPS_C = NUM(:,4); % FIPS_C


%%
load data.mat
for i=1:size(data,1);
    temp = ( data{i,3} );
    if temp=='19430'
        data{i,3} = '19380';
    end
end
%T4 = array2table(data, 'VariableNames' , {'cityFIPS', 'CountryFIPS' , 'CBSA' ,  'flightD'});

MetroID = []; for i=1:size(data,1); MetroID = [MetroID; str2num( data{i,3} )]; end
CityID  = []; for i=1:size(data,1); CityID = [CityID; str2num( data{i,1} )]; end
CityAndMetro = [CityID MetroID];
CityAndMetro_U = unique(CityAndMetro, 'rows');
tempMG = unique(MetroID);

%%
filename = 'schoolopenings/schoolopenings-2009.xlsx';
% [~, Humidity_Time, Humidity_CBSA, Humidity_mean] = textread(filename,  '%s%f%f%f', 'delimiter',',', 'headerlines', 1);
[ndata, text, alldata] = xlsread(filename);

% in xlsx of 38581 denoting 18 Aug 2009
% in matlab of 734003 denoting 18 Aug 2009
TimeGap = 734003-38581;
ndata = ndata+TimeGap;
tempM = [];
for i=1:size(ndata,2)
    temp = ndata(:,i);
    tempM = [tempM; mean(temp(~isnan(temp)))];
end

% Hawaii is set as 09/7/30 in State of Hawaii – Department of Education
% 2009-2010 OFFICIAL SCHOOL CALENDAR due to: https://www.nctq.org/dmsView/23-09_7375

% cityFIPS, metropolitan ID, meanOpenTime
% cityFIPS, stateID, meanOpenTime
[~, tempb] = ismember(StateFullName, text);
OpenTime500 = tempM(tempb);
temp2 = [FIPS_C, OpenTime500];
temp2 = unique(temp2, 'rows');

%  metropolitan ID, meanOpenTime,
% [CityAndMetro_U(:,2) temp2(:,2) ];
meanOpenTime = [];
for i=1:length(tempMG)
    meanOpenTime = [meanOpenTime; median(temp2(  CityAndMetro_U(:,2) == tempMG(i),2) )];
end

meanCloseTime_2019Summer = meanOpenTime-60;
meanOpenTime_2019Autumn = meanOpenTime;
% meanCloseTime_2019Winter = meanOpenTime;


% tNum_Begin = 20090401;
% tNum_End = 20091231;
DateNum_Begin = datenum( num2str(tNum_Begin),'yyyymmdd');
DateNum_End = datenum( num2str(tNum_End),'yyyymmdd');

VacRatioM = ones(DateNum_End-DateNum_Begin+1, length(tempMG));

for i=1:length(meanCloseTime_2019Summer)
    temp = (meanCloseTime_2019Summer(i)-DateNum_Begin+1):(meanOpenTime_2019Autumn(i)-DateNum_Begin+1);
    VacRatioM( round(temp), i) = VacRatio;
end

for i=1:size(VacRatioM,1)
    temp = weekday(DateNum_Begin-1+i);
    if temp==1 | temp==7
        VacRatioM(i,:) = VacRatio;
    end
end

% November 25-27, 2009 Thanksgiving Break 
DateNum_Begin2 = datenum( num2str(20091125),'yyyymmdd');
DateNum_End2 = datenum( num2str(20091127),'yyyymmdd');
VacRatioM( (DateNum_Begin2-DateNum_Begin+1):(DateNum_End2-DateNum_Begin+1),:) = VacRatio;


% Sept 7, 2009 Labor day Break 
DateNum_Begin2 = datenum( num2str(20090907),'yyyymmdd');
DateNum_End2 = datenum( num2str(20090907),'yyyymmdd');
VacRatioM( (DateNum_Begin2-DateNum_Begin+1):(DateNum_End2-DateNum_Begin+1),:) = VacRatio;


% Dec.21, 2009 to 3 Jan 2010 Chrismas Break 
DateNum_Begin2 = datenum( num2str(20091221),'yyyymmdd');
DateNum_End2 = datenum( num2str(20100103),'yyyymmdd');
VacRatioM( (DateNum_Begin2-DateNum_Begin+1):(DateNum_End2-DateNum_Begin+1),:) = VacRatio;





