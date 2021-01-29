%% Codes:
% Get list of cities
% Generate work flows between counties and flight flowsbetween airports

%% Metropolitan of 500 cities
% Get list of cities
MyFileInfo = strcat('list2_Sep_2018.xlsx');
%list2_Sep_2018.xls: CBSA Code	CBSA Title	Metropolitan/Micropolitan Statistical Area	Principal City Name	FIPS State Code	FIPS Place Code
[NUM2,TXT2,RAW2] = xlsread( strcat('DataSource/', MyFileInfo)); %
TXT2 = TXT2(4:(end-4), :);

Metropolitan = TXT2(:,2);
for i=1:length(Metropolitan)
    Metropolitan{i} = char(strrep(Metropolitan{i}, " ", ''));
end
[Metropolitan b] = unique(Metropolitan);
CBSA_Code =  TXT2(b,1);

Metro_Sht = [];
for i=1:length(Metropolitan)
    temp = Metropolitan{i};
    temp = strsplit(temp, ',');
    Metro_Sht{i,1} = temp{1};
end
% CityName2 = TXT2(:,4);
% temp = Metropolitan(ismember(FIPS_C2, FIPS_C));


%% Names of 500 cities
MyFileInfo = strcat('500_Cities_List.xlsx');
%500_Cities_List.xlsx: StateAbbr	StateDesc	CityName	UniqueID	PopulationCount	GeoLocation	CityFIPS
[NUM,TXT,RAW] = xlsread( strcat('DataSource/', MyFileInfo)); %
State500 = TXT(2:end,2);

CityName = TXT(2:end,3);
% StateName= lower(TXT(2:end,2));
StateName= (TXT(2:end,1));
TXT_Mas6 = TXT(2:end,[6]);
Num_Mas6 = []; % for airport
for i=1:length(TXT_Mas6)
    temp = TXT_Mas6{i};
    temp = strsplit(temp,{'(', ')', ','});
    Num_Mas6 = [Num_Mas6; str2num(temp{2}) str2num(temp{3})];
end

FIPS_C = NUM(:,4); % FIPS_C
clear CityName  i  MyFileInfo  NUM  RAW temp  TXT  TXT_Mas6  TXTR

%%
MyFileInfo = strcat('data_Politan.txt');
[data_Politan] = textread( strcat('WebScraler/', MyFileInfo), '%s'); %
for i = 1:length(data_Politan)
    data_Politan{i} = strcat(data_Politan{i}, ",", StateName{i});
    data_Politan{i} = char( strrep(data_Politan{i}, '"', '') );
    if data_Politan{i} == "Philadelphia-Camden-Wilmington,DE"
        data_Politan{i} = "Philadelphia-Camden-Wilmington,PA-NJ-DE-MD";
    end
    if data_Politan{i} == "Chicago-Naperville-Elgin,IN"
        data_Politan{i} = "Chicago-Naperville-Elgin,IL-IN-WI";
    end
    if data_Politan{i} == "KansasCity,KS" 
        data_Politan{i} = "KansasCity,MO-KS";
    end
    if data_Politan{i} == "Providence-Warwick,MA" 
        data_Politan{i} = "Providence-Warwick,RI-MA";
    end
    if data_Politan{i} == "Lansing-EastLansing.,MI" 
        data_Politan{i} = "Lansing-EastLansing,MI" ;
    end
    if data_Politan{i} == "Philadelphia-Camden-Wilmington,NJ" 
        data_Politan{i} = "Philadelphia-Camden-Wilmington,PA-NJ-DE-MD" ;
    end
    if data_Politan{i} == "NewYork-Newark-JerseyCity,NJ" 
        data_Politan{i} = "NewYork-Newark-JerseyCity,NY-NJ-PA" ;
    end
    if data_Politan{i} == "Dayton,OH"
        data_Politan{i} = "Dayton-Kettering,OH";
    end
    if data_Politan{i} == "Charlotte-Concord-Gastonia,SC"
        data_Politan{i} = "Charlotte-Concord-Gastonia,NC-SC";
    end
    
    if data_Politan{i} == "Washington-Arlington-Alexandria,VA"
        data_Politan{i} = "Washington-Arlington-Alexandria,DC-VA-MD-WV";
    end
    if data_Politan{i} == "Portland-Vancouver-Hillsboro,WA"
        data_Politan{i} = "Portland-Vancouver-Hillsboro,OR-WA";
    end
    if data_Politan{i} == "Chicago-Naperville-Elgin,WI"
        data_Politan{i} = "Chicago-Naperville-Elgin,IL-IN-WI";
    end
end

temps = [];
for i=1:length(data_Politan)
    temp = find( contains(Metropolitan, data_Politan{i}) == 1 );
    if length(temp)>0
        temps = [temps temp];
    else
        temps = [temps -1];
    end
end
% find(temps==0,1)
% sum(temps)
temps_Name = temps;
temps_CBSA = CBSA_Code(temps); % 217 unique politans

MetroID = [];
for i=1:length(temps_CBSA)
    MetroID = [MetroID; str2num(temps_CBSA{i})];
end
% temps_CBSA temps_Name

%%
% AIRPORT_SEQ_ID	AIRPORT_ID	AIRPORT	DISPLAY_AIRPORT_NAME	DISPLAY_AIRPORT_CITY_NAME_FULL	AIRPORT_WAC	AIRPORT_COUNTRY_NAME	AIRPORT_COUNTRY_CODE_ISO	AIRPORT_STATE_NAME	AIRPORT_STATE_CODE	AIRPORT_STATE_FIPS	CITY_MARKET_ID	DISPLAY_CITY_MARKET_NAME_FULL	CITY_MARKET_WAC	LAT_DEGREES	LAT_HEMISPHERE	LAT_MINUTES	LAT_SECONDS	LATITUDE	LON_DEGREES	LON_HEMISPHERE	LON_MINUTES	LON_SECONDS	LONGITUDE	AIRPORT_START_DATE	AIRPORT_THRU_DATE	AIRPORT_IS_CLOSED	AIRPORT_IS_LATEST
MyFileInfo = strcat('518711423_T_MASTER_CORD.xlsx');
[NUM_Mas, TXT_Mas, RAW_Mas] = xlsread( strcat('DataSource/', MyFileInfo)); %

NUM_Ind = find(NUM_Mas(:,end)==1);
NUM_Mas = NUM_Mas(NUM_Ind,:);
TXT_Mas = TXT_Mas(NUM_Ind,:);

NUM_Mas1924 = NUM_Mas(:,[19 24]);
% Num_Mas6
d1kms2 = [];
for i=1:length(NUM_Mas1924)
    i/length(NUM_Mas1924)
    for j=1:length(Num_Mas6)
        [d1km d2km]=lldistkm(NUM_Mas1924(i,:), Num_Mas6(j,:)); %(km)
        d1kms2(i,j) = d1km;
    end
end

airports = [];
for i=1:size(d1kms2,2)
    airports =[airports; find( d1kms2(:,i)==min( d1kms2(:,i) ), 1)];
end
flightID = NUM_Mas(airports, 2);
% (1.2) [FIPS_C flightID] as cityID, flightID
% (1.6) [FIPS_C flightID MetroID] as cityID, flightID, MetroID
% 385 airports shared by 500 cities
%% YEAR,	QUARTER,	ORIGIN_AIRPORT_ID,	ORIGIN_AIRPORT_SEQ_ID,	ORIGIN_CITY_MARKET_ID,	DEST_AIRPORT_ID,
% DEST_AIRPORT_SEQ_ID,	DEST_CITY_MARKET_ID,	PASSENGERS
% MyFileInfo = strcat('900381120_T_DB1B_MARKET.xlsx');
% [NUM_fl, TXT_fl, RAW_fl] = xlsread( strcat('DataSource/', MyFileInfo)); %

fid = fopen('DataSource/900381120_T_DB1B_MARKET.csv');
NUM_fl = textscan(fid, '%f,%f,%f,%f,%f,%f,%f,%f,%f,', 'Headerlines', 1);
fclose(fid);
NUM_fl = cell2mat(NUM_fl);

NUM_fl2 = NUM_fl(ismember(NUM_fl(:,3), flightID) & ismember(NUM_fl(:,6), flightID), [3,6,9] );
% (1.5) NUM_fl2 as [flightO, flightD, flight travel Num] in this YEAR & QUARTER
Metro_GroundNum = [];
for i=1:size(NUM_fl2, 1)
    if mod(i,10000)==1
        i/size(NUM_fl2, 1)
    end
    flightO = NUM_fl2(i, 1);
    cityO = FIPS_C(find(flightID==flightO,1));
    metroO = MetroID(find(flightID==flightO,1));
    
    flightD = NUM_fl2(i, 2);
    cityD = FIPS_C(find(flightID==flightD,1));
    metroD = MetroID(find(flightID==flightD,1));
    
    Metro_GroundNum = [Metro_GroundNum; metroO metroD NUM_fl2(i,3)];
end
% (3) Metro_GroundNum as [MetroO, MetroD, flight travel Num]
%%
% %%
% % (1) [FIPS_C, MetroID] as cityID and metroID from web
% %% Data from R codes for 500 cities
% MyFileInfo = strcat('R_Metro_Dem.xlsx');
% %500_Cities_List.xlsx: StateAbbr	StateDesc	CityName	UniqueID	PopulationCount	GeoLocation	CityFIPS
% [NUMR,TXTR,RAWRR] = xlsread( MyFileInfo ); %
% % NUMR(:,[2,7]); % CBSA and STPLFIPS as metro/micropolitans vs. cities
% % NUMR(:,[2,5]); % CBSA vs. FIPS as CBSA and Country
% CCC = NUMR(:,[2, 5, 7]); % CBSA and Country and cities, and further flights
% % citiesR = NUMR(:,[7]);
% CCC_New = [];
% for i=1:size(CCC, 1)
%     if MetroID(find(FIPS_C==CCC(i,3))) == CCC(i,1)
%         CCC_New = [CCC_New; CCC(i,:)  flightID( find(FIPS_C==CCC(i,3)) ) ];
%     end
% end
% 
% 
% T = array2table(CCC_New, 'VariableNames' , {'CBSA' , 'Country' , 'cities', 'flights'});
% writetable(T, strcat('CCC_New.csv'), 'WriteVariableNames', true);

%%
County_states = lower(State500);
[NUM2,TXT2,RAW2] = xlsread('R_countiesALL.xlsx'); %
TXT2 = lower(TXT2);
TXT2N = TXT2(:,2);

MyFileInfo = strcat('data_County.xlsx');
[NUM3,TXT3,RAW3] = xlsread(strcat('WebScraler/', MyFileInfo)); %
TXT32 = TXT3;
TXT3 = lower(TXT3);

i_states = 0;
for i = 1:size(TXT3,1)
    for j = 1:size(TXT3,2)
        [i j]
        temp = TXT3{i, j};
        if size(temp, 1)>0
            i_states = i_states+1;
            if County_states{i} == "north carolin"
                County_states{i} = "north carolina";
            end
            if County_states{i} == "south carolin"
                County_states{i} = "south carolina";
            end
            temp2 = strcat( char(temp), ", ", char(County_states{i}) );
            temp2 = strrep(temp2, '"', '');
            if temp2=="district of columbia federal district, district of c"
                temp2 = lower("District of Columbia, District of Columbia");
            end
            if temp2=="doã±a ana county, new mexico"
                temp2 = lower("Do<U+00F1>a Ana County, New Mexico")
            end
            temp3 = find(ismember(TXT2N, temp2)==1, 1);
            if length(temp3)==0
                1
            else
                TXT32{i,j} = temp3;
            end
        else
            TXT32{i,j} = 0;
        end
    end
end

data = [];
for i = 1:size(TXT32,1)
    for j = 1:size(TXT32,2)
        if TXT32{i,j}>0
            temp = {num2str(FIPS_C(i)), num2str(RAW2{TXT32{i,j}, 2}), num2str(MetroID(i)), num2str(flightID(i))};
            if FIPS_C(i)<10^6 & FIPS_C(i)~=15003
                temp{1} = strcat('0', num2str(FIPS_C(i)) );
            end
            if RAW2{TXT32{i,j}, 2}<10000
                temp{2} = strcat('0', num2str(RAW2{TXT32{i,j}, 2}) );
            end
            data = [data; temp];
        end
    end
end

T4 = array2table(data, 'VariableNames' , {'cityFIPS', 'CountryFIPS' , 'CBSA' ,  'flightD'});
writetable(T4, strcat('T4.csv'), 'WriteVariableNames', true);


length(unique(data(:,1)))
length(unique(data(:,2)))
length(unique(data(:,3)))
length(unique(data(:,4)))
save data.mat data % cityid, countryid, CBSA, flightid

%% 
% CBSA Code	Metro Division Code	CSA Code	CBSA Title	Level of CBSA	Status, 1=metro 2=micro	Metropolitan Division Title	CSA Title	Component Name	State	FIPS	County Status
% MyFileInfo = strcat('list3.xlsx');
% [NUM_Country, TXT_Country, RAW_Country] = xlsread( strcat('DataSource/', MyFileInfo)); %
% MetroInd = ismember(TXT_Country(:, 5), 'Metropolitan Statistical Area');
% TXT_Country = TXT_Country(MetroInd,:);
% RAW_Country = RAW_Country(MetroInd,:);
% CBSA_Country = TXT_Country(:,1);
% FIPS_Country = TXT_Country(:,11);

% County	State	SSA State county code	FIPS State county code	CBSA (Blanks are Rural)	CBSA NAME (Blanks are Rural)
MyFileInfo = strcat('FY 2017 FR Cty CBSA Xwalk and CBSA Con Cty.xlsx');
[NUM_Country, TXT_Country, RAW_Country] = xlsread( strcat('DataSource/', MyFileInfo)); %
CBSA_Country = TXT_Country(2:end, 5);
FIPS_Country = TXT_Country(2:end, 4);

% (4) FIPS_Country, CBSA_Country
% [FIPS_Country CBSA_Country]


%%
MyFileInfo = strcat('table3_shorten.txt');
% State_FIPS_Code,County_FIPS_Code,Minor_Civil Division FIPS Code,State Name,County Name,Minor Civil Division Name,State FIPS Code,County FIPS Code,Minor Civil Division FIPS Code,State Name,County Name,Minor Civil Division Name,Workers in Commuting Flow,Margin of Error
[State_FIPS_Code_L, County_FIPS_Code_L, Minor_Civil_Division_FIPS_Code_L, ...
    State_Name_L, County_Name_L, Minor_Civil_Division_Name_L, ...
    State_FIPS_Code, County_FIPS_Code,Minor_Civil_Division_FIPS_Code, ...
    State_Name,  County_Name,Minor_Civil_Division_Name, ...
    Workers_in_Commuting_Flow, Margin_of_Error, a] ...
    = textread( strcat('DataSource/', MyFileInfo), ...
    '%s%s%s%s%s%s%s%s%s%s%s%s%f%s%s', 'delimiter',',');%, 1, 'A1:N100'); %

% (5) countryO, countryD, GroundTravelNum
for i=1:length(State_FIPS_Code_L)
    if mod(i,10000)==1
        i/length(State_FIPS_Code_L)
    end
    countryO{i} = strcat(State_FIPS_Code_L{i}, County_FIPS_Code_L{i});
    countryD{i} = strcat(State_FIPS_Code{i}, County_FIPS_Code{i});
    GroundTravelNum{i} = Workers_in_Commuting_Flow(i);
end

clear State_FIPS_Code_L  County_FIPS_Code_L  Minor_Civil_Division_FIPS_Code_L  ...
    State_Name_L  County_Name_L  Minor_Civil_Division_Name_L  ...
    State_FIPS_Code  County_FIPS_Code Minor_Civil_Division_FIPS_Code  ...
    State_Name   County_Name Minor_Civil_Division_Name  ...
    Workers_in_Commuting_Flow  Margin_of_Error  a;

GroundTravelNums = []; temp_MetroGNum=[];
for i=1:length(GroundTravelNum)
    temp_MetroGNum = [temp_MetroGNum; GroundTravelNum{i}];
    
    if mod(i,10000)==1
        i/length(GroundTravelNum)
        GroundTravelNums = [GroundTravelNums; temp_MetroGNum];
        temp_MetroGNum = [];
    end
end
GroundTravelNums = [GroundTravelNums; temp_MetroGNum];

%%
for i=1:length(countryO)
    if mod(i,10000)==1
        i/length(countryO)
    end
    
    temp1 = find(ismember(FIPS_Country, countryO{i})==1);
    temp2 = find(ismember(FIPS_Country, countryD{i}(2:end))==1);
    
    if length(temp1)>0 & length(temp2)>0
        MetroGO{i} = CBSA_Country{temp1};
        MetroGD{i} = CBSA_Country{ temp2 };
        GroundTN{i} = GroundTravelNum{i};
%     else
%         1
    end
end

MetroGNum = []; temp_MetroGNum = [];
for i=1:length(MetroGO)
    if length(MetroGO{i})>0 & length(MetroGD{i})>0 
        temp_MetroGNum = [temp_MetroGNum; str2num(MetroGO{i}) str2num(MetroGD{i}) (GroundTN{i}) ];
    end
    
    if mod(i,10000)==1
        i/length(MetroGO)
        MetroGNum = [MetroGNum; temp_MetroGNum];
        temp_MetroGNum = [];
    end
end
MetroGNum = [MetroGNum; temp_MetroGNum];

% MetroGNum and Metro_GroundNum
% save loadCDCData_Web.mat -v7.3
% save('loadCDCData_Web_short.mat', 'MetroGNum', 'Metro_GroundNum', '-v7.3');