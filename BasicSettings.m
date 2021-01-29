
setTpeirod_hourlyPerD;

%%
addpath('mossong_2008')
load data.mat 
for i=1:size(data,1);
    temp = ( data{i,3} );
    if temp=='19430'
        data{i,3} = '19380';
    end
end
MetroID = []; for i=1:size(data,1); MetroID = [MetroID; str2num( data{i,3} )]; end
CityID  = []; for i=1:size(data,1); CityID = [CityID; str2num( data{i,1} )]; end
CityAndMetro = [CityID MetroID];
CityAndMetro_U = unique(CityAndMetro, 'rows');
tempMG = unique(MetroID);

%%
ageG = zeros(5,1);
ageG2 = zeros(5,1);

[Pop_Metro Pop_Metro17]= getPopFromRv4_17(tempMG, CityAndMetro_U, ageG2); 
Pop_MetroAll = sum(sum(sum(Pop_Metro,2),3),4);
Contact_Day = getMobility(tempMG);
Contact_Day(logical( eye(size(Contact_Day,1)) ))=0;

ageNum = size(Pop_Metro,2);
Px = ones(ageNum,1)/ageNum;
for i=1:size(Contact_Day,1)
    for j=1:size(Contact_Day,2)
        for a=1:ageNum
            xi(i,j,a) = Contact_Day(i,j)*Px(a)/ sum(sum(Pop_Metro(i,a,:,:)));
        end
    end
end

%%
vaccin.vaccin_PriorG = [1,2]; vaccin.vaccin_Prior = 10^5;

%%
% 217 locations     2 years    12 months   31 days   24 hours
Humidity_mean_data = readCsvs(tempMG);
Humidity_mean_day_data = mean(Humidity_mean_data, 5);

save BasicSettings.mat