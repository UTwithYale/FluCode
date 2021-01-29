addpath('mossong_2008')
load tmpCaseRateUS.mat
load BasicSettings.mat;
shiftWeek = 0;
vaccin.vaccin_Num = 0;
vaccin.drug_Num = 0;
hourlyPerD = 5;
% EPI, VX, SC, AV, flag_realVX
Strategy_List = [1,0,0,0,0;1,1,0,0,0;1,2,0,0,0;1,0,1,0,0;1,0,2,0,0;1,0,3,0,0;1,0,4,0,0;1,0,0,1,0;1,1,1,1,0;1,1,2,1,0;1,2,3,1,0;1,2,4,1,0;2,0,0,0,0;2,1,0,0,0;2,2,0,0,0;2,0,1,0,0;2,0,2,0,0;2,0,3,0,0;2,0,4,0,0;2,0,0,1,0;2,1,1,1,0;2,1,2,1,0;2,2,3,1,0;2,2,4,1,0;3,0,0,0,1;3,3,0,0,1;3,0,0,1,1;3,3,0,1,1;3,0,3,0,1;3,0,4,0,1;3,3,3,1,1;3,3,4,1,1];
clear R0ps

for i = 1
    temp = [0.2, 2, 5.2];
    parfor scalej= 1:20
        temp2 = temp;
        temp2(1) = (0.10+0.01*scalej); %EPIs(EPI);
        j=1;
        EPI = Strategy_List(i,j); j=j+1;
        VX = Strategy_List(i,j); j=j+1;
        SC = Strategy_List(i,j); j=j+1;
        AV = Strategy_List(i,j); j=j+1;
        flag_realVX = Strategy_List(i,j); j=j+1;  % whether using real vaccine supply in US
        
        % set transmission rate
        iX = temp2(1); betasRatio = temp2(2); iI0 = temp2(3);
        
        % set rand seeds
        rng(1);
        
        for jjj=1:100
            randBin(100, 0.2, 1); randPois(100, 1);
        end
        
        % 26420 is for houston
        LocationNumTest = find(tempMG==26420);
        
        Tpeirod = 30;
        R0p5=[];
        humidityYesOrNot = 1; HoustonFlag=1; sto1div5=0;
        for fjj = 1:2
            fitness = model_SECIR_V20_sto_Jul([iX iI0], xi, Pop_Metro, ...
                vaccin, Humidity_mean_data, ageG2, ageG, ...
                Pop_MetroAll, Humidity_mean_day_data, [], LocationNumTest, ...
                [], [], ...
                shiftWeek, Pop_Metro17, betasRatio, [], VX, SC, AV, flag_realVX, EPI, Tpeirod, sto1div5, humidityYesOrNot, HoustonFlag);
            
            fitness12{scalej} = fitness;
            
            temp3 = 100000*(sum(sum(sum(fitness{1},5),4),3))';
            I_day = sum(temp3(LocationNumTest, :),1);
            I_dayS100 = I_day./sum(Pop_MetroAll(LocationNumTest));
            
            temp4 = squeeze(sum(reshape(I_dayS100, hourlyPerD, length(I_dayS100)/hourlyPerD),1));
            indexD = 7;%datenum('2009-10-30')-datenum('2009-10-15')+1;
            a = 1+log(temp4([indexD:(indexD+6)])/temp4([indexD:(indexD+6)]-1))* 3.628;
            
            R0p5 = [R0p5; a];
        end
        R0ps{scalej, 1} =  [ mean(R0p5), iX];
    end
    R0psAll{i} =  R0ps;
end

save getR0dym_houston3.mat


