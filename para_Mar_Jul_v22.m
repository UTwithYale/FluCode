
%% A. In Sleet of WHO_NREVSS_Clinical_Labs2009 of WHO_NREVSS_Clinical_Labs.xls
% Sum between periods from # to #
% 2009 Week 15: April 6, 2009	April 12, 2009
% 2009 Week 32	August 3, 2009	August 9, 2009
P_immunity = 997.6614/10000;

% [49094, 89668, 49499, 6715]; AGE 0-4, AGE 5-24, 25-64, 65
% 0-4, 5-17, 18-24, 25-64, 65+
temp = [49094, 89668*(17-5-1)/(24-4), 89668*(24-17)/(24-4), 49499, 6715];
temp = temp./sum(temp);
% P_immunity = P_immunity*temp;
if flag_realVX==1
    P_immunity = P_immunity*temp;
else
    P_immunity = zeros(1,5);
end

Humidity_mean_day_data191 = Humidity_mean_day_data(:, 1, 9, 1);

%%
if length(ageG) == 5
    if flag_realVX==1
        temp = [1.96, 1.96, 1,    1,    0.17];
    else
        temp = [1.96, 1.96, 1.96, 1.96, 1.96];
    end
    phi = mangle_ZD16_Mark(Pop_Metro17)/hourlyPerD;
    for i=1:5
        phi(i,:) = phi(i,:)*temp(i);
    end
    
    PMatrix_noSchool = [0.54, 0.77, 0.91,0.93,0.990000000000000;0.72,0.54,0.83,0.88,0.990000000000000;0.89,0.72,0.92,0.97,0.990000000000000;0.76,0.58,0.93,0.97,0.990000000000000;0.97,0.95,0.98,0.97,0.980000000000000];
    PMatrix_noWorknoSchool = [0.54,0.77,0.91000000000,0.93,0.990000000000000;0.72,0.53,0.75,0.81,0.990000000000000;0.89,0.65,0.48,0.48,0.950000000000000;0.76,0.51,0.52,0.68,0.970000000000000;0.97,0.95,0.95,0.95,0.980000000000000];
    phiVacation = phi.*PMatrix_noWorknoSchool;
    phiSchool = phi.*PMatrix_noSchool;
end


%% Para: PertVacci Iy the probability of taking vaccine and efficacy.
% 0-4, 5-17, 18-24, 25-64, 65+
if length(ageG) == 5
    PertVacciNew = [62, 62, 62*2/7+74*5/7, 74,  62]*0.01;
end
PertVacci = PertVacciNew;
% PertVacci = ones(5, 1);
%% relative probability to be infected for indivItual unvaccinated or vaccinated
pv = [1 1];

%% Parameter settings
vaccin_Num = 0;
% P_vaSucceed as probability of vaccination successful, Vaccine effectiveness
P_drugSucceed = 1;

%  PertDrug Iy the probability of taking vaccine.
PertDrug = 0.1;
% PertDrug = 0;
drug_Num = vaccin.drug_Num; %0% 10^6; % 10^5;
delayVac = 14;
% pai as probability of vaccination successful
% π: probability of treatment
pai = 0;

% Now set vaccine group v as 2 for vaccianted group
LocationNum = size(Pop_Metro,1);
ageNum = size(Pop_Metro,2); riskNum = size(Pop_Metro,3); vacNum = size(Pop_Metro,4);


% initial compartments
S = zeros(Tpeirod*hourlyPerD, LocationNum, ageNum, riskNum, vacNum); % (t,i,a,r,v)
L = S; C = S; Ia = S; Iy = S; It = S; R = S; DL=S; D = S; YT = S;
% Vaccin Num
V = zeros(Tpeirod*hourlyPerD, LocationNum, ageNum, riskNum); % 7M storage
% Drug Num
Dg = zeros(Tpeirod*hourlyPerD, LocationNum, ageNum, riskNum, vacNum); % 7M storage

S(1,:,:,:,:) = 1*Pop_Metro;
% totalIniInfectionPerc = (Pop_MetroAll(i) )*I0; %

I0 = poissrnd(I0*1000)/1000;
totalIniInfectionPerc = (Pop_MetroAll ).*I0; %
numInfectedF = totalIniInfectionPerc';

t=1; a=2; r=2; v=1;
if length(LocationNumTest)==1 && flag_realVX==0
    % Place of Emergence: U.S.A. (100 cases in Wichita, KS)
    % LocationNumTest should be the
    i = LocationNumTest;
    numInfectedF = 100;
    numInfectedF = poissrnd(numInfectedF);
    S(t,i,a,r,v) = S(t,i,a,r,v)-numInfectedF;
    Iy(t,i,a,r,v) = Iy(t,i,a,r,v)+numInfectedF;
else
    i=1:LocationNum;
    S(t,i,a,r,v) = S(t,i,a,r,v)-numInfectedF(i);
    Iy(t,i,a,r,v) = Iy(t,i,a,r,v)+numInfectedF(i);
    
    S(find(S<0))=0;
    Iy(find(Iy<0))=0;
end

% A. Considering the first wave for preimmunity
% v=1:2 % 1 as unvacination; 2 as vacination
% P_immunity
% temp = P_immunity*(S(t,i,a,r, 1)+S(t,i,a,r, 2));
% S(t,i,a,r, 1) = S(t,i,a,r, 1)-temp;
% S(t,i,a,r, 2) = S(t,i,a,r, 2)+temp;
for ia=1:size(Pop_Metro,2)
    S(t,:,ia,:, 1) = S(t,:,ia,:, 1)-P_immunity(ia)*(S(t,:,ia,:,1)+S(t,:,ia,:,2));
    S(t,:,ia,:, 2) = S(t,:,ia,:, 2)+P_immunity(ia)*(S(t,:,ia,:,1)+S(t,:,ia,:,2));
end

% probability for drug over age groups
% P_drug = 0;
if AV == 0
    P_d = 0*ones(1,ageNum);
end
if AV == 1
    % Percentage of persons with influenza symptoms who seek healthcare and are diagnosed with flu (%):
    P_d = [67, 52, 37, 43, 56]*0.01;
end


% probability for symptom over age groups
% τ: symptomatic proportion (probability that an exposed indivItual will progress to the symptomatic compartment rather than the asymptomatic compartment)

if sto1div5==0
    tau = 0.55*ones(1,ageNum);
else
    tau = 0.55*ones(1,ageNum);
end
% 0-4, 5-17, 18-24, 25-64, 65+
if length(ageG) == 5
    ea = 1-[0, 9.5, 18.5, 6.8, 15.7]*0;
end

if length(ageG) == 4
    DeathHighR = [.0296, .0296, .0414, .0414, ...
        (.0414*3+.2082*2)/5, .2082, .2082, .2082, .2082, .2082, ...
        .2082, 0.2191, 0.2191, 0.2191, 2.88, 2.88, 2.88]*0.01;
    va = [];
    for k=1:length(ageG2)
        tempG = ageG2{k};
        va = [va; mean(DeathHighR(tempG))];
    end
end
% 0-4, 5-17, 18-24, 25-64, 65+
if length(ageG) == 5
    va = [0.0296, 0.0414, 0.2082, 0.2082, 2.88]*0.01;
    vaL = [0.0032, 0.0027, 0.0106, 0.0106, .416]*0.01;
    
    if EPI==2
        va(1:4) = va(1:4)*2;
        vaL(1:4) = vaL(1:4)*2;
    end
end

DeathLowR = [.0032, .0032, .0027, .0027, ...
    (.0027*3+.0106*2)/5, .0106, .0106, .0106, .0106, .0106, ...
    .0106, 0.0182, 0.0182, 0.0182, .416, .416, .416]*0.01;


%% relative probability of infection for Ia, Iy, It compartment
%
% ωa,X: relative infectiousness of infectious indivItuals of age a in compartment X (X indicates either asymptomatic (A), symptomatic (Y) or treated (T))
% ωa,A=0.5 and ωa,Y=1 (Elveback et al. 1976; Yang et al. 2009)
% persons who are symptomatic are 0.5X as likely to transmit as those who are symptomatic, as 0.5
omega_c = 0.68; omega_a = 0.36; omega_y=1; omega_t= 1;%0.5;

% 1/days from E to C
% κ: exposed rate (average duration between E and C Iy 1/κ)
% 1/0.98 with the mean incubation period Iy 0.98 day (Yang et al. 2009)
k_E = 1/0.98/hourlyPerD;

% 1/days from C to I, σ
sigma_C = 1/0.5/hourlyPerD*ones(ageNum,1);
gamma_a2r = 1/3/hourlyPerD;
gamma_t2r=[1/(5-17.6/24), 1/(5-17.6/24), 1/(5-25.2/24), 1/(5-25.2/24),1/(5-25.2/24)]/hourlyPerD; % gamma_t2r=1/4.3/hourlyPerD;

if sto1div5==0
    gamma_y2r=1/5/hourlyPerD;
else
    pd4 = makedist('Triangular','a',3,'b',5,'c',7);
    %     pd4 = makedist('Uniform','lower', 3, 'upper', 7);
    tempPd4  =random(pd4,1,1);
    gamma_y2r=1/tempPd4/hourlyPerD;
    gamma_t2r=[1/(tempPd4-17.6/24), 1/(tempPd4-17.6/24), ...
        1/(tempPd4-25.2/24), 1/(tempPd4-25.2/24),1/(tempPd4-25.2/24)]/hourlyPerD;
end

Pop_US = 306.8*10^6;
Pop_217 = sum(Pop_MetroAll);

%%
flagSchoolClose = 0;
flagSchoolOpenAgain = 0;
vac_num_used = 0;

everTransferLocation = zeros(LocationNum,1);
WeekNum = 0; WeekNum_last = 0;

WeekNum_Start = 10^5;
if VX==1 || SC==1 || SC==2
    WeekNum_Start = weeknum(datenum( '20090901','yyyymmdd')+100);
end

if VX==2 || SC==3 || SC==4
    WeekNum_Start = weeknum(datenum( '20090901','yyyymmdd')+150);
end

if VX==3
    Vac_shift_day = DateNum_Begin - datenum( '20090901','yyyymmdd');
    Vac_date_num(:,1) = Vac_date_num(:,1)-Vac_shift_day;
end

% allocate memory for lambdas matrix
lambdas = zeros(ageNum, riskNum, vacNum);
tempt_Past = -1; % record of last time, used to decide whether t is another day and thus to update phi2
DayNum_last = 0;

for t = 2 : (Tpeirod*hourlyPerD)
    WeekNum = floor( (t-1) /hourlyPerD)/7 + weeknum('2009-09-01') -1;
    DayNum = floor( (t-1) /hourlyPerD)+1;
    tempt = ceil( (t-1+0.1)/hourlyPerD)+1;
    if VX>0
        if flag_realVX==0
            if WeekNum~=WeekNum_last && WeekNum<=(10+WeekNum_Start) && WeekNum>=WeekNum_Start
                vaccin_Num = vaccin_Num+25*10^6 * Pop_217/Pop_US;
            end
        else
            if DayNum~=DayNum_last
                tempVN = Vac_date_num(find(Vac_date_num(:,1)== floor( (t-1) /hourlyPerD) ),2);
                if length(tempVN)==0;  tempVN=0; end
                vaccin_Num = vaccin_Num + tempVN;
            end
        end
        %%
        applyVaccinCDC;
        for i = 1:LocationNum
            for a = 1:ageNum
                for r = 1:riskNum
                    for v = 1
                        Move_pop_VaccinCDC;
                    end
                end
            end
        end
    end
    %% humItity-specific transmIysion rate
    if sto1div5==0
        tempTimeVector = datevec( DateNum_Begin+ceil( (t-1+0.1)/hourlyPerD)+1);
    else
        tempTimeVector = datevec( DateNum_Begin+ceil( (t-1+0.1)/hourlyPerD)+1 - randi(14)+7);
    end
    tempT = [tempTimeVector(2), tempTimeVector(3)];
    temp_HumItity = Humidity_mean_day_data(:, tempTimeVector(1)-2008, tempT(1), tempT(2));
    
    if temp_HumItity==0
        temp_HumItity = Humidity_mean_day_data(:, tempTimeVector(1)-2009, tempT(1), tempT(2));
    end
    
    if humidityYesOrNot==0
        temp_HumItity = Humidity_mean_day_data191;
    end
    betaH = (b1*exp(-180*temp_HumItity)+b2);
    
    temp_Si = [];
    temp_DLi = [];
    temp_Li = [];
    temp_Ci = [];
    temp_Iai = [];
    temp_Iyi = [];
    temp_Iti = [];
    temp_Ri = [];
    temp_Di = [];
    
    
    Pop_Metroij = sum(sum(Pop_Metro(:,:,:,:),4),3);
    % tempt is daily scale, only when a new day, we update xi and phi2
    if tempt_Past<tempt
        tempt_Past = tempt;
        
        for i = 1:LocationNum
            %% Transmission
            phi2 = phi;
            beta = betaH;
            tempt = ceil( (t-1+0.1)/hourlyPerD)+1;
            if VacationRatioM(tempt,i)<1
                for phiii = 1:5
                    phi2(phiii,phiii) = phi(phiii,phiii)*0.45;
                end
            end
            
            %% Closure Trigger (CT)
            if SC==0
                flagSchoolClose = 0;
            end
            if flagSchoolClose==0
                if t>7*hourlyPerD
                    % Closure Trigger (CT)
                    if SC==1 || SC ==3
                        if sum( DL( (( t-7*hourlyPerD):(t-1) ),:,:,:,:) , 'all')> Pop_217*0.01
                            flagSchoolClose = 1;
                        end
                    end
                    if SC==2 || SC ==4
                        if flag_realVX == 0
                            %             “Close at Day X” (October 28, 2009: first Monday of October)
                            schoolClosure_Day = hourlyPerD * (datenum('2009-10-28')-DateNum_Begin);
                            if t>=schoolClosure_Day
                                flagSchoolClose = 1;
                            end
                        else
                            % For SC4, schools do not reopen after the summer break
                            schoolClosure_Day = hourlyPerD * (datenum('2009-09-01')-DateNum_Begin);
                            if t>=schoolClosure_Day
                                flagSchoolClose = 1;
                            end
                        end
                    end
                end
            end
            
            %% School closure (SC) scenarios
            % SC1: Illness prevalence trigger, early (VX1) vaccine availability
            % SC2: Close at Day X, early (VX1) vaccine availability
            % SC3: Illness prevalence trigger, normal (VX2) vaccine availability
            % SC4: Close at Day X, normal (VX2) vaccine availability
            % SCs = 0:1:4; % 0 for no SC
            %     US pop        306.8 million (2009)
            if flagSchoolClose==1 && flagSchoolOpenAgain==0
                tempFlag = 0;
                if flag_realVX == 0
                    
                    if SC==1 || SC==2
                        tempFlag = (t>hourlyPerD * (datenum('2009-09-01')-DateNum_Begin+100+7+35-1));
                    end
                    if SC==4 || SC==3
                        tempFlag = (t>hourlyPerD * (datenum('2009-09-01')-DateNum_Begin+150+7+35-1));
                    end
                    
                else
                    if SC==3
                        tempFlag = t>hourlyPerD * (datenum('2019-12-02')-DateNum_Begin);
                    end
                    if SC==4
                        tempFlag = t>hourlyPerD * (datenum('2019-12-02')-DateNum_Begin);
                    end
                    
                end
                if tempFlag
                    flagSchoolOpenAgain = 1;
                end
            end
            
            if VacationRatioM(tempt,i)==1 && flagSchoolClose==1 && flagSchoolOpenAgain==0
                if VacationRatioM(tempt,i)==1
                    phi2 = phiVacation;
                else
                    phi2 = phiSchool;
                end
            end
            
            xis{i} = xi;
            phi2s{i} = phi2;
        end
    end
    
    Ets = zeros(LocationNum, 1);
    for i = 1:LocationNum
        Ets(i) = sum(L(t-1,i, :,:,:), 'all')+sum(C(t-1,i, :,:,:), 'all')+sum(It(t-1,i, :,:,:), 'all')+sum(Iy(t-1,i, :,:,:), 'all')+sum(Ia(t-1,i, :,:,:), 'all');
    end
    LocationNum_j = find(Ets>0);
    
    everTransferLocation(find(Ets>=100)) = 1;
    
    for i = 1:LocationNum
        % relative probability of infection for Ia, Iy, It compartment
        Pop_Metroi = sum(sum(Pop_Metro(i,:,:,:),4),3);
        
        temp_S = zeros(ageNum, riskNum, vacNum);
        temp_E=zeros(ageNum, riskNum, vacNum);
        
        temp_DL = zeros(ageNum, riskNum, vacNum);
        temp_L=zeros(ageNum, riskNum, vacNum);
        temp_C=zeros(ageNum, riskNum, vacNum);
        temp_Ia=zeros(ageNum, riskNum, vacNum);
        temp_Iy=zeros(ageNum, riskNum, vacNum);
        temp_It=zeros(ageNum, riskNum, vacNum);
        temp_R=zeros(ageNum, riskNum, vacNum);
        temp_D=zeros(ageNum, riskNum, vacNum);
        temp_YT = zeros(ageNum, riskNum, vacNum);
        Et = Ets(i); %; sum(sum(sum(L(t-1,i, :,:,:))))+sum(sum(sum(It(t-1,i, :,:,:))))+sum(sum(sum(Iy(t-1,i, :,:,:))))+sum(sum(sum(Ia(t-1,i, :,:,:))));
        
        phi2_i = phi2s{i};
        xi = xis{i};
        
        if Et < 1 && mod(t, hourlyPerD)==1 && everTransferLocation(i) ==0
            
            % between nodes
            lambdas = zeros(ageNum, riskNum, vacNum);
            lambdas(:,:,:) = 0;
            for a = 1:ageNum
                for r = 1:riskNum
                    for v = 1%:vacNum
                        tempSb = 0; temp1Sb = 0 ;
                        for jL = 1:length(LocationNum_j)
                            j=LocationNum_j(jL);
                            Pop_Metroj = Pop_Metroij(j,:);
                            phi2_j = phi2s{j};
                            temp_SE=0; temp1S=0;
                            if i~=j && xi(i,j,a) > 0
                                tempP52 = C(t-1, j, :, :, 1); tempP52 = reshape(tempP52, size(tempP52, 3), size(tempP52, 4));
                                tempA52 = Ia(t-1, j, :, :, 1); tempA52 = reshape(tempA52, size(tempA52, 3), size(tempA52, 4));
                                tempY52 = Iy(t-1, j, :, :, 1); tempY52 = reshape(tempY52, size(tempY52, 3), size(tempY52, 4));
                                tempT52 = It(t-1, j, :, :, 1); tempT52 = reshape(tempT52, size(tempT52, 3), size(tempT52, 4));
                                temp_SE = sum( phi2_j(a,:)./Pop_Metroj * beta(j) * (omega_a * tempA52 +  omega_y* tempY52 + omega_c' .* tempP52  + omega_t * tempT52 )  ) ;
                                temp1S  = sum( phi2_i(a,:)./Pop_Metroj * beta(i) * (omega_a * tempA52 + omega_c' .* tempP52 )) ;
                            end
                            tempSb = tempSb +xi(i,j,a)*temp_SE;  % indivIHuals of age a travel from i to j wIHh a daily probabilIHy of xi(i,j,a)
                            temp1Sb= temp1Sb+xi(j,i,a)*temp1S; % /hourlyPerD
                            tempSb(isnan(tempSb)) = 0;
                            temp1Sb(isnan(temp1Sb)) = 0;
                        end
                        lambdas(a,r,v) = (tempSb+temp1Sb);
                    end
                end
            end
            temp = S( t-1, i, :, :, :);
            temp(find(temp<0))=0;
            lambdas(find(lambdas>1))=1; stoFlag=1;
            dEa = permute( randBin( temp, lambdas*hourlyPerD , stoFlag ), [3,4,5,1,2] );
            temp_S = -dEa ;
            temp_E = +dEa ;
        end
        % within nodes
        j = i;
        for a = 1:ageNum
            for r = 1:riskNum
                for v = 1%:vacNum
                    temp10 = 0;
                    for b = 1:ageNum
                        for r2 = 1:riskNum
                            for v2 = 1%:vacNum
                                temp10 = temp10+ ...
                                    beta(i) * phi2(a,b) *  S(t-1,i,a,r,v)/Pop_Metroi(b) * (omega_c *C(t-1,j,b,r2,v2)+...
                                    omega_a *Ia(t-1,j,b,r2,v2)+...
                                    omega_y *Iy(t-1,j,b,r2,v2)+...
                                    omega_t *It(t-1,j,b,r2,v2) );
                            end
                        end
                    end
                    temp10s = temp10;
                    temp1 = sum(temp10s) * ea(a);
                    if isnan(temp1)
                        temp1=0;
                    end
                    
                    if r==1 || r==2 % preg and high as high
                        tempVa = va;
                    else
                        tempVa = vaL;
                    end
                    
                    temp2 = k_E*L(t-1,i,a,r,v);
                    %                     temp2c = sigma_C(a)*(1-tempVa(a))*C(t-1,i,a,r,v);
                    temp2c = sigma_C(a)*C(t-1,i,a,r,v);
                    % temp3 = gamma_a2r*(1-tempVa(a))*Ia(t-1,i,a,r,v);
                    temp3 = gamma_a2r*Ia(t-1,i,a,r,v);
                    temp4 = gamma_y2r*(1-tempVa(a))*Iy(t-1,i,a,r,v);
                    temp5 = gamma_t2r(a)*(1-(1-0.36)* tempVa(a))*It(t-1,i,a,r,v);
                    
                    
                    %                     temp2d = sigma_C(a)*tempVa(a)*C(t-1,i,a,r,v);
                    temp2d = 0;%sigma_C(a)*tempVa(a)*C(t-1,i,a,r,v);
                    
                    temp3d = 0;% gamma_a2r*tempVa(a)*Ia(t-1,i,a,r,v);
                    temp4d = gamma_y2r*tempVa(a)*Iy(t-1,i,a,r,v);
                    
                    % Reduction in mortality (critically ill children)
                    % aOR: 0.36 (0.16–0.84)
                    % Louie 2013: https://pediatrics.aappublications.org/content/132/6/e1539
                    temp5d =  gamma_t2r(a)*(1-0.36)* tempVa(a)*It(t-1,i,a,r,v);
                    
                    if t>0
                        if temp1>1; temp1 = randPois(temp1, 1); end % randBin( temp, lambdas , stoFlag )
                        if isinf(temp1); temp1 = 0; end
                        
                        if temp2c>1; temp2c = randPois(temp2c, 1); end
                        if isinf(temp2c); temp2c = 0; end
                        
                        
                        if temp3>1; temp3 = randPois(temp3, 1); end
                        if isinf(temp3); temp3 = 0; end
                        
                        if temp4>1; temp4 = randPois(temp4, 1); end
                        if isinf(temp4); temp4 = 0; end
                        
                        if temp5>1; temp5 = randPois(temp5, 1); end
                        if isinf(temp5); temp5 = 0; end
                        
                        if temp2d>1; temp2d = randPois(temp2d, 1); end
                        if isinf(temp2d); temp2d = 0; end
                        
                        if temp3d>1; temp3d = randPois(temp3d, 1); end
                        if isinf(temp3d); temp3d = 0; end
                        
                        if temp4d>1; temp4d = randPois(temp4d, 1); end
                        if isinf(temp4d); temp4d = 0; end
                        
                        if temp5d>1; temp5d = randPois(temp5d, 1); end
                        if isinf(temp5d); temp5d = 0; end
                    end
                    
                    if temp_S(a,r,v)<0
                        1;
                    end
                    dS = -temp1;
                    temp_S(a,r,v) = S(t-1,i,a,r,v)+dS+ temp_S(a,r,v);
                    if temp_S(a,r,v)<0
                        temp1 = S(t-1,i,a,r,v)+ temp_S(a,r,v);
                        temp_S(a,r,v)=0;
                    end
                    if temp_S(a,r,v)<0
                        1;
                    end
                    
                    dL = temp1 - temp2;
                    temp_L(a,r,v) = L(t-1,i,a,r,v)+dL + temp_E(a,r,v);
                    if temp_L(a,r,v)<0
                        temp2 = L(t-1,i,a,r,v)+temp1+ temp_E(a,r,v);
                        temp_L(a,r,v)=0;
                    end
                    
                    %                     temp_DL(a,r,v) = tau(a)*(1-P_d(a))*temp2c;
                    temp_DL(a,r,v) = tau(a)*temp2c;
                    if temp_DL(a,r,v)<0
                        temp2c = 0;
                        temp_DL(a,r,v)=0;
                    end
                    
                    dC = temp2 - temp2c- temp2d;
                    temp_C(a,r,v) = C(t-1,i,a,r,v)+dC;
                    if temp_C(a,r,v)<0
                        temp2c = C(t-1,i,a,r,v)+temp2-temp2d;
                        temp_C(a,r,v)=0;
                    end
                    
                    
                    dIa= (1-tau(a))*temp2c - temp3- temp3d;
                    temp_Ia(a,r,v) = Ia(t-1,i,a,r,v)+dIa;
                    if temp_Ia(a,r,v)<0
                        temp3 = Ia(t-1,i,a,r,v)+(1-tau(a))*temp2c-temp3d;
                        temp_Ia(a,r,v)=0;
                    end
                    
                    
                    dIy= tau(a)*(1-P_d(a))*temp2c - temp4- temp4d;
                    temp_Iy(a,r,v) = Iy(t-1,i,a,r,v)+dIy;
                    if temp_Iy(a,r,v)<0
                        temp4 = Iy(t-1,i,a,r,v)+tau(a)*(1-P_d(a))*temp2c-temp4d;
                        temp_Iy(a,r,v)=0;
                    end
                    
                    
                    temp_YT(a,r,v) =  tau(a)*P_d(a)*temp2c;
                    
                    dIt= temp_YT(a,r,v) - temp5 - temp5d;
                    temp_It(a,r,v) = It(t-1,i,a,r,v)+dIt;
                    if temp_It(a,r,v)<0
                        temp5 = It(t-1,i,a,r,v)+temp_YT(a,r,v)- temp5d;
                        temp_It(a,r,v)=0;
                    end
                    
                    
                    dR = temp3 + temp4 + temp5;
                    temp_R(a,r,v) = R(t-1,i,a,r,v)+dR;
                    
                    dD = temp2d + temp3d + temp4d + temp5d;
                    temp_D(a,r,v) = D(t-1,i,a,r,v)+dD;
                    
                end
            end
        end
        %             end
        temp_Si{i} = temp_S;
        temp_DLi{i} = temp_DL;
        temp_Li{i} = temp_L;
        temp_Ci{i} = temp_C;
        temp_Iai{i} = temp_Ia;
        temp_Iyi{i} = temp_Iy;
        temp_Iti{i} = temp_It;
        temp_Ri{i} = temp_R;
        temp_Di{i} = temp_D;
        temp_YTi{i} = temp_YT;
        
    end
    
    
    for i = 1:LocationNum
        temp_S = temp_Si{i};
        temp_DL = temp_DLi{i};
        temp_L = temp_Li{i};
        temp_C = temp_Ci{i};
        temp_Ia = temp_Iai{i};
        temp_Iy = temp_Iyi{i};
        temp_It = temp_Iti{i};
        temp_R = temp_Ri{i};
        temp_D = temp_Di{i};
        temp_YT = temp_YTi{i};
        
        for a = 1:ageNum
            for r = 1:riskNum
                for v = 1%:vacNum
                    S(t,i,a,r,v) = temp_S(a,r,v);
                    DL(t,i,a,r,v) = temp_DL(a,r,v);
                    L(t,i,a,r,v) = temp_L(a,r,v);
                    C(t,i,a,r,v) = temp_C(a,r,v);
                    Ia(t,i,a,r,v) =  temp_Ia(a,r,v) ;
                    Iy(t,i,a,r,v) = temp_Iy(a,r,v);
                    It(t,i,a,r,v) = temp_It(a,r,v);
                    R(t,i,a,r,v) = temp_R(a,r,v);
                    D(t,i,a,r,v) = temp_D(a,r,v);
                    YT(t,i,a,r,v) = temp_YT(a,r,v);
                end
            end
        end
        
    end
    
    WeekNum_last = WeekNum;
    DayNum_last = DayNum;
end