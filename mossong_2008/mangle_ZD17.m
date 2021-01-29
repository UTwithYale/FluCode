% #!/usr/bin/python
% import os
% import numpy
% import pickle
% ageG{1} = 1; % 0-4
% ageG{2} = 2:4; % 5-9, 10-14, 15-19
% ageG{3} = 5:13; % 20-64
% ageG{4} = 14:15; % 65-69, 70+
% phis = mangle_ZD(ageG)
function phis = mangle_ZD17()

fileList = dir('mossong_2008/*_contact.dat');
for i=1:length(fileList)
    filename = strcat('mossong_2008/', fileList(i).name);
    temp = textread(filename);
    ConM{i} = temp;
end

fileList = dir('mossong_2008/*_population.dat');
for i=1:length(fileList)
    filename = strcat('mossong_2008/', fileList(i).name);
    temp = textread(filename);
    PopM{i} = temp;
end

%% Divide
% 0–4, 5-9, 
% ageG{1} = 1; % 0-4
% ageG{2} = 2:4; % 5-9, 10-14, 15-19
% ageG{3} = 5:13; % 20-64
% ageG{4} = 14:15; % 65-69, 70+

% 0-4, 5-17, 18-24, 25-64, 65+
% 0-4, 
% 5-9, 10-14, 15-17
% 18-19, 20-24, 
% 8: 25-29, 30-34, 35-39, 40-44, 45-49, 50-54, 55-69,60-64 
% 65-69, 70+

phis = 0;
for i= 1:length(PopM)
% for i= [5 5]
    tempCon = ConM{i};
    tempPop = PopM{i};
    
    tempPop4 = [];
    tempCon4 = [];
    
    tempPop4 = [tempPop(1) % 0-4
    sum(tempPop(2)+tempPop(3)+tempPop(4)*3/5) % 5-9, 10-14, 15-17
    sum(tempPop(4)*2/5+ tempPop(5)) %18-19, 20-24, 
    sum(tempPop(6:13)) % 8: 25-29, 30-34, 35-39, 40-44, 45-49, 50-54, 55-69,60-64 
    sum(tempPop(14:15)) % 65-69, 70+
    ];
    
    tempCon4 = [tempCon(:,1) ... % 0-4
    tempCon(:,2)+tempCon(:,3)+tempCon(:,4)*3/5  ... % 5-9, 10-14, 15-17
    tempCon(:,4)*2/5+ tempCon(:,5) ... %18-19, 20-24, 
    sum(tempCon(:,6:13),2)  ...% 8: 25-29, 30-34, 35-39, 40-44, 45-49, 50-54, 55-69,60-64 
    sum(tempCon(:,14:15),2) % 65-69, 70+
    ];
    
    tempCon5 = [tempCon4(1,:)  % 0-4
    tempCon4(2,:)*tempPop(2)/tempPop4(2)+tempCon4(3,:)*tempPop(3)/tempPop4(2)+tempCon4(4,:)*3/5*tempPop(4)/tempPop4(2)   % 5-9, 10-14, 15-17
    tempCon4(4,:)*2/5*tempPop(4)/tempPop4(3)+ tempCon4(5,:)*tempPop(5)/tempPop4(3)  %18-19, 20-24, 
    sum( tempCon4(6:13,:).*tempPop(6:13)/tempPop4(4) ,1)  % 8: 25-29, 30-34, 35-39, 40-44, 45-49, 50-54, 55-69,60-64 
    sum( tempCon4(14:15,:).*tempPop(14:15)/tempPop4(5),1) % 65-69, 70+
    ];

    phi0 = tempCon5./tempPop4;% * sum(tempPop4);
    phi0 = (phi0+phi0')/2;
    phis = phis+phi0;
end
phis = phis/length(PopM);
% 1;
% phis = phis/2;
% phis = phis/length(PopM);
