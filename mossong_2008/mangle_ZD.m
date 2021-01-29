% #!/usr/bin/python
% import os
% import numpy
% import pickle
% ageG{1} = 1; % 0-4
% ageG{2} = 2:4; % 5-9, 10-14, 15-19
% ageG{3} = 5:13; % 20-64
% ageG{4} = 14:15; % 65-69, 70+
% phis = mangle_ZD(ageG)
function phis = mangle_ZD(ageG)

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

%% Traditional
% phis = 0;
% for i=1:length(PopM)
%     tempCon = ConM{i};
%     tempPop = PopM{i};
%     phi0 = tempCon./tempPop;
%     phi0 = (phi0+phi0')/2;
%     phis = phis+phi0;
% end
% phi = phis/length(PopM);


%% Divide
% 0–4, 5-9, 
% ageG{1} = 1; % 0-4
% ageG{2} = 2:4; % 5-9, 10-14, 15-19
% ageG{3} = 5:13; % 20-64
% ageG{4} = 14:15; % 65-69, 70+
phis = 0;
for i=1:length(PopM)
    tempCon = ConM{i};
    tempPop = PopM{i};
    
    tempPop4 = [];
    tempCon4 = [];
    for k=1:length(ageG)
        tempG = ageG{k};
        tempPop4 = [tempPop4; sum(tempPop(tempG))];
        tempCon4 = [tempCon4 sum(tempCon(:, tempG),2)];
    end
    
    tempCon5= [];
    for j=1:length(ageG)
        tempG = ageG{j};
        if length(tempG)>1
            tempCon5 = [tempCon5; sum(tempCon4(tempG,:).*tempPop(tempG))/sum(tempPop(tempG))];
        else
            tempCon5 = [tempCon5; (tempCon4(tempG,:).*tempPop(tempG))/sum(tempPop(tempG))];
        end
    end
    phi0 = tempCon5./tempPop4;
    phi0 = (phi0+phi0')/2;
    phis = phis+phi0;
end
phis = phis/length(PopM);
%%
% # imbed in 17x17
% Phi = numpy.empty((17, 17))
% Phi[1 : 16, 1 : 16] = phi
% Phi[0, 1 : 16] = phi[0, :]
% Phi[1 : 16, 0] = phi[:, 0]
% Phi[-1, 1 : 16] = phi[-1, :]
% Phi[1 : 16, -1] = phi[:, -1]
% Phi[0, 0] = phi[0, 0]
% Phi[0, 16] = phi[0, 14]
% Phi[16, 0] = phi[14, 0]
% Phi[16, 16] = phi[14, 14]
