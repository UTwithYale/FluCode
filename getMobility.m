function Contact_Day = getMobility(tempMG)

load loadCDCData_Web_short.mat;

% MetroGNum
ContactGM = zeros(length(tempMG), length(tempMG));
for i=1:length(MetroGNum)
    temp1 = find(tempMG==MetroGNum(i,1));
    temp2 = find(tempMG==MetroGNum(i,2));
    if length(temp1)>0 & length(temp2)>0
        ContactGM(temp1, temp2) = ContactGM(temp1, temp2)+MetroGNum(i,3);
    end
end
% Contact_Day = ContactM;
% find(MetroGNum(:,1)==11260)

% ContactFM = zeros(length(tempMG), length(tempMG));
% for i=1:length(Metro_GroundNum)
%     temp1 = find(tempMG==Metro_GroundNum(i,1));
%     temp2 = find(tempMG==Metro_GroundNum(i,2));
%     ContactFM(temp1, temp2) = ContactFM(temp1, temp2)+Metro_GroundNum(i,3);
% end
% ceil(ContactFM*10/(30*3))

[NUM,TXT,RAW]=xlsread('Daily_Flows_2009Q4.xlsx');
ContactFM_Remy = NUM(2:end, 2:end);
Contact_Day = ContactGM+ceil(ContactFM_Remy);

% csvwrite('Workflow_CBSA.csv', ContactGM)
% csvwrite('Flightflow_CBSA.csv', ContactFM)