function [I100, CityAndMetro_U, tempMG, StateID, Pop_MetroAll, shiftWeek] = loadWeekLagFiles(iX, iTest)

fileName = strcat('weeklagfiles/Iasd_Single_os_217_weeklag0_all_iX_', num2str(iX), '_I0', num2str(iTest),'_*') ;
listing = dir(fileName);
load( strcat('weeklagfiles/', listing(1).name  ) )
I100 = Iasd_Single_os(:,1);