% Vacations consItered now:
% weekends, school opeing, Thanksgiving Break for November 25-27, 2009
% tNum_End = 20100331;
VacationRatioM = loadSchoolOpen(tNum_Begin, tNum_End, theta);
VacationRatioM = [VacationRatioM; VacationRatioM(2:end, :)];
VacationRatioM = VacationRatioM( ( datenum( num2str(tNum_Begin_sim),'yyyymmdd') +7*shiftWeek - datenum( num2str(tNum_Begin),'yyyymmdd')+1):end ,:);
% VacationRatioM: timeshift since the first day of tNum_Begin  x 217 locations
