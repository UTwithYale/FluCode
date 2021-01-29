function Iasd_Single7 = aggDayToWeek(Iasd_Single, hourlyPerD)
Iasd_Single7 = [];
for i=1: (7*hourlyPerD):length(Iasd_Single)
    Iasd_Single7=[Iasd_Single7; sum( Iasd_Single(i:(i+7*hourlyPerD-1)) )];
end