function I_day = aggHourlyToDay(hourlyPerD, I_o)
I_day = [];
for i =1: (hourlyPerD):(length(I_o)-hourlyPerD)
    I_day=[I_day; sum( I_o(i:(i+hourlyPerD-1)) )];
end