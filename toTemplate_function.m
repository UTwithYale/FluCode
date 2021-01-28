function toTemplate_function( INTERVAL_PER_DAY, TOTAL_TIME, NUM_SIM, metric, start_week, incidence, filename_out, first_weeks_flat_scen, scenario_type)
% first_weeks_flat_scen: sets the first weeks equal to zero if start_week not 1
% scenario_type: HP or RL
Ts = [];
AgeGroupNames = {'Overall', '0-4 yr', '5-17 yr', '18-49 yr',  '50-64 yr', '65+ yr'};

%% Bins depend on metric
switch metric
    case "SymIllness"
        bins = [0.05:0.1:19.95];
        bins_display = [0.0:0.1:19.9];
        cml_multiplier = 2;
        
    case "Hosp"
        switch scenario_type
            case "RL"
                bins = [0.05:0.1:59.95];
                bins_display = [0.0:0.1:59.9];
                cml_multiplier = 2;
            case "HP"
                bins = [0.125:0.25:149.875];
                bins_display = [0.0:0.25:149.75];
                cml_multiplier = 4;
        end
                
    case "Deaths"
        switch scenario_type
            case "RL"
                bins = [0.005:0.01:5.995];
                bins_display = [0.0:0.01:5.99];
                cml_multiplier = 2;
            case "HP"
                bins = [0.0125:0.025:14.9875];
                bins_display = [0.0:0.025:14.975];
                cml_multiplier = 4;
        end
        
    case "AntiviralTX"
        bins = [0.005:0.01:3.995];
        bins_display = [0.0:0.01:3.99];
        cml_multiplier = 2;
end

bin_cmls = cml_multiplier * bins;
bin_cmls_display = cml_multiplier * bins_display;
nb_bins = length(bins);

% Extra metrics
extra_bins = {'sm.mean';'sm.median';'sm.perc2p5';'sm.perc5';'sm.perc25';'sm.perc75';'sm.perc95';'sm.perc97p5';'sm.peak'};
n_extra_bins = length(extra_bins);

% Bin names to display
bin_name = cell(nb_bins+n_extra_bins,1);
bin_cml_name = cell(nb_bins+n_extra_bins,1);
for i = 1:nb_bins
    lb = sprintf('%g',bins_display(i));
    lb_cml = sprintf('%g',bin_cmls_display(i));
    if i < nb_bins
        ub = sprintf('%g',bins_display(i+1));
        ub_cml = sprintf('%g',bin_cmls_display(i+1));
    else
        ub = 'Inf';
        ub_cml = 'Inf';
    end
    
    bin_i = strcat('[', lb, '-', ub, ')');
    bin_i_cml = strcat('[', lb_cml, '-', ub_cml, ')');
    
    bin_name{i} = bin_i;
    bin_cml_name{i} = bin_i_cml;
end

for i = 1:n_extra_bins
    bin_name{nb_bins+i} = extra_bins{i};
    bin_cml_name{nb_bins+i} = extra_bins{i};
end


%% Parse
for agei=1:length(AgeGroupNames) % agei = 1
    % template
    OutputM = zeros(nb_bins, 52);
    
    IY_hours = [];
    for i_sim=1:NUM_SIM % number of simulations i_sim = 1
        temp = squeeze(incidence(:,agei,i_sim)); % row vector 
        IY_hours = [IY_hours; temp']; 
    end
    
    %% Get weekly data
    Week_begin = start_week;
    
    % If start_week is not 1 and we keep the weeks before that equal to
    % 0, do the adjustment here (effectively discard last weeks)
    if (start_week > 1) && (first_weeks_flat_scen > 0)
        Week_end = 52;% TOTAL_TIME/7 - 1; % 55 weeks bugs
    else
        Week_end = Week_begin + 52 -1;% TOTAL_TIME/7 - 1; % 55 weeks bugs
    end
    
    weekly_IY = zeros(NUM_SIM,52); % would be 53 if 53 MMWR that year
    for Wi = Week_begin:Week_end
        if Wi > 52
            Wi_report = Wi-52;
        else
            Wi_report = Wi;
        end
        
        week_start = (Wi-Week_begin)*7*INTERVAL_PER_DAY+1;
        week_end = (Wi-Week_begin+1)*7*INTERVAL_PER_DAY;
        temp = sum(IY_hours(:,week_start:week_end),2); % aggregates simulation time steps into weeks
        weekly_IY(:,Wi_report) = temp;
    end

    Week_end = Week_begin + 52 -1;% TOTAL_TIME/7 - 1; % 55 weeks bugs
    for Wi = Week_begin:Week_end
        if Wi > 52
            Wi_report = Wi-52;
        else
            Wi_report = Wi;
        end
        %Wi % might need to adjust for first partial week % MMWR weeks start on Sundays
        temp = weekly_IY(:, Wi_report); % aggregates simulation time steps into weeks
        [N,X]  = hist(temp, bins);
        OutputM(:,Wi_report) = N./sum(N);
    end 
    
    tempM = [mean(weekly_IY);quantile(weekly_IY, [0.5, 0.025, 0.05, 0.25, 0.75, 0.95, 0.975])];
    OutputM = [OutputM; tempM];
	
    % Based on weekly data
    Peak_M = [];
    Peak_Week = [];
    Peak_Magnitude = [];
    for i = 1:NUM_SIM
        [M,I] = max(weekly_IY(i,:));
        Peak_M = [Peak_M;  M];
        Peak_Week = [Peak_Week; I];
    end
    [N,X]  = hist(Peak_M, bins);
    Peak_Magnitude = N'./sum(N);
    
    tempM = [mean(Peak_M), quantile(Peak_M, [0.5, 0.025, 0.05, 0.25, 0.75, 0.95, 0.975])];
    Peak_Magnitude = [Peak_Magnitude; tempM'; NaN];
    
    %% Get probability peak is in a given week
    week_bins = [1: 1 : 52];
    [N,X]  = hist(Peak_Week, week_bins);
    Peak_Week_Dist = N'./sum(N);
    OutputM = [OutputM; Peak_Week_Dist'];
    
    
    %% Cumulative
    Cml_sum = sum(IY_hours');
    [N,X]  = hist(Cml_sum, bin_cmls);
    Cml = N'./sum(N);
    
    tempM = [mean(Cml_sum), quantile(Cml_sum, [0.5, 0.025, 0.05, 0.25, 0.75, 0.95, 0.975])];
    Cml = [Cml; tempM'; NaN];
    
    %%
    LocationCell = cell(length(Peak_Magnitude),1);
    LocationCell(1:length(LocationCell)) = {'US National'};
    Agegroup(1:length(LocationCell), 1) = AgeGroupNames(agei);
    bin_cml_name = bin_cml_name(1:length(LocationCell));
    bin_name = bin_name(1:length(LocationCell));
    temp_Data = [LocationCell, Agegroup, bin_cml_name, num2cell(Cml) bin_name,num2cell(Peak_Magnitude), num2cell(OutputM) ];%, Peak_Magnitude];
    VariableNames = ["Location", "Agegroup", "Bin_cml", "Cml", "Bin", "Peak.Magnitude", "Week1", "Week2", "Week3", "Week4", "Week5", "Week6", "Week7", "Week8", "Week9", "Week10", "Week11", "Week12", "Week13", "Week14", "Week15", "Week16", "Week17", "Week18", "Week19", "Week20", "Week21", "Week22", "Week23", "Week24", "Week25", "Week26", "Week27", "Week28", "Week29", "Week30", "Week31", "Week32", "Week33", "Week34", "Week35", "Week36", "Week37", "Week38", "Week39", "Week40", "Week41", "Week42", "Week43", "Week44", "Week45", "Week46", "Week47", "Week48", "Week49", "Week50", "Week51", "Week52"];
    T = cell2table(temp_Data,'VariableNames', VariableNames);
    Ts = [Ts; T];
end

writetable(Ts, filename_out)
