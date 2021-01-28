%% Formatting parameters
clc;
clear;
scaling_factor_epi = 1.0;
scaling_factor_rl = 1.0;

% Change folder
% cd 'C:\Users\remyp\Research\CDC Flu Pandemic\Results\2020 XX XX'

% Scenario names
n_scenarios = 32;
scenario_names = strings(1,n_scenarios);
start_week_list = zeros(1,n_scenarios);
first_week_offset_list = zeros(1,n_scenarios); % for MMWR weeks
first_weeks_flat = zeros(1,n_scenarios);

for i = 1:24
    scenario_names(i) = "HP" + sprintf('%02d',i);
    start_week_list(i) = 1; %36 in 2019
    first_week_offset_list(i) = 0;%2; % September 1st is a Tuesday
	scaling_factor(i) = scaling_factor_epi;
end

for i = 1:8
    scenario_names(24+i) = "RL" + sprintf('%02d',i);
    start_week_list(24+i) = 11; % September 14 is week 37, week 27 needed
    first_week_offset_list(24+i) = 0;%1; % September 14 is a Monday
    first_weeks_flat(24+i) = 1; % first 10 weeks equal to zero
	scaling_factor(24+i) = scaling_factor_rl;
end

% Outcomes
outcome_names = ["SymIllness", "Hosp", "Deaths", "AntiviralTX"];

% 100 indicates results are per 100
outcome_denominator = [100, 100000, 100000, 100];

% Model name
model_name = "UTA";

% Number of days in report
nb_days = 364; % 52  weeks
nb_groups = 6; % 5 age groups and total
nb_metrics = 4; % cases, hosp, deaths, antivirals


%% Parsing
% Get list of data files in the directory
dirData = dir('mats/*.mat');

% Number of simulations = number of files
nb_sim = length(dirData);

% Get population in each age group from first file
file_details = dirData(1);
file_path = fullfile(file_details.folder, file_details.name);
data_1 = load(file_path);


% To get the right age groups
groups_17 = ["0-0.5"; "0.5-4"; "5-9"; "10-14"; "15-19"; "20-24"; ...
    "25-29"; "30-34"; "35-39"; "40-44"; "45-49"; "50-54"; ...
    "55-59"; "60-64"; "65-69"; "70-74"; "75+"];
data_groups = ["0-4"; "5-17"; "18-24"; "25-64"; "65+"];
target_groups = ["0-4"; "5-17"; "18-49"; "50-64"; "65+"];
Pop_Metro17 = squeeze(sum(sum(sum(data_1.Pop_Metro17,4),3),1));

% To go from data groups to target groups
w_group = sum(Pop_Metro17(7:11)) / sum(Pop_Metro17(7:14)); % 25-49 over 25-64
weights_groups = [1, 0, 0, 0, 0; ...
                  0, 1, 0, 0, 0; ...
                  0, 0, 1, w_group, 0; ...
                  0, 0, 0, 1-w_group, 0; ...
                  0, 0, 0, 0, 1];
pop_vector_data = squeeze(sum(sum(sum(data_1.Pop_Metro,4),3),1));

% check the population still corresponds to data_groups listed above

pop_vector = zeros(6,1); % 6 age groups including total
pop_vector(1) = sum(Pop_Metro17); % Total
pop_vector(2) = sum(Pop_Metro17(1:2));
pop_vector(3) = sum(Pop_Metro17(3:4)) + 0.6*Pop_Metro17(5);
pop_vector(4) = 0.4*Pop_Metro17(5) + sum(Pop_Metro17(6:11));
pop_vector(5) = sum(Pop_Metro17(12:14));
pop_vector(6) = sum(Pop_Metro17(15:17));
pop_proportion = pop_vector(2:6) / sum(pop_vector(2:6)); % for demo run


% Loop through scenarios
for i_scen = 1:n_scenarios % i_scen = 1   i_scen = 25
    data_scenario_counts = zeros(nb_metrics,nb_days,nb_groups,nb_sim);
    data_scenario_norm = zeros(nb_metrics,nb_days,nb_groups,nb_sim);
    first_week_offset = first_week_offset_list(i_scen);
    
    % Get all simulations
    for i_file = 1:nb_sim % i_file = 1
        file_details = dirData(i_file);
        file_path = fullfile(file_details.folder, file_details.name);
        data_i = load(file_path);
        res = data_i.res4s; % 1x32 cell
    
        scen_data = res{1,i_scen}; % 1x5 cell
        sym_ill = scen_data{1,1}; % 370x5 matrix
        hosp = scen_data{1,2};
        deaths = scen_data{1,3};
        treated = scen_data{1,4};
        
        % Add simulation to data
        day1 = 1+first_week_offset; % MMWR start on Sundays
        data_end = nb_days-first_week_offset;
	
		for i_age = 2:nb_groups
			sf = scaling_factor(i_scen);
			w_i = weights_groups(i_age-1,:);
			data_scenario_counts(1,day1:end,i_age,i_file) = sym_ill(1:data_end,:) * w_i' * sf;
			data_scenario_counts(2,day1:end,i_age,i_file) = hosp(1:data_end,:) * w_i' * sf;
			data_scenario_counts(3,day1:end,i_age,i_file) = deaths(1:data_end,:) * w_i' * sf;
			data_scenario_counts(4,day1:end,i_age,i_file) = treated(1:data_end,:) * w_i' * sf;
		end
    end
    
    % Get total counts
    data_scenario_counts(:,:,1,:) = sum(data_scenario_counts(:,:,2:nb_groups,:),3);
    
    % Normalize based on population in age group and denominator of metric
    for i_age = 1:nb_groups
        pop_age_i = pop_vector(i_age);
        
        for i_outcome = 1:nb_metrics
            denominator_i = outcome_denominator(i_outcome);
            data_scenario_norm(i_outcome,:,i_age,:) = ...
                data_scenario_counts(i_outcome,:,i_age,:) * denominator_i / pop_age_i;
        end
    end
    
    % To get when weekly incidence goes above 1
    ppt_numbers = false;
    if ppt_numbers
        daily_in = median(data_scenario_norm(1,:,1,:),4);
        n_days = length(daily_in);
        weekly_in = zeros(n_days,1);
        weekly_in(1:6) = daily_in(1:6);
        for i_d = 7:n_days
            weekly_in(i_d) = sum(daily_in(i_d-6:i_d));
        end
        t0 = datetime(2009,9,1,0,0,0,'TimeZone','America/New_York');
        t0 + days(82) % Nov 21  for EPI1 ~1% 7-day incidence, schools closed

        % 25% of population vaccinated by ~January 7 with early
        jan7 = floor(days(datetime(2010,1,7,0,0,0,'TimeZone','America/New_York')-t0))+1;
        nb_cases = sum(daily_in(1:jan7));
        prop_cases = nb_cases/sum(daily_in);
    end
    
    % Call formatting function, with each metric separate
    start_week = start_week_list(i_scen);
    first_weeks_flat_scen = first_weeks_flat(i_scen);
    
    % Proportions
    for i_outcome = 1:nb_metrics % i_outcome = 1
        outcome_current = outcome_names(i_outcome);
        scenario_name_i = scenario_names(i_scen);
        out_filename = strcat("ExcelOutput/",scenario_name_i,"_",outcome_current,"_",model_name,".csv");
        incidence = squeeze(data_scenario_norm(i_outcome,:,:,:)); % daily
        scen_name_char = char(scenario_name_i);
        scenario_type = scen_name_char(1:2);
        toTemplate_function(1, nb_days, nb_sim, outcome_current, start_week, incidence, out_filename, first_weeks_flat_scen, scenario_type)
    end
    
    % Counts
    % for i_outcome = 1:nb_metrics % i_outcome = 1
        % outcome_current = outcome_names(i_outcome);
        % scenario_name_i = scenario_names(i_scen);
        % out_filename = strcat("ExcelOutputCounts/",scenario_name_i,"_",outcome_current,"_",model_name,"_Counts.csv");
        % incidence = squeeze(data_scenario_counts(i_outcome,:,:,:)); % daily
        % scen_name_char = char(scenario_name_i);
        % scenario_type = scen_name_char(1:2);
        % toTemplate_function(1, nb_days, nb_sim, outcome_current, start_week, incidence, out_filename, first_weeks_flat_scen, scenario_type)
    % end
end


