%% Skript to analyze Simulation results 

%% Structure
%
% 1) Set necessary variables and files settings
% 2) Load simulation data 
% 3) Analyze the simulation data
% 4) Plot and save important stuff
%
%% 

% 1) Set necessary variables and files settings

epsilons	= [0, 0.01, 0.02, 0.03, 0.04]; % N
div_factors = [1.005, 1.1];			 % m

N	= length(epsilons);
M	= length(div_factors);

% prallocate

[m1_t, u, x, p_xt, n_xt, theta_t] = deal(cell(N,M));

base_folder = '../Python_Files/csv_files_kk2/';

m1_u_file	= 'Time_Dependent_Values.csv';
n_file		= 'Full_time_Distributions.csv';
p_file		= 'Full_time_transformed_Distribution.csv';
theta_file	= 'Full_time_theta.csv';

% 2) Load simulation data 

for n = 1:N		% epsilon
	for m = 1:M % div_factor
		this_folder_name = sprintf('Epsilon_%0.3f_DivFac_%0.3f',epsilons(n),div_factors(m));

		% m1, u 
		data_file	= fullfile(base_folder,this_folder_name,m1_u_file);
		datatable	= readtable(data_file,'Delimiter',',');
		m1_cell		= datatable{:,2};
		m1_cell		= strrep(m1_cell,'j','i');
		m1_cell		= strrep(m1_cell,'(','');
		m1_cell		= strrep(m1_cell,')','');
		m1_t{n,m}	= str2double(m1_cell);
		u{n,m}		= table2array(datatable(2:end,1));
		
		% x, n_xt
		
		data_file	= fullfile(base_folder,this_folder_name,n_file);
		datatable	= readtable(data_file,'Delimiter',',','ReadVariableNames',0,'ReadRowNames',1);
		dataarray	= table2array(datatable)';
		x{n,m}		= dataarray(:,1);
		n_xt{n,m}	= dataarray(:,2:end);
		
		% x, p_xt,	
		data_file	= fullfile(base_folder,this_folder_name,p_file);
		datatable	= readtable(data_file,'Delimiter',',','ReadVariableNames',0,'ReadRowNames',1);
		dataarray	= table2array(datatable)';
		p_xt{n,m}	= dataarray(:,2:end);
		
		% theta 
		data_file	= fullfile(base_folder,this_folder_name,theta_file);
		datatable	= readtable(data_file,'Delimiter',',','ReadVariableNames',0,'ReadRowNames',1);
		theta_t{n,m}= table2array(datatable)';
		
	end
end

%% 3) Analyze the simulation data

a = 2;
b = 2;
t = 0:length(u{1});

grect = [20 20 1000 600];
fh2= figure('Color','w','Position',grect);

for n = 1:N
	for m = 1:M
		subplot(a,b,1)
		plot(t,[0;u{n,m}])
		hold on
		subplot(a,b,2)
		plot(t,abs(m1_t{n,m}))
		hold on
		subplot(a,b,4)
		plot(t,angle(m1_t{n,m}))
		hold on
		
	end
end
%% plot

grect = [20 20 1000 600];
fh2= figure('Color','w','Position',grect);

a = 1;
b = 3;
take = 50;
for n = 1:N
	for m = 1:M
		subplot(a,b,1)
		p1(n,m) = plot(x{n,m},p_xt{n,m}(:,1));
		ylim([0,0.5])
		hold on
		subplot(a,b,2)
		p2(n,m) = plot(x{n,m},n_xt{n,m}(:,1));
		ylim([0,0.5])
		hold on
		subplot(a,b,3)
		p3(n,m) = polarplot(theta_t{n,m}(1:take:end,1),ones(size(theta_t{n,m}(1:take:end,1)))+n/10,'o');
		hold on
		
	end
end

for i = 2:3:length(t)
for n = 1:N
	for m = 1:M
		p1(n,m).YData = p_xt{n,m}(:,i);
	
		p2(n,m).YData = n_xt{n,m}(:,i);

		p3(n,m).XData = theta_t{n,m}(1:take:end,i);
		pause(0.01)
	end
end	
end


