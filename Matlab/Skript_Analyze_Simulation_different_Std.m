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

Std			= [0.2,0.3];
epsilons	= [ 0.02, 0.03]; % N
div_factors = [1.005, 10000];			 % m

K	= length(Std);
N	= length(epsilons);
M	= length(div_factors);

% prallocate

[m1_t, u, x, p_xt, n_xt, theta_t] = deal(cell(N,M,K));

for k = 1:K

base_folder = sprintf('../Python_Files/csv_files_MD_%0.1f/',Std(k));

m1_u_file	= 'Time_Dependent_Values.csv';
n_file		= 'Full_time_Distributions.csv';
p_file		= 'Full_time_transformed_Distribution.csv';
theta_file	= 'Full_time_theta.csv';


% 2) Load simulation data 

for n = 1:N		% epsilon
	for m = 1:M % div_factor
		this_folder_name = sprintf('Epsilon_%0.3f_DivFac_%0.3f',epsilons(n),div_factors(m));
		this_folder_name = strrep(this_folder_name,'Inf','inf');
% 		this_folder_name = sprintf('Epsilon_%0.3f_DivFac_%s',epsilons(n),div_names{m});

		% m1, u 
		data_file	= fullfile(base_folder,this_folder_name,m1_u_file);
		datatable	= readtable(data_file,'Delimiter',',');
		m1_cell		= datatable{:,2};
		m1_cell		= strrep(m1_cell,'j','i');
		m1_cell		= strrep(m1_cell,'(','');
		m1_cell		= strrep(m1_cell,')','');
		m1_t{n,m,k}	= str2double(m1_cell);
		u{n,m,k}		= table2array(datatable(2:end,1));
		
		% x, n_xt
		
		data_file	= fullfile(base_folder,this_folder_name,n_file);
		datatable	= readtable(data_file,'Delimiter',',','ReadVariableNames',0,'ReadRowNames',1);
		dataarray	= table2array(datatable)';
		x{n,m,k}		= dataarray(:,1);
		n_xt{n,m,k}	= dataarray(:,2:end);
		
		% x, p_xt,	
		data_file	= fullfile(base_folder,this_folder_name,p_file);
		datatable	= readtable(data_file,'Delimiter',',','ReadVariableNames',0,'ReadRowNames',1);
		dataarray	= table2array(datatable)';
		p_xt{n,m,k}	= dataarray(:,2:end);
		
		% theta 
		data_file	= fullfile(base_folder,this_folder_name,theta_file);
		datatable	= readtable(data_file,'Delimiter',',','ReadVariableNames',0,'ReadRowNames',1);
		theta_t{n,m,k}= table2array(datatable)';
		
	end
end
end
%% 3) Analyze the simulation data

a = 1;
b = 2;
t = 0:length(u{1});

grect = [20 20 1000 600];
fh2= figure('Color','w','Position',grect);
legstr = cell(N*M,1);
for k = 1:K
for n = 1:N
	for m = 1:M
% 		subplot(a,b,1)
% 		plot(t,[0;u{n,m}])
% 		xlabel('time in h')
% 		ylabel('input')
% 		hold on
		
		subplot(a,b,k);
		plot(t,abs(m1_t{n,m,k}))
		ylabel('length first moment')
		xlabel('time in h')
		title(sprintf('Std = %0.1f',Std(k)));
		
		hold on
% 		subplot(a,b,3)
% 		plot(t,angle(m1_t{n,m}))
% 		hold on
% 		ylabel('angle first moment')
% 		xlabel('time in h')
% 		
		legstr{sub2ind([M,N],m,n)} = sprintf('e = %0.2f, d = %0.3f',epsilons(n),div_factors(m));
		
	end
end
% 	subplot(a,b,2)
legend(legstr);
end

%% plot

grect = [20 20 1000 600];
fh3= figure('Color','w','Position',grect);

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

for i = 2:1:length(t)
for n = 1:N
	for m = 1:M
		p1(n,m).YData = p_xt{n,m}(:,i);
	
		p2(n,m).YData = n_xt{n,m}(:,i);

		p3(n,m).XData = theta_t{n,m}(1:take:end,i);
		pause(0.01)
	end
end	
end


