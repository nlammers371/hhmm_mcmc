% Script to make illustrative inference figures for 2 state case
clear 
close all

addpath(genpath('../utilities'))

% make path to write figures
figPath = '../../fig/mcmc_validation/';
mkdir(figPath)

% specify path to read inference results from
resultsPath = 'C:\Users\nlamm\Dropbox (Personal)\hhmm_MCMC_data\mcmc_validation_basic_v3\';

% specify type of inference to show
nSteps = 7;
n_traces = 25;
n_reps = 2;
nStates = 3;
discrete_data_flag_vec = [0 1];
master_struct = struct;

% specify other key parameters
burn_in = 4e3;

for discrete_data_flag = 0:1
    readString = ['K' sprintf('%01d',nStates) '_W' sprintf('%03d',10*round(nSteps,1)) '_nt' sprintf('%03d',n_traces) ...
                      '_nrep' sprintf('%1d',n_reps) '_disc' sprintf('%1d',discrete_data_flag) ];

   

    % get list of matching inference files
    inf_files = dir([resultsPath '*' readString '*']);

    mcmc_results = struct;
    for i = 1:length(inf_files)
        load([inf_files(i).folder filesep inf_files(i).name],'mcmcInfo');
        fnames = fieldnames(mcmcInfo);
        for f = 1:length(fnames)
            mcmc_results(i).(fnames{f}) = mcmcInfo.(fnames{f});
        end
        clear mcmcInfo
    end
    master_struct(discrete_data_flag+1).mcmc_results = mcmc_results;
end

%% Make figures 
close all
n_chains_total = sum(mcmc_results(1).n_chains);
n_chains_keep = 125;
n_best = 1; % number of chains to use
rng(242);

trueParams = mcmc_results(1).trueParams;
tres = trueParams.tres / 60;
r_true_vec = trueParams.v / tres;

% obtain ML estimates 
for k = 1:length(master_struct)
    mcmc_results = master_struct(k).mcmc_results;
    
    v_inf_array = [];
    v_sort_indices = NaN(length(mcmc_results),nStates);
    ml_i_vec = NaN(1,length(mcmc_results));
    for i = 1:length(mcmc_results)        
        use_ids = randsample(1:n_chains_total, n_chains_keep, false);
        logL_vec = mean(mcmc_results(i).logL_vec(burn_in:end,use_ids),1);
        [~,si] = sort(logL_vec,'descend');
        ml_i_vec(i) = use_ids(si(1));
        for j = 1:n_best                                
            v_inf_array_best = mcmc_results(i).v_inf_array(burn_in:end,:,use_ids(si(j)));
            [~,si_v] = sort(nanmean(v_inf_array_best,1));
            v_sort_indices(i,:) = si_v;
            v_inf_array = [v_inf_array ; v_inf_array_best(:,si_v)];
        end
    end
    
    master_struct(k).v_inf_array = v_inf_array;
    master_struct(k).v_sort_indices = v_sort_indices;
    master_struct(k).ml_i_vec = ml_i_vec;
    
    master_struct(k).r_mean_vec = mean(v_inf_array,1) / tres;
    master_struct(k).r_ste_vec = std(v_inf_array,[],1) / tres;
end
% Emission rates

bkg_color = [228,221,209]/255;
markerSize = 100;

v_fig = figure('Position', [360   198   280   420]);
hold on
cmap = brewermap([],'Set2');
s(1) = scatter(1:3,r_true_vec,markerSize,'o','MarkerFaceColor',brighten(cmap(6,:),-0.5),'MarkerEdgeColor','k');
% results for rate simulations
errorbar((1:3)-0.1,master_struct(1).r_mean_vec,master_struct(1).r_ste_vec,'.','Color','k','LineWidth',1.5)
s(2) = scatter((1:3)-0.1,master_struct(1).r_mean_vec,markerSize*0.75,'s','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k');
% results for prob simulations
errorbar((1:3)+0.1,master_struct(2).r_mean_vec,master_struct(2).r_ste_vec,'.','Color','k','LineWidth',1.5)
s(3) = scatter((1:3)+0.1,master_struct(2).r_mean_vec,markerSize*0.75,'d','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k');


% xlabel('inference step')
ylabel('emission rate (au/min)')

set(gca,'Fontsize',14)
ylim([-0.5 15])
xlim([0.75 3.25])
% legend(s,'ground truth','estimates (continuous)','estimates (discrete)','Location','northwest')
set(gca,'Color',bkg_color) 
box on
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

v_fig.InvertHardcopy = 'off';
set(gcf,'color','w');
set(gca,'xtick',1:3)

saveas(v_fig,[figPath 'v_validation_fig.png'])
saveas(v_fig,[figPath 'v_validation_fig.pdf'])
%% Transition rates

% obtain ML estimates for R
for k = 1:length(master_struct)
    mcmc_results = master_struct(k).mcmc_results;
    
    % get A estimates first
    A_inf_array = [];        
    for i = 1:length(mcmc_results)                        
        ml_i = master_struct(k).ml_i_vec(i);
        si_v = master_struct(k).v_sort_indices(i,:);
        
        A_inf_array_best = mcmc_results(i).A_inf_array(si_v,si_v,burn_in:end,ml_i);        
        A_inf_array = cat(3, A_inf_array , A_inf_array_best);        
    end
    
    master_struct(k).A_inf_array = A_inf_array;
   
    % now attempt to convert to transition rates
    R_inf_array = NaN(size(A_inf_array));
    R_inf_array_alt = NaN(size(A_inf_array));
    wb = waitbar(0,'Conducting R fits...');
    for i = 1:size(R_inf_array,3)
        waitbar(i/size(R_inf_array,3),wb);
        
        A = A_inf_array(:,:,i);
        R_est = logm(A)/tres;
%         if ~isreal(R_est) || sum(R_est(:)<-0.01) > nStates           
%             out = prob_to_rate_fit_sym(A, tres*60, 'gen', 1/60, 10/60);                        
%             R_est = out.R_out*60;     
%         end
        % NL: using ad hoc adjustment based off of patterns in inference
        % results. Assigning negative fluo from 3->1 to 3->2
        R_est_alt = R_est;
        if R_est_alt(1,3) < 0
            R_est_alt(3,2) = R_est_alt(3,2) + R_est_alt(1,3);
            R_est_alt(1,3) = 0;
        end
        R_inf_array(:,:,i) = R_est;
        R_inf_array_alt(:,:,i) = R_est_alt;
    end
    master_struct(k).R_inf_array = R_inf_array;
    
    R_mean_array = mean(R_inf_array,3);
    R_mean_array(eye(3)==1) = 0;
    R_mean_array(eye(3)==1) = -sum(R_mean_array);
    master_struct(k).R_mean = R_mean_array;
    master_struct(k).R_ste = std(R_inf_array,[],3);
    
    R_mean_array_alt = mean(R_inf_array_alt,3);
    R_mean_array_alt(eye(3)==1) = 0;
    R_mean_array_alt(eye(3)==1) = -sum(R_mean_array_alt);
    master_struct(k).R_mean_alt = R_mean_array_alt;
    master_struct(k).R_ste_alt = std(R_inf_array_alt,[],3);
    
    delete(wb);
end

%%
R_true_array = trueParams.R * 60;
lin_inds = [2 3 4 6 7 8];

R_fig = figure('Position', [360 198 560 420]);
cmap = brewermap([],'Set2');
hold on
x_vec = 1:6;
s(1) = scatter(x_vec,R_true_array(lin_inds),markerSize,'o','MarkerFaceColor',brighten(cmap(6,:),-0.5),'MarkerEdgeColor','k');
% results for rate simulations
errorbar(x_vec-0.1,master_struct(1).R_mean(lin_inds),master_struct(1).R_ste(lin_inds),'.','Color','k','LineWidth',1.5)
% results for prob simulations
errorbar(x_vec+0.1,master_struct(2).R_mean(lin_inds),master_struct(2).R_ste(lin_inds),'.','Color','k','LineWidth',1.5)
% results for rate simulations
s(3) = scatter(x_vec+0.1,master_struct(2).R_mean(lin_inds),markerSize*0.75,'d','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k');
% results for prob simulations
s(2) = scatter(x_vec-0.1,master_struct(1).R_mean(lin_inds),markerSize*0.75,'s','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k');

ylabel('transition rate (events per minute)')

set(gca,'Fontsize',14)
ylim([-1 5])
xlim([0.75 6.25])
legend(s,'ground truth','estimates (continuous data)','estimates (discrete data)','Location','northwest')
set(gca,'Color',[228,221,209]/255) 
box on
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

R_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(R_fig,[figPath 'R_validation_fig.png'])
saveas(R_fig,[figPath 'R_validation_fig.pdf'])

%% Occupancy and dwell times
close all
for k = 1:length(master_struct)
    mcmc_results = master_struct(k).mcmc_results;
     
    % SS vector
    pi0 = master_struct(k).A_inf_array;
    pi0_inf_array = [];        
    for i = 1:length(mcmc_results)                        
        ml_i = master_struct(k).ml_i_vec(i);
        si_v = master_struct(k).v_sort_indices(i,:);
        
        pi0_inf_array_best = mcmc_results(i).pi0_inf_array(burn_in:end,si_v,ml_i);        
        pi0_inf_array = cat(1, pi0_inf_array , pi0_inf_array_best);        
    end
    
    master_struct(k).ss_array = pi0_inf_array;
    master_struct(k).ss_mean = mean(pi0_inf_array,1);
    master_struct(k).ss_ste = std(pi0_inf_array,[],1);        
    
    % Dwell times
    R_inf_array = master_struct(k).R_inf_array;
    tau_array = NaN(size(R_inf_array,3),3);
    for i = 1:size(R_inf_array,3)
        tau_array(i,:) = -1./diag(R_inf_array(:,:,i));
    end
    master_struct(k).dwell_mean = mean(tau_array,1);
    master_struct(k).tau_array = tau_array;
    master_struct(k).dwell_ste = std(tau_array,[],1);  
end


ss_true_vec = trueParams.pi0;

ss_fig = figure('Position', [360   198   280   420]);
cmap = brewermap([],'Set2');
hold on

s(1) = scatter(1:3,ss_true_vec,markerSize,'o','MarkerFaceColor',brighten(cmap(6,:),-0.5),'MarkerEdgeColor','k');
% results for rate simulations
errorbar((1:3)-0.1,master_struct(1).ss_mean,master_struct(1).ss_ste,'.','Color','k','LineWidth',1.5)
% results for prob simulations
errorbar((1:3)+0.1,master_struct(2).ss_mean,master_struct(2).ss_ste,'.','Color','k','LineWidth',1.5)
% results for rate simulations
s(3) = scatter((1:3)+0.1,master_struct(2).ss_mean,markerSize*0.75,'d','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k');
% results for prob simulations
s(2) = scatter((1:3)-0.1,master_struct(1).ss_mean,markerSize*0.75,'s','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k');

% xlabel('inference step')
ylabel('state probability')

set(gca,'Fontsize',14)
set(gca,'xtick',1:3)
ylim([0 0.6])
xlim([0.75 3.25])
% legend(s,'ground truth','estimates (continuous data)','estimates (discrete data)','Location','northwest')
set(gca,'Color',[228,221,209]/255) 
box on
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

ss_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(ss_fig,[figPath 'ss_validation_fig.png'])
saveas(ss_fig,[figPath 'ss_validation_fig.pdf'])

dwell_vec_true = -1./(diag(trueParams.R)*60);

dwell_time_fig = figure('Position', [360   198   280   420]);
cmap = brewermap([],'Set2');
hold on

s(1) = scatter(1:3,dwell_vec_true,markerSize,'o','MarkerFaceColor',brighten(cmap(6,:),-0.5),'MarkerEdgeColor','k');
% results for rate simulations
errorbar((1:3)-0.1,master_struct(1).dwell_mean,master_struct(1).dwell_ste,'.','Color','k','LineWidth',1.5)
% results for prob simulations
errorbar((1:3)+0.1,master_struct(2).dwell_mean,master_struct(2).dwell_ste,'.','Color','k','LineWidth',1.5)
% results for rate simulations
s(3) = scatter((1:3)+0.1,master_struct(2).dwell_mean,markerSize*0.75,'d','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k');
% results for prob simulations
s(2) = scatter((1:3)-0.1,master_struct(1).dwell_mean,markerSize*0.75,'s','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k');

% xlabel('inference step')
ylabel('state dwell time (minutes)')

set(gca,'Fontsize',14)
set(gca,'xtick',1:3)
ylim([0 1])
xlim([0.75 3.25])
% legend(s,'ground truth','estimates (continuous data)','estimates (discrete data)','Location','northeast')
set(gca,'Color',[228,221,209]/255) 
box on
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

ss_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(ss_fig,[figPath 'dwell_time_validation_fig.png'])
saveas(ss_fig,[figPath 'dwell_time_validation_fig.pdf'])

%% Validation for effective 2 state system
close all
for k = 1:length(master_struct)
    mcmc_results = master_struct(k).mcmc_results;
     
    % SS vector
    ss_array = master_struct(k).ss_array;
    R_array = master_struct(k).R_inf_array;
    r_inf_array = master_struct(k).v_inf_array / tres;
    
    effective_transition_rates = NaN(size(ss_array,1),2);
    r2_array = NaN(size(ss_array,1),2);
    for i = 1:size(ss_array,1)
        ss_vec = ss_array(i,:);
        % transition rates
        freq = -1./R_array(1,1,i);
        dur = freq * (1/ss_vec(1) - 1);
        effective_transition_rates(i,:) = [freq dur];
        % emission rates
        r = r_inf_array(i,:);

        r2 = (r(2) * ss_vec(2) + r(3) * ss_vec(3)) / (ss_vec(2)+ss_vec(3));
        r2_array(i,:) = [r(1) r2];
    end
    master_struct(k).r_eff_mean = mean(r2_array,1);
    master_struct(k).r_eff_ste = nanstd(r2_array,[],1);
    master_struct(k).R_eff_mean = mean(effective_transition_rates,1);
    master_struct(k).R_eff_ste = nanstd(effective_transition_rates,[],1);
end

R_true = trueParams.R*60;
freq_true = -1./R_true(1,1);
dur_true = freq_true * (1/ss_true_vec(1) - 1);
r2_true = (r_true_vec(2) * ss_true_vec(2) + r_true_vec(3) * ss_true_vec(3)) / (ss_true_vec(2)+ss_true_vec(3));


R_eff_fig = figure('Position', [360   198   280   420]);
hold on

s(1) = scatter(1:2,[freq_true dur_true],markerSize,'o','MarkerFaceColor',brighten(cmap(6,:),-0.5),'MarkerEdgeColor','k');
% results for rate simulations
errorbar((1:2)-0.05,master_struct(1).R_eff_mean,master_struct(1).R_eff_ste,'.','Color','k','LineWidth',1.5)
% results for prob simulations
errorbar((1:2)+0.05,master_struct(2).R_eff_mean,master_struct(2).R_eff_ste,'.','Color','k','LineWidth',1.5)
% results for rate simulations
s(3) = scatter((1:2)+0.05,master_struct(2).R_eff_mean,markerSize*0.75,'d','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k');
% results for prob simulations
s(2) = scatter((1:2)-0.05,master_struct(1).R_eff_mean,markerSize*0.75,'s','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k');

ylabel('transition rate (events per minute)')

set(gca,'Fontsize',14)
ylim([0 1.2])
xlim([0.75 2.25])
% legend(s,'ground truth','estimates (continuous data)','estimates (discrete data)','Location','northwest')
set(gca,'Color',[228,221,209]/255) 
box on
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

R_eff_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(R_eff_fig,[figPath 'R_eff_validation_fig.png'])
saveas(R_eff_fig,[figPath 'R_eff_validation_fig.pdf'])

r_eff_fig = figure('Position', [360 198 280 420]);
hold on

s(1) = scatter(1:2,[0 r2_true],markerSize,'o','MarkerFaceColor',brighten(cmap(6,:),-0.5),'MarkerEdgeColor','k');
% results for rate simulations
errorbar((1:2)-0.05,master_struct(1).r_eff_mean,master_struct(1).r_eff_ste,'.','Color','k','LineWidth',1.5)
% results for prob simulations
errorbar((1:2)+0.05,master_struct(2).r_eff_mean,master_struct(2).r_eff_ste,'.','Color','k','LineWidth',1.5)
% results for rate simulations
s(3) = scatter((1:2)+0.05,master_struct(2).r_eff_mean,markerSize*0.75,'d','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k');
% results for prob simulations
s(2) = scatter((1:2)-0.05,master_struct(1).r_eff_mean,markerSize*0.75,'s','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k');

ylabel('emission rate (au/min)')

set(gca,'Fontsize',14)
ylim([-0.25 8])
xlim([0.75 2.25])
% legend(s,'ground truth','estimates (continuous data)','estimates (discrete data)','Location','northwest')
set(gca,'Color',[228,221,209]/255) 
box on
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

r_eff_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(r_eff_fig,[figPath 'r2_eff_validation_fig.png'])
saveas(r_eff_fig,[figPath 'r2_eff_validation_fig.pdf'])