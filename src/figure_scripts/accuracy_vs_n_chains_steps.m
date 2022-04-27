% Script to make illustrative inference figures for 2 state case
clear 
close all

addpath(genpath('../utilities'))

% make path to write figures
figPath = '../../fig/mcmc_validation_low_noise/';
mkdir(figPath)

% specify path to read inference results from
resultsPath = 'C:\Users\nlamm\Dropbox (Personal)\hhmm_MCMC_data\mcmc_validation_basic_v2\';

% specify type of inference to show
nSteps = 7;
n_traces = 25;
n_reps = 2;
nStates = 3;
discrete_data_flag_vec = [0 1];
master_struct = struct;
for discrete_data_flag = 0:1
    readString = ['K' sprintf('%01d',nStates) '_W' sprintf('%03d',10*round(nSteps,1)) '_nt' sprintf('%03d',n_traces) ...
                      '_nrep' sprintf('%1d',n_reps) '_disc' sprintf('%1d',discrete_data_flag) ];

    % specify other key parameters
    burn_in = 4e3;

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

%% 


avg_window = 250;
start_vec = 1:250:4501;

close all
n_chains_total = sum(mcmc_results(1).n_chains);
n_chains_keep = 125;
n_best = 1; % number of chains to use
rng(242);

% trueParams = mcmc_results(1).trueParams;
% tres = trueParams.tres / 60;
% r_true_vec = trueParams.v / tres;

% obtain ML estimates 
for k = 1:length(master_struct)
    mcmc_results = master_struct(k).mcmc_results;
    
    master_struct(k).v_mean_array = NaN(n_chains_total,3,length(master_struct),length(start_vec));
    master_struct(k).logL_array = NaN(n_chains_total,length(master_struct),length(start_vec));
    master_struct(k).A_mean_array = NaN(3,3,n_chains_total,length(master_struct),length(start_vec));
       
    for s = 1:length(start_vec)
        start_i = start_vec(s);
        stop_i = start_vec(s)+avg_window-1;
        for i = 1:length(mcmc_results)        
            master_struct(k).A_mean_array(:,:,:,i,s) = permute(mean(mcmc_results(i).A_inf_array(:,:,start_i:stop_i,:),3),[1 2 4 3]);
            master_struct(k).v_mean_array(:,:,i,s) = permute(mean(mcmc_results(i).v_inf_array(start_i:stop_i,:,:),1),[3 2 1]);
            master_struct(k).logL_array(:,i,s) = mean(mcmc_results(i).logL_vec(start_i:stop_i,:),1);   
        end
    end
    
end


% Use bootstrap resampling to estimate model error as a function of the
% number of mcmc steps and the number of unique inference chains
trueParams = mcmc_results(1).trueParams;
v_true = trueParams.v;
v_true_norm = v_true / v_true(end);
A_true = trueParams.A;

nBoots = 100;
n_chain_vec = floor(linspace(1,125,16));

for k = 1:length(master_struct)
  
    A_deviation_array = NaN(length(start_vec),length(n_chain_vec),nBoots);
    v_deviation_array = NaN(length(start_vec),length(n_chain_vec),nBoots);
    total_deviation_array = NaN(length(start_vec),length(n_chain_vec),nBoots);
    
    mcmc_results = master_struct(k).mcmc_results;
    wb = waitbar(0,'Calculating model deviations...');    
    for s = 1:length(start_vec)
        waitbar(s/length(start_vec),wb);
        for c = 1:length(n_chain_vec)
            n_chains = n_chain_vec(c);
            for n = 1:nBoots
                A_array = NaN(3,3,length(mcmc_results));
                v_array = NaN(length(mcmc_results),3);
                for m = 1:length(mcmc_results)
                    % randomly select ids to use
                    use_ids = randsample(1:n_chains_total,n_chains,true);
                    [~,si] = sort(master_struct(k).logL_array(use_ids,m,s),'descend');
                    v_array(m,:) = master_struct(k).v_mean_array(use_ids(si(1)),:,m,s);
                    A_array(:,:,m) = master_struct(k).A_mean_array(:,:,use_ids(si(1)),m,s);
                   
                        
                end
                A_mean = mean(A_array,3);
                v_mean = mean(v_array,1);
                
                A_dev = sqrt(mean((A_true(:)-A_mean(:)).^2))./mean(A_true(:));
                A_deviation_array(s,c,n) = A_dev;
                v_dev = sqrt(mean((v_true(:)-v_mean(:)).^2))./mean(v_true(:));
                v_deviation_array(s,c,n) = v_dev;
                total_deviation_array(s,c,n) = sqrt([A_dev.^2 + v_dev.^2])/2;
%                 deviation_array(s,c,n) = 
            end
        end
    end
    delete(wb);
    master_struct(k).A_deviation_array = A_deviation_array;
    master_struct(k).v_deviation_array = v_deviation_array;
    master_struct(k).total_deviation_array = total_deviation_array;
end
                  
%% Make deviation plots
close all
bkg_color = [228,221,209]/255;
markerSize = 100;

for k = 1:2
    
    A_dev_mean = nanmean(master_struct(k).A_deviation_array,3);
    A_dev_ste = nanstd(master_struct(k).A_deviation_array,[],3);
    
    v_dev_mean = nanmean(master_struct(k).v_deviation_array,3);
    v_dev_ste = nanstd(master_struct(k).v_deviation_array,[],3);
    
    tot_dev_mean = nanmean(master_struct(k).total_deviation_array,3);
    tot_dev_ste = nanstd(master_struct(k).total_deviation_array,[],3);
  
    % A
    A_dev_fig = figure;
    hold on
    cmap = flipud(brewermap(length(n_chain_vec),'Spectral'));
    colormap(cmap);
    for c = 1:length(n_chain_vec)
        plot(start_vec+avg_window,A_dev_mean(:,c),'Color',cmap(c,:),'LineWidth',1.5);
    end

    h = colorbar;
    h.YTick = linspace(0,1,6);
    h.YTickLabel = n_chain_vec(1:3:end);
    % xlabel('inference step')
    ylabel('normalized A deviation')
    ylabel(h,'number of MCMC chains')
    xlabel('number of MCMC steps')
    
    set(gca,'Fontsize',14)
   
    set(gca,'Color',bkg_color) 
    box on
    grid on
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';
    xlim([0 5000])
    A_dev_fig.InvertHardcopy = 'off';
    set(gcf,'color','w');
%     set(gca,'xtick',1:3)
    
    if k == 1
        suffix = '_continuous';
    else
        suffix = '_discrete';
    end
    saveas(A_dev_fig,[figPath 'A_dev_' suffix '_nr' num2str(n_reps) '.png'])
    saveas(A_dev_fig,[figPath 'A_dev_' suffix '_nr' num2str(n_reps) '.pdf'])
    
    % V
    v_dev_fig = figure;
    hold on
    cmap = flipud(brewermap(length(n_chain_vec),'Spectral'));
    colormap(cmap);
    for c = 1:length(n_chain_vec)
        plot(start_vec+avg_window,v_dev_mean(:,c),'Color',cmap(c,:),'LineWidth',1.5);
    end

    h = colorbar;
    h.YTick = linspace(0,1,6);
    h.YTickLabel = n_chain_vec(1:3:end);
    % xlabel('inference step')
    ylabel('normalized v deviation')
    ylabel(h,'number of MCMC chains')
    xlabel('number of MCMC steps')
    
    set(gca,'Fontsize',14)
   
    set(gca,'Color',bkg_color) 
    box on
    grid on
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';

    v_dev_fig.InvertHardcopy = 'off';
    set(gcf,'color','w');
%     set(gca,'xtick',1:3)
    xlim([0 5000])
    if k == 1
        suffix = '_continuous';
    else
        suffix = '_discrete';
    end
    saveas(v_dev_fig,[figPath 'v_dev_' suffix '_nr' num2str(n_reps) '.png'])
    saveas(v_dev_fig,[figPath 'v_dev_' suffix '_nr' num2str(n_reps) '.pdf'])
    
    
    %%%%%%%%%%%%
    % Total
    tot_dev_fig = figure;
    hold on
    cmap = flipud(brewermap(length(n_chain_vec),'Spectral'));
    colormap(cmap);
    for c = 1:length(n_chain_vec)
        plot(start_vec+avg_window,tot_dev_mean(:,c),'Color',cmap(c,:),'LineWidth',1.5);
    end

    h = colorbar;
    h.YTick = linspace(0,1,6);
    h.YTickLabel = n_chain_vec(1:3:end);
    % xlabel('inference step')
    ylabel('normalized model error')
    ylabel(h,'number of MCMC chains')
    xlabel('number of MCMC steps')
    
    set(gca,'Fontsize',14)
   
    set(gca,'Color',bkg_color) 
    box on
    grid on
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.XAxis(1).Color = 'k';

    tot_dev_fig.InvertHardcopy = 'off';
    set(gcf,'color','w');
%     set(gca,'xtick',1:3)
    xlim([0 5000])
    if k == 1
        suffix = '_continuous';
    else
        suffix = '_discrete';
    end
    saveas(tot_dev_fig,[figPath 'tot_dev_' suffix '_nr' num2str(n_reps) '.png'])
    saveas(tot_dev_fig,[figPath 'tot_dev_' suffix '_nr' num2str(n_reps) '.pdf'])
end    