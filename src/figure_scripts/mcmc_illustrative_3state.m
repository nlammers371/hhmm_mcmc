% Script to make illustrative inference figures for 2 state case
clear 
close all

addpath(genpath('../utilities'))

% make path to write figures
figPath = '../../fig/mcmc_validation/';
mkdir(figPath)

% specify path to read inference results from
resultsPath = 'C:\Users\nlamm\Dropbox (Personal)\hhmm_MCMC_data\mcmc_validation_basic_v2\';

% specify type of inference to show
discrete_data_flag = 1;
nSteps = 7;
n_traces = 25;
n_reps = 2;
nStates = 3;

readString = ['K' sprintf('%01d',nStates) '_W' sprintf('%03d',10*round(nSteps,1)) '_nt' sprintf('%03d',n_traces) ...
                  '_nrep' sprintf('%1d',n_reps) '_disc' sprintf('%1d',discrete_data_flag) ];
        
% specify other key parameters
burn_in = 4e3;

% get list of matching inference files
inf_files = dir([resultsPath '*' readString '*']);

mcmc_results = struct;
for i = 1%:length(inf_files)
    load([inf_files(i).folder filesep inf_files(i).name],'mcmcInfo');
%     fnames = fieldnames(mcmcInfo);
%     for f = 1:length(fnames)
%         mcmc_results(i).(fnames{f}) = mcmcInfo.(fnames{f});
%     end
%     clear mcmcInfo
end

%% Make figures 
close all
n_mcmc_steps = mcmcInfo.n_mcmc_steps(1);

bkg_color = [228,221,209]/255;
[logL_sorted, si] = sort(mean(mcmcInfo.logL_vec(burn_in:end,:),1),'descend');
err_sig = 1;
lb_prct = 100*normcdf(-err_sig,0,1);
ub_prct = 100*normcdf(err_sig,0,1);
trueParams = mcmcInfo.trueParams;

% v and A inference over time
v_fig = figure;
cmap = brewermap([],'Set2');
hold on

plot(mcmcInfo.v_inf_array(:,1,si(1)),'Color',[cmap(3,:) 0.4],'LineWidth',4)
plot(mcmcInfo.v_inf_array(:,2,si(1)),'Color',[cmap(5,:) 0.4],'LineWidth',4)
plot(mcmcInfo.v_inf_array(:,3,si(1)),'Color',[cmap(2,:) 0.4],'LineWidth',4)

plot(repelem(trueParams.v(3),n_mcmc_steps),'-.','Color',brighten(cmap(2,:),-0.5),'LineWidth',2)
plot(repelem(trueParams.v(2),n_mcmc_steps),'-.','Color',brighten(cmap(5,:),-0.5),'LineWidth',2)
plot(repelem(trueParams.v(1),n_mcmc_steps),'-.','Color',brighten(cmap(3,:),-0.5),'LineWidth',2)

xlabel('inference step')
ylabel('emission rate (v)')

set(gca,'Fontsize',14)
ylim([-0.1 4.5])
set(gca,'Color',bkg_color) 
box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

v_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(v_fig,[figPath 'v3_fig.png'])
saveas(v_fig,[figPath 'v3_fig.pdf'])

%%
a_fig = figure;
cmap = brewermap([],'Set2');
hold on
lineAlpha = 0.5;
% out of state 1
plot(reshape(mcmcInfo.A_inf_array(3,1,:,si(1)),1,[]),'Color',[cmap(3,:) lineAlpha],'LineWidth',2)
plot(reshape(mcmcInfo.A_inf_array(2,1,:,si(1)),1,[]),'Color',[cmap(3,:) lineAlpha],'LineWidth',2)

% out of state 3
plot(reshape(mcmcInfo.A_inf_array(2,3,:,si(1)),1,[]),'Color',[cmap(2,:) lineAlpha],'LineWidth',2)
plot(reshape(mcmcInfo.A_inf_array(1,3,:,si(1)),1,[]),'Color',[cmap(2,:) lineAlpha],'LineWidth',2)

% out of state 2
plot(reshape(mcmcInfo.A_inf_array(3,2,:,si(1)),1,[]),'Color',[cmap(5,:) lineAlpha],'LineWidth',2)
plot(reshape(mcmcInfo.A_inf_array(1,2,:,si(1)),1,[]),'Color',[cmap(5,:) lineAlpha],'LineWidth',2)

plot(repelem(trueParams.A(3,1),n_mcmc_steps),'-.','Color',brighten(cmap(3,:),-0.5),'LineWidth',2)
plot(repelem(trueParams.A(2,1),n_mcmc_steps),'-.','Color',brighten(cmap(3,:),-0.5),'LineWidth',2)

plot(repelem(trueParams.A(3,2),n_mcmc_steps),'-.','Color',brighten(cmap(5,:),-0.5),'LineWidth',2)
plot(repelem(trueParams.A(1,2),n_mcmc_steps),'-.','Color',brighten(cmap(5,:),-0.5),'LineWidth',2)

plot(repelem(trueParams.A(2,3),n_mcmc_steps),'-.','Color',brighten(cmap(2,:),-0.5),'LineWidth',2)
plot(repelem(trueParams.A(1,3),n_mcmc_steps),'-.','Color',brighten(cmap(2,:),-0.5),'LineWidth',2)

xlabel('inference step')
ylabel('transition probability (a)')

set(gca,'Fontsize',14)
ylim([0 0.75])
set(gca,'Color',[228,221,209]/255) 
box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

a_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(a_fig,[figPath 'a3_fig.png'])
saveas(a_fig,[figPath 'a3_fig.pdf'])

%% Make trace plots
close all
n_keep = 25;
sample_fluo_raw = reshape(mcmcInfo.sample_fluo(:,:,1:25:end),120,125);
sample_fluo_array = sample_fluo_raw(:,si(1:n_keep));
% plot_frames = [1:10 10:10:1e3];

time_axis = (0:mcmcInfo.seq_length-1)*trueParams.tres/60;

fluo_fit_fig = figure;
hold on
cmap2 = (brewermap(n_keep+2,'Reds'));

plot(time_axis, mcmcInfo.observed_fluo(:,1),'-k','LineWidth',3);

for p = 1:n_keep
    plot(time_axis, sample_fluo_array(:,p),'Color',[cmap2(p,:) 0.75], 'LineWidth',1.5);
end    
  


xlabel('time (minutes)')
ylabel('fluorescence (au)')

set(gca,'Fontsize',14)
ylim([-0.5 20])
xlim([0 40])
% set(gca,'Color',[228,221,209]/255) 
box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

fluo_fit_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(fluo_fit_fig,[figPath 'trace_fit_fig_3state.png'])
saveas(fluo_fit_fig,[figPath 'trace_fit_fig_3state.pdf'])

%% get average predictions and bounds
burn_in = 1;
close all

fluo_pd_mean = nanmean(sample_fluo_array,2);
fluo_pd_ste = nanstd(sample_fluo_array,[],2);
fluo_ub = fluo_pd_mean + fluo_pd_ste;
fluo_lb = fluo_pd_mean - fluo_pd_ste;

fluo_fit_fig2 = figure;
hold on
cmap3 = (brewermap(n_keep,'set2'));

fill([time_axis fliplr(time_axis)], [fluo_ub' fliplr(fluo_lb')], brighten(cmap3(5,:),0.35),'FaceAlpha',1,'EdgeAlpha',0,'EdgeColor','k')
p1 = plot(time_axis, mcmcInfo.observed_fluo(:,1),'-k','LineWidth',2);
p2 = plot(time_axis, fluo_pd_mean,'-','Color',brighten(cmap3(5,:),-0.5),'LineWidth',2);

  
xlabel('time (minutes)')
ylabel('fluorescence (au)')

set(gca,'Fontsize',14)
ylim([-1.5 20])
xlim([0 40])
legend([p1 p2], 'simulation', 'inference','Location','northwest')
% set(gca,'Color',[228,221,209]/255) 
box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

fluo_fit_fig2.InvertHardcopy = 'off';
set(gcf,'color','w');


saveas(fluo_fit_fig2,[figPath 'trace_fit_fig2_3state.png'])
saveas(fluo_fit_fig2,[figPath 'trace_fit_fig2_3state.pdf'])


%% Look at naive state decoding
sample_states_raw = reshape(mcmcInfo.sample_chains(:,:,1:25:end),120,125)-1;
sample_state_array = sample_states_raw(:,si(1:n_keep));

state_pd_mean = nanmean(sample_state_array,2);
state_pd_ml = sample_state_array(:,1);
% state_pd_ste = nanstd(sample_state_array,[],2);
% fluo_ub = state_pd_mean + state_pd_ste;
% fluo_lb = state_pd_mean - state_pd_ste;

true_state_vec = trueParams.true_states(:,1)-1;

state_fit_fig = figure;
hold on
cmap3 = brewermap([],'set2');
cmap4 = brewermap([],'set1');

% fill([time_axis fliplr(time_axis)], [state_95' fliplr(state_05')], brighten(cmap3(2,:),0.35),'FaceAlpha',1,'EdgeAlpha',0.5,'EdgeColor','k')

sh = stairs(time_axis, true_state_vec,'-k','LineWidth',2);
bottom = 0;
x = [sh.XData(1),repelem(sh.XData(2:end),2)];
y = [repelem(sh.YData(1:end-1),2),sh.YData(end)];
f = fill([x,fliplr(x)],[y,bottom*ones(size(y))], cmap3(6,:),'FaceAlpha',0.65,'EdgeAlpha',0);
delete(sh);

s2 = stairs(time_axis, state_pd_ml,'-','Color',cmap4(1,:),'LineWidth',1.5);
s1 = stairs(time_axis, state_pd_mean,'-.','Color','k','LineWidth',1.5);


xlabel('time (minutes)')
ylabel('promoter state')
set(gca,'ytick',[0 1 2])
set(gca,'Fontsize',14)
ylim([-0.05 2.3])
xlim([0 40])
legend([f,s1,s2],'simulation','fit (posterior)','fit (ML)','Location','northwest')
% set(gca,'Color',[228,221,209]/255) 
box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

state_fit_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(state_fit_fig,[figPath 'state_fit_fig_3state.png'])
saveas(state_fit_fig,[figPath 'state_fit_fig_3state.pdf'])

%% Make heatma
%% Make bivariate plots (just kon)
% cmap_gra_short = brewermap(8,'Greys');
burn_in = 4000;
p_dim = 9;
param_min_vec = [-0.1 3.5 0.15   0.4];% -5  5  0.75];
param_max_vec = [ 0.1 4.5 0.35   0.8];%  0  6.5  0.77];
true_val_vec = [trueParams.v' trueParams.A(2,1) trueParams.A(3,1)];
% tickCell = {-10:2:-4, 3:5, 3:4};%, -5:2:0, 5:0.5:6.5,[.75 .76 0.77]};
labelCell = {'v_1', 'v_2', 'v_3', 'a_{21}', 'a_{31}', 'a_{12}', 'a_{32}', 'a_{13}', 'a_{23}' };%, 'H (k_{off})', 'k_{off} max','r'};
close all
% colormap_cell = {cmap_bu,cmap_gra_short,cmap_bu,cmap_gr,cmap_gr,cmap_rd};

% concatneate results into a master array
master_param_array = [mcmcInfo.v_inf_array(burn_in:end,1:3,si(1))...
            reshape(mcmcInfo.A_inf_array(2,1,burn_in:end,si(1)),[],1) reshape(mcmcInfo.A_inf_array(3,1,burn_in:end,si(1)),[],1) ...
            reshape(mcmcInfo.A_inf_array(1,2,burn_in:end,si(1)),[],1) reshape(mcmcInfo.A_inf_array(3,2,burn_in:end,si(1)),[],1) ...
            reshape(mcmcInfo.A_inf_array(1,3,burn_in:end,si(1)),[],1) reshape(mcmcInfo.A_inf_array(2,3,burn_in:end,si(1)),[],1)];
          
grid_size = 25;

hist_fig = figure;
cmap_gra = brewermap(128,'Greys');
cmap_s2 = brewermap([],'Set2');
cmap_hist = cmap_gra;%[cmap_s2(3,:) ; cmap_s2(2,:) ; cmap_s2(2,:) ; cmap_s2(3,:)];
cmap_gra(1,:) = bkg_color;

for i = 1:p_dim
    param1_i = i;
    param1_min = floor(10*prctile(master_param_array(:,param1_i),1))/10;
    param1_max = ceil(10*prctile(master_param_array(:,param1_i),99))/10;
    x_axis = linspace(param1_min,param1_max,grid_size);
    for j = i+1:p_dim
              
        param2_i = j;       
        param2_min = floor(10*prctile(master_param_array(:,param2_i),0.1))/10;
        param2_max = ceil(10*prctile(master_param_array(:,param2_i),99.9))/10;
       
        y_axis = linspace(param2_min,param2_max,grid_size);
        [x_grid, y_grid] = meshgrid(x_axis,y_axis);
        [point_density, point_grid,bw] = ksdensity([master_param_array(:,param1_i) master_param_array(:,param2_i)],[x_grid(:),y_grid(:)]);
        dim = sqrt(length(point_density));
      
        % make tile
        s_ind = (j-1)*p_dim + i;
        subplot(p_dim,p_dim,s_ind);
        colormap(cmap_gra);
        
        % make contour plot
        [~,h] = contourf(x_grid, y_grid,reshape(point_density,dim,dim));                                
%         xlim(
        % turn off tick labels for axes not on outer edges
        if mod(s_ind,p_dim)~=1
            set(subplot(p_dim,p_dim,s_ind), 'YTickLabels', '')
        else
%             set(subplot(p_dim,p_dim,s_ind), 'YTick',tickCell{j})
            ylabel(labelCell{j})
        end
        
        if s_ind < p_dim*(p_dim-1)+1
            set(subplot(p_dim,p_dim,s_ind), 'XTickLabels', '')
        else
%             set(subplot(p_dim,p_dim,s_ind), 'XTick',tickCell{i})
            xlabel(labelCell{i})
        end
%         set(subplot(p_dim,p_dim,s_ind),'Color',bkg_color)
    end
    
    % make histogram plot
    s_ind = (i-1)*p_dim + i;
    subplot(p_dim,p_dim,s_ind);
    histogram(master_param_array(:,i),x_axis,'FaceColor',cmap_hist(i,:),'EdgeAlpha',0,'Normalization','probability')
    hold on
    % add info about optimal values, uncertainties, etc.
    param_true = true_val_vec(i);%mean(master_param_array(burn_in:end,i),1);
    param_mean = mean(master_param_array(:,i),1);
    param_ub = prctile(master_param_array(:,i),ub_prct,1);
    param_lb = prctile(master_param_array(:,i),lb_prct,1);
                
    y_lim = get(gca,'YLim');
    x_lim = get(gca,'XLim');
    
    fill([[param_ub, param_ub] fliplr([param_lb, param_lb])], [y_lim(1) y_lim(2) y_lim(2) y_lim(1)], 'k', 'FaceAlpha',0.15,'EdgeAlpha',0)
    plot([param_mean param_mean], [y_lim(1) y_lim(2)], '-.k', 'LineWidth',1)
    plot([param_true param_true], [y_lim(1) y_lim(2)], '-', 'LineWidth',1.5,'Color',cmap(5,:))
    
    % add mean value
    txt = ['$ \overline{' labelCell{i} '}$=' num2str(round(param_mean,1)) ];
    text(0.97*x_lim(1)+0.03*x_lim(2),0.85*y_lim(2),txt,'Interpreter','latex');
    
    if i ~= p_dim
        set(subplot(p_dim,p_dim,s_ind), 'XTickLabels', '')    
    else
%         set(subplot(p_dim,p_dim,s_ind), 'XTick',tickCell{i})
        xlabel(labelCell{i})
    end
    set(subplot(p_dim,p_dim,s_ind), 'YTickLabels', '')
    set(subplot(p_dim,p_dim,s_ind),'Color',bkg_color)
    
end    
hist_fig.InvertHardcopy = 'off';
hist_fig.Renderer = 'Painters';
set(gcf,'color','w');
saveas(hist_fig,[figPath 'mcmc_results_3state.png'])
saveas(hist_fig,[figPath 'mcmc_results_3state.pdf'])