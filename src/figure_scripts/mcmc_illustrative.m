% Script to make illustrative inference figures for 2 state case
clear 
close all

addpath(genpath('../utilities'))

% make path to write figures
figPath = '../../fig/illustrative_mcmc/';
mkdir(figPath)

% initialize info structure
trueParams = setParamsBasic2state;
trueParams.sigma = 2;

%%%%%%%%%%%%%%%%%%%%% Simulated data %%%%%%%%%%%%%%%%
% basic inference params 
rng(458)

n_mcmc_steps = 1e3;
mcmcInfo.n_mcmc_steps = n_mcmc_steps; % number of MCMC steps (need to add convergence criteria)
mcmcInfo.burn_in = 100;
n_chains = 10;
mcmcInfo.n_chains = n_chains; % number of parallel MCMC chains to run
mcmcInfo.n_reps = 1; % number of chain state resampling passes per inference step

% characteristics of simulated data
mcmcInfo.n_traces = 1;
mcmcInfo.seq_length = 240; % length of simulated traces in time steps



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set MCMC options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nSteps = 5;
trueParams.nSteps = nSteps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trueParams.n_traces = mcmcInfo.n_traces;
trueParams.seq_length = mcmcInfo.seq_length;
trueParams = generateSimulatedData(trueParams);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize inference arrays and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mcmcInfo.nStates = trueParams.nStates;    
mcmcInfo.alpha_frac = trueParams.alpha_frac;
mcmcInfo.observed_fluo = trueParams.observed_fluo;

mcmcInfo = setMCMCOptions(mcmcInfo);

if ~mcmcInfo.inferNStepsFlag
    mcmcInfo.nSteps = trueParams.nSteps;
end    
mcmcInfo = initializeInferenceArrays(mcmcInfo);
mcmcInfo = initializeVariablesBasicRandom(mcmcInfo);

mcmcInfo.save_trace_results = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conduct inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
mcmcInfo = inferenceWrapper(mcmcInfo);
toc

%% Make figures 
close all
bkg_color = [228,221,209]/255;
[logL_sorted, si] = sort(mean(mcmcInfo.logL_vec(end-100:end,:),1),'descend');
err_sig = 1;
lb_prct = 100*normcdf(-err_sig,0,1);
ub_prct = 100*normcdf(err_sig,0,1);

% v and A inference over time
v_fig = figure;
cmap = brewermap([],'Set2');
hold on

plot(mcmcInfo.v_inf_array(:,1,si(1)),'Color',[cmap(3,:) 0.4],'LineWidth',3)
plot(mcmcInfo.v_inf_array(:,2,si(1)),'Color',[cmap(2,:) 0.4],'LineWidth',3)

plot(repelem(trueParams.v(2),n_mcmc_steps),'-.','Color',brighten(cmap(2,:),-0.5),'LineWidth',2)
plot(repelem(trueParams.v(1),n_mcmc_steps),'-.','Color',brighten(cmap(3,:),-0.5),'LineWidth',2)

xlabel('inference step')
ylabel('emission rate (v)')

set(gca,'Fontsize',14)
ylim([-0.1 5])
set(gca,'Color',bkg_color) 
box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

v_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(v_fig,[figPath 'v2_fig.png'])
saveas(v_fig,[figPath 'v2_fig.pdf'])

a_fig = figure;
cmap = brewermap([],'Set2');
hold on

plot(reshape(mcmcInfo.A_inf_array(1,2,:,si(1)),1,[]),'Color',[cmap(3,:) 0.4],'LineWidth',3)
plot(reshape(mcmcInfo.A_inf_array(2,1,:,si(1)),1,[]),'Color',[cmap(2,:) 0.4],'LineWidth',3)

plot(repelem(trueParams.A(2,1),n_mcmc_steps),'-.','Color',brighten(cmap(2,:),-0.5),'LineWidth',2)
plot(repelem(trueParams.A(1,2),n_mcmc_steps),'-.','Color',brighten(cmap(3,:),-0.5),'LineWidth',2)

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

saveas(a_fig,[figPath 'a2_fig.png'])
saveas(a_fig,[figPath 'a2_fig.pdf'])

%% Sigma plot
sigma_fig = figure;
cmap = brewermap([],'Set2');
hold on

plot(mcmcInfo.sigma_inf_array(:,si(1)),'Color',[cmap(8,:) 0.6],'LineWidth',3)
plot(repelem(trueParams.sigma,n_mcmc_steps),'-.','Color',brighten(cmap(8,:),-0.5),'LineWidth',2)


xlabel('inference step')
ylabel('fluorescence noise (\sigma)')

set(gca,'Fontsize',14)
ylim([1 3])
set(gca,'Color',[228,221,209]/255) 
box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

sigma_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(sigma_fig,[figPath 'sigma_fig.png'])
saveas(sigma_fig,[figPath 'sigma_fig.pdf'])

%% Make trace plots
close all

sample_fluo_array = permute(mcmcInfo.sample_fluo_inf_array(:,si(1),:,:),[1 3 4 2]);
% plot_frames = [1:10 10:10:1e3];
% 
time_axis = (0:mcmcInfo.seq_length-1)*trueParams.tres/60;
% 
% fluo_fit_fig = figure;
% hold on
% cmap2 = (brewermap(length(plot_frames)+2,'Reds'));
% 
% plot(time_axis, mcmcInfo.observed_fluo(:,1),'-k','LineWidth',3);
% 
% for p = 1:length(plot_frames)
%     plot(time_axis, sample_fluo_array(:,1,plot_frames(p)),'Color',[cmap2(p,:) 0.75], 'LineWidth',1.5);
% end    
%   
% 
% 
% xlabel('time (minutes)')
% ylabel('fluorescence (au)')
% 
% set(gca,'Fontsize',14)
% ylim([-0.5 20])
% xlim([0 40])
% % set(gca,'Color',[228,221,209]/255) 
% box on
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.XAxis(1).Color = 'k';
% 
% fluo_fit_fig.InvertHardcopy = 'off';
% set(gcf,'color','w');
% 
% saveas(fluo_fit_fig,[figPath 'trace_fit_fig.png'])
% saveas(fluo_fit_fig,[figPath 'trace_fit_fig.pdf'])

% get average predictions and bounds
burn_in = 250;
close all

fluo_pd_mean = nanmean(sample_fluo_array(:,:,burn_in:end),3);
fluo_95 = prctile(sample_fluo_array(:,:,burn_in:end),95,3);
fluo_05 = prctile(sample_fluo_array(:,:,burn_in:end),5,3);

fluo_fit_fig2 = figure('Position',[100 100 768 256]);
hold on
cmap3 = (brewermap([],'set2'));

fill([time_axis fliplr(time_axis)], [fluo_95' fliplr(fluo_05')], brighten(cmap3(5,:),0.35),'FaceAlpha',0.5,'EdgeAlpha',0.5,'EdgeColor','k')
p1 = plot(time_axis, mcmcInfo.observed_fluo(:,1),'-k','LineWidth',2.5);
p2 = plot(time_axis, fluo_pd_mean,'-','Color',brighten(cmap3(5,:),-0.5),'LineWidth',1.5);

  
xlabel('time (minutes)')
ylabel('fluorescence (au)')

set(gca,'Fontsize',14)
ylim([-1.5 20])
xlim([0 40])
legend([p1 p2], 'simulation', 'fit')
% set(gca,'Color',[228,221,209]/255) 
box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

fluo_fit_fig2.InvertHardcopy = 'off';
set(gcf,'color','w');


saveas(fluo_fit_fig2,[figPath 'trace_fit_fig2.png'])
saveas(fluo_fit_fig2,[figPath 'trace_fit_fig2.pdf'])


%% Look at naive state decoding
close all
sample_state_array = permute(mcmcInfo.sample_states_inf_array(:,si(1),:,:),[1 3 4 2])-1;
% gind the most likely iteration within selected chain 
[~,si_chain] = sort(mcmcInfo.logL_vec(:,si(1)),'descend');
state_pd_mean = nanmean(sample_state_array(:,:,burn_in:end),3);
state_95 = prctile(sample_state_array(:,:,burn_in:end),95,3);
state_05 = prctile(sample_state_array(:,:,burn_in:end),5,3);

state_pd_ml = sample_state_array(:,:,si_chain(1));

true_state_vec = trueParams.true_states-1;

state_fit_fig = figure('Position',[100 100 768 256]);
cmap4 = brewermap([],'set1');
hold on

% fill([time_axis fliplr(time_axis)], [state_95' fliplr(state_05')], brighten(cmap3(2,:),0.35),'FaceAlpha',1,'EdgeAlpha',0.5,'EdgeColor','k')

sh = stairs(time_axis, state_pd_ml,'-k','LineWidth',2);
bottom = 0;
x = [sh.XData(1),repelem(sh.XData(2:end),2)];
y = [repelem(sh.YData(1:end-1),2),sh.YData(end)];
f = fill([x,fliplr(x)],[y,bottom*ones(size(y))], cmap_s2(5,:),'FaceAlpha',0.65,'EdgeAlpha',0);
delete(sh);

s2 = stairs(time_axis, true_state_vec,'-','Color','k','LineWidth',1.5);
% s1 = stairs(time_axis, state_pd_mean,'-.','Color','k','LineWidth',1.5);

xlabel('time (minutes)')
ylabel('promoter state')

set(gca,'Fontsize',14)
ylim([-0.05 1.1])
xlim([0 40])
legend([s2,f],'simulation','fit','Location','southeast')
% legend([f,s1,s2],'simulation','fit (posterior)','fit (ML)','Location','northwest')
% set(gca,'Color',[228,221,209]/255) 
box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

state_fit_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(state_fit_fig,[figPath 'state_fit_fig.png'])
saveas(state_fit_fig,[figPath 'state_fit_fig.pdf'])

%% Make bivariate plots (just kon
% cmap_gra_short = brewermap(8,'Greys');
p_dim = 5;
param_min_vec = [-0.1 3.25 0.15 0.35 1.5];% -5  5  0.75];
param_max_vec = [ 0.1 4.5 0.35 0.75 2.5];%  0  6.5  0.77];
true_val_vec = [trueParams.v' trueParams.A(2,1) trueParams.A(1,2) trueParams.sigma];
% tickCell = {-10:2:-4, 3:5, 3:4};%, -5:2:0, 5:0.5:6.5,[.75 .76 0.77]};
labelCell = {'v_1', 'v_2', 'a_{21}', 'a_{12}','\sigma' };%, 'H (k_{off})', 'k_{off} max','r'};
close all
% colormap_cell = {cmap_bu,cmap_gra_short,cmap_bu,cmap_gr,cmap_gr,cmap_rd};

% concatneate results into a master array
master_param_array = [mcmcInfo.v_inf_array(:,1,si(1)) mcmcInfo.v_inf_array(:,2,si(1))...
            reshape(mcmcInfo.A_inf_array(2,1,:,si(1)),[],1) reshape(mcmcInfo.A_inf_array(1,2,:,si(1)),[],1) mcmcInfo.sigma_inf_array(:,si(1))];
          
grid_size = 75;

hist_fig = figure;
cmap_gra = brewermap(128,'Greys');
cmap_s2 = brewermap([],'Set2');
cmap_hist = [cmap_s2(3,:) ; cmap_s2(2,:) ; cmap_s2(2,:) ; cmap_s2(3,:) ; cmap_s2(8,:)];
cmap_gra(1,:) = bkg_color;

for i = 1:p_dim
    param1_i = i;
    param1_min = param_min_vec(i);%floor(10*prctile(master_param_array(:,param1_i),1))/10;
    param1_max = param_max_vec(i);%ceil(10*prctile(master_param_array(:,param1_i),99))/10;
    x_axis = linspace(param1_min,param1_max,grid_size);
    for j = i+1:p_dim
              
        param2_i = j;       
        param2_min = param_min_vec(j);%floor(10*prctile(master_param_array(:,param2_i),0.1))/10;
        param2_max = param_max_vec(j);%ceil(10*prctile(master_param_array(:,param2_i),99.9))/10;
       
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
saveas(hist_fig,[figPath 'mcmc_results_kon.png'])
saveas(hist_fig,[figPath 'mcmc_results_kon.pdf'])

%% Make heatmap of promoter states over time
plot_array = permute(sample_state_array,[3 1 2]);
close all

hm_fig = figure;
cmap_cm = brewermap([],'YlOrRd');
colormap(cmap_cm);
p = pcolor(imgaussfilt(plot_array,1));
% p = pcolor(flipud(plot_array));
p.EdgeAlpha = 0.01;
p.FaceAlpha = 0.5;

ylabel('MCMC inference step')
xlabel('time')
set(gca,'Fontsize',14)
set(gca,'YDir','reverse')

h = colorbar;
ylabel(h,'promoter state')

saveas(hm_fig,[figPath 'promoter_state_hm.png'])
saveas(hm_fig,[figPath 'promoter_state_hm.pdf'])

%% Make illustrative sampling plots
red = brighten([232 177 157]/256,-0.25);
green = brighten([220 236 203]/256,-0.25);
blue = brighten([169 191 227]/256,-0.25);

close all

promoter_state_hm = figure('Position',[100 100 1024 128]);
cmap_s2 = brewermap([],'Set2');
colormap([blue ; green]);
pcolor(repmat(state_pd_ml(1:120)',2,1))

% set(gca,'xticklabels',[round(20/3),round(40/3), 20, round(80/3), round(100/3), 40])
xlabel('time steps')
h = colorbar;
ylabel(h,'promoter state')
set(gca,'Fontsize',14)
set(gca,'ytick',[]);

saveas(promoter_state_hm,[figPath 'promoter_state_series_hm.png'])
saveas(promoter_state_hm,[figPath 'promoter_state_series_hm.pdf'])

%%
close all
toggle_ind = 79;
e_vec_1 = v_true(true_state_vec+1);
e_vec_1(toggle_ind) = v_true(1);
e_vec_2 = v_true(true_state_vec+1);
e_vec_2(toggle_ind) = v_true(2);

v_true = trueParams.v;

e_fig = figure('Position',[100 100 1024 128]);
hold on

% fill([time_axis fliplr(time_axis)], [state_95' fliplr(state_05')], brighten(cmap3(2,:),0.35),'FaceAlpha',1,'EdgeAlpha',0.5,'EdgeColor','k')
stairs(1:length(state_pd_ml), e_vec_2,'-','Color',green,'LineWidth',2);
stairs(1:length(state_pd_ml), e_vec_1,'-.','Color',blue,'LineWidth',2);
% bottom = 0;
% x = [sh.XData(1),repelem(sh.XData(2:end),2)];
% y = [repelem(sh.YData(1:end-1),2),sh.YData(end)];
% f = fill([x,fliplr(x)],[y,bottom*ones(size(y))], green,'FaceAlpha',1,'EdgeAlpha',0);
% delete(sh);

% s2 = stairs(1:length(state_pd_ml), state_pd_ml,'-','Color','k','LineWidth',1.5);
% s1 = stairs(time_axis, state_pd_mean,'-.','Color','k','LineWidth',1.5);

xlabel('time steps')
ylabel('emission rate (e)')
% ylabel(vertcat({'emission'},{'rate (e)'}))

set(gca,'Fontsize',14)
ylim([-0.2 4.3])
xlim([0 120])
% legend([f,s2],'simulation','fit')
% legend([f,s1,s2],'simulation','fit (posterior)','fit (ML)','Location','northwest')
% set(gca,'Color',[228,221,209]/255) 
box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

e_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(e_fig,[figPath 'e_fig.png'])
saveas(e_fig,[figPath 'e_fig.pdf'])

%% Now look at what differences this leads to in fluo
coeff_MS2 = mcmcInfo.coeff_MS2;     
f_vec_1 = fluo_conv_fun(e_vec_1,coeff_MS2);
f_vec_2 = fluo_conv_fun(e_vec_2,coeff_MS2);

close all

f_fig = figure('Position',[100 100 1024 412]);
hold on

% fill([time_axis fliplr(time_axis)], [state_95' fliplr(state_05')], brighten(cmap3(2,:),0.35),'FaceAlpha',1,'EdgeAlpha',0.5,'EdgeColor','k')
p1 = plot(1:length(state_pd_ml), mcmcInfo.observed_fluo(:,1),'-k','LineWidth',2);
plot(1:length(state_pd_ml), f_vec_2,'-','Color',green,'LineWidth',3);
plot(1:length(state_pd_ml), f_vec_1,'-.','Color',blue,'LineWidth',3);


xlabel('time steps')
ylabel('emission rate (e)')
% ylabel(vertcat({'emission'},{'rate (e)'}))

set(gca,'Fontsize',14)
ylim([-3 20])
xlim([0 120])

box on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';

f_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(f_fig,[figPath 'f_fig.png'])
saveas(f_fig,[figPath 'f_fig.pdf'])



