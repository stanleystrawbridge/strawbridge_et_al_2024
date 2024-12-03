classdef model_inference < handle
    
properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    data_directory      
    results_directory
    chimera_exclusion 
    chimera_exclusion_clustering = struct('groups',{[]});        
    chimera_crowding 
    embryo_litter_mate = struct('name',{cell(1)},...
        'stage',{cell(1)},...
        'stage_group',{0},...
        'time',{0},...
        'regime',{0},...
        'cell_type',{{'DN','DP','Epi','PrE','TE'}},...
        'cells',{[0 0 0 0 0]},...
        'embyro_ratio',{[]},...
        'icm_ratio',{[]});
    time = linspace(0,48,501);
    embryo_abc_priors = struct('alpha', [.055 0.07],... % Proliferation rate constant
     'beta', [0.01 1],...    % B -> T, basal rate constant
      'rho', [.3 0.38],...     % B -> C, lineage bias
     'zeta', [5e-4 1e-3],... % C -> E, rate constant [1e-5 1e-2]
      'eta', [9e-5 3e-4],... % C -> P, rate constant
        'l', [-0.5 2.5],...    % PrE signal feedback to ICM (??)
        'm', [-0.5 2.5]);      % FGF4 signal feedback from EPI and ICM to ICM
    embryo_abc_parameter_space
    embryo_abc_IC = [8 0 0 0 0];
    best_fit_parameters_embryo    
    chimera_abc_priors = struct('alpha_D', [.01 .035],... % Donor Cell proliferation rate
     'a', [.08 .2],...% spatial crowding constant
     'n', [-0.5 2.5]);   % Donor cell feedback to Blastomeres (crowding)
    chimera_abc_parameter_space
    chimera_abc_IC = [8 0 0 0 0  0 15;
                      8 0 0 0 0  0 10;
                      8 0 0 0 0  0  0;
                      8 0 0 0 0 10  0;
                      8 0 0 0 0 15  0;];
    best_fit_parameters_chimera

end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
function obj = model_inference() %=========================================
   
    close all

    obj.load_data
    obj.write_properties                   

end %======================================================================

    
function load_data(obj) %==================================================
    
    %LOAD_DATA imports data from DATA folder into the objects:
    %embryo_litter_mate
    
    % ADD data path and 'saiz-et-al_2016-master'
    obj.data_directory = fullfile(pwd,'data');

    % get the Raw Data from  Saiz et al 2015 
    model_inference.import_litter_mate_data(genpath(obj.data_directory))  

    % make a save path
    obj.results_directory = fullfile(fileparts(obj.data_directory),...
        'results');
    if ~exist(obj.results_directory,'dir')
        mkdir(obj.results_directory)
    end

end %======================================================================
        

function write_properties(obj) %===========================================           
    
    %WRITE_PROPERTIES write raw data from DATA folder into the
    %objects: embryo_litter_mate, chimera_exclusion, embryo_crowding
    % Keep the data path      
    addpath(obj.data_directory);
    %EMBRYO LITTERMATES----------------------------------------------------
    % 'obj.embryo_litter_mate' structure
    % *.name            str Unique Embryo Name
    % *.stage           str binned grouping for embryo size, 
    %                       [<32,32_64,64_90,90_120,>120] cells
    % *.stageGroup      int numerical values assigned to different 
    %                       binned grouping for the purpose of 
    %                       plotting. (1:5)
    % *.time            int Indicated time of treatment (hours)
    % *.cells           Nx5 The number of cells for each cell type
    %                       [Double Negative,Double Positive,
    %                       Epiblast,Primitive, Endoderm,
    %                       Trophectoderm]
    % *.embryo_ratio    Nx5 The proportion of the cell type to the 
    %                       entire embryo [Double Negative,Double 
    %                       Positive,Epiblast,Primitive Endoderm,
    %                       Trophectoderm]
    % *.icm_ratio       Nx4 The proportion of the cells of the
    %                       within the ICM relative to the ICM
    %                       [Double Negative,Double Positive,
    %                       Epiblast,Primitive Endoderm]
    
    load FGFallpooledtrans.mat FGFallpooledtrans           

    % Embryo_ID
    x = FGFallpooledtrans;
    clear FGFallpooledtrans
    
    x.Tt_length(isnan(x.Tt_length)) = 0;

    % Regime
    % idx  Value Explanation
    % 1  = 0,    littermate, Collected, fixed, and stained at collection
    % 2  = 1,          8 cell  + 48 hours 
    % 3  = 3,     64- 92 cells + 24 hours
    % 4  = 4,     90-120 cells + 20 hours
    % 5  = 5,     32- 64 cells + 30 hours
    % 6  = 6,    120-150 cells + 15 hours
    % 7  = 8,          8 cell  + 24 hours treatment + 24 hours Control
    % 8  = 9,          8 cell  + 30 hours treatment + 18 hours Control
    x.Regime(isnan(x.Regime)) = 0;
    regime = unique(x.Regime);
    regimeSet = [0:9];
    regimeIndex = (ismember(x.Regime, regimeSet));

    % Condition that embryos were in
    % 1 = 'AZD_1'
    % 2 = 'Control'
    % 3 = 'FGF42PD'
    % 4 = 'FGF4_1000'
    % 5 = 'Littermate' Collected, fixed, and stained at collection
    % 6 = 'PD03_1'
    % 7 = 'SU_10'
    % 8 = 'SU_20'
    condition = unique(x.Treatment);
    conditionSet = [1 2 3 4 6 7 8];
    conditionIndex = ~ismember(x.Treatment,condition(conditionSet)); 

    % X point
    % 1 = 'NA' 
    % 2 = 'ep' end point
    % 3 = 'sp' 'start point'
    % 4 = 'xp' exchange point'
    xPoint = unique(x.Xpoint);
    xPointNum = [1 2 3 4];
    xPointIndex = ismember(x.Xpoint,xPoint(xPointNum));
    
    % Experimenter
    % 1 = KW
    % 2 = NS
    experimenter = unique(x.Experimenter);
    experimenterSet = [1 2];
    experimenterIndex = ...
        ismember(x.Experimenter,experimenter(experimenterSet));

   % get the index of control embryos
    ctlIndex = ...
        [regimeIndex & conditionIndex & xPointIndex & experimenterIndex];

    embryoID = unique(x.Embryo_ID(ctlIndex))'; % embryo names
    nEmbryo = numel(embryoID); % number of embryos

    % get cell types [DN,DP,EPI,PRE,TE]
    cellType = unique(x.Identity);
    nCellType = numel(cellType);

    % get staging IDS ['<32','32_64','64_90','90_120','120-150','>150'];
    % also numericall ordered 1:5 for plotting purposes
    stages = {'<32';'32_64';'64_90';'90_120';'120_150';'>150'};
    stageSet = [ 32,64,90,120,150,inf];
    
    % step through all embryos
    for i = 1:nEmbryo
        % Set name of embryo i
        obj.embryo_litter_mate.name{i} = embryoID(i);

        % get the indext for all cells withing embryo i
        embryoIndex = x.Embryo_ID==embryoID(i);
    
        % Step through all cell types
        for j = 1:nCellType                   
            % Save number of cells for celltype j in embryp i
            obj.embryo_litter_mate.cells(i,j) = ...
                sum(x.Identitykm(embryoIndex)==cellType(j));                    
        end
        
        % stage ID
        cellNum =  sum(obj.embryo_litter_mate.cells(i,:));
        stageIndex = find(cellNum<stageSet,1,'first');
        
        obj.embryo_litter_mate.stage(i) = stages(stageIndex);    
        % numerical stage 
        obj.embryo_litter_mate.stage_group(i) = stageIndex;
        
        % Indicated time of treatment
        obj.embryo_litter_mate.time(i) = unique(x.Tt_length(embryoIndex));
        
        % Treatment Regime
        obj.embryo_litter_mate.regime(i) = unique(x.Regime(embryoIndex));
    end
    
    clear x
    
    %Calculate the ratios of the tissues-----------------------------------    
    % Whole Embryo Tissue Ratio [Double Negative,Double
    % Positive,Epiblast,Primitive Endoderm,Trophectoderm]
    obj.embryo_litter_mate.embyro_ratio = ...
        obj.embryo_litter_mate.cells./...
        repmat(sum(obj.embryo_litter_mate.cells,2),1,5); 
    % ICM Tissue Ratio
    % [Double Negative,Double Positive,Epiblast,Primitive Endoderm]
    obj.embryo_litter_mate.icm_ratio =...
        obj.embryo_litter_mate.cells(:,1:4)./...
       repmat(sum(obj.embryo_litter_mate.cells(:,1:4),2),1,4);           

    %CHIMERA_EXCLUSION-----------------------------------------------------
    fileName = fullfile('data','chimera_exclusion_data_file.txt');
    obj.chimera_exclusion = load(fileName);
    clear fileName
    
    %CHIMERA_CROWDING------------------------------------------------------
    fileName = fullfile('data','chimera_crowding_data_file.txt');
    obj.chimera_crowding = load(fileName);
    clear fileName

end %======================================================================
        
       
function analysis_chimera_exclusion_clustering(obj) %======================
    
    close all
    save_folder_name = fullfile(obj.results_directory,...
        'chimera_exclusion_clustering');
    if ~exist(save_folder_name,'dir')
        mkdir(save_folder_name)
    end    
    
    % Get index of different conditions
    noninjected_index = obj.chimera_exclusion(:,1)==2;
    wt_injected_index = obj.chimera_exclusion(:,1)==3;
    null_injected_index = obj.chimera_exclusion(:,1)==1;

    % get the number of host derived cells from each embryo for each
    % injection type
    noninjected   = obj.chimera_exclusion(noninjected_index,[4 5]);
    wt_injected   = obj.chimera_exclusion(wt_injected_index,[4 5]);
    null_injected = obj.chimera_exclusion(null_injected_index,[4 5]);
    
    n_test_groups = 6;
    clust = zeros(size(null_injected,1),n_test_groups);
    clust_data = log10(null_injected+10);  
    
    for i=1:6
        clust(:,i) = kmeans(clust_data,i,'emptyaction','singleton',...
            'replicate',1000);
    end

    eva = evalclusters(null_injected,clust,'silhouette');

    obj.chimera_exclusion_clustering.groups = clust(:,eva.OptimalK);

    mu = [mean(null_injected(obj.chimera_exclusion_clustering.groups==1,1)); ...
        mean(null_injected(obj.chimera_exclusion_clustering.groups==2,1))];
    whoIsBigger = diff(mu);

    if whoIsBigger > 0
        obj.chimera_exclusion_clustering.groups =...
            obj.chimera_exclusion_clustering.groups==1;    
    else
        obj.chimera_exclusion_clustering.groups =...
            obj.chimera_exclusion_clustering.groups==2;
   end        

    fig = figure(1);
    clf
    hold on

    set(gcf,'Position',[50 50 300 420])

    x_ticks = 2:n_test_groups;
    plot(x_ticks, eva.CriterionValues(2:end),'k','LineWidth',4)
    scatter( eva.OptimalK, max(eva.CriterionValues),...
        250,'filled','r')
    
    xlabel('Number of Populations')
    ylabel('Silhouette Index')
    set(gca,'XTick',x_ticks)
    xlim([2 n_test_groups])
    ylim([-0.1 1])
    grid on

    set(gca,'FontName','Arial','FontWeight','bold','FontSize',12)

    saveFileName = ['silhouette'];
    fullFileName = fullfile(save_folder_name,saveFileName);
    print(fig, fullFileName,'-dtiffn')
    print(fig, fullFileName,'-vector', '-dsvg')
    savefig(fig, fullFileName)
    close(fig)

    % Plot scatter plots and clustering -----------------------------------
    fig = figure(1);
    clf
    hold on

    ax_lims = [0 40 0 60];

    scatter(noninjected(:,1),noninjected(:,2),100,...
        [0.25 0.25 0.25],'filled','MarkerFaceAlpha',0.5);    
    scatter(wt_injected(:,1),wt_injected(:,2),100,...
        [0.70 0.14 0.14],'filled','MarkerFaceAlpha',0.5);

    axis equal
    axis(ax_lims)

    legend('Noninjected','10+/+')    
    
    ylabel 'Donor Derived Epiblast'
    xlabel 'Host Derived Epiblast'
    
    grid on

    set(gca,'FontName','Arial','FontWeight','bold','FontSize',12)

    saveFileName = ['scatter_epiblast_composition'];
    fullFileName = fullfile(save_folder_name,saveFileName);
    print(fig, fullFileName,'-dtiffn')
    print(fig, fullFileName,'-vector', '-dsvg')
    savefig(fig, fullFileName)
    close(fig)

    fig = figure(1);
    clf
    hold on

    ax_lims = [0 40 0 60];

    h_lo = scatter(null_injected(obj.chimera_exclusion_clustering.groups==0,1),...
        null_injected(obj.chimera_exclusion_clustering.groups==0,2),...
        100, [0.5 0.25 0],'filled','MarkerFaceAlpha',0.5,'Marker','square');
    h_lo = scatter(null_injected(obj.chimera_exclusion_clustering.groups==1,1),...
        null_injected(obj.chimera_exclusion_clustering.groups==1,2),...
        100, [1 0.549 0],'filled','MarkerFaceAlpha',0.5,'Marker','^');    
    axis equal
    axis(ax_lims)

    legend('10-/- Lo','10-/- Hi')
    
    ylabel 'Donor Derived Epiblast'
    xlabel 'Host Derived Epiblast'
    
    grid on

    set(gca,'FontName','Arial','FontWeight','bold','FontSize',12)

    saveFileName = ['scatter_epiblast_composition_clustering'];
    fullFileName = fullfile(save_folder_name,saveFileName);
    print(fig, fullFileName,'-dtiffn')
    print(fig, fullFileName,'-vector', '-dsvg')
    savefig(fig, fullFileName)
    close(fig)


    obj.chimera_exclusion(null_injected_index) = ...
        ~obj.chimera_exclusion_clustering.groups;

end %======================================================================
 

function analysis_chimera_exclusion(obj) %=================================
    
    close all
    addpath(fullfile('funs','plotSpread','plotSpread'))
                
    % Embryo/Chimera Composition Data Grops:
    % 1             2               3           4     
    % 10-/- hi      10-/- lo        Non         10+/+ 
    
    % Columns of 'data'
    % 1     number          Group
    % 2     number          Cell Injected
    % 3     Cell number     PrE
    % 4     Cell number     Host Epiblast
    % 5     Cell number     Donor Epiblast
    % 6     Cell number     Whole Epiblast
    % 7     Cell number     ICM
    % 8     Epi Fraction    Host/Total Epiblast
    % 9     ICM Fraction    PrE
    % 10    ICM Fraction    Host Epiblast
    % 11    ICM Fraction    Donor Epiblast
    % 12    ICM Fraction    Whole Epiblast   
    


    data = [obj.chimera_exclusion(:,1)+1, ...
        obj.chimera_exclusion(:,2:end)];
    
    analysis_folder_name = fullfile(obj.results_directory,...
        'chimera_exclusion');
    if ~exist(analysis_folder_name,'dir')
        mkdir(analysis_folder_name)
    end
                 
    % Set up the names for labelling/saving
    quantities = {'pre_#',...           1
                  'epi_host_#',...      2
                  'epi_donor_#',...     3
                  'epi_total_#',...     4
                  'ICM_#',...           5        
                  'epi_host_total_%',...6
                  'pre_%',...           7
                  'epi_host_%',...      8
                  'epi_donor_%',...     9
                  'epi_total_%'}; %     10
              
    % resort the data order to get the correct numbers
    data_order = [5, 1, 4, 2, 7, 10, 8, 6];
    cell_type = {'ICM',...               1
                 'PRE',...               2
                 'Epiblast',...          3  
                 'Host Epiblast',...     4     
                 'PRE',...               5
                 'Epiblast',...          6
                 'Host Epiblast',...     7
                 'Host/Total Epiblast'};%8
    % Give the tissue colors the 
    tissue_color =[0.4 0.1 0.4; % Cell #        ICM
                   1.0 0.0 1.0; % Cell #        PrE 
                   0.0 1.0 1.0; % Cell #        Epi 
                   0.0 0.5 0.5; % Cell #        Epi Host
                   1.0 0.0 1.0; % Cell % of ICM PrE
                   0.0 1.0 1.0; % Cell % of ICM Epi
                   0.0 0.5 0.5; % Cell % of ICM Epi Host
                   0.4 0.4 0.4];% Cell % of ICM Epi host/Total

    y_lims_num = [ 0  50 125; % Cell #        ICM
                   0  25  50; % Cell #        PrE 
                   0  25 100; % Cell #        Epi 
                   0  25 100];% Cell #        Epi Host
               
    % Name of Conditions
    group_names = {'10-/- Hi', '10-/- Lo', 'Non','10+/+'};

    % Set the number to round up to for the y-axis
    round_num = 25;

    % BOX_SWARM------------------------------------------------------------
    % Step through each tissue type and plot Box_swarm for each figure
    % Creat a folder to save the box_swarm plots in
    
    save_folder_name = fullfile(analysis_folder_name,'box_swarm');
    if ~exist(save_folder_name,'dir')
        mkdir(save_folder_name)
    end

    save_folder_name_stats = fullfile(analysis_folder_name,'stats');
    if ~exist(save_folder_name_stats,'dir')
        mkdir(save_folder_name_stats)
    end
    
    % Step Through each tissue type
    for i = 1:numel(data_order)

        % Get the index of for cells in the current tissue
        tissue_index = data_order(i)+2;

        % Perform Anova ---------------------------------------------------
        close all hidden
        [~,~,st,~] = anovan(data(:,tissue_index),data(:,1),[]);
        close all hidden            
        [c,~] = multcompare(st,'alpha',0.05,'display','off');   
        starz = sum(c(:,end)< [0.05 0.01 0.001 0.0001],2);
        p_value =  [c(:,1:2) c(:,end) starz];

        [m,sd,n,gname] = grpstats(data(:,tissue_index),data(:,1),...
             {@mean 'std' 'numel' 'gname'});
                                   
        % save ANOVA to a text file
        saveFileName = ['stats_' quantities{data_order(i)} '.txt'];
        fullFileName = fullfile(save_folder_name_stats,saveFileName);
        fileID = fopen(fullFileName,'w');
        
        fprintf(fileID,'%-12s  \n\n', ['% ', quantities{data_order(i)}]);
        fprintf(fileID,'%-12s  \n\n',['% SUMMARY STATS']);
        formatSpec0 = '%12s %-25s \n';
        fprintf(fileID,formatSpec0,'%    Group #',	' Condition');
        fprintf(fileID,formatSpec0,'%          1',...
        	' 10 null ES cell injected (Hi contribution)');
        fprintf(fileID,formatSpec0,'%          2',...
        	' 10 null (Lo contribution)');
        fprintf(fileID,formatSpec0,'%          3',	' noninjected');
        fprintf(fileID,formatSpec0,'%          4',...
        	' 10 WT ES cell injected');
        
        fprintf(fileID,'\n\n');
        
        fprintf(fileID,'%12s %12s %12s %12s \n',...
            '%    Group #','Mean','StdDev','N' );
       
        formatSpec = '%12u %12.2f %12.4f %12u \n';
        fprintf(fileID,formatSpec, [(str2double(gname)), m, sd, n]');
        
        fprintf(fileID,'%-12s  \n\n','% ANOVA');
        fprintf(fileID,formatSpec0,'%    Columns',' Description');         
        fprintf(fileID,formatSpec0,'%        1,2',' Groups Compared, Null Hypothesis H_0 is "The two groups are drawn from the same distribution?"');
	    fprintf(fileID,formatSpec0,'%          3',' p_value, level of significance is 0.05.');
        fprintf(fileID,formatSpec0,'%          4',' Accept or reject H_0? If p_vale > 0.05, acceptreject. Accept = 0; Reject < 0 idicates number of stars.');
        
        fprintf(fileID,'\n\n');
        
        fprintf(fileID,'%9s %9s %12s %9s \n','% Group 1','Group 2','p-value','# Stars' );
        
        formatSpec = '%9u %9u %12.4e %9u \n';
        fprintf(fileID,formatSpec, p_value');
        fclose(fileID);
                     
        % Ploty Figures ---------------------------------------------------
        fig = figure(10+tissue_index);

        % Plot a box plot with a superimposed swarm plot (plotSpread)
        hold on 
        h = plotSpread(data(:,tissue_index),'distributionIdx',data(:,1),...
                'binWidth',.2,'distributionMarkers','o');

        boxplot(data(:,tissue_index),data(:,1),'Colors','k',...
            'labels',group_names,'Symbol','','Whisker',inf)
        H = gca;

        % Step over each box and scatter plot and set the color
        for j = 1: numel(H.Children)
            switch H.Children(j).Type
                case 'line'
                    set(H.Children(j),'markerSize',7.5,'linewidth',.5,...
                        'MarkerFaceColor',tissue_color(i,:),...
                        'MarkerEdgeColor','k')
                case 'hggroup'
                    set(findall(H.Children(j).Children,'Type','Line'),...
                        'linewidth',1.5)
                    box_data = ...
                        get(findall(H.Children(j).Children,'Tag','Box'),...
                        {'XData','Ydata'});
                    for k = 1:length(box_data)
                        fill(box_data{k,1},box_data{k,2},tissue_color(i,:),...
                            'FaceAlpha',.25)
                        uistack(H.Children(1),'bottom')
                    end
                otherwise
                    0;
            end
        end

        % Set Axis Limits
        if any(ismember(quantities{data_order(i)},'#'))
            index = find(data_order==tissue_index-2);
            y_max_limits = ceil(max(data(:,tissue_index))/round_num)...
                *round_num;
            ylabel 'Cell #'
            y_ticks = y_lims_num(index,1):y_lims_num(index,2):y_lims_num(index,3);
            set(gca,'YTick', y_ticks)
           axis([0.5,4.5,0,y_max_limits])
        else
            ylabel 'Cell %'
            set(gca,'YTick',0:.25:1)
            axis([0.5,4.5,0,1])
        end

        % Add figure labels
        title(cell_type(i))
        set(gca,'FontSize',15,'FontWeight','bold', 'XTickLabelRotation',90,'YGrid','on')
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.54, .1875, 0.5]) 
        
        % Save Box_swarm
        saveFileName = ['box_swarm_' quantities{data_order(i)}];
        fullFileName = fullfile(save_folder_name,saveFileName);
        print(fig, fullFileName,'-dtiffn')
        print(fig, fullFileName, '-vector', '-dsvg')
        savefig(fig, fullFileName)

    end

    % STACKED_BARS---------------------------------------------------------
    % Stacked Bar Plots by Condition
    % Step through each Condition and make a stacked bar plot showing the
    % composition of ICM by number and ratio for each condition

    % Creat a folder to save the bar plots in
    save_folder_name = fullfile(analysis_folder_name,'bars');
    if ~exist(save_folder_name,'dir')
        mkdir(save_folder_name)
    end  
    % Creat a folder to save the bar plots in
    save_folder_name_means = fullfile(analysis_folder_name,'mean_bars');
    if ~exist(save_folder_name_means,'dir')
        mkdir(save_folder_name_means)
    end  
    % Set up the consitions to step through
    conditions = {'10nullHi','10nullLo','Non','10wt'};
    n_conditions = numel(conditions);

    % Plot both cell number and fraction for tissues
    compsition_type = {'Number','Fraction'};
    n_compsition_type = numel(compsition_type);

    % set up the sets of data to plot
    tissue_sets = {'ICM'};
    n_tissue_sets = numel(tissue_sets);

    % Plot the cell number and ratio for the tissues of the the ICM
    % {PrEand Epi}
    tissue_sets_index = {[1 2 3],[7 8 9]};

    % Set up the sorting such that the order of the embryos are
    % ascending fraction of PrE in the ICM
    [sort_data,sort_index] = sortrows(data(:,[1 9 4]));

    % Step through the different tissues sets and plot stacked bar
    % plots for the individual embryos and the mean embryos
    for i = 1:n_tissue_sets

        for j = 1:n_compsition_type

            index = n_tissue_sets*(i-1) + j;
            current_set = tissue_sets_index{index}+2;

            if strcmp(tissue_sets(i),'ICM')
                face_colors = {[1 0 1]; [0 1 1];[0.4 0.4 0.4]};
                legend_names = {'PrE','Host','Donor'};
                y_max_lim =  110;
                y_step = 25;
            else
                error('Wrong Composition Type')
             end

            if strcmp(compsition_type(j),'Number')
                y_ticks =  0:y_step:y_max_lim;
            elseif strcmp(compsition_type(j),'Fraction')
                y_max_lim = 1.05;
                y_ticks = 0:.25:1;
            else
                error('Wrong Composition Type')
            end

            % Calulate the group stats to be used for the mean bar
            % plots
            [m,sd,n,sem] = grpstats(data(:,current_set),data(:,1),...
                {@mean 'std' 'numel' 'sem'});

            % Plot mean stacked bar values
            fig = figure(1);
            clf
            h2 = bar(m,'stacked','barwidth',1);
            hold on
            errorbar(cumsum(m,2),sem,'.k','linewidth',2);
            set(h2,{'FaceColor'},face_colors)
            axis([0.5,numel(group_names)+.5,0, y_max_lim])

            set(gcf, 'Units', 'Normalized', 'OuterPosition', ...
                [0, 0.04, .125, 0.5])

            set(gca,'FontSize',10,'FontWeight','bold','YTick',y_ticks)
            set(gca,'XTickLabels','','XTickLabelRotation',90,'Position',...
                [0.28 0.21 0.66 0.72])
            set(gca,'box','off')
                
            set(0,'DefaultAxesColor','none')
                
            % Save Mean Bar plots
            saveFileNameMean = [compsition_type{j} '_' tissue_sets{i}];
            fullFileName = fullfile(save_folder_name_means,saveFileNameMean);
            print(fig, fullFileName,'-dtiffn')
            print(fig, fullFileName, '-vector', '-dsvg')
            savefig(fig, fullFileName)

            % Step through the individual inject types and make stacked
            % bar plots for the individual embryos
            for k = 1:n_conditions

                embryo_order = sort_index(sort_data(:,1)==k);
                current_data = data(embryo_order ,current_set);

                fig = figure(k+1);
                clf
                h1 = bar(current_data,'stacked','barwidth',1);
                set(h1,{'FaceColor'},face_colors);
                set(gca,'XTickLabels',[])

                ylabel(compsition_type{j})
                title(group_names{k})
                axis([0.5,size(current_data,1)+.5,0, y_max_lim])
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.125*(k), 0.04, .125, 0.5])
                legend(legend_names,'Location','SouthOutside','Orientation','Horizontal')
                set(gca,'FontSize',12,'FontWeight','bold','YTick',y_ticks,...
                    'XTickLabelRotation',90,'Position', [0.28 0.21 0.66 0.72])
                
                % Save Bar plots
                saveFileName = [saveFileNameMean '_' conditions{k}];
                fullFileName = fullfile(save_folder_name,saveFileName);
                print(fig, fullFileName,'-dtiffn')
                print(fig, fullFileName, '-vector', '-dsvg')
                savefig(fig, fullFileName)
            end
        end
    end

end %======================================================================

      
function analysis_chimera_crowding(obj)  %=================================
    
    close all
    addpath(fullfile('funs','plotSpread','plotSpread'))
    
    %  Groups:
    % 1         2           3           4           5
    % 15-/-     10-/-       Non         10+/+       15+/+
    
    % Columns of obj.chimera_crowding
    % 1     number          Group
    % 2     Cell number     TE 
    % 3     Cell number     PrE
    % 4     Cell number     Epiblast
    % 5     Cell number     Embryo
    % 6     Cell number     ICM
    % 7     Embryo Fraction TE
    % 8     Embryo Fraction ICM
    % 9     Embryo Fraction PrE
    % 10    Embryo Fraction Epiblast
    % 11    ICM Fraction    PrE
    % 12    ICM Fraction    Epiblast       
    
    data = obj.chimera_crowding;
    
    % Set up the names for labelling/saving
    quantities = {'te_#',...        1
                  'pre_#',...       2
                  'epi_#',...       3
                  'embryo_#',...    4
                  'ICM_#',...       5
                  'te_%',...        6
                  'icm_%',...       7
                  'pre_embryo_%',...8
                  'epi_embryo_%',...9
                  'pre_icm_%',...   10
                  'epi_icm_%'};%    11
    
    % resort the data order to get the correct numbers
    data_order = [4 1 5 3 2 6 7 11 10]; 
    cell_type = {   'All',...
                    'TE',...
                    'ICM',...
                    'Epiblast',...
                    'PrE',...
                    'TE',...
                    'ICM',...
                    'Epiblast',...
                    'PrE'};
    
    % Give the tissue colors the 
    tissue_color =[0.4 0.4 0.4; % Cell #            Whole Embryo
                   0.1 0.1 1.0; % Cell #            TE 
                   0.4 0.1 0.4; % Cell #            ICM 
                   0.0 1.0 1.0; % Cell #            Epiblast
                   1.0 0.0 1.0; % Cell #            PrE 
                   0.1 0.1 1.0; % Cell % of Embryo  TE
                   0.4 0.1 0.4; % Cell % of Embryo  ICM
                   0.0 1.0 1.0; % Cell % of ICM     Epiblast
                   1.0 0.0 1.0]; % Cell % of ICM     PrE
    
    y_lims_num = [  0  50 300; % Cell #            Whole Embryo
                    0  50 225; % Cell #            TE 
                    0  50 125; % Cell #            ICM 
                    0  25 100; % Cell #            Epiblast
                    0  25  50];% Cell #            PrE 
    
    % Name of Conditions
    group_names = {'15-/-','10-/-','Non','10+/+','15+/+'};
    
    % Set the number to round up to for the y-axis
    round_num = 25;    
    
    % BOX_SWARM------------------------------------------------------------
    % Step through each tissue type and plot Box_swarm for each figure
    % Creat a folder to save the box_swarm plots in
    save_folder_name = fullfile(obj.results_directory,...
        'chimera_crowding','box_swarm');
    if ~exist(save_folder_name,'dir')
        mkdir(save_folder_name)
    end
    
    save_folder_name_stats = fullfile(obj.results_directory,...
        'chimera_crowding','stats');
    if ~exist(save_folder_name_stats,'dir')
        mkdir(save_folder_name_stats)
    end
    
    % Step Through each tissue type
    for i = 1:numel(data_order)
    
        % Get the index of for cells in the current tissue
        tissue_index = data_order(i)+1;
       
        % Perform Anova ---------------------------------------------------
        close all hidden
        [~,~,st,~] = anovan(data(:,tissue_index),data(:,1),[]);    
        close all hidden
        [c,~] = multcompare(st,'alpha',0.05,'CType','lsd','display','off');      
        starz = sum(c(:,end)< [0.05 0.01 0.001 0.0001],2);
        p_value =  [c(:,1:2) c(:,end) starz];            
    
        [m,sd,n,gname] = grpstats(data(:,tissue_index),data(:,1),...
             {@mean 'std' 'numel' 'gname'});
         
        % save ANOVA to a text file
        saveFileName = ['stats_' quantities{data_order(i)} '.txt'];
        fullFileName = fullfile(save_folder_name_stats,saveFileName);
        fileID = fopen(fullFileName,'w');
        
        fprintf(fileID,'%-12s  \n\n', ['% ', quantities{data_order(i)}]);
        fprintf(fileID,'%-12s  \n\n',['% SUMMARY STATS']);
        formatSpec0 = '%12s %-25s \n';
        fprintf(fileID,formatSpec0,'%    Group #',	' Condition');
        fprintf(fileID,formatSpec0,'%          1',	' 15 null ES cell injected');
        fprintf(fileID,formatSpec0,'%          2',	' 10 null');
        fprintf(fileID,formatSpec0,'%          3',	' noninjected');
        fprintf(fileID,formatSpec0,'%          4',	' 10 WT ES cell injected');
        fprintf(fileID,formatSpec0,'%          5',	' 15 WT');
        
        fprintf(fileID,'\n\n');
        
        fprintf(fileID,'%12s %12s %12s %12s \n','%    Group #','Mean','StdDev','N' );
       
        formatSpec = '%12u %12.2f %12.4f %12u \n';
        fprintf(fileID,formatSpec, [(str2double(gname)), m, sd, n]');
        
        fprintf(fileID,'%-12s  \n\n','% ANOVA');
        fprintf(fileID,formatSpec0,'%    Columns',' Description');         
        fprintf(fileID,formatSpec0,'%        1,2',' Groups Compared, Null Hypothesis H_0 is "The two groups are drawn from the same distribution?"');
        fprintf(fileID,formatSpec0,'%          3',' p_value, level of significance is 0.05.');
        fprintf(fileID,formatSpec0,'%          4',' Accept or reject H_0? If p_vale > 0.05, acceptreject. Accept = 0; Reject < 0 idicates number of stars.');
        
        fprintf(fileID,'\n\n');
        
        fprintf(fileID,'%9s %9s %12s %9s \n','% Group 1','Group 2',...
            'p-value','# Stars' );
        
        formatSpec = '%9u %9u %12.4e %9u \n';
        fprintf(fileID,formatSpec, p_value');
        fclose(fileID);
                  
        % Ploty Figures ---------------------------------------------------
        fig = figure(10+tissue_index);
    
        % Plot a box plot with a superimposed swarm plot (plotSpread)
         hold on 
         h = plotSpread(data(:,tissue_index),'distributionIdx',data(:,1),...
                'binWidth',.5,'distributionMarkers','o');
        boxplot(data(:,tissue_index),data(:,1),'Colors','k',...
            'labels',group_names,'Symbol','','Whisker',inf)
        H = gca;
    
        % Step over each box and scatter plot and set the color
        for j = 1: numel(H.Children)
            switch H.Children(j).Type
                case 'line'
                    set(H.Children(j),'markerSize',7.5,'linewidth',.5,...
                        'MarkerFaceColor',tissue_color(i,:),...
                        'MarkerEdgeColor','k')
                case 'hggroup'
                    set(findall(H.Children(j).Children,'Type','Line'),...
                        'linewidth',1.5)
                    box_data = ...
                        get(findall(H.Children(j).Children,'Tag','Box'),...
                        {'XData','Ydata'});
                    for k = 1:length(box_data)
                        fill(box_data{k,1},box_data{k,2},...
                            tissue_color(i,:),'FaceAlpha',.25)
                        uistack(H.Children(1),'bottom')
                    end
                otherwise
                    0;
            end
        end
        
        % Set Axis Limits
        if any(ismember(quantities{data_order(i)},'#'))
            index = find(data_order==tissue_index-1);
            if index == 10
               index = 6;
            end                
            y_max_limits =...
                ceil(max(data(:,tissue_index))/round_num)*round_num;
            ylabel 'Cell #'
            y_ticks = ...
                y_lims_num(index,1):y_lims_num(index,2):y_lims_num(index,3);
            set(gca,'YTick', y_ticks)
           axis([0.5,5.5,0,y_max_limits])
        else
            ylabel 'Cell %'
            set(gca,'YTick',0:.25:1)
            axis([0.5,5.5,0,1.1])
        end
    
        % Add figure labels
        title(cell_type(i))
        set(gca,'FontSize',15,'FontWeight','bold',...
            'XTickLabelRotation',90,'YGrid','on')
        set(gcf, 'Units', 'Normalized',...
            'OuterPosition', [0, 0.54, .1875, 0.5])           
       
        % Save Box_swarm
        saveFileName = ['box_swarm_' quantities{data_order(i)}];
        fullFileName = fullfile(save_folder_name,saveFileName);
        print(fig, fullFileName,'-dtiffn')
        print(fig, fullFileName, '-vector', '-dsvg')
        savefig(fig, fullFileName)
    end        

    % STACKED_BARS---------------------------------------------------------
    % Stacked Bar Plots by Condition
    % Step through each Condition and make a stacked bar plot showing the
    % composition of each embryo by number and ratio and the ICM by number
    % and ratio for each condition
    
    % Creat a folder to save the bar plots in
    save_folder_name  = fullfile(obj.results_directory,...
        'chimera_crowding','bars');
    if ~exist(save_folder_name,'dir')
        mkdir(save_folder_name)
    end  
    % Creat a folder to save the bar plots in
    save_folder_name_means  = fullfile(obj.results_directory,...
        'chimera_crowding','mean_bars');
    if ~exist(save_folder_name_means,'dir')
        mkdir(save_folder_name_means)
    end  
    % Set up the conditions to step through
    conditions = {'15null','10null','Non','10wt','15wt'};
    n_conditions = numel(conditions);
    
    % Plot both cell number and fraction for tissues
    compsition_type = {'Number','Fraction'};
    n_compsition_type = numel(compsition_type);
    
    % set up the sets of data to plot
    tissue_sets = {'Embryo','ICM'};
    n_tissue_sets = numel(tissue_sets);
    
    % Plot the cell number and ratio for the tissues of the whole
    % embryo {TE,PrE,and Epi} and the ICM {PrEand Epi}
    tissue_sets_index = {[1 2 3],[6 8 9],[2 3],[10 11]};
    
    % Set up the sorting such that the order of the embryos are
    % ascending fraction of PrE in the ICM
    [sort_data,sort_index] = sortrows(data(:,[1 11 4]));
    
    % Step through the different tissues sets and plot stacked bar
    % plots for the individual embryos and the mean embryos
    for i = 1:n_tissue_sets
    
        for j = 1:n_compsition_type
    
            index = n_tissue_sets*(i-1) + j;
            current_set = tissue_sets_index{index}+1;
    
             if strcmp(tissue_sets(i),'Embryo')
                face_colors = {[0 0 1]; [1 0 1]; [0 1 1]};
                legend_names = {'TE','PrE','EPI'};
                y_max_lim =  275;
                 y_step = 50;
            elseif strcmp(tissue_sets(i),'ICM')
                face_colors = {[1 0 1]; [0 1 1]};
                legend_names = {'PrE','EPI'};
                y_max_lim =  110;
                y_step = 25;
            else
                error('Wrong Composition Type')
             end
    
            if strcmp(compsition_type(j),'Number')
                y_ticks =  0:y_step:y_max_lim;
            elseif strcmp(compsition_type(j),'Fraction')
                y_max_lim = 1.05;
                y_ticks = 0:.25:1;
            else
                error('Wrong Composition Type')
            end
    
            % Calulate the group stats to be used for the mean bar
            % plots
            [m,sd,n,sem] = grpstats(data(:,current_set),data(:,1),...
             {@mean 'std' 'numel' 'sem'});        
            
            % Plot mean stacked bar values
            fig = figure(1);
            clf
            h2 = bar(m,'stacked','barwidth',1);
            hold on
            errorbar(cumsum(m,2),sem,'.k','linewidth',2);
            set(h2,{'FaceColor'},face_colors)
            axis([0.5,numel(group_names)+.5,0, y_max_lim])
    
            set(gcf, 'Units', 'Normalized', 'OuterPosition',...
                [0, 0.04, 0.125, 0.5])
    
            set(gca,'FontSize',10,'FontWeight','bold','YTick',y_ticks)
            set(gca,'XTickLabels','','Position', [0.28 0.21 0.66 0.72]);
            set(gca,'box','off')
            
            set(0,'DefaultAxesColor','none')
            
            % Save Mean Bar plots
            saveFileNameMean = [compsition_type{j} '_' tissue_sets{i}];
            fullFileName = fullfile(save_folder_name_means,saveFileNameMean);
            print(fig, fullFileName,'-dtiffn')
            print(fig, fullFileName, '-vector', '-dsvg')
            savefig(fig, fullFileName)
    
            % Step through the individual inject types and make stacked
            % bar plots for the individual embryos
            for k = 1:n_conditions
    
                embryo_order = sort_index(sort_data(:,1)==k);
                current_data = data(embryo_order ,current_set);
    
                fig = figure(k+1);
                clf
                h1 = bar(current_data,'stacked','barwidth',1);
                set(h1,{'FaceColor'},face_colors);
                set(gca,'XTickLabels',[])
    
                ylabel(compsition_type{j})
                title(group_names{k})
                axis([0.5,size(current_data,1)+.5,0, y_max_lim])
                set(gcf, 'Units', 'Normalized', 'OuterPosition',...
                    [0.125*(k), 0.04, .125, 0.5])
                legend(legend_names,'Location','SouthOutside',...
                    'Orientation','Horizontal')
                set(gca,'FontSize',12,'FontWeight','bold','YTick',y_ticks,...
                    'XTickLabelRotation',90,'Position', ...
                    [0.28 0.21 0.66 0.72])
                
                % Save Bar plots
                saveFileName = [saveFileNameMean '_' conditions{k}];
                fullFileName = fullfile(save_folder_name,saveFileName);
                print(fig, fullFileName,'-dtiffn')
                print(fig, fullFileName, '-vector', '-dsvg')
                savefig(fig, fullFileName)
            end
        end
    end
        
end %======================================================================
        

function embryo_inference(obj) %===========================================
    
    close all
    
    save_folder_name = fullfile(obj.results_directory,'abc_mcmc_embryos');
    if ~exist(save_folder_name,'dir')
        mkdir(save_folder_name)
    end

    %----------------------------------------------------------------------
    % Initialize parameters, time, cell state, and ABC parameters    
    
    % Set an initial guess for the parameter set
    theta_0 = [0.062 0.18 0.36 0.012 0.00014 0 2; 
               0.061 0.2 0.33 0.0007 0.00014 2 2]; 
    % theta_0 = [0.062564 0.21311 0.35441 8.4736e-05 0.00014127 2];

    % For the ABC MCMC how many parameter sets do you want to keep
    num_accepted = 500;
    
    % After how many itereations do you want to stop trying
    num_samples = 1e6;
    
    % set the lower upper-bound for accepted sum_squared_error for the
    % inference
    
    % error threshold determined by running err_types(2) which finds a
    % minima over the parameter space. Then setting err to that value and
    % rerunning abc   
    err = [17300, 13500];     
    pre_feedback = {'no','yes'};
    err_types = {'const','sse'};
    err_type = err_types(1);
    
    %----------------------------------------------------------------------
    % Load Training data
    data = [obj.embryo_litter_mate.cells(:,2),...% Double Positive
        sum(obj.embryo_litter_mate.cells(:,[1 3]),2),...% combined DN/EPI
        obj.embryo_litter_mate.cells(:,4),...% PrE
        obj.embryo_litter_mate.cells(:,5),...%  TE
        sum(obj.embryo_litter_mate.cells,2)];% whole Embryo

    %----------------------------------------------------------------------
    % Perform ABC
    obj.best_fit_parameters_embryo = struct('pre_feedback','','stats',[],...
        'accepted',[],'sse',[], 'sampled', []);

    for i = 1:2
        
        obj.embryo_abc_parameter_space{i}(:,:) = [obj.embryo_abc_priors.alpha;
            obj.embryo_abc_priors.beta;
            obj.embryo_abc_priors.rho;
            obj.embryo_abc_priors.zeta;
            obj.embryo_abc_priors.eta;
            obj.embryo_abc_priors.l
            obj.embryo_abc_priors.m]';

        if strcmp(pre_feedback{i},'no')
            obj.embryo_abc_parameter_space{i}(:,4) = [0.008 0.02]';            
        end

        [best_fit_parameters] = ...
            model_inference.abc_mcmc_embryo(obj.embryo_abc_parameter_space{i},...
            obj.time,obj.embryo_abc_IC,data,theta_0(i,:),err(i),err_type,num_samples,...
            num_accepted,pre_feedback{i});
        
        obj.best_fit_parameters_embryo(i) = best_fit_parameters;

    end
    
end %======================================================================
        

function plotEmbryoParamerters(obj) %======================================
   
    % Plot matrix of parameters -------------------------------------------    
    parameter_name = {'\alpha','\beta','\rho','\zeta','\eta','l','m','error'};

    for k = 1:numel(obj.best_fit_parameters_embryo)
       
        accepted = [obj.best_fit_parameters_embryo(k).accepted,...
            obj.best_fit_parameters_embryo(k).sse(:,1)];
        
        if k == 1
            accepted(:,5) = accepted(:,5).*(1e4); 
        else
            accepted(:,[4 5]) = accepted(:,[4 5]).*(1e4);
        end
        n_parameters = 5;%size(accepted,2);
        
        H = figure(1);
        clf 
        hold on
        
        x_lim = [obj.embryo_abc_parameter_space{k}, ...
            [min(obj.best_fit_parameters_embryo(k).sse(:,1));...
            max(obj.best_fit_parameters_embryo(k).sse(:,1))] ]';

        if k == 1
            x_lim(5,:) = x_lim(5,:).*(1e4); 
        else
             x_lim([4 5],:) = x_lim([4 5],:).*(1e4);            
        end

        for i = 1:n_parameters        
    
            % plot univariate distribution 
            subplot(n_parameters,n_parameters,(i-1)*n_parameters + i)    
           
            hold on
            % x_ticks = min(x_lim):10^floor(log10(diff(x_lim))):max(x_lim)
    
            med = median(accepted(:,i));
            
            h = histogram(accepted(:,i),'Normalization','pdf',...
                'EdgeColor','none','FaceColor','k','FaceAlpha',0.6);
            ay = ylim;
            plot([med med], ay,'k','linewidth',2)        
            
            xlim(x_lim(i,:))
            ylim([0 max(h.Values)])
    
            axis square
            grid on
            
            title(parameter_name{i})
            ylabel('pdf')
    
            for j = 1:i-1
    
                % plot bivariate distribution
                subplot(n_parameters,n_parameters,(i-1)*n_parameters + j)
                scatter(accepted(:,j),accepted(:,i),5,'filled',...
                    'MarkerFaceColor','k','MarkerFaceAlpha',0.5,...
                    'MarkerEdgeColor','none');
                axis([x_lim(j,:) x_lim(i,:)])
                axis square
                grid on
    
            end
        end
    
        save_name = 'pre_feedback=';
        
        if strcmp(obj.best_fit_parameters_embryo(k).pre_feedback,'yes')
            save_name = [save_name,'yes'];
        else
            save_name = [save_name,'no'];
        end

        fullFileName = fullfile(obj.results_directory,'abc_mcmc_embryos',...
            ['parameter_posteriors_',save_name]);
        print(H, [fullFileName, '.tif'],'-dtiffn')
        print(H, [fullFileName, '.svg'],'-vector','-dsvg')  

    end

end %======================================================================


function plotEmbryoSimulation(obj) %=======================================

    IC = [8 0 0 0 0];
        
    for i = 1:numel(obj.best_fit_parameters_embryo)
        parameters = obj.best_fit_parameters_embryo(i).stats(1,:);
        [xn(:,:,i),~] = model_inference.numerical_solution_pre_feedback(...
            parameters,obj.time,IC);
    end

    t = obj.time';
    tInt = [t(1) t(end)];
    
    % get the total number if cells present at time t
    sumXn = sum(xn,2);
    sumEmbryo = sum(obj.embryo_litter_mate.cells,2);
    % normalize cell number by total number of cell
    normXn = xn(:,[1 4 3 5 2],:)./repmat(sumXn,1,5,1);  
     
    cell_colors = [0.6000 0.6000 0.6000;
                   0.0000 0.0000 1.0000;
                   1.0000 1.0000 1.0000;
                   1.0000 0.0000 1.0000;
                   0.0000 1.0000 1.0000];
    
    embryo_keep_idx = sum(obj.embryo_litter_mate.cells,2)<151;
    
    % Plotting! ===========================================================
         
    %----------------------------------------------------------------------
    % Tissue Fraction Fill Plots!
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % Embryo tissue Types
    % Recapitulate Saiz (2016) Nat Comm Figure

    for i = 1:2
        H = figure(1);
        clf
        hold on
        
        fill([t; t(end)],1-[sum(normXn(:,[3 4 5],i),2); 0],cell_colors(3,:))
        fill([t; t(end)],1-[sum(normXn(:,[ 4 5],i),2); 0],cell_colors(4,:))
        fill([t; t(end)],1-[normXn(:,5,i); 0],cell_colors(2,:))
        fill([t; t(end); t(1)],[sum(normXn(:,[1 2],i),2);0; 0],cell_colors(5,:))
        fill([t; t(1)],[normXn(:,1,i); 0],cell_colors(1,:))
           
        xlim(tInt)
        ylim([0 1])
        
        % Labels!
        xlabel 'Time (hours)'
        ylabel 'Cell Fraction'
        set(gca,'FontSize',15,'FontWeight','Bold',...
            'XTick',[0:8:48],'YTick',[0:0.25:1],...
            'YTickLabel',num2str([0:0.25:1]','%.2f'))
        title('Embryo Composition','FontSize',17)    
        
        % Saving figure!!!
        save_name = ['fill_embryo_pre_feedback=', ...
            obj.best_fit_parameters_embryo(i).pre_feedback];
        fullFileName = fullfile(pwd,'results','abc_mcmc_embryos',save_name);
        print(H, [fullFileName, '.tif'],'-dtiffn')
        print(H, [fullFileName, '.svg'], '-vector', '-dsvg')    

    

    %----------------------------------------------------------------------
    % ICM tissue Types
    
    % Recapitulate Saiz (2016) Nat Comm Figure
        
        H = figure(1);
        clf
        hold on
        norm3 = [xn(:,3:5,i)]./repmat(nansum(xn(:,3:end,i),2),1,3);
            
        fill([t(1) t(end) t(end) t(1)],[0 0 1 1],cell_colors(3,:))
        fill([t; t(end)],[0;norm3([2:end],2);0],cell_colors(5,:))
        fill([t; t(end)],[1;1-norm3([2:end],3);1],cell_colors(4,:))
        
        xlim(tInt)
        ylim([0 1])
        xlabel 'Time (hours)'
        ylabel 'Cell Fraction'
        
        set(gca,'FontSize',15,'FontWeight','Bold',...
            'XTick',[0:8:48],'YTick',[0:0.25:1],...
            'YTickLabel',num2str([0:0.25:1]','%.2f'))
        title ('Tissue Fraction of ICM','FontSize',17)
    
        % Saving figure!!!
        save_name = ['fill_icm_pre_feedback=', ...
            obj.best_fit_parameters_embryo(i).pre_feedback];
        fullFileName = fullfile(pwd,'results','abc_mcmc_embryos',save_name);
        print(H, [fullFileName, '.tif'],'-dtiffn')
        print(H, [fullFileName, '.svg'], '-vector', '-dsvg')    
       
    end

    %----------------------------------------------------------------------
    % ScatterPlots!
    %----------------------------------------------------------------------
    
    line_color = [0.6 0.6 0.6];

    %----------------------------------------------------------------------
    % {ICM,PrE,Epi} v TE
    H = figure(1);
    clf
    % Plot the empirical data and the model#
    scatter(obj.embryo_litter_mate.cells(embryo_keep_idx,5),...
        sum(obj.embryo_litter_mate.cells(embryo_keep_idx,1:4),2),[],...
        sumEmbryo(embryo_keep_idx),'filled')
    hold on

    plot(xn(:,2,1),sum(xn(:,3:5,1),2),':','Color',...
        line_color,'linewidth',3)
    plot(xn(:,2,2),sum(xn(:,3:5,2),2),'-','Color',...
        0.5*line_color,'linewidth',3)
    % Set color sca-le
    c = colorbar;
    c.Location = 'EastOutside';
    c.Position = [.9 .125 .03 .75];
    color_map = [linspace(.2, .8,1000)+.2;linspace(.2, .8,1000);linspace(.2, .8,1000)]';
	colormap(color_map)
    clim([0 150])

    % Labels!
    ylabel '# ICM + EPI + PrE'
    xlabel '# TE'
    title  'Whole Embyro'
    legend('Empirical','No PrE Feedback','PrE Feedback','Location','northwest')
    grid on 

    set(gca,'FontSize',12,'FontWeight','Bold','FontName','Arial')

    % Saving figure!!!
    fullFileName = fullfile(pwd,'results','abc_mcmc_embryos',...
        'scatter_icm_v_TE');
    print(H, [fullFileName, '.tif'],'-dtiffn')
    print(H, [fullFileName, '.svg'], '-vector','-dsvg')    
       
    axis([0 120 0 60])

    %----------------------------------------------------------------------
    % DP v Epi v PrE
    pause(0.1)
    H = figure(1);
    clf
    % Plot the empirical data and the model#
    scatter3(obj.embryo_litter_mate.cells(embryo_keep_idx,2),...
        sum(obj.embryo_litter_mate.cells(embryo_keep_idx,[1 3]),2),...
        obj.embryo_litter_mate.cells(embryo_keep_idx,4),...
        [],sumEmbryo(embryo_keep_idx),'filled')
    hold on
    plot3(xn(:,3,1),xn(:,4,1),xn(:,5,1),':','Color',...
        line_color,'linewidth',3)
    plot3(xn(:,3,2),xn(:,4,2),xn(:,5,2),'-','Color',...
        0.5*line_color,'linewidth',3)
    axis equal
    axis([0 30 0 30 0 30])

    % Set color sca-le
    c = colorbar ;
    c.Location = 'EastOutside';
    c.Position = [.9 .125 .03 .75];
    color_map = [linspace(.2, .8,1000)+.2;linspace(.2, .8,1000);linspace(.2, .8,1000)]';
	colormap(color_map)
    clim([0 150])

    % Labels!
    xlabel '# ICM'
    ylabel '# EPI'
    zlabel '# PrE'
    
    title  'Total ICM'

    legend('Empirical','No PrE Feedback','PrE Feedback','Location','northeast')
    grid on 

    set(gca,'FontSize',12,'FontWeight','Bold','FontName','Arial')    
    
    view(-38,11)

    % Saving figure!!!
    fullFileName = fullfile(pwd,'results','abc_mcmc_embryos',...
        'scatter_icm_epi_pre');
    print(H, [fullFileName, '.tif'],'-dtiffn')
    print(H, [fullFileName, '.svg'], '-vector', '-dsvg')   
    pause(0.1)

    % Plot lineage size as a function of embryo size and ICM size ---------
    cell_colors = [...
        0.0000 0.0000 1.0000;
        0.6000 0.6000 0.6000;
        0.0000 1.0000 1.0000;
        1.0000 0.0000 1.0000];

    n_accepted = 500;
    Xn = zeros(numel(obj.time),5,n_accepted,2);

    for i = 1:numel(obj.best_fit_parameters_embryo)        
        for j = 1:n_accepted
        [Xn(:,:,j,i),tn] = ...
            model_inference.numerical_solution_pre_feedback(...
            obj.best_fit_parameters_embryo(i).accepted(j,:),...
            obj.time,[8 0 0 0 0]);    
        end
    end

    reference = {'# Whole Embryo','# Whole ICM'};
    cell_types = {'# TE','# ICM', '# EPI','# PrE'};
    y_max = [110 30 30 30];
    emperical_cell_type_idx = {5, 2,[1 3], 4};

    for i = 1:numel(reference)
        
        if strcmp(reference{i},'# Whole Embryo')
            empirical_x_data = sum(obj.embryo_litter_mate.cells,2);
            simulation_x_data = sum(Xn,2);     
            fullFileName = fullfile(pwd,'results','abc_mcmc_embryos',...
                'whole_embryo_');
            c_max = 150;
             mu_x = mean(simulation_x_data,3);
            x_lims = [0 max(mu_x(:))];
            x_ticks = x_lims(1):25:x_lims(2);
        else
            empirical_x_data = sum(obj.embryo_litter_mate.cells(:,1:4),2);
            simulation_x_data = sum(Xn(:,3:end,:,:),2);     
            fullFileName = fullfile(pwd,'results','abc_mcmc_embryos',...
                'whole_icm_');
            c_max = 50;
             mu_x = mean(simulation_x_data,3);
            x_lims = [0 max(mu_x(:))];
            x_ticks = x_lims(1):10:x_lims(2);
        end

        for j = 1:size(cell_types,2)

            fig1 = figure(1);
            clf
            empirical_y_data = ...
                sum(obj.embryo_litter_mate.cells(:,emperical_cell_type_idx{j}),2);

            simulation_y_data = Xn(:,j+1,:,:);

            scatter(empirical_x_data(embryo_keep_idx),...
                empirical_y_data(embryo_keep_idx),[],...
                empirical_x_data(embryo_keep_idx),'filled')
            hold on

           
            mu_y = mean(simulation_y_data,3);
            sigma_y_plus = mu_y + std(simulation_x_data,[],3);
            sigma_y_minus =  mu_y - std(simulation_x_data,[],3);

            plot( mu_x(:,:,:,1),mu_y(:,:,:,1),':','Color',...
                cell_colors(j,:),'linewidth',4)
            plot( mu_x(:,:,:,2),mu_y(:,:,:,2),'-','Color',...
                0.5.*cell_colors(j,:),'linewidth',4)

            fill([mu_x(:,:,:,1);flip(mu_x(:,:,:,1))],...
                [sigma_y_plus(:,:,:,1);flip(sigma_y_minus(:,:,:,1))],...
                cell_colors(j,:),'FaceAlpha',0.25,'LineStyle','none')
           fill([mu_x(:,:,:,2);flip(mu_x(:,:,:,2))],...
               [sigma_y_plus(:,:,:,2);flip(sigma_y_minus(:,:,:,2))],...
                0.5.*cell_colors(j,:),'FaceAlpha',0.25,'LineStyle','none')
                        
            % Set color sca-le
            c = colorbar;
            c.Location = 'EastOutside';
            c.Position = [.9 .125 .03 .75];
            color_map = [linspace(.2, .8,1000)+.2;...       
                         linspace(.2, .8,1000);
                         linspace(.2, .8,1000)]';
            colormap(color_map)
            clim([0 c_max])

            % Labels!
            xlabel(reference{i})
            ylabel(cell_types{j})            
            grid on 
            
            xlim(x_lims)
            ylim([0 y_max(j)])
            set(gca,'FontSize',12,'FontWeight','Bold','FontName','Arial')
            set(gca,'XTick',x_ticks)
            saveFileName = [fullFileName, extractAfter(cell_types{j},'# ')];
            print(fig1, [saveFileName, '.tif'],'-dtiffn')
            print(fig1, [saveFileName, '.svg'], '-vector', '-dsvg')  
                    
        end
    end

end %======================================================================


function chimera_inference(obj) %==========================================
    
    close all
    
    save_folder_name = fullfile(obj.results_directory,'abc_mcmc_chimera');
    if ~exist(save_folder_name,'dir')
        mkdir(save_folder_name)
    end
    
    obj.chimera_abc_parameter_space = [obj.chimera_abc_priors.alpha_D;
        obj.chimera_abc_priors.a;
        obj.chimera_abc_priors.n]';

    obj.best_fit_parameters_chimera = struct('model','','stats',[],...
        'accepted',[],'sse',[]);
    model_names = {'Same Growth Rate', 'Different Growth Rate',...
        'Different Growth Rate and Spatial Crowding'};

    for i = 1:numel(model_names)
        obj.best_fit_parameters_chimera(i).model = model_names{i};
    end

    % for the Same Growth Rate model, the growth rate of donor cells uses
    % the the inferred growth rate from the embryo. 
    obj.best_fit_parameters_chimera(1).stats = ...
        [obj.best_fit_parameters_embryo.stats(:,[1:end, 1])];
    obj.best_fit_parameters_chimera(1).accepted =  ...
        [obj.best_fit_parameters_embryo.accepted(:,[1:end, 1]),...
        zeros(size(obj.best_fit_parameters_embryo.accepted,1),1)];
    obj.best_fit_parameters_chimera(1).sse = nan;
     obj.best_fit_parameters_chimera(1).sampled = nan;

    % Set an initial guess for the parameter set
    theta_0 = {[0.018], [0.028 0.11 1]};
    parameters_to_infer = {[1] , [1 2 3]};
    embryo_parameters = obj.best_fit_parameters_embryo.stats(1,:);

    % For the ABC MCMC how many parameter sets do you want to keep
    num_accepted = 500;
    
    % After how many itereations do you want to stop trying
    num_samples = 1e6;
    
    % set the lower upper-bound for accepted sum_squared_error for the
    % inference
    
    % error threshold determined by running err_types(2) which finds a
    % minima over the parameter space. Then setting err to that value and
    % rerunning abc   
    err = [1.238e5, 6.2e4]; 
    % err = [1e6 1e6]; 
  
    err_types = {'const','sse'};
    err_type = err_types(1);
    % err_type = err_types(2)

    %----------------------------------------------------------------------
    % Load Training data
    data = obj.chimera_crowding;

    % perform inference for two other models
    for i = 2:3

        %------------------------------------------------------------------
        % Perform ABC 
        [best_fit_parameters] = model_inference.abc_mcmc_chimera(embryo_parameters,...
            obj.chimera_abc_parameter_space(:,parameters_to_infer{i-1}),...
            obj.time,obj.chimera_abc_IC, data, theta_0{i-1},...
            err(i-1),err_type,num_samples,num_accepted);

        % Store results
        obj.best_fit_parameters_chimera(i).stats = ...
            best_fit_parameters.stats;
        obj.best_fit_parameters_chimera(i).accepted = ...
            best_fit_parameters.accepted;
        obj.best_fit_parameters_chimera(i).sse = ...
            best_fit_parameters.sse;
        obj.best_fit_parameters_chimera(i).sampled = ...
            best_fit_parameters.sampled;

    end

end %======================================================================


function plotChimeraParameters(obj) %======================================

    % Plot matrix of parameters -------------------------------------------
    parameter_name = {'\alpha_D','a','n'};
    model_parameter_number = {1, [1 2]};

    for k = 1:numel(model_parameter_number)
        accepted = [obj.best_fit_parameters_chimera(k+1).accepted,...
            obj.best_fit_parameters_chimera(k+1).sse(:,1)];
        n_parameters = 5;%size(accepted,2);
        H = figure(1);
        clf 
        hold on
        x_lim = [...
            obj.chimera_abc_parameter_space(:,model_parameter_number{k}),...
            [min(obj.best_fit_parameters_chimera(k+1).sse(:,1));...
            max(obj.best_fit_parameters_chimera(k+1).sse(:,1))] ]';
        for i = 1:numel(model_parameter_number{k})
    
            % plot univariate distribution 
            subplot(n_parameters,n_parameters,(i-1)*n_parameters + i)    
           
            hold on
            % x_ticks = min(x_lim):10^floor(log10(diff(x_lim))):max(x_lim)
    
            med = median(accepted(:,i));
            
            h = histogram(accepted(:,i),'Normalization','pdf',...
                'EdgeColor','none','FaceColor','k','FaceAlpha',0.6);
            ay = ylim;
            plot([med med], ay,'k','linewidth',2)        
            
            xlim(x_lim(i,:))
            ylim([0 max(h.Values)])
    
            axis square
            grid on
            
            title(parameter_name{i})
            ylabel('pdf')
    
            for j = 1:i-1
    
                % plot bivariate distribution
                subplot(n_parameters,n_parameters,(i-1)*n_parameters + j)
                scatter(accepted(:,j),accepted(:,i),5,'filled',...
                    'MarkerFaceColor','k','MarkerFaceAlpha',0.5,...
                    'MarkerEdgeColor','none');
                axis([x_lim(j,:) x_lim(i,:)])
                axis square
                grid on
    
            end
        end        
    
        if strcmp(obj.best_fit_parameters_chimera(k+1).model,...
                'Different Growth Rate')
            save_name = 'differentialGrowth';
        else
            save_name = 'differentialGrowthAndCrowding';
        end

        fullFileName = fullfile(obj.results_directory,'abc_mcmc_chimera',...
            ['parameter_posteriors_', save_name]);
        print(H, [fullFileName, '.tif'],'-dtiffn')
        print(H, [fullFileName, '.svg'], '-vector', '-dsvg')  
    end
       

end %======================================================================


function plotChimeraSimulation(obj) %======================================
       
    % simulate chimeras for different models 
    chimera_endpoint = zeros(5,6,3);

    for i = 1:size(obj.chimera_abc_IC,1)
        for j = 1:numel(obj.best_fit_parameters_chimera)   
            
            if j == 1
                parameter_set = obj.best_fit_parameters_chimera(j).stats(1,:);
            else
                parameter_set = ...
                    [obj.best_fit_parameters_embryo.stats(1,:)...
                    obj.best_fit_parameters_chimera(j).stats(1,:)];
            end
               
            [xn,~] = model_inference.numerical_solution_chimeras(...
                parameter_set,obj.time,obj.chimera_abc_IC(i,:));
            embryo_endpoint = sum(xn(end,:));
            icm_endpoint = sum(xn(end,3:end));
            endpoint = [xn(end,2) sum(xn(end,[4 6 7])) xn(end,5)];
            endpoint = [endpoint, endpoint(1)./embryo_endpoint,...
                endpoint(2:end)./icm_endpoint];

            chimera_endpoint(i,:,j) = endpoint;
        end
    end

    tissue_color = [0.0 0.0 1.0;  % TE 
                 0.0 1.0 1.0;  % Epiblast
                 1.0 0.0 1.0]; % PrE 
    
    y_lims_num =  [0  50 225; % Cell #            TE 
                0  50 350; % Cell #            Epiblast
                0  10  55]; % Cell #            PrE                 
    
    group_names = {'15-/-','10-/-','Non','10+/+','15+/+'};
    group_index = obj.chimera_crowding(:,1); 
    
    tissue = {'TE', 'EPI', 'PrE', 'TE', 'EPI', 'PrE'}; 
    tissue_index = [2 4 3 7 12 11];
    tissue_color_index = [1 2 3 1 2 3];
    
    % Step Through each tissue type
    for i = 1:numel(tissue_index)
    
        c = tissue_color(tissue_color_index(i),:);
        
        if i <= 3
            y_label = 'Number';
        else
            y_label = 'Fraction';
        end
        
        fig = figure(1);
        clf
        hold on
        
        % plot boxes of empirical data ------------------------------------
        boxplot(obj.chimera_crowding(:,tissue_index(i)),group_index,'Colors','k',...
            'labels',group_names,'Symbol','','Whisker',inf)
        h = gca;
        
        for j = 1:numel(h.Children.Children)            
            h.Children(1).Children(j).LineWidth = 1.5;
            if strcmp(h.Children(1).Children(j).Tag,'Box')                                                                      
                x_data = h.Children(1).Children(j).XData;
                y_data = h.Children(1).Children(j).YData;                                  
                h_fill = fill(x_data,y_data,c,'FaceAlpha',.3);
                uistack(h_fill,'bottom');                                      
            end
        end
    
        % plot lines of simulations ---------------------------------------
         plot([1:5],chimera_endpoint(:,i,1),'s:','Color',0.6*ones(3,1),...
            'markersize',5,'linewidth',3,'MarkerFaceColor',0.6*ones(3,1))
        plot([1:5],chimera_endpoint(:,i,2),'^--','Color',0.3*ones(3,1),...
            'markersize',5,'linewidth',3,'MarkerFaceColor',0.3*ones(3,1))
        plot([1:5],chimera_endpoint(:,i,3),'o-','Color',0*ones(3,1),...
            'markersize',5,'linewidth',3,'MarkerFaceColor',0*ones(3,1))

        set(gca,'XTick',1:5,'XTickLAbel',group_names,'XTickLabelRotation',90)
        
        title(tissue{i})
        ylabel(y_label)
        
        if i <=3
            axis([0.5 5.5 y_lims_num(i,[1 3])])
            set(gca,'YTick',y_lims_num(i,1):y_lims_num(i,2):y_lims_num(i,3))
            save_name = 'number';
        else
            axis([0.5 5.5 0 1])
            set(gca,'YTick',0:0.25:1,'YTickLabel',num2str([0:0.25:1]','%0.2f'))
            save_name = 'fraction';
        end

        grid on

        set(gca,'FontSize',12,'FontWeight','bold','FontName','Arial')
        set(gcf, 'Units', 'Normalized', 'InnerPosition', [0.6292    0.2741    0.1354    0.3954])

        save_name = [save_name,'_',tissue{i}];
        fullFileName = fullfile('results','abc_mcmc_chimera',save_name);
        print(fig, fullFileName,'-dtiffn')
        print(fig, fullFileName, '-vector', '-dsvg')

    end


end %======================================================================


function validateHostEpiblastExclusion(obj) %==============================
       
     close all
    
    save_folder_name = fullfile(obj.results_directory,'model_validation');
    if ~exist(save_folder_name,'dir')
        mkdir(save_folder_name)
    end

    % simulate chimeras for different models 
    chimera_endpoint = zeros(3,8,3);

    for i = 1:3
        for j = 1:numel(obj.best_fit_parameters_chimera)   
            
            if j == 1
                parameter_set = obj.best_fit_parameters_chimera(j).stats(1,:);
            else
                parameter_set = ...
                    [obj.best_fit_parameters_embryo.stats(1,:)...
                    obj.best_fit_parameters_chimera(j).stats(1,:)];
            end
               
            [xn,~] = model_inference.numerical_solution_chimeras(...
                parameter_set,obj.time,obj.chimera_abc_IC(i+1,:));
            embryo_endpoint = sum(xn(end,:));
            icm_endpoint = sum(xn(end,3:end));
            endpoint = [xn(end,2) sum(xn(end,4)) sum(xn(end,[6 7])) xn(end,5)];
            endpoint = [endpoint, endpoint(1)./embryo_endpoint,...
                endpoint(2:end)./icm_endpoint];

            chimera_endpoint(i,:,j) = endpoint;
            
        end
    end

    data = [obj.chimera_exclusion_clustering.groups...
        obj.chimera_exclusion(:,2:end)];

    host_epiblast = data(:,[1 4 10]);
    host_epiblast(host_epiblast(:,1)==2,:) = [];
    host_epiblast(host_epiblast(:,1)~=1,1) = ...
        host_epiblast(host_epiblast(:,1)~=1,1) - 1;

    simulated_host_epiblast = chimera_endpoint(:,[2 6],:);

    c = [0.0 0.5 0.5];
    group_names = {'10-/-', 'Non','10+/+'};
    y_labels = {'Number', 'Fraction'};
    
    for i = 1:2

         y_label = y_labels{i};

        fig = figure(1);
        clf
        hold on
        
        boxplot(host_epiblast(:,i+1),host_epiblast(:,1),'Colors','k',...
            'labels',group_names,'Symbol','','Whisker',inf)
        h = gca;

       for j = 1:numel(h.Children.Children)            
            h.Children(1).Children(j).LineWidth = 1.5;
            if strcmp(h.Children(1).Children(j).Tag,'Box')                                                                      
                x_data = h.Children(1).Children(j).XData;
                y_data = h.Children(1).Children(j).YData;                                  
                h_fill = fill(x_data,y_data,c,'FaceAlpha',.3);
                uistack(h_fill,'bottom');                                      
            end
        end

        % plot lines of simulations ------------------------------------
         plot([1:3],simulated_host_epiblast(:,i,1),'s:','Color',0.6*ones(3,1),...
            'markersize',5,'linewidth',3,'MarkerFaceColor',0.6*ones(3,1))
        plot([1:3],simulated_host_epiblast(:,i,2),'^--','Color',0.3*ones(3,1),...
            'markersize',5,'linewidth',3,'MarkerFaceColor',0.3*ones(3,1))
        plot([1:3],simulated_host_epiblast(:,i,3),'o-','Color',0*ones(3,1),...
            'markersize',5,'linewidth',3,'MarkerFaceColor',0*ones(3,1))

        title('Host EPI')
        ylabel(y_label)

        if i == 1
            axis([0.5 3.5 0 35])
            set(gca,'YTick',0:10:35)
            save_name = 'number';
        else
            axis([0.5 3.5 0 1])
            set(gca,'YTick',0:0.25:1,'YTickLabel',num2str([0:0.25:1]','%0.2f'))
            save_name = 'fraction';
        end

        grid on

        set(gca,'FontSize',12,'FontWeight','bold','FontName','Arial')
        set(gcf, 'Units', 'Normalized', 'InnerPosition', [0.6292    0.2741    0.1354    0.3954])

        save_name = [save_name,'_hostEpiblast'];
        fullFileName = fullfile(save_folder_name,save_name);
        print(fig, fullFileName,'-dtiffn')
        print(fig, fullFileName, '-vector', '-dsvg')

    end

end %======================================================================


end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
methods(Static) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
function import_litter_mate_data(path_name) %==============================

    addpath(path_name)

    %  Setup the Import Options
    opts = delimitedTextImportOptions('NumVariables', 36);

    % Specify range and delimiter
    opts.DataLines = [2, Inf];
    opts.Delimiter = ',';

    % Specify column names and types
    opts.VariableNames = {'Embryo_ID', 'Experiment', 'Experimenter', ...
        'Regime', 'Treatment', 'Tt_length', 'Cellcount', 'Cell_ID', ...
        'Identity', 'Size', 'X', 'Y', 'Z', 'CH1Avg', 'CH1Sum', 'CH2Avg', ...
        'CH2Sum', 'CH3Avg', 'CH3Sum', 'CH4Avg', 'CH4Sum', 'CH5Avg', ...
        'CH5Sum', 'TE_ICM', 'Xpoint', 'Stage', 'Markers', 'CH1logCor',...
        'CH2logCor', 'CH3logCor', 'CH4logCor', 'CH5logCor', ...
        'CH1ebLogCor', 'CH4ebLogCor', 'CH5ebLogCor', 'Identitykm'};
    opts.VariableTypes = {'categorical', 'categorical', 'categorical', ...
        'double', 'categorical', 'double', 'double', 'double', ....
        'categorical', 'double', 'double', 'double', 'double',...
        'double', 'double', 'double', 'double', 'double', 'double',...
        'double', 'double', 'double', 'double', 'categorical', ...
        'categorical', 'categorical', 'categorical', 'double', 'double',...
        'double', 'double', 'double', 'double', 'double', 'double',...
        'categorical'};
    opts = setvaropts(opts, [4, 6], 'TrimNonNumeric', true);
    opts = setvaropts(opts, [4, 6], 'ThousandsSeparator', ',');
    opts = setvaropts(opts, [1, 2, 3, 5, 9, 24, 25, 26, 27, 36],...
        'EmptyFieldRule', 'auto');
    opts.ExtraColumnsRule = 'ignore';
    opts.EmptyLineRule = 'read';
    % Import the data
    FGFallpooledtrans = ...
        readtable('\data\saiz-et-al_2016-master\FGF_all_pooled_trans.csv',...
        opts);

    %  Clear temporary variables
    clear opts
    rmpath(path_name)

    % Save table in data folder
    save('data\FGFallpooledtrans','FGFallpooledtrans')

end %======================================================================


function plot_chimera_exclusion_clustering(data, varargin) %===============
    x_lims = varargin{1};
    binWidth = varargin{2};
    histEdgeColor = varargin{3};
    faceColor = varargin{4};
    subplotNum = varargin{5};
    titleLabel = varargin{6};
    yLabelName = varargin{7};
    normStyle = varargin{8};           
    
    subplot(subplotNum)
    histogram(data,'FaceColor',faceColor,...
        'binwidth',binWidth ,'DisplayStyle','Bar',...
        'EdgeColor',histEdgeColor,'linewidth',.1,...
        'Normalization',normStyle)
    hold on
    histogram(data,'binwidth',binWidth ,'DisplayStyle','Stairs',...
        'EdgeColor','k','lineWidth',3,'Normalization',normStyle)
    xlim(x_lims)
    title(titleLabel)
    ylabel(yLabelName)
    xlabel('Number of Host Epiblast')
    set(gca, 'FontSize',15)                       
end %======================================================================


function plot_GMM(GMModel,null_injected,v,ax_lims,titleName,subplotNum)  %=
    
    subplot(subplotNum)
    gmPDF = @(x1,x2)reshape(pdf(GMModel,[x1(:) x2(:)]),size(x1));
    h_con = fcontour(gmPDF,ax_lims,'Fill','on');
      hold on;  
    h_gs = gscatter(null_injected(:,1),null_injected(:,2),...
        ~v,[1 0 1;0 1 1],[],30);
    
    colormap([repmat(linspace(1,0,100)',1,3)])
    caxis([0,6e-3])
    axis(ax_lims)
    legend('GMM','10-/- Hi','10-/- Lo')
    title(titleName)
    ylabel 'Donor Derived Epiblast'
    xlabel 'Host Derived Epiblast'
    set(gca, 'FontSize',15)
    axis equal

end %======================================================================


function best_fit_parameters = abc_mcmc_embryo(parameter_space,time,IC,data,varargin) 
   
    %ABC_MCMC 
    % 
    %
    % Optional Arguments
    % THETA_0 initial guess for 
    % ERR for lower upper-bound on accepted error
    % ERR_TYPE is set to constant
    % NUM_SAMPLES Number of attempts at getting the number accepted
    % NUM_ACCEPTED number of accepted values desired
    
    % Record the number of parameters
    num_parameters = size(parameter_space,2);
    
    pre_feedback = varargin{end};
    %----------------------------------------------------------------------
    % Default Parameters set up
    
    numvarargs = length(varargin);
    
    assert(numvarargs < 7,'requires at most 5 optional inputs');
    theta_sampled = [];
    % Sample THETA_0 from the prior distribution (assumed to be uniform)
    min_matrix = parameter_space(1,:);
    max_matrix = parameter_space(2,:);
    theta_0 = min_matrix + (max_matrix-min_matrix).*rand(1,num_parameters);
    
    % Set defaults for optional inputs
    %optargs =  {theta_0, err, err_type, num_samples, num_accepted}
    optargs = {theta_0 2e4 'const' 1e6 500};
    
    % Get the index of nonempty VARARGIN
    varargin_index = find(~cellfun(@isempty,varargin));
    
    % Now put the specified values from varargins into the valuesToUse cell
    % array, and overwrite the defaults.
    [optargs{varargin_index}] = varargin{varargin_index};
    
    % Place optional args in memorable variable names
    [theta_0, err, err_type, num_samples, num_accepted] = optargs{:};
  
    % Add one because the first guess will be  deleted at the end
    num_accepted = num_accepted + 1;

    %----------------------------------------------------------------------
    % Set up storage vectors
    
    % Storage for accepted parameter sets
    accepted = zeros(num_accepted, num_parameters); 
    % initialize with THETA_0 however, this will be overwritten with the
    % first accepted parameter set
    accepted(1,:) = theta_0 ;
    
    % Storage for error and err threshold
    error_parameters = NaN(num_accepted, 2);
    % initialize with SUM_SQUARED_ERROR by calculating sse for THETA_0.
    % However, this will be overwritten with the first accepted parameter
    % set.
    if strcmp(pre_feedback,'no')
        theta_0(end-1) = 0;
    end

    sse = model_inference.sum_squared_error_embryo(theta_0,time,IC,data,'yes');
    error_parameters(1,:) = [sse err];
            
    %----------------------------------------------------------------------
    % Initalize counters and data
    sample = 1;
    accepted_counter = 1;
    
    %----------------------------------------------------------------------
    % Set up proposal_variance such that it sweeps 1/4 of the interal with
    % a 99.9% CI(ish). I.e. the normally distributed variable will look
    % within 1/4 the distace of the the parameter space interval away from
    % the current point with 99.9% confidence.
    proposal_variance = diff(parameter_space)/(2*3.29);
    
    % Start looking over the parameter space until you hit the desired
    % number of acceptances or until
    while accepted_counter < num_accepted && sample < num_samples
        
        % Propose theta^* according to a proposal distribution
        theta_star = accepted(accepted_counter,:) ...
            + proposal_variance.*randn(1,num_parameters);
        theta_star(end-1:end) = randi(3,1,2)-1; 
        
        % Ensure that theta_star is within the parameter_space
        while any(parameter_space(1,:) >= theta_star) ...
                || any(parameter_space(2,:) <= theta_star) 
            theta_star = accepted(accepted_counter,:) ...
                + proposal_variance.*randn(1,num_parameters);
            theta_star(end-1:end) = randi(3,1,2)-1; 
        end

        if strcmp(pre_feedback,'no')
            theta_star(end-1) = 0;
        end

        theta_sampled = [theta_sampled;theta_star];
        % Simulate a dataset x^* using theta^* and compute d(x_0, x^*) =
        % sse.
        sse = model_inference.sum_squared_error_embryo(theta_star,...
            time,IC,data,'yes');
        
        % If sse < epsilon, record theta_{i} set theta_{i+1} to theta^*
        if sse < err
                    
            % Update accepted_counter
            accepted_counter = accepted_counter + 1;
            
            % Record theta^*, sse, and current err
            accepted(accepted_counter,:) = theta_star;
            error_parameters(accepted_counter,:) = [sse err]; 
            
            % Determine if err update is required
            switch err_type{:}
                case 'const'
                    0;
                case 'sse'
                    err = min(err,sse);
            end
            
            % Show the user things are still running
            clc
            disp(['accepted counter = ' num2str(accepted_counter-1)])
            disp(['sample number = ' num2str(sample)])
            disp(['err = ' num2str(err)])
            disp(['sse = ' num2str(sse)])         
            disp(['theta = [' num2str(theta_star) ']'])
            
        end    
    
        % Update sample counter
        sample = sample + 1;
    
    end
    
    % Show the user the final results 
    disp(['Accepted: ' num2str(accepted_counter-1)])
    disp(['Steps: ' num2str(sample)])
    
    % eliminate any remaining 0s and NaNs from the storage vectors
    error_parameters(~any(accepted,2),:) = [];   
    accepted(~any(accepted,2),:) = [];
    accepted_parameters = accepted(2:end,:);

    % Calculate the mean and Std for the first round
    stats = [median(accepted_parameters);std(accepted_parameters)];
    stats(1,end) = round(stats(1,end));

    best_fit_parameters.pre_feedback = pre_feedback;
    best_fit_parameters.stats = stats;
    best_fit_parameters.accepted = accepted_parameters;
    best_fit_parameters.sse = error_parameters(2:end,:);
    best_fit_parameters.sampled = theta_sampled;
    
end %======================================================================


function [ sse ] = sum_squared_error_embryo(parameters,time,IC,data,feedback)
   
    % Calculate the analytical solution for the times observations occured
    if strcmp(feedback,'yes')
        [xn,~] = model_inference.numerical_solution_pre_feedback(...
            parameters,time,IC);
    elseif strcmp(feedback,'no')
        [xn,~] = model_inference.numerical_solution_sans_pre_feedback(...
            parameters,time,IC);
    end

    length_x = size(xn,1);

    obs = data(:,[4 1 2 3]);
    del_index = sum(obs,2)>150;
    obs(del_index,:) = [];
    obs = [[8 0 0 0 0] ;zeros(size(obs,1),1) obs];

    % Initalize error vector
     sse = 0;
   
    % IMPOSE PENALTY 
    % Max embryo size
     sse  = (max(sum(xn,2)) - max(sum(obs,2))).^2;
    % Max TE size
     sse = sse + (max(xn(:,2)) - max(obs(:,2))).^2;
    % Max ICM size
     sse = sse + (max(xn(:,3)) - max(obs(:,3))).^2;
    % Max Epiblast
    sse = sse + (max(xn(:,4)) - max(obs(:,4))).^2;
    % Max PrE
    sse = sse + (max(xn(:,5)) - max(obs(:,5))).^2;
    % t_final should have not B or ICM
    sse = sse +  10*xn(end,1) +  10*xn(end,3);

    % Get the number of data points
    num_obs = size(obs,1);

    % Step through data points
    for k=1:num_obs
        % Find the minimal distance between this observation and the model
        % solution        
        squared_distances = sum((xn - repmat([obs(k,:)],length_x,1)).^2,2);

        sse = sse + min(squared_distances);
    end

end %======================================================================


function [xn,tn] = numerical_solution_pre_feedback(parameters,time,IC) %===
%NUMERICAL_SOLUTION solves a system of coupled ODEs, using ode45, that
%describe the time evolution of cell species in the ealry embryo from the
%E2.5 8-cell stage morula to the E4.5 late blastocyst by taking in
%arguments 1X6 row vector PARAMETERS, Nx1 column vector TIME, and 1X5 row
%vector IC. Then Returns the numerical solution XN over the time interval
%TN.
%
% PARAMETERS are {alpha,beta,rho, zeta, and eta} which are {effective
% growth rate, Blastomere specification rate, Blastomere bias to the ICM,
% ICM to Epiblast specification rate, and ICM to Primitive endoderm
% specification rate} defined in the trasition rates,
%
% Proliferation
% X -> X+1, lambda*X.
% X -> X-1,    mui*X.
% s.t. alpha = lambda - mu; Doubling Time Constant
%
% Blastomere specification
% B -> T, (1-rho)*beta*B.
% B -> C,    rho *beta*B.
%
% ICM specification
% C -> E,        zeta*P*C.
% C -> P, eta*(C + E)*C.
%
% TIME is a vector over which ODE45 will solve the system of equations.
% This should be from 0:48 of varying step size (E.G. TIME =
% LINSPACE(0,48,501)).
%
% IC are the initial conditions for the stat variables x(1:7) = {B,T,C,E,P}
% which are the number of {Blastomere, Trophectoderm, ICM, Epiblast,
% Primitive Endoderm, FGF-WT Donor ES cells, and FGF-null Donor ES cells}
% cells present at TIME = 0. This will be [8 0 0 0 0] where the last two
% entries are subsets of the integers

% time vector, make sure time is a column vector
[a,b] = size(time);
if a<b
   time = time' ;
end

% Write parameters for human eyes
alpha = parameters(1); % Proliferation Rate
beta = parameters(2); % B -> T, basal rate
rho = parameters(3); % B -> C
zeta = parameters(4); % C -> E 
eta = parameters(5); % C -> P
l = parameters(6); % PrE signal feedback to ICM (??)
m = parameters(7); % FGF4 signal feedback from EPI and ICM to ICM

% Set up ODEs function
if l == 0 % No PrE Feedback    
    f = @(t,x) [(alpha - beta).*x(1);...dB/dt
     alpha.*x(2) + (1-rho).*beta.*x(1);...dT/dt                
    (alpha - (zeta + eta.*(x(3)+x(4)).^m)).*x(3) + rho.*beta*x(1);...% dC/dt
     alpha.*x(4) + zeta.*x(3);...% dE/dt
     alpha.*x(5) + eta.*((x(3) + x(4)).^m).*x(3)];% dP/dt  
else % PrE feedback
    f = @(t,x) [(alpha - beta).*x(1);...dB/dt
     alpha.*x(2) + (1-rho).*beta.*x(1);...dT/dt                
    (alpha - (zeta.*x(5).^l + eta.*(x(3)+x(4)).^m)).*x(3) + rho.*beta*x(1);...% dC/dt
     alpha.*x(4) + zeta.*(x(5).^l).*x(3);...% dE/dt
     alpha.*x(5) + eta.*((x(3) + x(4)).^m).*x(3)];% dP/dt  
end   
         
% Numerically solve system of ODEs
[tn,xn] =  ode45(f,time,IC);

end %======================================================================


function best_fit_parameters = abc_mcmc_chimera(embryo_parameters,parameter_space,time,IC,data,varargin) 
   
    num_parameters = size(parameter_space,2);
    theta_sampled = [];
    % Sample THETA_0 from the prior distribution (assumed to be uniform)
    min_matrix = parameter_space(1,:);
    max_matrix = parameter_space(2,:);    
    theta_0 = min_matrix + (max_matrix-min_matrix).*rand(1,num_parameters);

    % Set defaults for optional inputs
    %optargs =  {theta_0, err, err_type, num_samples, num_accepted}
    optargs = {theta_0 2e4 'const' 1e6 500};
    
    % Get the index of nonempty VARARGIN
    varargin_index = find(~cellfun(@isempty,varargin));
    
    % Now put the specified values from varargins into the valuesToUse cell
    % array, and overwrite the defaults.
    [optargs{varargin_index}] = varargin{varargin_index};
    
    % Place optional args in memorable variable names
    [theta_0, err, err_type, num_samples, num_accepted] = optargs{:};

    % Add one because the first guess will be  deleted at the end
    num_accepted = num_accepted + 1;
    
    %----------------------------------------------------------------------
    % Set up storage vectors
    
    % Storage for accepted parameter sets
    accepted = zeros(num_accepted, num_parameters); 
    % initialize with THETA_0 however, this will be overwritten with the
    % first accepted parameter set
    accepted(1,:) = theta_0;
    
    % Storage for error and err threshold
    error_parameters = NaN(num_accepted, 2);
    % initialize with SUM_SQUARED_ERROR by calculating sse for THETA_0.
    % However, this will be overwritten with the first accepted parameter
    % set.
    sse = model_inference.sum_squared_error_chimera(...
        [embryo_parameters,theta_0],time,IC,data);
    error_parameters(1,:) = [sse err];

    %----------------------------------------------------------------------
    % Initalize counters and data
    sample = 1;
    accepted_counter = 1;

    %----------------------------------------------------------------------
    % Set up proposal_variance such that it sweeps 1/4 of the interal with
    % a 99.9% CI(ish). I.e. the normally distributed variable will look
    % within 1/4 the distace of the the parameter space interval away from
    % the current point with 99.9% confidence.

    if numel(theta_0) == 1
        proposal_variance = diff(parameter_space)/(0.25*3.29);
    else
        proposal_variance = diff(parameter_space)/(0.5*3.29);
    end
        
    % Start looking over the parameter space until you hit the desired
    % number of acceptances or until
    while accepted_counter < num_accepted && sample < num_samples

         % Propose theta^* according to a proposal distribution
        theta_star = accepted(accepted_counter,:) ...
            + proposal_variance.*randn(1,num_parameters);
        if numel(theta_star) == 3
            theta_star(end) = randi(3,1)-1;
        end

        % Ensure that theta_star is within the parameter_space
        while any(parameter_space(1,:) >= theta_star) ...
                || any(parameter_space(2,:) <= theta_star) 
            theta_star = accepted(accepted_counter,:) ...
                + proposal_variance.*randn(1,num_parameters);
            if numel(theta_star) == 3
                theta_star(end) = randi(3,1)-1;
            end
        end
        theta_sampled = [theta_sampled;theta_star];
        % Simulate a dataset x^* using theta^* and compute d(x_0, x^*) =
        % sse.
        sse = model_inference.sum_squared_error_chimera(...
            [embryo_parameters, theta_star],time,IC,data);
        
        % If sse < epsilon, record theta_{i} set theta_{i+1} to theta^*
        if sse < err
                    
            % Update accepted_counter
            accepted_counter = accepted_counter + 1;
            
            % Record theta^*, sse, and current err
            accepted(accepted_counter,:) = theta_star;
            error_parameters(accepted_counter,:) = [sse err]; 
            
            % Determine if err update is required
            switch err_type{:}
                case 'const'
                    0;
                case 'sse'
                    err = min(err,sse);
            end
            
            % Show the user things are still running
            clc
            disp(['accepted counter = ' num2str(accepted_counter-1)])
            disp(['sample number = ' num2str(sample)])
            disp(['err = ' num2str(err)])
            disp(['sse = ' num2str(sse)])         
            disp(['theta = [' num2str(theta_star) ']'])
            
        end    

        % Update sample counter
        sample = sample + 1;

    end
    
     % Show the user the final results 
    disp(['Accepted: ' num2str(accepted_counter-1)])
    disp(['Steps: ' num2str(sample)])
    
    % eliminate any remaining 0s and NaNs from the storage vectors
    error_parameters(~any(accepted,2),:) = [];   
    accepted(~any(accepted,2),:) = [];
    accepted_parameters = accepted(2:end,:);

    % Calculate the mean and Std for the first round
    stats = [median(accepted_parameters);std(accepted_parameters)];    

    best_fit_parameters.stats = stats;
    best_fit_parameters.accepted = accepted_parameters;
    best_fit_parameters.sse = error_parameters(2:end,:);
    best_fit_parameters.sampled = theta_sampled;

end %======================================================================


function [ sse ] = sum_squared_error_chimera(parameters,time,IC,data)

    % Initalize error vector
    sse = 0;
    
    % Step through each injection type and calculte the sse of the
    % simulation endpoint and the empirical data
    for i = [1 2 4 5]
    
        % get the empirical data for that injection type
        observed_index = data(:,1) == i;
        current_data = data(observed_index,2:4);
        obs_length = size(current_data,1);
        
        % Simulate
        [xn,~] = model_inference.numerical_solution_chimeras(parameters,time,IC(i,:));
    
        % collect the [TE PrE EPI+D^{+,-}]
        xn_comparison = [xn(end,2) xn(end,5) sum(xn(end,[4 6 7]),2)];
        xn_comparison = repmat(xn_comparison,obs_length,1);
    
        squared_error =  (current_data - xn_comparison).^2;
        sse = sse  +  sum(squared_error(:));
           
    end

end %======================================================================


function [xn,tn] = numerical_solution_chimeras(parameters,time,IC)
%NUMERICAL_SOLUTION_chimeras solves a system of coupled ODEs, using ode45,
%that describe the time evolution of cell species of an ealry embryo
%injected with varying numbers of fgf4-WT D^+ and fgf4-null D^- doner ES
%cells then allowed to develop from the E2.5 8-cell stage morula to the
%E4.5 late blastocyst by taking in arguments 1X6 row vector PARAMETERS, Nx1
%column vector TIME, and 1X7 row vector IC. Then Returns the numerical
%solution XN over the time interval TN.
%
% PARAMETERS are {alpha,beta,rho, zeta, eta, a, and alpha_D} which are {doubling time
% constant, Blastomere specification rate constant, Blastomere bias to the
% ICM, ICM to Epiblast specification rate constant, and ICM to Primitive
% endoderm specification rate constant, spatial crowding term} defined in the trasition rates,
%
% Proliferation
% X -> X+1, lambda*X.
% X -> X-1,    mui*X.
% s.t. alpha = lambda - mu; Doubling Time Constant
%
% Blastomere specification
% B -> T, (1-rho)*beta*B.
% B -> C,    rho *beta*B.
%
% ICM specification
% C -> E,        zeta*P*C.
% C -> P, eta*(C + E + D)*C.
%
% TIME is a vector over which ODE45 will solve the system of equations.
% This should be from 0:48 of varying step size (E.G. TIME =
% LINSPACE(0,48,501)).
%
% IC are the initial conditions for the stat variables x(1:7) =
% {B,T,C,E,P,D^+,^D-} which are the number of {Blastomere, Trophectoderm,
% ICM, Epiblast, Primitive Endoderm, FGF-WT Donor ES cells, and FGF-null
% Donor ES cells} cells present at TIME = 0. This will be [8 0 0 0 0 n^+
% n^-] where the last two entries are the different number of donor clls

% time vector, make sure time is a column vector
[p,q] = size(time);
if p<q
   time = time' ;
end

% Write parameters for human eyes

% fixed values determined by embryo ABC
alpha = parameters(1); % Proliferation Rate
beta = parameters(2); % B -> T, basal rate
rho = parameters(3); % B -> C, bias
zeta = parameters(4); % C -> E, rate
eta = parameters(5); % C -> P, rate
l = parameters(6); % PrE signal feedback to ICM parameter (??)
m = parameters(7); % EPI signal feedback to ICM  parameter (FGF4)

% Current Parameters being fit by ABC-MCMC
alpha_D = parameters(8); % Donor cell proliferation rate
if numel(parameters) == 8
    a = 0;
    n = 0; % Donor cell feedback to Blastomeres (crowding)
elseif numel(parameters) == 10
    a = parameters(9); % spatial crowding constant
    n = parameters(10); % Donor cell feedback to Blastomeres (crowding)
end

% Set up ODEs function
if a == 0 % FGF4 Signalling
    f = @(t,x) [(alpha - beta).*x(1);...dB/dt
     alpha.*x(2) + (1-rho).*beta.*x(1);...dT/dt                
    (alpha - (zeta.*(x(5).^l) + eta.*(x(3) + x(4) + x(6)).^m)).*x(3) + rho.*beta*x(1);...% dC/dt
     alpha.*x(4) + zeta.*(x(5).^l).*x(3);...% dE/dt
     alpha.*x(5) + eta.*((x(3) + x(4) + x(6)).^m).*x(3);...% dP/dt         
     alpha_D.*x(6);...% dD^+/dt
     alpha_D.*x(7)];% dD^-/dt    
else % FGF4 Signalling and Spatial Crowding    
    f = @(t,x) [(alpha - beta).*x(1);...dB/dt
     alpha.*x(2) + (1-rho./(1+a.*((x(6) + x(7)).^n))).*beta.*x(1);...dT/dt                
    (alpha - (zeta.*(x(5).^l) + eta.*(x(3) + x(4) + x(6)).^m)).*x(3) + (rho./(1+a.*((x(6) + x(7)).^n))).*beta*x(1);...% dC/dt
     alpha.*x(4) + zeta.*(x(5).^l).*x(3);...% dE/dt
     alpha.*x(5) + eta.*((x(3) + x(4) + x(6)).^m).*x(3);...% dP/dt         
     alpha_D.*x(6);...% dD^+/dt
     alpha_D.*x(7)];% dD^-/dt    
end

% Numerically solve system of ODEs
[tn,xn] =  ode45(f,time,IC);

end %======================================================================

end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% 
% original priors
% embryo_abc_priors = struct('alpha', [.05 0.08],... % Proliferation rate constant
%      'beta', [.01 10],...    % B -> T, basal rate constant
%       'rho', [eps 1],...     % B -> C, lineage bias
%      'zeta', [1e-4 1e-2],... % C -> E, rate constant [1e-5 1e-2]
%       'eta', [1e-5 1e-3],... % C -> P, rate constant
%         'l', [-0.5 2.5],...    % PrE signal feedback to ICM (??)
%         'm', [-0.5 2.5]); 
% chimera_abc_priors = struct('alpha_D', [1e-3 1e-1],... % Donor Cell proliferation rate
%      'a', [1e-2 10],...% spatial crowding constant
%      'n', [-0.5 2.5]);   % Donor cell feedback to Blastomeres (crowding)


% x_axis_scale = {'linear','linear','linear','linear','linear','linear','linear'};
% num_parameters = size(obj.embryo_abc_parameter_space,2);
% 
% for i = 1:num_parameters
% 
%     H = figure(1);
%     clf
%     hold on
% 
%     accepted_parameters = ...
%             obj.best_fit_parameters_embryo.accepted(:,i);
%     x_lim = obj.embryo_abc_parameter_space(:,i)';
%     x_ticks = min(x_lim):10^floor(log10(diff(x_lim))):max(x_lim);        
%     x_tick_labels = string(x_ticks');
% 
%     if strcmp(parameter_name{i},'\rho')
%         x_tick_labels = num2str(x_ticks','%0.1f');
%     elseif strcmp(parameter_name{i},'l') || strcmp(parameter_name{i},'m')
%         x_lim = ceil(obj.embryo_abc_parameter_space(:,i)');
%         x_ticks = min(x_lim)-0.:10^floor(log10(diff(x_lim))):max(x_lim);        
%         x_tick_labels = string(x_ticks');
%         x_lim = [x_lim(1)-0.5, x_lim(2)+1];
%     end
% 
%     if strcmp(x_axis_scale{i},'log')
%         accepted_parameters = log10(accepted_parameters);
%         x_lim = log10(x_lim);
%         x_ticks = min(x_lim):max(x_lim);
%         x_tick_labels = strcat('10^{',string(x_ticks'),'}');
%     end
% 
%     med = median(accepted_parameters);
% 
%     h = histogram(accepted_parameters,'Normalization','pdf',...
%         'EdgeColor','none','FaceColor','k','FaceAlpha',0.6);
% 
%     ay = ylim;
%     plot([med med], ay,'k','linewidth',2)
% 
%     xlim(x_lim);
%     ylim([0 max(h.Values)])
% 
%     grid on
% 
%     title(parameter_name{i})
%     ylabel('pdf')
% 
%     set(gca,'XTick',x_ticks,'XTickLabel',x_tick_labels)
% 
%     set(gca,'FontSize',12,'FontWeight','Bold','FontName','Arial')        
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05 0.05 0.1250 0.5]) 
% 
%     % Saving figure!!!
%     save_name_ending = parameter_name{i};
%     if contains(save_name_ending,'\')
%         save_name_ending = save_name_ending(2:end);
%     end
%     fullFileName = fullfile(obj.results_directory,'abc_mcmc_embryos',...
%         ['parameter_',save_name_ending]);
%     print(H, [fullFileName, '.tif'],'-dtiffn')
%     print(H, [fullFileName, '.svg'],'-dsvg')    
% 
% end       

% parameter_name = {'\alpha_D','a','n'};
% x_axis_scale = {'linear','linear','linear'};
% model_parameter_number = {[1 3] [1 2 3]};
% 
% for i = 1:2
% 
%     num_parameters = numel(model_parameter_number{i});
% 
%     for j = 1:num_parameters
% 
%         H = figure(1);
%         clf
%         hold on
% 
%         if strcmp(x_axis_scale{model_parameter_number{i}(j)},'log')
%             accepted_parameters = ...
%                 log10(obj.best_fit_parameters_chimera(i+1).accepted(:,j));
% 
%             x_lim = log10(obj.chimera_abc_parameter_space(:,j)');
%             x_ticks = min(x_lim):max(x_lim);
%             x_tick_labels = strcat('10^{',string(x_ticks'),'}');
%         else
%             accepted_parameters = ...
%                 obj.best_fit_parameters_chimera(i+1).accepted(:,j);
%             x_lim = obj.chimera_abc_parameter_space(:,model_parameter_number{i}(j))';
%             x_ticks = ceil(min(x_lim)):floor(max(x_lim));        
%             x_tick_labels = string(x_ticks');
%         end
% 
%         med = median(accepted_parameters);
% 
%         h = histogram(accepted_parameters,'Normalization','pdf',...
%             'EdgeColor','none','FaceColor','k','FaceAlpha',0.6);
% 
%         ay = ylim;
%         plot([med med], ay,'k','linewidth',2)
% 
%         xlim(x_lim);
%         ylim([0 1.0*max(h.Values)])
% 
%         grid on
% 
%         title(parameter_name{model_parameter_number{i}(j)})
%         ylabel('pdf')
% 
%         set(gca,'XTick',x_ticks,'XTickLabel',x_tick_labels)
% 
%         set(gca,'FontSize',12,'FontWeight','Bold','FontName','Arial')        
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05 0.05 0.1250 0.5])
% 
%         % Saving figure!!!
%         save_name_parameter = parameter_name{model_parameter_number{i}(j)};
%         if contains(save_name_parameter,'\')
%             save_name_parameter = save_name_parameter(2:end);
%         end
% 
%         if strcmp(obj.best_fit_parameters_chimera(i+1).model,'Different Growth Rate')
%             save_name = ['diffGrowth_',save_name_parameter];
%         else
%             save_name = ['diffGrowthAndCrowding_',save_name_parameter];
%         end
% 
%         fullFileName = fullfile(obj.results_directory,...
%             'abc_mcmc_chimera',save_name );
%         print(H, [fullFileName, '.tif'],'-dtiffn')
%         print(H, [fullFileName, '.svg'],'-dsvg')  
% 
%     end
% 
% end
% 
%  %TWO DIMENSIONAL CLUSTERING--------------------------------------------
% 
%     % get the number of host derived cells from each embryo for each
%     % injection type
%     noninjected   = obj.chimera_exclusion(noninjected_index,[4 5]);
%     wt_injected   = obj.chimera_exclusion(wt_injected_index,[4 5]);
%     null_injected = obj.chimera_exclusion(null_injected_index,[4 5]);
% 
%     clust_data = null_injected ...
%         - repmat(mean(null_injected),size(null_injected,1),1);
%     clust_data = clust_data ...
%         ./ repmat(std(clust_data),size(null_injected,1),1);
%     n_test_groups = 4;
%     obj.chimera_exclusion_clustering.BIC = zeros(1,n_test_groups);
% 
%     % clust_data = clust_data(:,1)
% 
%     for i = 1:n_test_groups
%         gmm = fitgmdist(clust_data,i, 'Replicates', 1000);
%         obj.chimera_exclusion_clustering.BIC(i) = gmm.BIC;
%     end
% 
%     gmm = fitgmdist(clust_data,...
%         find(obj.chimera_exclusion_clustering.BIC...
%         ==min(obj.chimera_exclusion_clustering.BIC)),...
%         'Replicates', 1000);
%     obj.chimera_exclusion_clustering.groups = gmm.cluster(clust_data);
% 
%     whoIsBigger = diff(gmm.mu(:,1));
%     if whoIsBigger > 0
%         obj.chimera_exclusion_clustering.groups = ...
%             obj.chimera_exclusion_clustering.groups==1;    
%     else
%         obj.chimera_exclusion_clustering.groups = ...
%             obj.chimera_exclusion_clustering.groups==2;
%     end
% 
%     % Plot BIC ------------------------------------------------------------
%     fig = figure(1);
%     clf
%     hold on
% 
%     x_ticks = 1:n_test_groups;
%     plot(x_ticks, obj.chimera_exclusion_clustering.BIC,'k','LineWidth',4)
%     scatter(gmm.NumComponents, min(obj.chimera_exclusion_clustering.BIC),...
%         250,'filled','r')
% 
%     xlabel('Number of Populations')
%     ylabel('BIC')
%     set(gca,'XTick',x_ticks)
%     xlim([1 n_test_groups])
%     grid on
% 
%     set(gca,'FontName','Arial','FontWeight','bold','FontSize',12)
% 
%     saveFileName = ['bic'];
%     fullFileName = fullfile(save_folder_name,saveFileName);
%     print(fig, fullFileName,'-dtiffn')
%     print(fig, fullFileName,'-vector', '-dsvg')
%     savefig(fig, fullFileName)
%     close(fig)
% 
%     % Plot scatter plots and clustering -----------------------------------
%     fig = figure(1);
%     clf
%     hold on
% 
%     ax_lims = [0 45 0 60];
% 
%     scatter(noninjected(:,1),noninjected(:,2),100,...
%         [0.25 0.25 0.25],'filled','MarkerFaceAlpha',0.5);    
%     scatter(wt_injected(:,1),wt_injected(:,2),100,...
%         [0.70 0.14 0.14],'filled','MarkerFaceAlpha',0.5);
% 
%     axis equal
%     axis(ax_lims)
% 
%     legend('Noninjected','10+/+')
% 
%     title 'Epiblast Composition'
%     ylabel 'Donor Derived'
%     xlabel 'Host Derived'
% 
%     grid on
% 
%     set(gca,'FontName','Arial','FontWeight','bold','FontSize',12)
% 
%     saveFileName = ['scatter_epiblast_composition'];
%     fullFileName = fullfile(save_folder_name,saveFileName);
%     print(fig, fullFileName,'-dtiffn')
%     print(fig, fullFileName,'-vector', '-dsvg')
%     savefig(fig, fullFileName)
%     close(fig)
% 
%     fig = figure(1);
%     clf
%     hold on
% 
%     ax_lims = [0 45 0 60];
% 
%     h_lo = scatter(null_injected(obj.chimera_exclusion_clustering.groups==0,1),...
%         null_injected(obj.chimera_exclusion_clustering.groups==0,2),...
%         100, 'k','filled','MarkerFaceAlpha',0.5);
%     h_lo = scatter(null_injected(obj.chimera_exclusion_clustering.groups==1,1),...
%         null_injected(obj.chimera_exclusion_clustering.groups==1,2),...
%         100,[1 0.549 0],'filled','MarkerFaceAlpha',0.5);    
%     axis equal
%     axis(ax_lims)
% 
%     legend('10-/- Lo','10-/- Hi')
% 
%     title 'Epiblast Composition'
%     ylabel 'Donor Derived'
%     xlabel 'Host Derived'
% 
%     grid on
% 
%     set(gca,'FontName','Arial','FontWeight','bold','FontSize',12)
% 
%     saveFileName = ['scatter_epiblast_composition_clustering'];
%     fullFileName = fullfile(save_folder_name,saveFileName);
%     print(fig, fullFileName,'-dtiffn')
%     print(fig, fullFileName,'-vector', '-dsvg')
%     savefig(fig, fullFileName)
%     close(fig)
% 
% 
% 
% 
%     %ONE DIMENSIONAL CLUSTERING--------------------------------------------    
%     % Get index of different conditions
%     noninjected_index = obj.chimera_exclusion(:,1)==2;
%     wt_injected_index = obj.chimera_exclusion(:,1)==3;
%     null_injected_index = obj.chimera_exclusion(:,1)==1;
% 
%     % get the number of host derived cells from each embryo for each
%     % injection type
%     noninjected   = obj.chimera_exclusion(noninjected_index,4);
%     wt_injected   = obj.chimera_exclusion(wt_injected_index,4);
%     null_injected = obj.chimera_exclusion(null_injected_index,4);
% 
%     % Plot the 1D host contribution to the epiblast            
%     % set up axis attributes
%     x_lims = [0 45];
%     binWidth = 3;
%     histEdgeColor = [.2 .2 .2];
% 
%     % Plot raw data
%     fig = figure(1); clf; hold on
%     model_inference.plot_chimera_exclusion_clustering(noninjected,...
%         x_lims,binWidth,histEdgeColor,'r',321,'Noninjected',...
%         'Counts','count')
%     model_inference.plot_chimera_exclusion_clustering(wt_injected,...
%         x_lims,binWidth,histEdgeColor,'b',323,'10+/+','Counts','count')
%     model_inference.plot_chimera_exclusion_clustering(null_injected,...
%         x_lims,binWidth,histEdgeColor,[0.5 0.5 0.5],325,'10-/-',...
%         'Counts','count')
% 
% 
% set(gcf,'Position',[50 50 1800 900])
%     ax_lims = [0 45 0 60];
% 
%     % Perform clustering on the null data set to get the two
%     % populations                       
%     %KMEANS
%     [index, stats] = kmeans(null_injected,2,'Replicates',1000);
%     whoIsBigger = diff(stats(:,1));
%     if whoIsBigger > 0
%         v = index==1;    
%     else
%         v = index==2;
%     end                        
% 
%     % save clustering
%     obj.chimera_exclusion_clustering.one_D = v;
% 
%     %GAUSIAN MIXTURE MODEL
%     x = linspace(0,50,1000)';
%     GMModel_1 = fitgmdist(null_injected,1, 'Replicates', 100);
%     GMModel_2 = fitgmdist(null_injected,2, 'Replicates', 100);
% 
%     % Get Results
%     obj.chimera_exclusion_clustering.one_D_BIC =  ...
%         {['Components = 1: ' num2str(GMModel_1.BIC)];
%          ['Components = 2: ' num2str(GMModel_2.BIC)]};             
% 
%     % Plot Results of clustering
%     figure(1)
%     model_inference.plot_chimera_exclusion_clustering(null_injected,...
%         x_lims,binWidth,histEdgeColor,[.5 .5 .5],322,...
%         '10-/- Mixed population','Density','pdf')
%     hold on
%     h_1C = plot(x,pdf(GMModel_1,x),':','Color',[.25 .25 .25],'linewidth',4);
%     h_2C =plot(x,pdf(GMModel_2,x),'Color',[.25 .25 .25],'linewidth',4);           
%     legend([h_1C,h_2C],{'1 Component GMM','2 Component GMM'})
% 
%     model_inference.plot_chimera_exclusion_clustering(null_injected(v),...
%         x_lims,binWidth,histEdgeColor,[1 0 1],324,...
%         '10-/- High Contribution','Density','pdf')
%     model_inference.plot_chimera_exclusion_clustering(null_injected(~v),...
%         x_lims,binWidth,histEdgeColor,[0 1 1],326,...
%         '10-/- Low Contribution','Density','pdf')
% 
%      set(gcf,'Position',[50 50 900 900])
% 
%      % save and close figure
%     saveFileName = ['one_d_clustering'];
%     fullFileName = fullfile(save_folder_name,saveFileName);
%     print(fig, fullFileName,'-dtiffn')
%     print(fig, fullFileName,'-vector', '-dsvg')
%     savefig(fig, fullFileName)
%     close(fig)
% 
%     %TWO DIMENSIONAL CLUSTERING--------------------------------------------
%     fig = figure(3); clf
%     set(gcf,'Position',[50 50 1800 900])
% 
%     % get the number of host derived cells from each embryo for each
%     % injection type
%     noninjected   = obj.chimera_exclusion(noninjected_index,[4 5]);
%     wt_injected   = obj.chimera_exclusion(wt_injected_index,[4 5]);
%     null_injected = obj.chimera_exclusion(null_injected_index,[4 5]);
% 
%     %KMEANS null_injected 
%     [index, stats] = kmeans([null_injected(:,1:2)*[1.5 0 ;0 1]],2,...
%         'Replicates',1000);
% 
%     % save clustering
%     whoIsBigger = diff(stats(:,1));
%     if whoIsBigger > 0
%         v = index==1;    
%     else
%         v = index==2;
%     end
%     % Get index of different conditions
%     noninjected_index = obj.chimera_exclusion(:,1)==2;
%     wt_injected_index = obj.chimera_exclusion(:,1)==3;
%     null_injected_index = obj.chimera_exclusion(:,1)==1;
% 
%     % get index of embryos with Hi contribution of donor cells to
%     % the Epiblast
%     null_index_number = find(null_injected_index);            
%     index_modifier = [find([noninjected_index | wt_injected_index]); ...
%                             null_index_number(~v)];
% 
%     obj.chimera_exclusion_clustering.groups = obj.chimera_exclusion(:,1);
%     obj.chimera_exclusion_clustering.groups(index_modifier,1) = ...
%         obj.chimera_exclusion_clustering.groups(index_modifier,1)+1;
% 
%     % save clustering
%     obj.chimera_exclusion_clustering.two_D = v;
%     %GAUSIAN MIXTURE MODEL
% 
% 
%     GMModel_1 = fitgmdist(null_injected,1, 'Replicates', 100);
%     GMModel_2 = fitgmdist(null_injected,2, 'Replicates', 100);
%     GMModel_3 = fitgmdist(null_injected,2, 'Replicates', 100);
% 
% 
%     % Get Results
%     obj.chimera_exclusion_clustering.two_D_BIC =  ...
%         {['Components = 1: ' num2str(GMModel_1.BIC)];
%          ['Components = 2: ' num2str(GMModel_2.BIC)]};
% 
%     % Plot Results of clustering
%     subplot(242)
%     h_non = scatter(noninjected(:,1),noninjected(:,2),[], 'r','linewidth',3);
%     hold on
%     h_pos = scatter(wt_injected(:,1),wt_injected(:,2),[], 'b','linewidth',3);
%     axis equal
%     axis(ax_lims)
%     legend('Noninjected','10+/+')
%     title 'Epiblast Composition'
%     ylabel 'Donor Derived Epiblast'
%     xlabel 'Host Derived Epiblast'
%     set(gca, 'FontSize',15)
% 
%     subplot(246)
%     h_gs = scatter(null_injected(:,1),null_injected(:,2),60,[0.5 0.5 0.5],...
%         'Filled');
%     axis equal
%     axis(ax_lims)
%     legend('10-/-')
%     set(gca, 'FontSize',15)
%     title 'Epiblast Composition'
%     ylabel 'Donor Derived Epiblast'
%     xlabel 'Host Derived Epiblast'
% 
%     model_inference.plot_GMM(GMModel_1,null_injected,v,ax_lims,...
%         '1 Component GMM',143);
%     model_inference.plot_GMM(GMModel_2,null_injected,v,ax_lims,...
%         '2 Component GMM',144);
% 
%     % save and close figure
%     saveFileName = ['two_d_clustering'];
%     fullFileName = fullfile(save_folder_name,saveFileName);
%     print(fig, fullFileName,'-dtiffn')
%     print(fig, fullFileName, '-vector', '-dsvg')
%     savefig(fig, fullFileName)
%     close(fig)