%% PAPER TITLE
clear
% Must be in '.code\' directory

addpath('data','funs');

% Construct Data obj
obj = model_inference;

%% Generates pannels in Figure S2
% Figures can be found in  '.code\results\chimera_exclusion_clustering'
 obj.analysis_chimera_exclusion_clustering
 
%% Generates pannels in Figure 2
% Figures can be found in  '.code\results\chimera_exclusion'
obj.analysis_chimera_exclusion

%% Generates pannels in Figure 3
% Figures can be found in  '.code\results\chimera_crowding'
obj.analysis_chimera_crowding

%% Generates pannels in Figure 1 and S1 
% Figures can be found in  '.code\results\abc_mcmc_emrbyos'
 obj.embryo_inference
 
obj.plotEmbryoParamerters
obj.plotEmbryoSimulation

%% Generates pannels in Figures 4 and S4 
% Figures can be found in  '.code\results\abc_mcmc_chimeras'
obj.chimera_inference

obj.plotChimeraParameters
obj.plotChimeraSimulation

%% Generates pannels in Figure 4
% Figures can be found in  '.code\results\abc_mcmc_chimeras'
  
obj.validateHostEpiblastExclusion
           
%% Generates pannels in Figure S3
% Figures can be found in  '.code\results\chimera_crowding' 
three_tissue_classification
