%DSI_PIPELIN Script for Differential Symmetry Index (DSI) computation
%   DSI_PIPELIN generates a set of synthetic connectivity networks, obatins
%   the Jaccard index and approximates a bilinear model in order to obtain
%   the DSI values from a set of fNIRS connectivity networks.
% ------------------------  Setup variables  --------------------------
%    'pathbct'            - The path for brain connectivity toolbox
%    'comp_cbm'           - Compute a set of brain complex measures
%                           (1=Yes/0=No)
%    'nodes'              - Number of nodes in connectivity networks
%    'repetitions'        - Repetitions for simulation
%    'filename_info'      - Structure of information for files access:  
%      .path                - Path of the adjacency matrices files
%      .file_prefix         - Prefix of input ajacency matrices files
%      .subject_prefix      - Prefix of the subjects
%      .conditions_prefix   - Prefix of blocks or conditions
%      .header              - Header on the file? (1=Yes/0=No)
%    'n_subj'                - Number of subjects to process
% ------------------------  Output variables  -----------------------
%    'tb_rnd_jaccard'        - Jaccard symmetry of simulation of random networks
%     .denshbo   Connection density in HbO net
%     .denshbr   Connection density in HbR net
%     .mean      Mean value of Jaccard index across repetitions
%     .std       Standard deviation of Jaccard index
%     .max       Max value Jaccard index across repetitions
%     .min       Min value Jaccard index across repetitions
%     .edgeshbo  Number of edges in HbO net
%     .edgeshbr  Number of edges in HbR net
%
%    'tbresults_sym'  - Table with DSI and Jaccard symmetry results:
%      .subj          - Subject ID
%      .nodes         - Nodes
%      .condition     - Condition
%      .denshbo       - HbO connection density
%      .denshhb       - HHb connection density
%      .jaccardInd    - Jaccard index value
%      .dsi           - DSI value
%
%    'tbresults_cbm'  - Table with complex brain measures results
%      .subj          - Subject ID
%      .nodes         - Nodes
%      .condition     - Condition
%      .hbx           - Hb signal (HbO/HbR)
%      .denshbx       - 'HbX connection density,...
%      .tra           - Transitivity
%      .cpl           - Characteristic path length
%      .gef           - Global efficiency
%      .mod           - Modularity
%

clearvars;
pathbct='C:\BCT';
addpath(pathbct);
% Compute complex brain measures (1=Yes/0=No)
comp_cbm = 1;
    
nodes = 16;
repetitions = 20;
tb_rnd_jaccard = ji_rand_nets(nodes,repetitions);
save('tb_rnd_jaccard.mat','tb_rnd_jaccard');

%'am-SUBJ_ID-CONDITION-HBX.csv'
filename_info.path = 'c:/dsi-data-example/';
filename_info.file_prefix = 'am';
filename_info.subject_prefix = 'S';
filename_info.conditions_prefix = {'BL','OG','SocPM','NonSocPM','OGc'};
filename_info.sys_prefix = {'Sys','NoSys'};
filename_info.header = 1; % Header on the file? (1=Yes/0=No)
n_subj = 7;
[tbresults_sym,tbresults_cbm] = dsi_nirs_nets(nodes,filename_info,n_subj,tb_rnd_jaccard,comp_cbm);

% Plot of results
group_by = tbresults_sym.condition;
figure(1);
clf(1);
colorsgroups = colormap(hsv(length(unique(group_by))));

sp1=scatterhist(1:2:(2*length(tbresults_sym.dsi)),tbresults_sym.dsi, ...
    'NBins',11,...
    'Marker','.',...
    'MarkerSize',15,...
    'Group',group_by,...
    'Color',colorsgroups,...
    'PlotGroup','off');
xlabel('fNIRS Networks'); ylabel('Differential Symmetry Index'); 
l1=line([-10 155],[0.2 0.2],'Color',[.38 .38 .38],'LineWidth',2,'LineStyle','--');

