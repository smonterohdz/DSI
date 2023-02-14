function [tbresults_sym,tbresults_cbm] = dsi_nirs_nets(n_nodes,filename_info,n_subj,tb_rnd_jaccard, comp_cbm)
%JI_NIRS_NETS Jaccard index of fnirs networks
%   JI_NIRS_NETS(n_nodes, path, file_prefix, n_nets) computes Jaccard index 
%   and measures of interest from a simulation of random connectivity networks.
%   Input: 
%    'n_nodes'               - Desired number of nodes in simulated HbO and HHb networks
%    'filename_info'           - A structure containing:
%      .path               - Path of the adjacency matrices files
%      .file_prefix        - Prefix of input ajacency matrices files
%      .subject_prefix     - Prefix of the subjects
%      .conditions_prefix  - Prefix of blocks or conditions
%      .sys_prefix         - Prefix of existence of systemic variables
%    'n_subj'                - Number of subjects to process
%    'tb_rnd_jaccard'        - Jaccard symmetry of simulation of random networks
%   Output: 
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
%   Montero-Hernandez - 2018 May 
%
% Computing symmetry, complex brain measures 
%



%%%%%%%%%%%%%%%%%%%%%%%%
% Bilinear model
Xrnd = [tb_rnd_jaccard.denshbo, tb_rnd_jaccard.denshbr];
SymJacRnd = tb_rnd_jaccard.mean;
% Calculating the coefficients of the Bilinear Model [C][a]=[Z]
% where [C] = [[1], X, Y, X*Y], [Z] = f([Q]), and [Q] = [X,Y]
N = size(Xrnd,1);
C = [ones(N,1) Xrnd(:,1) Xrnd(:,2) (Xrnd(:,1).*Xrnd(:,2))];
a = C\SymJacRnd;
% Bilinear interpolation model
fbil = @(X,Y,a) a(1) + a(2)*X + a(3)*Y + a(4)*(X.*Y);

% Evaluating random Nets with Bilinear model
%fJac  = fbil(Xrnd(:,1),Xrnd(:,2),a);
%fnirs = fbil(Xnirs(:,1),Xnirs(:,2),a);

%%%%%%%%%%%%%%%%%%%%%%%
% Filename format 'am-SUBJ_ID-CONDITION-HBX.csv'
path = filename_info.path;
file = filename_info.file_prefix;
subject = filename_info.subject_prefix;
conditions = filename_info.conditions_prefix;
%sys_prefix = filename_info.sys_prefix;
nconditions = length(conditions);
header = filename_info.header;


% dataframes for symmetry
nsamp_symm = n_subj*nconditions;
dfresults_sym = cell(nsamp_symm,7);

% dataframes for individual HbX Complex Brain Measures
if comp_cbm == 1
    nsamp_cbm = nsamp_symm*2;
    dfresults_cbm = cell(nsamp_cbm,9);
end

dfisymm = 1;
dficbm = 1;
modrep=10;
nedges = (n_nodes*(n_nodes-1))/2;

%for sys = 1:length(sys_prefix)
    for s=1:n_subj
        fprintf('S%d\n',s);
        for i=1:nconditions
            
            % Filename format 'am-SID-CONDITION-HBX.csv'
            if exist([path,file,'-',sprintf('%s%d-',subject,s),conditions{i},'-hbo.csv'], 'file') == 2
                amhbo = csvread([path,file,'-',sprintf('%s%d-',subject,s),conditions{i},'-hbo.csv'],header,0);
            else
                disp(['Not found ',path,file,'-',sprintf('%s%d-',subject,s),conditions{i},'-hbo.csv']);
                continue;
            end
            if exist([path,file,'-',sprintf('%s%d-',subject,s),conditions{i},'-hbr.csv'], 'file') == 2
                amhbr = csvread([path,file,'-',sprintf('%s%d-',subject,s),conditions{i},'-hbr.csv'],header,0);
            else
                continue;
            end
            
            amhbo = double(amhbo~=0); amhbr = double(amhbr~=0);
            
            %-- SYMM
            symJaccard = sum(and(amhbo(:),amhbr(:)))/sum(or(amhbo(:),amhbr(:)));
            den_hbo = (sum(sum(amhbo))/2)/nedges;
            den_hbr = (sum(sum(amhbr))/2)/nedges;
            
            dfresults_sym(dfisymm,1) = {sprintf('%s%d',subject,s)};
            dfresults_sym(dfisymm,2) = {n_nodes};
            dfresults_sym(dfisymm,3) = {conditions{i}};
            dfresults_sym(dfisymm,4) = {den_hbo};
            dfresults_sym(dfisymm,5) = {den_hbr};
            dfresults_sym(dfisymm,6) = {symJaccard};
            symBlm = fbil(den_hbo,den_hbr,a);
            dfresults_sym(dfisymm,7) = {abs(symJaccard - symBlm)};
            %dfresults_sym(dfisymm,8) = {sys_prefix{sys}};
            dfisymm = dfisymm +1;
            if comp_cbm == 1
                %-- CBN for HBO
                dfresults_cbm(dficbm,1) = {sprintf('%s%d',subject,s)};
                dfresults_cbm(dficbm,2) = {n_nodes};
                dfresults_cbm(dficbm,3) = {conditions{i}};
                dfresults_cbm(dficbm,4) = {'HbO'};
                dfresults_cbm(dficbm,5) = {den_hbo};
                dfresults_cbm(dficbm,6) = {transitivity_bu(amhbo)};  %TRA
                Do = distance_bin(amhbo);
                [cplo,gefficiencyo] = charpath(Do);
                dfresults_cbm(dficbm,7) = {cplo};                      %CPL
                dfresults_cbm(dficbm,8) = {gefficiencyo};              %GEf
                mod=0;
                for m=1:modrep
                    [~, modt] = modularity_und(amhbo);
                    mod = mod + modt;
                end
                mod = mod/modrep;
                dfresults_cbm(dficbm,9) = {mod};                       %MOD
                %dfresults_cbm(dficbm,10) = {sys_prefix{sys}};
                dficbm = dficbm +1;
                
                %-- CBN for HBR
                dfresults_cbm(dficbm,1) = {sprintf('%s%d',subject,s)};
                dfresults_cbm(dficbm,2) = {n_nodes};
                dfresults_cbm(dficbm,3) = {conditions{i}};
                dfresults_cbm(dficbm,4) = {'HbR'};
                dfresults_cbm(dficbm,5) = {den_hbr};
                dfresults_cbm(dficbm,6) = {transitivity_bu(amhbr)};  %TRA
                Dr = distance_bin(amhbr);
                [cplr,gefficiencyr] = charpath(Dr);
                dfresults_cbm(dficbm,7) = {cplr};                      %CPL
                dfresults_cbm(dficbm,8) = {gefficiencyr};              %GEf
                mod=0;
                for m=1:modrep
                    [~, modt] = modularity_und(amhbr);
                    mod = mod + modt;
                end
                mod = mod/modrep;
                dfresults_cbm(dficbm,9) = {mod};                       %MOD
                %dfresults_cbm(dficbm,10) = {sys_prefix{sys}};
                dficbm = dficbm +1;
            end
        end
    end
%end

tbresults_sym = cell2table(dfresults_sym,...
    'VariableNames',{'subj','nodes','condition','denshbo',...
    'denshhb','jaccardInd','dsi'});
tbresults_sym.Properties.Description = 'DSI and Jaccard index values in fnirs networks.';
tbresults_sym.Properties.VariableDescriptions = {
    'Subject ID',...
    'Nodes',...
    'Condition',...
    'HbO connection density',...
    'HHb connection density',...
    'Jaccard index value',...
    'DSI value'};
idx_empty = cellfun('isempty', tbresults_sym.subj);
tbresults_sym(idx_empty,:) = [];
tbresults_sym.nodes = cell2mat(tbresults_sym.nodes);
tbresults_sym.denshbo = cell2mat(tbresults_sym.denshbo);
tbresults_sym.denshhb = cell2mat(tbresults_sym.denshhb);
tbresults_sym.jaccardInd = cell2mat(tbresults_sym.jaccardInd);
tbresults_sym.dsi = cell2mat(tbresults_sym.dsi);

tbresults_cbm = nan;
if comp_cbm == 1
    tbresults_cbm = cell2table(dfresults_cbm,...
        'VariableNames',{'subj','nodes','condition','hbx','denshbx',...
        'tra','cpl','gef','mod'});
    tbresults_cbm.Properties.Description = 'Complex brain measures for fnirs networks.';
    tbresults_cbm.Properties.VariableDescriptions = {
        'Subject ID',...
        'Nodes',...
        'Condition',...
        'Hb signal (HbO/HbR)',...
        'HbX connection density',...
        'Transitivity',...
        'Characteristic path length',...
        'Global efficiency',...
        'Modularity'};
end
idx_empty = cellfun('isempty', tbresults_cbm.subj);
tbresults_cbm(idx_empty,:) = [];
tbresults_cbm.denshbx = cell2mat(tbresults_cbm.denshbx);
tbresults_cbm.tra = cell2mat(tbresults_cbm.tra);
tbresults_cbm.cpl = cell2mat(tbresults_cbm.cpl);
tbresults_cbm.gef = cell2mat(tbresults_cbm.gef); 
tbresults_cbm.mod = cell2mat(tbresults_cbm.mod);
disp('Done!');
end
