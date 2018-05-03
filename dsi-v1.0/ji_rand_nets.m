function [tbresults] = ji_rand_nets(n_nodes,rep)
%JI_RAND_NETS Jaccard index of random networks
%   JI_RAND_NETS(n_nodes,rep) computes Jaccard index and measures of
%   interest from a simulation of random connectivity networks.
%   Input: 
%   'n_nodes'     - Desired number of nodes in simulated HbO and HHb networks
%   'rep'         - Desired number of repetitons of simulation
%   Output: 
%   'tbresults'   - Table with the results of simulation with columns:
%        .denshbo   Connection density in HbO net
%        .denshbr   Connection density in HbR net
%        .mean      Mean value of Jaccard index across repetitions
%        .std       Standard deviation of Jaccard index
%        .max       Max value Jaccard index across repetitions
%        .min       Min value Jaccard index across repetitions
%        .edgeshbo  Number of edges in HbO net
%        .edgeshbr  Number of edges in HbR net
%
%   Montero-Hernandez - 2018 May  
%
if n_nodes < 2
    error('n_nodes must be > 2');
end
if rep < 5 || rep > 21
    error('rep must be > 5 and <= 20');
end
% Start with a complete graph
am_rnd = ones(n_nodes);
% Delete self edges (main diagonal)
am_rnd(1:(n_nodes+1):(n_nodes*n_nodes))=0;
% Compute number of total edges (https://en.wikipedia.org/wiki/Complete_graph)
tot_edges = (n_nodes*(n_nodes-1))/2;
% Variable storing the results
dfresults = zeros((tot_edges-1).^2,8); %[nedges-1 rows X (mean,std,sparsity,maxji,minji)]
dfi = 1;
for nrndhbo = 1:(tot_edges-1)
    for nrndhbr = 1:(tot_edges-1)
        maxji = 0;
        minji = 2;
        dfji = zeros(rep,2);
        for ii = 1 :rep%((n*n)-n)
            hbo = am_rnd;
            hhb = am_rnd;
            [hbo] = del_edge(hbo,nrndhbo);
            [hhb] = del_edge(hhb,nrndhbr);
            ji = jaccard(hbo,hhb);
            dfji(ii,1) = ji;
            if ji> maxji
                maxji = ji;
            end
            if ji<minji
                minji = ji;
            end
        end
        
        dfresults(dfi,1) = (sum(sum(hbo))/2)/tot_edges;
        dfresults(dfi,2) = (sum(sum(hhb))/2)/tot_edges;
        dfresults(dfi,3) = mean(dfji(:,1));
        dfresults(dfi,4) = std(dfji(:,1));
        dfresults(dfi,5) = maxji;
        dfresults(dfi,6) = minji;
        dfresults(dfi,7) = tot_edges - nrndhbo;
        dfresults(dfi,8) = tot_edges - nrndhbr;
        if mod(dfi,100) == 0
            fprintf('(%i of %i) \t %.2f \t %.4f \n',dfi,(tot_edges-1).^2,dfresults(dfi,1),dfresults(dfi,2));
        end
        dfi = dfi +1;
    end
end
tbresults = array2table(dfresults,...
    'VariableNames',{'denshbo','denshbr','mean',...
    'std','max','min','edgeshbo','edgeshbr'});
tbresults.Properties.Description = 'Measures of Jaccard index in simulated networks.';
tbresults.Properties.VariableDescriptions = {
    'HbO connection density',...
    'HbR connection density',...
    'mean Jaccard value across repetitions',...
    'std dev Jaccard value',...
    'Max Jaccard value across repetitions',...
    'Min Jaccard value across repetitions',...
    'Edges in HbO net',...
    'Edges in HbR net'};
end

function [mat] = del_edge (mat,nrnd)
n = size(mat,1);
for ii = 1: nrnd
    f = true;
    while f == true
        rndlink = randi([1,(n*n)],1);
        if mat(rndlink)==1
            f = false;
        end
    end
    [i,j] = ind2sub([n,n],rndlink);
    mat(i,j) = 0;
    mat(j,i) = 0;
end
end
