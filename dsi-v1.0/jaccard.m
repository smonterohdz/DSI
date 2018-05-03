function [ index ] = jaccard( G,H )
%JACCARD Jaccard index between two networks
%   JACCARD(G,H) computes Jaccard index between G and H adjacency matrices
%   derived from two connectivity networks.
%   Input:
%   'G,H'     - Ajacency networks derived from two connectivty networks
%   represented by undirected graphs. Both G and H must have the same
%   dimensions. Rows and columns represent the number of nodes in
%   corresponding networks. Both matrices encode present and absent links
%   with '1' and '0' respectively. The main diagonal must contains '0'.
%   Output:
%   'index'   - The value of Jaccard index in the range [0,1] where 0 and 1
%   are completely assymetric and symmetric matrices respectively.
%
%   Montero-Hernandez - 2018 March 28
%
index = sum(and(G(:),H(:))) / sum(or(G(:),H(:)));

end

