function norm_expression = normalize_expression_loop(selectGenes, probeInformation, parcelExpression)
% normalize_expression    normalizes expression measures for selected genes using z-score normalization
%
%---INPUTS:
% selectGenes, cell containing gene names
% probeInformation, structure containing information about genes (columns
% in parcelExpression variable)
% parcelExpression, matrix containing regional gene expression measures
% (rows - regions, columns - genes)
%
%---OUTPUT:
% normalized gene expression values

% 1. an inefficient version using a loop
% define variables
numGenes = length(selectGenes); 
norm_expression = zeros(size(parcelExpression,1), numGenes); 
Vnorm = zeros(size(parcelExpression,1),1); 

% loop over selected genes
for g=1:numGenes
    
    % find corresponding column in the matrix
    indG = strcmp(probeInformation.GeneSymbol, selectGenes{g});
    % select values from the expression matrix
    expVal = parcelExpression(:, indG); 
    
    % loop over values and z-score them
    for i=1:size(parcelExpression,1)
    Vnorm(i) = (expVal(i) - nanmean(expVal))./nanstd(expVal); 
    end
    
    % assing normalised values to the matrix
    norm_expression(:,g) = Vnorm; 
    
end

end
