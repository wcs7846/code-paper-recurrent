function [foldSizes] = computeFoldSizes(vecsPerCat, numFolds)
%COMPUTEFOLDSIZES Compute the fold sizes for the provided data set.
%  [foldSizes] = computeFoldSizes(vecsPerCat, numFolds);
%
%   Part of performing n-fold cross-validation. See 'getFoldVectors'.
%
%  Parameters:
%    vecsPerCat - Column vector with the number of vectors in the data set for
%                 each category.
%    numFolds   - The number of folds you'll be using for cross validation.
%
%  Returns:
%    A matrix (with size numCategories x numFolds) containing the number of 
%    vectors to include in each of the 'numFolds' folds for each category.

% $Author: ChrisMcCormick $    $Date: 2013/07/31 22:00:00 $    $Revision: 1.0 $

% Get the number of categories present.
numCats = numel(vecsPerCat);

% For each category...
for (i = 1 : numCats)
    
    % Get the number of vectors for this category;
    numVecs = vecsPerCat(i, 1);
    
    % Verify that there are at least 'numFolds' samples.
    if (numVecs < numFolds)
        disp("ERROR! Each category must have at least 'numFolds' samples.");
        return;
    end
end
    
% foldSizes will be a matrix holding the number of vectors to place in each fold
% for each category. The number of folds may not divide evenly into the number
% of vectors, so we need to distribute the remainder.
foldSizes = zeros(numCats, numFolds);

% For each category...
for (i = 1 : numCats)
    
    % Get the the number of vectors for this category.
    numVecs = vecsPerCat(i, 1);
            
    % For each of the ten folds...
    for (fold = 1 : numFolds)
        
        % Divide the remaining number of vectors by the remaining number of folds.
        foldSize = ceil(numVecs / (numFolds - fold + 1));
        
        % Store the fold size.
        foldSizes(i, fold) = foldSize;
        
        % Update the number of remaining vectors for this category.
        numVecs = numVecs - foldSize;
    end
end

% Verify the fold sizes sum up correctly.
if (any(sum(foldSizes, 2) ~= vecsPerCat))
    disp("ERROR! The sum of fold sizes did not equal the number of category vectors.");
    return;
end

end
