function [X_sorted, y_sorted] = randSortAndGroup(X, y, categories)
%RANDSORTANDGROUP Randomly sort the vectors in X, then group them by category.
% [X_sorted, y_sorted] = randSortAndGroup(X, y, categories);
%
% Parameters:
%  X - The matrix of input vectors, one per row.
%  y - The category or class label of the associated input vector.
%  categories - A column vector listing the category values in use by this data set.
%
% Returns:
%   The inputs X and y but with the vectors grouped by category and in a random
%   order within their category.

% $Author: ChrisMcCormick $    $Date: 2013/07/31 22:00:00 $    $Revision: 1.0 $

% ======================================
%       Randomly Sort The Vectors
% ======================================

% Get the total number of input vectors.
[totalVecs,~] = size(X);

% Get a random order of the indeces.
randOrder = randperm(totalVecs)';

% Sort the vectors and categories with the random order.
randVecs = X(randOrder, :);
randCats = y(randOrder, :);

X_sorted = [];
y_sorted = [];

% =======================================
%      Group The Vectors By Category
% =======================================

% Re-group the vectors according to category.
Num_categories = numel(categories);
for (i = 1 : Num_categories)
    
    % Get the next category value.
    cat = categories(i);
    
    % Select all of the vectors for this category.
    catVecs = randVecs((randCats == cat), :);
    catCats = randCats((randCats == cat), :);
   
    % Append the vectors for this category.
    X_sorted = [X_sorted; catVecs];
    y_sorted = [y_sorted; catCats];
end

end
