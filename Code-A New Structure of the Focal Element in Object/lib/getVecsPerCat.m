function [vecsPerCat] = getVecsPerCat(X, y, categories)
%GETVECSPERCAT Get the number of vectors in X belonging to each category.
%  vecsPerCat = getVecsPerCat(X, y, categories);
%
%  The vectors may be supplied in any order (they don't need to be grouped by
%  category) and the category values do not have to be contiguous (you could 
%  have categories 1, 2 and 5, for example).
%
%  Parameters:
%    X          - The matrix of input vectors, one per row.
%    y          - The category or class label of the associated input vector.
%    categories - A column vector listing the category values in use by this 
%               data set.
%
%  Returns:
%    A column vector containing the number of vectors in X belonging to each 
%    category.

% $Author: ChrisMcCormick $    $Date: 2013/07/31 22:00:00 $    $Revision: 1.0 $

% Get the number of categories present.
numCats = length(categories);

% 'vecsPerCat' will store the number of input vectors belonging to each category.
vecsPerCat = zeros(numCats, 1);

% For each category...
for (i = 1 : numCats)
    
    % Get the ith category; store the category value in column 1.
    cat = categories(i);
    
    % Count the number of input vectors with that category.
    vecsPerCat(i, 1) = sum(y == cat);    
end

end