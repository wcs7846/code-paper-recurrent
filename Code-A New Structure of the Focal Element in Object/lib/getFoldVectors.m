function [X_train, y_train, X_val, y_val] = getFoldVectors(X_sorted, y_sorted, ...
                                                categories, vecsPerCat, ...
                                                foldSizes, roundNumber)
%GETFOLDVECTORS Selects the vectors to use for training and validation for the
% specified round number.
%  [X_train, y_train, X_val, y_val] = getFoldVectors(X_sorted, y_sorted, ...
%                                                categories, vecsPerCat, ...
%                                                foldSizes, roundNumber)
%
%  With n-fold cross-validation, the data set is divided into n folds, and then
%  training and validation are performed in 'n' separate rounds. In each round,
%  1 fold is used for validation while the remaining n - 1 folds are used for 
%  training.
%
%  Parameters:
%    X_sorted    - The matrix of input vectors, one per row, grouped by 
%                  category.
%    y_sorted    - The category or class label of the associated input vector.
%    categories  - A column vector listing the category values in use by this
%                  data set.
%    vecsPerCat  - Column vector listing the number of vectors in the dataset 
%                  for each category.
%    foldSizes   - The precomputed fold sizes for each category.
%    roundNumber - The current round of validation (a value between 1 and the 
%                  number of folds).
%
%  Returns:
%    Vectors to use for training in this round (with their associated 
%    categories) and vectors to use for validation in this round (also with 
%    their associated categories).   

% $Author: ChrisMcCormick $    $Date: 2013/07/31 22:00:00 $    $Revision: 1.0 $

X_train = [];
y_train = [];
X_val = [];
y_val = [];

% ================================
%         Verify Sorting
% ================================

% Verify the vectors are properly sorted.
catStart = 1;

% For each category...
for (i = 1 : numel(categories))
    
    % Compute the index of the last vector of this category.
    catEnd = catStart + vecsPerCat(i) - 1;

    % Verify that all of the vectors in the range catStart : catEnd have 
    % the expected category.
    if (any(y_sorted(catStart : catEnd) ~= categories(i)))
        disp("Input vectors are not properly sorted!");
        return;
    end
    
    % Set the starting index of the next category.
    catStart = catEnd + 1;
end

% ==================================
%          Select Vectors
% ==================================

% Get the number of folds from the foldSizes matrix.
[~,numFolds] = size(foldSizes);

catStart = 1;

% For each category...
for (catIndex = 1 : numel(categories))

    % Get the list of fold sizes for this category as a column vector.
    catFoldSizes = foldSizes(catIndex, :)';
    
    % Set the starting index of the first fold for this category.
    foldStart = catStart;
    
    % For each fold...
    for (foldIndex = 1 : numFolds)
        
        % Compute the index of the last vector in this fold.
        foldEnd = foldStart + catFoldSizes(foldIndex) - 1;
        
        % Select all of the vectors in this fold.
        foldVectors = X_sorted(foldStart : foldEnd, :);
        foldCats = y_sorted(foldStart : foldEnd, :);
        
        % If this fold is to be used for validation in this round...
        if (foldIndex == roundNumber)
            % Append the vectors to the validation set.
            X_val = [X_val; foldVectors];
            y_val = [y_val; foldCats];
        % Otherwise, use the fold for training.
        else
            % Append the vectors to the training set.
            X_train = [X_train; foldVectors];
            y_train = [y_train; foldCats];
        end
        
        % Update the starting index of the next fold.
        foldStart = foldEnd + 1;
    end
    
    % Set the starting index of the next category.
    catStart = catStart + vecsPerCat(catIndex);   
end

end
