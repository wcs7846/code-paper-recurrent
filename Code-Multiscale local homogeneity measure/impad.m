% IMPAD - adds zeros to the boundary of an image
%
% Usage:  paddedim = impad(im, b, v)
%
% Arguments:     im - Image to be padded (greyscale or colour)
%                 b - Width of padding boundary to be added
%                 v - This may be:
%                     1) A numeric value to be assigned to the padded area.
%                     2) 'replicate' which results in the edges of the image 
%                        being replicated outwards to form the padding. 
%                     4) 'taper' which results in the edges of the image 
%                        being tapered towards the cyclically opposite edge. 
%                     5) If v is omitted it defaults to 0.
%
% Returns: paddedim - Padded image of size rows+2*b x cols+2*b
%
% Tapering is perhaps the best for frequency domain filtering.
%
% See also: IMTRIM, IMSETBORDER

% Peter Kovesi
% Centre for Exploration Targeting
% The University of Western Australia
% http://peterkovesi.com/matlabfns/
% June    2010
% January 2011  Added optional padding value.
% October 2017  Added replicate and taper options.

function pim = impad(im, b, v)

    if b == 0
        pim = im;
        return;
    end

    if ~exist('v', 'var'), v = 0; end
    
    [rows, cols, channels] = size(im);    
    
    if isnumeric(im) && isnumeric(v);   % numeric padding
        pim = v*ones(rows+2*b, cols+2*b, channels, class(im));
        
    elseif isnumeric(im) && strcmp(v, 'replicate')
        pim = zeros(rows+2*b, cols+2*b, channels, class(im));
    
    elseif isnumeric(im) && strcmp(v, 'taper')
        pim = zeros(rows+2*b, cols+2*b, channels, 'double');        
        
    elseif islogical(im)                % Logical image
        if v == true
            pim = true(rows+2*b, cols+2*b, channels);
        else
            pim = false(rows+2*b, cols+2*b, channels);
        end        
    end
    
    pim(1+b:rows+b, 1+b:cols+b, :) = im;
    
    if strcmp(v, 'replicate') 
        % Replicate image edge values outwards
        for r = 1:b
            pim(r, 1+b:cols+b, :) = im(1,:,:);          % top
            pim(rows+b+r, 1+b:cols+b,:) = im(end,:,:);  % bottom
            
            pim(1+b:rows+b, r, :) = im(:,1,:);          % left
            pim(1+b:rows+b, cols+b+r, :) = im(:,end,:); % right
        end
        
        % Replicate corners
        for ch = 1:channels
            pim(1:b, 1:b, ch) = im(1,1,ch);             % top left
            pim(1:b, cols+b+1:end, ch) = im(1,cols,ch); % top right
            
            pim(rows+b+1:end, 1:b, ch) = im(rows,1,ch); % bottom left
            pim(rows+b+1:end, cols+b+1:end, ch) = im(rows,cols,ch); % bottom right
        end
        
    elseif strcmp(v, 'taper')
        delta = 1/(2*b+1);  % Fractional change that forms the tapered steps 
        bot_top = pim(rows+b,:) - pim(1+b,:);
        for n = 1:b
            pim(n,:) =  pim(1+b,:) + (b-n+1)*delta*bot_top;     % top
            pim(rows+b+n,:) =  pim(rows+b,:) - n*delta*bot_top; % bottom
        end

        right_left = pim(:,cols+b) - pim(:,1+b);        
        for n = 1:b
            pim(:,n) =  pim(:,1+b) + (b-n+1)*delta*right_left;     % left
            pim(:,cols+b+n) =  pim(:,cols+b) - n*delta*right_left; % right
        end
    end