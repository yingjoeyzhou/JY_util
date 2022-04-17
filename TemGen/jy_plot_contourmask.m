function h = jy_plot_contourmask( x, y, mask )
% This function adds a contour mask to the image (e.g., TFR).
% 
% INPUT:
%   x   : a 1-by-N vector for variable on the x-axis.
%   y   : a 1-by-M vector for variable on the y-axis.
%   mask: this is an N-by-M matrix with zeros and ones.
%
% OUTPUT:
%   h: handle of the contour plot.
% 
% ref: fieldtrip sub-function ft_plot_matrix. 
% 
% JY (March 2020)
% 



if size(mask,1)~=numel(x) | size(mask,2)~=numel(y)
    error('Check size of the inputs!')
end

% [X,Y] = meshgrid(x, y);
[X,Y] = meshgrid(y, x); %Y and X are of size( numel(y), numel(x) ).

X = interp2(X, 2); % change to 4 for round corners
Y = interp2(Y, 2); % change to 4 for round corners

contourlines = mask==1;
contourlines = interp2(contourlines, 2, 'nearest');  % change to 4 and remove 'nearest' for round corners

dx = mean(diff(X(1, :))); % remove for round corners
dy = mean(diff(Y(:, 1))); % remove for round corners

hold on;
[c h] = contour(X+dx/2,Y+dy/2,contourlines,1);



end