function [Pmono, direction] = monotone_1d(PP, varargin)
% Project vector to monotone space.
%
% Args:
%   PP: Vector of data values.
%   direction: Optional String 'up' or 'down'.
%
% Returns:
%   Pmono: Monotone Convex response variable.
%   direction: String, the direction chosen by the function.

% Check if user supplied known direction of monotonicity.
if nargin > 1
    direction = varargin{1};
else
    direction = 'unknown';
end

% Do projection.
if strcmp(direction, 'up')
    % Project to monotone ascending.
    Pmono = monotoneUp_1d(PP);
    
elseif strcmp(direction, 'down')
    % Project to monotone descending.
    Pmono = monotoneDown_1d(PP);
    
elseif strcmp(direction, 'unknown')
    % Try ascending, and if result is flat, try descending.
    Pmono = monotoneUp_1d(PP); 
    direction = 'up';
    if all(Pmono == Pmono(1))
        Pmono = monotoneDown_1d(PP); 
        direction = 'down';
    end
else
    error('Direction argument not recognized.')
end


end
    
    
    
 

   
