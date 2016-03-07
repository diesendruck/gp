function r = invgamrand(a, b, varargin)
%INVGAMRAND Random matrices from inverse gamma distribution
%
%   R = INVGAMRAND(A,B) returns a matrix of random numbers chosen   
%   from the inverse gamma distribution with parameters A and B.
%   The size of R is the common size of A and B if both are matrices.
%   If either parameter is a scalar, the size of R is the size of the other
%   parameter. Alternatively, R = GAMRND(A,B,M,N) returns an M by N matrix. 
% 
%   Note: Parameterization as in (Neal, 1996).
%      A is mean of the distribution
%      B is degrees of freedom
%   
%	See also GAMRAND
%
% Copyright (c) 1999 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

error('No mex-file for this archtitecture. See Matlab help and convert.m in ./linuxCsource or ./winCsource for help.')
