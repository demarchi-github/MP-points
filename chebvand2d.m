%--------------------------------------------------------------------------
% chebvand2d.
%--------------------------------------------------------------------------

function V = chebvand2d(deg,mymesh,rect)

%--------------------------------------------------------------------------
% OBJECT:
%-----------
%
% This procedure computes the Chebyshev-Vandermonde matrix at a 2d mesh
% in the product Chebyshev basis of a rectangle containing the mesh points.
%
%--------------------------------------------------------------------------
% INPUTS:
%----------
%
% deg : POLYNOMIAL DEGREE.
%
% mymesh: 2-columns array of mesh points coordinates (x,y).
%           
% rect : 4-components vector such that the rectangle
%       [rect(1),rect(2)] x [rect(3),rect(4)] contains the mesh
%
%----------
% OUTPUTS:
%----------
%
% V : Chebyshev-Vandermonde matrix at a 2d mesh in the product Chebyshev 
%     basis of a rectangle containing the mesh points.
%     It is a "dim(mymesh) x dim(polynomial space)" matrix.
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% Copyright (C) 2007-2010 Alvise Sommariva, Marco Vianello.
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
%%
%% Author:  Alvise Sommariva <alvise@euler.math.unipd.it>
%%          Marco Vianello   <marcov@euler.math.unipd.it>
%%
%% Date: May 17, 2010.
%--------------------------------------------------------------------------

if nargin < 3
    rect=[-1 1 -1 1];
end

% Couples with length less or equal to deg.
j=(0:1:deg);
[j1,j2]=meshgrid(j);
good=find(j1(:)+j2(:)<deg+1);
couples=[j1(good) j2(good)];

% Mapping the mesh in the square [-1,1]^2.
a=rect(1);b=rect(2);c=rect(3);d=rect(4);
map=[(2*mymesh(:,1)-b-a)/(b-a) (2*mymesh(:,2)-d-c)/(d-c)];

% Chebyshev-Vandermonde matrix at the mesh.
V1=cos(couples(:,1)*acos(map(:,1)'));
V2=cos(couples(:,2)*acos(map(:,2)'));
V=(V1.*V2)';

%ee=0.01;
%V=[diag((1+ee-map(:,1)).^(-2))*V];



