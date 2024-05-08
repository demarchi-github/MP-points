function [ tg1, tg2, ntg ] = target2 ( a1, b1, a2, b2, ntg1, ntgmax )

%*****************************************************************************80
%
%% TARGET returns the target points on the rectangle.
%
%  Discussion:
%
%    Target points (uniform grid) on the rectangle [A1,B1] x [A2,B2].
%    The number of target points is NTG = NTG1^2.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%  
%  Modified:
%
%    13 February 2014
%
%  Author:
%
%    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
%    Marco Vianello.
%    This MATLAB version by John Burkardt.
%
%  Reference:
%
%    Marco Caliari, Stefano De Marchi, Marco Vianello,
%    Algorithm 886:
%    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
%    ACM Transactions on Mathematical Software,
%    Volume 35, Number 3, October 2008, Article 21, 11 pages.
%
%  Parameters:
%
%    Input, real A1, B1, A2, B2, the coordinates of the extreme
%    corners of the rectangle.
%
%    Input, integer NTG1, a parameter determining the number 
%    of target points
%
%    Input, integer NTGMAX, the maximum number of target points.
%
%    Output, real TG1(NTG), TG2(NTG), the X and Y coordinates
%    of the target points.
%
%    Output, integer NTG, the number of target points computed.
%
 
  
  ntg = 0;
  for i = 1 : ntg1
    for j = 1 : ntg1
      ntg = ntg + 1;
      tg1(ntg) = a1 + ( j - 1 ) * ( b1 - a1 ) / ( ntg1 - 1 );
      tg2(ntg) = a2 + ( i - 1 ) * ( b2 - a2 ) / ( ntg1 - 1 );
    end
  end

  return
end