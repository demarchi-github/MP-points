% ----------------------------------------------------------------------
% Compute the Lebesgue constant for Morrow-Patterson points evaluating it
% on a grid of equally spaced or Chebyshev-Lobatto. 
% The code produces some plots, in particular,r the growth of the constant
% on varying the polynomial degree from 2 to 30 (on the evens), comparing
% the results with those of the best-known growth, (0.7n+1)^2.
%
%-----------------------------------------------------------------------
%
% Author: Stefano De Marchi, University of Padova
% Last update: May 2024
%
%-----------------------------------------------------------------------
clear all
close all

n_vett=2:2:30;    % vector degree

ntg1mx = 100;      % the maximum value of the parameter determining 
                  % the number of target points.
ntgmax = ntg1mx ^ 2;  % the maximum number of target points
ntg1=100;              % parameter determining the number of target points 

xyrange=[-1,1,-1,1];

% Compute the target points (equally spaced)
%[ tg1, tg2, ntg ] = target2 ( -1,1,-1,1, ntg1, ntgmax );

% Compute the target points (Chebyshev-Lobatto grid, whih is a Weakly admissible mesh)
 cheb_grid=-cos([0:ntg1]*pi/ntg1);
 [ tg1, tg2 ] = meshgrid(cheb_grid,cheb_grid)

%tg = [tg1',tg2'];
tg = [tg1(:),tg2(:)];

% Compute the Lebesgue constant
for index = 1:length(n_vett)
    deg = n_vett(index);
    V_leb = chebvand2d(deg,tg);
      
    
    [xmp, ymp,wmp] =  morr_patt (deg);
    mp = [xmp',ymp'];
    V_mp = chebvand2d(deg,mp);
   
    leb_mp(index) = norm(V_mp'\V_leb',1);
    leb_mp_true(index)=(n_vett(index)+1)*(n_vett(index)+2)/2;
    leb_approx(index)=(0.7*n_vett(index)+1)^2;
end

t=linspace(0,pi,1000);
plot(-cos((n_vett(index)+3)*t), -cos((n_vett(index)+2)*t),'r',xmp',ymp','o')
figure

lambda=sum(abs(V_leb'));
disp('Maximum of the Lebesgue function') 
Max_Leb=max(lambda)
Mt=length(cheb_grid);
%Mt=ntg1;  %in case of equally spaced points
funLeb=reshape(lambda,Mt,Mt);
xtg=reshape(tg1,Mt,Mt);
ytg=reshape(tg2,Mt,Mt);

mesh(xtg,ytg,funLeb)
title(['Lebesgue function for n='],num2str(n_vett(end)))
figure


semilogy(n_vett, leb_mp,'g-d','LineWidth',2)
hold on
%semilogy(n_vett, (0.7*n_vett+1).^2,'b-s','LineWidth',2);
semilogy(n_vett, leb_mp_true,'r-s','LineWidth',2);
hold off
xlabel('degree');
ylabel('Lebesgue constant');
legend('MP', '(n+2)*(n+1)/2' ) ;
%legend('MP', '(0.7*n+1).^2') ;

figure 
for i=1:length(leb_mp)
  err(i)=abs((leb_mp(i)-leb_mp_true(i))/leb_mp(i));
  err_approx(i)=abs((leb_mp(i)-leb_approx(i))/leb_mp(i));
end
semilogy(1:length(err),err,'b-o','LineWidth',2)
hold on
semilogy(1:length(err),err_approx,'r-s','LineWidth',2)
title('relative errors')
