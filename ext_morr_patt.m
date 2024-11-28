% Compute the Morrow-Patterson points

function [xmp, ymp, wmp] =  morr_patt (n)

 index = 1;
 
 alfa_n=cos(pi/(n+2)); beta_n=cos(pi/(n+3));
     for m = 1 : 1 : n + 1
         
        for k = 1 : 1 : n/2 +1
         
            xmp(index) = cos (m * pi / (n+2)) / alfa_n;
           if mod(m,2)==0
             ymp(index) = cos ((2*k-1)*pi/(n+3))/ beta_n;
           else
             ymp(index) = cos (2*k*pi/(n+3))/ beta_n;
           end
           wmp(index) = 8.0 /((n+2)*(n+3));
          % wmp(index) = 2.0 /((n+2)*(n+3));  %vecchi pesi
           index = index +1;
           
      end
 end
end