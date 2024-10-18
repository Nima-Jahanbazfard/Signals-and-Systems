% nimajahanbazfard 
clear all
clc
%define variables 
img = sqrt(-1);
width = [-0.25 0.25] ;
length = ( width(2) - width(1) )/ 2 ;
nop = 100 ; 
f = @(x) 1 ;
sum1 = 0 ;
sum2 = 0 ;
x = -1 : 0.01 : 1 ;
num = numel(x) ;
w0 = 2*pi;
% calculating integrals
a = -length ;
b = length ;
k = 100 ;
step = (b - a)/k ;
zoa0 =  zeros(1 , nop) ;
zoan =  zeros(1 , nop);
zobn =  zeros(1 , nop);

midpa0 = zeros(1 , nop) ;
midpan = zeros(1 , nop) ;
midpbn = zeros(1 , nop) ;
% calculating constants of fourie 
for j = 1 : nop
    
        f1 = @ (x) cos(j*x*(pi/length));
        f2 = @ (x) sin(j*x*(pi/length));
    
    for i = a : step : b
        % zozanaghe 
     zoa0(j) = zoa0(j) + (step/2)*( f(i) + f(i+1)  );
     zoan(j) = zoan(j) + (step/2)*( f1(i) + f1(i+1)  );
     zobn(j) = zobn(j) + (step/2)*( f2(i) + f2(i+1)  ) ;
     %noghte miyani 
     midpa0(j) = midpa0(j) +  step * ( f( (2*i+1)/2 ) );
     midpan(j) = midpan(j) +  step * ( f1( (2*i+1)/2 ) );
     midpbn(j) = midpbn(j) +  step * ( f2( (2*i+1)/2 ) );
    end 
end 
% ravesh simpson
for z= 1 : nop
    a0(z) = (  (zoa0(z) + 2*midpa0(z)) / 3   ) * (1/(2*length)) ;
    an(z) = (  (zoan(z) + 2*midpan(z)) / 3   ) * (1/length)     ;
    bn(z) = (  (zobn(z) + 2*midpbn(z)) / 3   ) * (1/length)     ;
    
end
 xt = 0;
   
    for x0 = 1 : num
        for n = 1 : nop
            
           sum1 = sum1 + an(n)*cos( n*x(x0)*(pi/length) );
           sum2 = sum2 + bn(n)*sin( n*x(x0)*(pi/length) );
           xt = a0(1) + sum1 + sum2;
      
        end         
    end
     ak = (an-img*bn)/2;
    
    
    xr = 0;
    
    r = [1 , 2 , 3 , 4 , 5 , 10 , 20 , 100];
    
    
     for x0 = 1 : num
        for n =1:8
          r0 = r(n);
          
             xr = xr + ak(r0)* exp(img*r0*w0*x0);
            
        end
     end
   %part a real series (an ,bn0
   an
   bn
   %part b 
   ak
   %part c 
   xr
   xt
   %part d
   d=xr-xt
                                  
    