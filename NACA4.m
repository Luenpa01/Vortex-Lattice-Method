function [zc, ze, zi] = NACA4(NACA, x)



f = str2double(NACA(1))/100; 

if f == 0 
    xf = 0.5;
else
    xf = str2double(NACA(2))/10;
end
t = str2double(NACA(3:4))/100; 



if x < xf
    zc = f/xf^2*(2*xf*x-x^2);
else
    zc = f/(1-xf)^2*((1-2*xf)+2*xf*x-x^2);
end



a = 0.2969;
b = -0.1260;
c = -0.3516;
d = 0.2843;
e = -0.1015;

zt = 5*t*(a*sqrt(x)+b*x+c*x^2+d*x^3+e*x^4);



ze = zc+zt;



zi = zc-zt;

end

