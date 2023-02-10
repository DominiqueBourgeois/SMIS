function [x2, y2]=apply_distortion(x,y,im_par)

xr=x-(im_par.n/2-0.5);
yr=y-(im_par.m/2-0.5);
deform=im_par.two_channel_deform;

shift_x=deform(1)*(1+deform(6)/100*randn(1,1));
shift_y=deform(2)*(1+deform(6)/100*randn(1,1));
str_x=deform(3)*(1+deform(6)/100*randn(1,1));
str_y=deform(4)*(1+deform(6)/100*randn(1,1));

theta=deform(5)*pi/180*(1+deform(6)/100*randn(1,1));
a=str_x*cos(theta);
b=str_y*sin(theta);
c=-str_x*sin(theta);
d=str_y*cos(theta);

x2=a*xr+b*yr+shift_x+(im_par.n/2-0.5);
y2=c*xr+d*yr+shift_y+(im_par.m/2-0.5);

end




