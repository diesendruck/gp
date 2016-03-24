function [ resid_sq ] = Resid_Sq_Gauss_Ker_Reg( h,x,y,xt,yt )

ys=gaussian_kern_reg(xt,x,y,h);
e=(yt-ys);
resid_sq = e*e';

end

