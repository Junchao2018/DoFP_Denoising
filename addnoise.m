function [xin,u]=addnoise(xorig,sig,state)
[m,n] = size(xorig);
randn('state',state)
u = sig*randn(m,n);
xin = xorig + u;
end

