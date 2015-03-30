close all
clear all
clc

%Global domain initialization (parameters)
B= [-1 1 -1 1];                     %left, right, bottom, top
N= 4;                               %Local poly order
K= [20 10];                         %Num elements along x,y
Ex= linspace(B(1),B(2),K(1)+1);     %Elem edges left-right
Ey= linspace(B(3),B(4),K(2)+1);     %Elem edges bottom-top
%Global domain initialization (calculated)
Np= N+1;                            %Number of interpolation points to achieve N
w= zeros(Np*K(2),Np*K(1));          %Global vorticity at interp points

