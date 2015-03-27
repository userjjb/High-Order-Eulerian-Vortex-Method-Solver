close all
clear all
clc

%Global domain initialization (parameters)
B= [1 -1 -1 1];                     %top, bottom, left, right
N= 4;                               %Local poly order
K= [10 20];                         %Num elements along y,x
Ex= linspace(B(3),B(4),K(2)+1);     %Elem edges left-right
Ey= linspace(B(2),B(1),K(1)+1);     %bottom-top
%Global domain initialization (calculated)
Np= N+1;                            %Number of interpolation points to achieve N
w= zeros(Np*K(1),Np*K(2));          %Global vorticity at interp points