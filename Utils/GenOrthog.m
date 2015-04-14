function [Qx, Qw] = GenOrthog(N,w,A,B)
%Generate the N roots 'Qx' of a generalized basis of orthogonal polynomials
%wrt weight function 'w' over the domain [A B]. These roots also form a set
%of quadrature abcissa Qx and associated quadrature weights 'Qw'. These are
%very useful for numerically integrating an integrand that is not well
%approximated by a polynomial interpolation. Instead express the integrand
%as f(x)*w(x) where 'f' is the polynomial-like part and 'w' is the
%troublesome part.

%Note that the below setup and first for() loop is kind of "hacky" to
%greatly improve calculation speed, see below "intuitive/nice" code
%commented out for a far more understandable/readable implementation
%
%For ortho poly n+1:
%cn_1(n+1) = <x*x_n,x_n>/<x_n,x_n>
%cn_2(n+1) = <x_n,x_n>/<x_n-1,x_n-1>
%Where inner products are wrt a weight w(x)
%Where we save I2(n) = <x_n,x_n> for reuse in cn_1 and cn_2
%We then use the recurrence relation
%p_n = (x-cn_1)*p_n-1 - cn_2*p_n-2
%to calculate the next ortho poly
%
%The hacky part is that rather than expressing in terms of the previous two
%polys we express in terms of one order poly lower still, by substituing
%the last term's fun algebraicly. The code is hard to read/follow, but the
%benefit is it's far quicker to eval a given fun; recursion is slow in
%Matlab. Ex: Eval p6 -> eval (p4 and p3) -> eval (p2 and p1) and p3
%Since we have "unrolled" p3, p2, and p1 this is much faster than the full
%recusive stack. It also means that a given poly traverses the stack faster
%by leapfrogging over adjacent lower polys
%Ex eval p8->eval (p7 and p6)->eval (p5 and p4) and (p4 and p3) ->
%eval [ (p3 and p2) and (p2 and p1) ] and [ (p3 and p2) and p3 ]
%If we were to right out the full form where only p1 was explicit and with
%no leapfrogging, it would be far more eval steps
p{1}=@(x) 1; %Start recurrence relation
I2(1)= integral(@(x) w(x),A,B,'RelTol',1e-14,'AbsTol',1e-16);
cn_1(2)= integral(@(x) x.*w(x),A,B,'RelTol',1e-14,'AbsTol',1e-16)/I2(1);
p{2}=@(x) x-cn_1(2); %Start recurr. relation (note hacky explicit definition)

I2(2)= integral(@(x) p{2}(x).^2.*w(x),A,B,'RelTol',1e-14,'AbsTol',1e-16);
cn_1(3)= integral(@(x) x.*p{2}(x).^2.*w(x),A,B,'RelTol',1e-14,'AbsTol',1e-16)/I2(2);
cn_2(3)= I2(2)/I2(1); %Save the cn_x terms outside loop for later use below
p{3}=@(x) (x-cn_1(3)).*(x-cn_1(2)) - cn_2(3); %Not req, unrolled for speed

for iter=4:N+1
    I2(iter-1)= integral(@(x) p{iter-1}(x).^2.*w(x),A,B,'RelTol',1e-14,'AbsTol',1e-16);
    cn_1(iter)= integral(@(x) x.*p{iter-1}(x).^2.*w(x),A,B,'RelTol',1e-14,'AbsTol',1e-16)/I2(iter-1);
    cn_2(iter)= I2(iter-1)/I2(iter-2);
    %Note this is the 1 layer deeper recurr. for p_n via p_n-2 and p_n-3
    p{iter}=@(x) ((x-cn_1(iter)).*(x-cn_1(iter-1))-cn_2(iter)).*p{iter-2}(x)-cn_2(iter-1)*(x-cn_1(iter)).*p{iter-3}(x);
end

%----Here's a much more readable version of the above
% p{1}=@(x) 1;
% I2(1)= integral(@(x) p{1}(x).^2.*w(x),A,B,'RelTol',1e-14,'AbsTol',1e-16);
% c2= integral(@(x) x.*p{1}(x).^2.*w(x),A,B,'RelTol',1e-14,'AbsTol',1e-16)/I2(1);
% p{2}=@(x) (x-c2).*p{1}(x);
% 
% for iter=3:N+1
%     I2(iter-1)= integral(@(x) p{iter-1}(x).^2.*w(x),A,B,'RelTol',1e-14,'AbsTol',1e-16);
%     c2= integral(@(x) x.*p{iter-1}(x).^2.*w(x),A,B,'RelTol',1e-14,'AbsTol',1e-16)/I2(iter-1);
%     c1= I2(iter-1)/I2(iter-2);
%     p{iter}=@(x) (x-c2).*p{iter-1}(x) - c1*p{iter-2}(x);
% end

%Here we find roots, we take advantage of the fact than exactly one root
%will appear for p_n+1 in between two roots (or domain bound) of p_n
orthRoot(2,1)= fzero(p{2},[A B]); %Special setup since p_1 has no roots
for poly=3:N+1
    for numroot=1:poly-1
        temp= [A, orthRoot(poly-1,:), B];
        %We have to sequentially look through each area bounded by two
        %roots since fzero finds only one root at a time
        orthRoot(poly,numroot)= fzero(p{poly},[temp(numroot), temp(numroot+1)]);
    end
end

%Here we calculate the associated quadrature weights with the nodes by
%integrating the associate Lagrange basis for the node
%This works because:
%\int f ~ \int interp_f = \int \sum y_i L_i(x) = \sum y_i \int L_i(x)
Qx=orthRoot(N,:)';
nn=elim(Qx(1:N)',Qx(1:N)',[1 3 2]);
Lag= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn(nv,:,:)),bsxfun(@minus,Qx(nv),nn(nv,:,:))),3);

for node=1:N
    Qw(1,node)=integral(@(x) Lag(x,node),A,B,'RelTol',1e-14,'AbsTol',1e-16);
end

end