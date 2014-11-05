%Evaluate a Nth order Lagrange interpolation for each basis 'n' in n_vect 
%at points 'x' in x_vect and returns a n by x matrix: LEval_n(x)
%Uses a Vandermonde for evaluation, this should be best for many x and
%fewer n (harder to vectorize for many n)
%While clever and elegant, this naive implementation only works for N<~30
%or so due to creeping precision errors
function LEval = LagrangeNaive(N,x_vect,n_vect)
    [Qx,~] = GLquad(N); %Calc interpolation points

    savecount =0;
    for n = n_vect %Scalability prob: not vectorized for n
        savecount =savecount+1;
        NOn =Qx([1:n-1, n+1:N]); %All interp points except the one for the n basis
        a =1; %Initialization param to get recurrence started
        for i=1:numel(NOn)
            %Calculate the coefficients in front of each term in the polynomial
            %(ordered by increasing power a_0 + a_1 x + a_2 x^2 +...)
            %Each step multiplies by the next term in (x-a)(x-b)(x-c)...
            %Essentially we shift the previous step's coefficients up one to
            %account for the multiplication by x and subtract the product of
            %the previous step and the constant in the new term to simulate
            %multiplying by (x-q) where 'q' is the constant from the next term
            %in the multiplicative sequence
            %Think of 'a' containing implicitly the poly with the order
            %determining which coefficient goes with what power of x
            %So if f=(x-2)(x-3) then the evolution of 'a' looks like
            % a= [1]                  f= 1
            % a= [0 1]+[1 0]*-2       f= 1(x-2)     = (0+x) + (-2+0x)
            %   a= [-2 1]               f= -2+x
            % a= [0 -2 1]+[-2 1 0]*-3 f= (x-2)(x-3) = (0-2x+x^2) + (6-3x+0)
            %   a= [6 -5 1]             f= 6-5x+1x^2
            a= [0 a] + [a.*(-NOn(i)) 0];
        end
        Coeff_n(savecount,:) =a./prod(Qx(n)-NOn); %Divide by the denominator of the Lagrange basis
    end
    %Calculate Vandermonde matrix and multiply by our poly coeffs to
    %evaluate our Lagrange basis at each point 'x'
    Vand=(repmat(x_vect',1,N).^repmat(0:N-1,length(x_vect),1))';
    LEval = Coeff_n*Vand;
end