function exampleGenField
%Generates correlated random field based on uniformly distributed random
%variables that is 100 by 100, mean 1 and correlation length 10. The
%covariance function used is the exponential of the distance between every
%point on the grid and every other point
N=100;
M=100;
mu=1;
l=10;
X=genFieldRealization(N,M,l,mu);
imagesc(X)
end


function X=genFieldRealization(N,M,l,mu)
%Generates a realization of an N by M random field
xi=rand(N*M,1);
z=(xi-1/2)/sqrt(12);
C=genCovarianceMatrix(N,M,l);
x=mu+1e-3*chol(C)*z;
X=reshape(x,N,M);
end

function C=genCovarianceMatrix(N,M,l)
%% Function generates the covariance matrix for the random field
C=zeros(N*M,N*M); %Pre-allocating matrix
C=diag(ones(N*M,1))/2; %Assigns the diagonal entries as 1/2 so when you 
%calculate only the upper triangular elements, you can use C+C' to populate
%the rest of the matrix
for i=1:(N*M);
    for j=(i+1):N*M
        %Using the fact that the element coordinates form a lattice of N
        %point in one direction and M in the other
        d2=ceil(i/N)^2-2*ceil(i/N)*ceil(j/N)+ceil(j/N)^2+(mod(i-1, M)-mod(j-1, M))^2;
        C(i,j)=exp(-sqrt(d2)/l); 
    end
end
C=(C'+C);
end



