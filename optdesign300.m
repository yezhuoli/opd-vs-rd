%% Optimization Function
% With a given matrix $z$ and number $k$, optimization function finds
% solution matrix $x$ and optimal value $opt$ in optimal design model in(at) $300s$.
%% Creat Optimization Function

function [x,opt]= optdesign300(z,k) 

[n,p]=size(z); %%%the size of matrix z.
zopt=[(ones(1,n))' z]; %%%add an all 1 column vector to z.

for i=1:n
    zz{i}= (zopt(i,:))'*zopt(i,:); 
    %%%ith matrix zz = (ith row of zopt)'*(ith row of zopt).
end

x=binvar(n,k); %%%define x as an n*k binary matrix.
d=sdpvar(1);  %%%define d as a symbolic decision variable.

%%%setup constraints of optimaization function here: 
J=cell(1,k);  
 for j=1:k
     H=0;
    for i=1:n
       H=H+(zz{i}*x(i,j));  
       %%%with a given j, H = the sum of [(ith matrix zz) * (xi)],
       %where i is from 1 to n.
    end
    J{j}=H-d*eye(p+1);  %%%with a given j, J{j}= H - d*(identity matrix).
 end

 constraint=[J{1}>=0]; 
 for i=2:k
    constraint=[constraint,J{i}>=0]; 
 end
 %%%constraints: J{i}>=0, for i from 1 to k.
 
 for i=1:n
    constraint=[constraint,sum(x(i,:))==1];  
 end
 %%%add the constraints: the sum of each row of x equals to 1.
 
 for i=1:k
    constraint=[constraint,sum(x(:,i))==n/k];  
 end
 %%%add the constraints: the sum of each column of x equals to n/k.
 
options=sdpsettings('solver','bnb','bnb.maxtime',300,'savesolveroutput',1);
 %%%set the max running time that maxtime=300s.

 
 diagnostics = optimize(constraint,-d,options); 
 %%%function for solving optimization problems with defined constraints,
 %%%objective and options.
 
 x=value(x); %%% x matrix 
 opt=value(d);  %%% optimal value 
 
end

%% Example
% In order to generate a different sequence of random numbers each time,
% 'shuffle' is used to initialize generator.
% Creat a 10-by-4 matrix of normal random numbers from the normal
% distribution with mean $0$ and standard deviation $1$. 
% Then set the number of treatments $k = 4$,
% and find solution and optimal value in optimal design model by using
% optimization function.

rng shuffle
z=normrnd(0,1,[10, 4]);
k=2;
[x,opt]=optdesign300(z,k)



