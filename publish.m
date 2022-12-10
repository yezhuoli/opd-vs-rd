%% Optimization Function
% With the given matrix $z$ and number $k$, optimization function finds
% solution matrix $x$ and optimal value $opt$ in optimal design model in
% $300s$.
%% Creat Optimization Function

function [x,opt]= optdesign300(z,k) 

%%
% Collect the number of samples $n$ and the number of side information $p$ from the size of a known matrix $z$.
[n,p]=size(z); 
%%
%% Add an all 1 column vector to $z$, and obtain an extended matrix $zopt$.
zopt=[(ones(1,n))' z]; 
%%
% Creat matrices $zz\left \{i  \right \}$, where $i$ is from $1$ to $n$, such that the $i_{th}$ matrix $zz\left \{i  \right \}=\left ( i_{th} \textit{ row of zopt}  \right ))^T \times \left ( i_{th} \textit{ row of zopt}  \right )$.
for i=1:n
    zz{i}= (zopt(i,:))'*zopt(i,:); 
end
%% 
% Define $x$ and $d$ as an $n*k$ binary matrix and a symbolic decision variable respectively.
x=binvar(n,k); 
d=sdpvar(1);  
%%
% Setup constraints of optimaization function: 
% With a given $j$, $H = \sum_{i=1}^{n} zz\left \{i  \right \} * x_{ij}$ , where $i$ starts from $1$ to $n$. Then $J\left \{ j \right \}= H - d*I$, where $I$ is an identity matrix and $j$ is from $1$ to $k$.
J=cell(1,k);  
 for j=1:k
     H=0;
    for i=1:n
       H=H+(zz{i}*x(i,j));  
    end
    J{j}=H-d*eye(p+1);  
 end
%%
% Constraint: all $J\left \{ j \right \} \ge 0$.
 constraint=[J{1}>=0]; 
 for i=2:k
    constraint=[constraint,J{i}>=0];
 end
 %%
 %% Add constraints: $\sum_{j=1}^{k} x_{ij}=1$.
 for i=1:n
    constraint=[constraint,sum(x(i,:))==1];  
 end
 %%
 %% Add constraints: $\sum_{i=1}^{n} x_{ij}=\frac{n}{k}$.
 for i=1:k
    constraint=[constraint,sum(x(:,i))==n/k];  
 end
 %% 
 % Set the max running time. In our design, it is $maxtime = 300s$.
 options=sdpsettings('solver','bnb','bnb.maxtime',300,'savesolveroutput',1); 
 %%
 % Apply $'optimize'$ function for solving optimization problems with defined constraints,objective and options.
 diagnostics = optimize(constraint,-d,options); 
 %%
 % Obtain solution matrix $x$ and optimal value $opt$.
 x=value(x);
 opt=value(d);
 
end

%% Example
% Generate a 10-by-4 vector of normal random numbers from the normal
% distribution with mean $0$ and standard deviation $1$. Then set the number of treatments $k = 4$,
% and find solution and optimal value in optimal design model by using
% optimization function.
z=normrnd(0,1,[10, 4]);
k=2;
[x,opt]=optdesign300(z,k)


