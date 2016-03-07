function [th,xi,squared_sse, prim_feas, time_vec] = FUNCvxVanilla(x,Y,eps1,eps2, Max_Iter, rho)


tol_thres = 0;

[n,d] = size(x);

X = zeros(n*d,1);       %converting into vector for convenience of .mex files
for i = 1:n,
    for k = 1:d,
        X((i-1)*d+k) = x(i,k);
    end
end  

% Some notation
Del = zeros(n*n*d,1); 
th = zeros(n,1); th_1 = zeros(n,1);
xi = zeros(n*d,1); xi_1 = zeros(n*d,1);
eta = zeros(n*n,1); eta_1 = zeros(n*n,1);
nu = zeros(n*n,1); nu_1 = zeros(n*n,1);
beta = zeros(n*n,1);


Del = a3(n,d,X);


% Defines sum of Delta_ij's
InvDel = zeros(d,d,n);
SumDel = zeros(d,d,n);
Inv = zeros(d,d);

if (d == 1)
    InvDel = a4(n,d,Del);
else
    for i=1:n,
        SumDel(:,:,i) = a5(n,d,(i-1),Del);
        InvDel(:,:,i) = inv(SumDel(:,:,i));
    end
end



%% initialization:

th = x*(x'*x)^(-1)*x'*Y; %%% Different starting point
for i = 1:n
   for k = 1:d
       xi((i-1)*d+k) = th(k);
   end
end
        
for i = 1:n
    for j = 1:n
        val = 0;
        for k =1:d
            val = val+Del((i-1)*n*d+(j-1)*d+k)*xi((j-1)*d+k);
        end
        eta((i-1)*n+j) = th(j)-th(i)+val-nu((i-1)*n+j)/rho;
        if (eta((i-1)*n+j)>0)
            eta((i-1)*n+j) = 0;
        end
    end
end


time_vec = zeros(Max_Iter,1);
time_sum = 0;
sqrt_sse = zeros(Max_Iter,1);
viol = zeros(Max_Iter,1);
prim_feas = zeros(Max_Iter,1);
mov_viol = 0;  mov_viol2 = 0;
grad_prim = zeros(Max_Iter,1);
 

m = 0; iter = 0;
nu = zeros(n*n,1);
k = (1>0);


% Starting the iterations
while (k==(1>0))
    
    tt = tic;
    mov_viol = mov_viol2;
    LL = (1>0);
    iter = iter+1;
    
    
        objold = norm(Y - th)^2;
        
        % Step 1: Updating the sub-gradients, i.e. the xi's "ksee's"
        xi_1 = d1(n,d,rho,nu,eta,th,Del,InvDel);
        
        % Step 2: Updating the function values, i.e. eventually the thetas.
        eta_T = d2(n,eta,Del,xi_1,d);
        
        
        % b1(n,a) = -D'*a
        tmp_th1= (Y - b1(n,nu) - rho*b1(n,eta_T)); 
        th_1 = tmp_th1/(2*n*rho+1) + (sum(sum(tmp_th1))*2*rho/(1+2*n*rho));

        
        % Step 3: Updating the residuals
        eta_1 = d3(n,rho,nu,th_1,Del,xi_1,d,tol_thres);
       

        
        beta = d4(n,eta_1,th_1,Del,xi_1,d);
        
        
        xi = xi_1; th = th_1; eta = eta_1;
    
        objnew = norm(Y - th_1)^2;  % Mo removed the ",2" in norm().
       
    
    % Step 4: Updating the dual variables
    nu_1 = nu + rho*beta;
    
    % Re-scale back to original scales
    
    nu = nu_1;

    sqrt_sse(iter) = norm(Y - th);
 
    grad_prim(iter) =norm( (th - Y) + b1(n,nu) );  %%% this is (th - Y - D'nu): checks grad cond wrt th
 
 
     time_vec(iter) = toc(tt);


    
    prim_feas(iter) = norm(beta)/sqrt(n);    


    %Termination KKT conditions
    if (iter>5)
        
                                                                %{
        tol_p = abs(sqrt_sse(iter) - sqrt_sse(iter-1))/normY;
        tol_f = prim_feas(iter);        
        
        if (tol_p > tol_f*mult);rho= rho*rho_mult; end
        if (tol_f > tol_p*mult);rho= rho/rho_mult; end
                                                                 %}
            
        
        if (prim_feas(iter) < eps1 && ((sqrt_sse(iter)-sqrt_sse(iter-1))/sqrt_sse(iter-1)) < eps2)
            k = (1<0);
        end
        
    end
    
    m = m + 1; 
    

    if (iter > Max_Iter)
        k = (1<0);
    end
    

end
    

                                                                
squared_sse = sqrt_sse.^2;
                                                                
                                                                
                                                                
                                                                
                                                                
