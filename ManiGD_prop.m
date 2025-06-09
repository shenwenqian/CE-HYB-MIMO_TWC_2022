function [F_RF] = ManiGD_prop(A,F_RF_ini,F_BB,H_obj)                 

    N = size(A,1);
    M = size(F_BB,1);
    F_RF = F_RF_ini;    
    miu_i = norm(A'*F_RF*F_BB*F_BB'*F_RF'*A-H_obj,'fro');

    f_RF = reshape(F_RF,[],1);                                             
    Grad_F = 2.*A*(2.*A'*F_RF*F_BB*F_BB'*F_RF'*A-H_obj-H_obj')*A'*F_RF*F_BB*F_BB';   
    Grad_Eu = reshape(Grad_F,[],1);                                        
    Grad_Ri = Grad_Eu-real(Grad_Eu.*conj(f_RF)).*f_RF;                    
    Dire = -Grad_Ri;                                                              
        
    eta = 10^(-6);                                                        
    eta_achi = Inf;                                                        
    Iter_out = 5000;
    Iter_in = 5000;
    iter_out = 1;
   
    while (iter_out <= Iter_out) && (eta_achi > eta)
    
        miu_i_last = miu_i;                                                
        Grad_Ri_last = Grad_Ri;                                            
        Dire_last = Dire;                                                  
        f_RF_last = f_RF;                                                  
         
        alpha_bar = 10^(0);                                                
        beta = 0.99;
        sigma = 0.025;
        Sign_Armijo = 0;                                                   
        iter_in = 1;
        while (iter_in <= Iter_in) && (Sign_Armijo ~= 1)            
            Stepsize = beta*alpha_bar;                                     
            Gradvec = Stepsize.*Dire_last;                                 
            f_RF_test = (f_RF+Gradvec)./abs(f_RF+Gradvec);                 
            F_RF_test = reshape(f_RF_test,N,M);
            miu_i_test = norm(A'*F_RF_test*F_BB*F_BB'*F_RF_test'*A-H_obj,'fro');
            ele_1 = miu_i_last-miu_i_test;
            ele_2 = real(-sigma.*Grad_Ri_last'*Gradvec);
            Sign_test = ele_1-ele_2;
            if Sign_test >= 0
                Sign_Armijo = 1;
            else
                beta = beta*beta;                                          
            end
            iter_in = iter_in+1;
        end
        f_RF = f_RF_test;
        F_RF = F_RF_test;
        
        Grad_F = 2.*A*(2.*A'*F_RF*F_BB*F_BB'*F_RF'*A-H_obj-H_obj')*A'*F_RF*F_BB*F_BB';
        Grad_Eu = reshape(Grad_F,[],1);
        Grad_Ri = Grad_Eu-real(Grad_Eu.*conj(f_RF)).*f_RF;
        Dire_trans = Dire_last-real(Dire_last.*conj(f_RF)).*f_RF;          
        beta_PR = Grad_Ri'*(Grad_Ri-Grad_Ri_last)/(Grad_Ri_last'*Grad_Ri_last);   
        Dire = -Grad_Ri+beta_PR.*Dire_trans;                               
        miu_i = norm(A'*F_RF*F_BB*F_BB'*F_RF'*A-H_obj,'fro');
        eta_achi = miu_i_last-miu_i;                                       
        if eta_achi <= eta
            f_RF = f_RF_last;
        end
        iter_out = iter_out+1;
        
    end

end  