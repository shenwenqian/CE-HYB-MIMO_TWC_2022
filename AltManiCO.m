function [F_HYB] = AltManiCO(A,N_RF,N_s,B)        
    
    Iter = 20;                                                            
    iter = 1;
    eta = 10^(-5);                                                         

    [N, G] = size(A);
    M = N_RF*B;
    T = N_s*B;
   
    F_RF = exp((1i*2*pi).*rand(N,M));                                      
    F_BB = zeros(M,T);
    for b = 1:B
        F_BB((b-1)*N_RF+1:b*N_RF,(b-1)*N_s+1:b*N_s) = randn([N_RF,N_s])+randn([N_RF,N_s])*1i;
    end     
    F_BB = F_BB./mean(diag(abs(A'*F_RF*F_BB*F_BB'*F_RF'*A)))^0.5;          
    F_HYB = F_RF*F_BB;                                                     
    
    miu_i = norm(A'*F_HYB*F_HYB'*A-eye(G),'fro');
    eta_achi = Inf;
    miu_i_min = Inf;
    
    while (iter <= Iter) && (eta_achi > eta)
    
        miu_i_last = miu_i;

        [F_RF] = ManiGD_prop(A,F_RF,F_BB,eye(G));
        miu_i = norm(A'*F_RF*F_BB*F_BB'*F_RF'*A-eye(G),'fro');
        
        if miu_i < miu_i_min
            miu_i_min = miu_i;
            F_RF_final = F_RF;
            F_BB_final = F_BB;
        end
        
        [F_BB] = Diagcvx(F_RF,A,N_RF,N_s,eye(G));       
        miu_i = norm(A'*F_RF*F_BB*F_BB'*F_RF'*A-eye(G),'fro');  

        if miu_i < miu_i_min
            miu_i_min = miu_i;
            F_RF_final = F_RF;
            F_BB_final = F_BB;
        end
        
        eta_achi = abs(miu_i-miu_i_last);
        iter = iter+1;
        
    end   
    
    F_RF = F_RF_final;
    F_BB = F_BB_final/norm(F_RF_final*F_BB_final,'fro')*T^(0.5);
    F_HYB = F_RF*F_BB;
    
end    