function [F, W] = Sens_cal_GCC(N_t,N_r,N_RF,B_t,B_r,A_t,A_r)         

    T_t = N_RF*B_t;
    T_r = N_RF*B_r;

    DFT_t = dftmtx(T_t);
    Index = [1:2:T_t T_t-1:-3:1 1:1:T_t];                                  
    F_RF = zeros(N_t,T_t);
    for row_index = 1:N_t
        F_RF(row_index,:) = DFT_t(:,Index(row_index)).';
    end
    DFT_r = dftmtx(T_r);
    Index = [1:2:T_r T_r-1:-3:1 1:1:T_r];
    W_RF = zeros(N_r,T_r);
    for row_index = 1:N_r
        W_RF(row_index,:) = DFT_r(:,Index(row_index)).';
    end

    F_BB = zeros(T_t,T_t);
    for b_t = 1:B_t
        F_RFi = F_RF(:,(b_t-1)*N_RF+1:b_t*N_RF);
        [U, Lambada] = eig(F_RFi.'*conj(A_t)*A_t.'*conj(F_RFi)); 
        F_BB((b_t-1)*N_RF+1:b_t*N_RF,(b_t-1)*N_RF+1:b_t*N_RF) = conj(U*(Lambada^(-0.5))');
    end
    F_BB = F_BB/norm(F_RF*F_BB,'fro')*sqrt(T_t);
    F = F_RF*F_BB;

    W_BB = zeros(T_r,T_r);
    for b_r = 1:B_r
        W_RFi = W_RF(:,(b_r-1)*N_RF+1:b_r*N_RF);
        [U, Lambada] = eig(W_RFi'*A_r*A_r'*W_RFi); 
        W_BB((b_r-1)*N_RF+1:b_r*N_RF,(b_r-1)*N_RF+1:b_r*N_RF) = U*(Lambada^(-0.5))';
    end
    W_BB = W_BB/norm(W_RF*W_BB,'fro')*sqrt(T_r);
    W = W_RF*W_BB;

end
