function [F, W] = Sens_cal_Ben(N_t,N_r,N_RF,N_s,B_t,B_r)              

    T_t = N_s*B_t;
    T_r = N_s*B_r;
    M_t = N_RF*B_t;
    M_r = N_RF*B_r;
    
    DFT_t = dftmtx(M_t);
    Index = [1:2:M_t M_t-1:-3:1 1:1:M_t];
    F_RF = zeros(N_t,M_t);
    for row_index = 1:N_t
        F_RF(row_index,:) = DFT_t(:,Index(row_index)).';
    end
    DFT_r = dftmtx(M_r);
    Index = [1:2:M_r M_r-1:-3:1 1:1:M_r];
    W_RF = zeros(N_r,M_r);
    for row_index = 1:N_r
        W_RF(row_index,:) = DFT_r(:,Index(row_index)).';
    end    

    F_BB = zeros(M_t,T_t);
    for b_t = 1:B_t
        U = dftmtx(N_RF);
        V = dftmtx(N_s);
        F_BB((b_t-1)*N_RF+1:b_t*N_RF,(b_t-1)*N_s+1:b_t*N_s) = U*([eye(N_s) zeros(N_s,N_RF-N_s)]).'*V';
    end
    F_BB = F_BB/norm(F_RF*F_BB,'fro')*sqrt(T_t);
    F = F_RF*F_BB;    

    W_BB = zeros(M_r,T_r);
    for b_r = 1:B_r
        U = dftmtx(N_RF);
        V = dftmtx(N_s);
        W_BB((b_r-1)*N_RF+1:b_r*N_RF,(b_r-1)*N_s+1:b_r*N_s) = U*([eye(N_s) zeros(N_s,N_RF-N_s)]).'*V';
    end
    W_BB = W_BB/norm(W_RF*W_BB,'fro')*sqrt(T_r);
    W = W_RF*W_BB;    

end
