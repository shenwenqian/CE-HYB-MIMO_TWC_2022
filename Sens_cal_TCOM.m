function [F, W] = Sens_cal_TCOM(N_t,N_r,N_RF,N_s)                     

    B_t = N_t/N_RF;
    B_r = N_r/N_RF;
    T_t = N_s*B_t;
    T_r = N_s*B_r;
    
    U = dftmtx(N_RF);
    V = dftmtx(N_s);

    DFT_t = dftmtx(N_t);
    Index = randperm(N_t);
    F_RF = zeros(N_t,N_t);
    for col_index = 1:N_t
        F_RF(:,col_index) = DFT_t(:,Index(col_index));
    end
    F_BB = zeros(N_t,T_t);
    for b_t = 1:B_t
        F_BB((b_t-1)*N_RF+1:b_t*N_RF,(b_t-1)*N_s+1:b_t*N_s) = U*([eye(N_s) zeros(N_s,N_RF-N_s)]).'*V';
    end
    F_BB = F_BB/norm(F_RF*F_BB,'fro')*sqrt(T_t);
    F = F_RF*F_BB;    

    DFT_r = dftmtx(N_r);
    Index = randperm(N_r);
    W_RF = zeros(N_r,N_r);
    for col_index = 1:N_r
        W_RF(:,col_index) = DFT_r(:,Index(col_index));
    end
    W_BB = zeros(N_r,T_r);
    for b_r = 1:B_r
        W_BB((b_r-1)*N_RF+1:b_r*N_RF,(b_r-1)*N_s+1:b_r*N_s) = U*([eye(N_s) zeros(N_s,N_RF-N_s)]).'*V';
    end
    W_BB = W_BB/norm(W_RF*W_BB,'fro')*sqrt(T_r);
    W = W_RF*W_BB;    

end
