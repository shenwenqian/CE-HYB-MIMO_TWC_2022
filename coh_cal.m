function [miu, miu_squ] = coh_cal(Sens)

    [M_BS, G_BS] = size(Sens);
    Sens_norm = zeros(M_BS,G_BS);
    for G_BS_Index = 1:G_BS
        Sens_norm(:,G_BS_Index) = Sens(:,G_BS_Index)./norm(Sens(:,G_BS_Index),'fro');
    end
    Gram_norm = Sens_norm'*Sens_norm;
    location = logical(true(G_BS,G_BS)-tril(true(G_BS,G_BS)));
    miu = max(abs(Gram_norm(location)));
    miu_squ = norm(Gram_norm - diag(diag(Gram_norm)),'fro')^2;
    
end
