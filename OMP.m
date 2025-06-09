function [z_est] = OMP(y_v,psi,L)                                      

    [M, G] = size(psi);
    z_index = [];
    Iter = L;
    eta = 0;
    y_res = zeros(M,Iter+1);
    y_res(:,1) = y_v;
    iter_index = 2;
    z_est = zeros(G,1);
    while ((iter_index-1 <= Iter) && (norm((y_res(:,iter_index-1)-y_res(:,iter_index)),2)^2 >= eta))
        product = psi'*y_res(:,iter_index-1);                             
        [~, pos] = max(abs(product)); 
        z_index = [z_index pos];
        z_est(z_index,:) = pinv(psi(:,z_index))*y_v;
        y_res(:,iter_index) = y_v-psi(:,z_index)*z_est(z_index,:);
        iter_index = iter_index+1;
    end

end