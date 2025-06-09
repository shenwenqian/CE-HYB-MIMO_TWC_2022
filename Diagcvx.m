function [F_BB] = Diagcvx(F_RF,A,N_RF,N_s,H_obj)
    
    M = size(F_RF,2);
    B = M/N_RF;
    T = N_s*B;
    F_BB = zeros(M,T);
    Q_BB = zeros(M,M);                                                     
    
    switch B
        
        case 16
            cvx_begin quiet
                variable Q_BB1(N_RF,N_RF) complex semidefinite
                variable Q_BB2(N_RF,N_RF) complex semidefinite
                variable Q_BB3(N_RF,N_RF) complex semidefinite
                variable Q_BB4(N_RF,N_RF) complex semidefinite
                variable Q_BB5(N_RF,N_RF) complex semidefinite
                variable Q_BB6(N_RF,N_RF) complex semidefinite
                variable Q_BB7(N_RF,N_RF) complex semidefinite
                variable Q_BB8(N_RF,N_RF) complex semidefinite
                variable Q_BB9(N_RF,N_RF) complex semidefinite
                variable Q_BB10(N_RF,N_RF) complex semidefinite
                variable Q_BB11(N_RF,N_RF) complex semidefinite
                variable Q_BB12(N_RF,N_RF) complex semidefinite
                variable Q_BB13(N_RF,N_RF) complex semidefinite
                variable Q_BB14(N_RF,N_RF) complex semidefinite
                variable Q_BB15(N_RF,N_RF) complex semidefinite
                variable Q_BB16(N_RF,N_RF) complex semidefinite
                minimize(norm(A'*F_RF*blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4,Q_BB5,Q_BB6,Q_BB7,Q_BB8,Q_BB9,Q_BB10,Q_BB11,Q_BB12,Q_BB13,Q_BB14,Q_BB15,Q_BB16)*F_RF'*A-H_obj,'fro'));
            cvx_end
            Q_BB = blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4,Q_BB5,Q_BB6,Q_BB7,Q_BB8,Q_BB9,Q_BB10,Q_BB11,Q_BB12,Q_BB13,Q_BB14,Q_BB15,Q_BB16);
            
        case 15
            cvx_begin quiet
                variable Q_BB1(N_RF,N_RF) complex semidefinite
                variable Q_BB2(N_RF,N_RF) complex semidefinite
                variable Q_BB3(N_RF,N_RF) complex semidefinite
                variable Q_BB4(N_RF,N_RF) complex semidefinite
                variable Q_BB5(N_RF,N_RF) complex semidefinite
                variable Q_BB6(N_RF,N_RF) complex semidefinite
                variable Q_BB7(N_RF,N_RF) complex semidefinite
                variable Q_BB8(N_RF,N_RF) complex semidefinite
                variable Q_BB9(N_RF,N_RF) complex semidefinite
                variable Q_BB10(N_RF,N_RF) complex semidefinite
                variable Q_BB11(N_RF,N_RF) complex semidefinite
                variable Q_BB12(N_RF,N_RF) complex semidefinite
                variable Q_BB13(N_RF,N_RF) complex semidefinite
                variable Q_BB14(N_RF,N_RF) complex semidefinite
                variable Q_BB15(N_RF,N_RF) complex semidefinite
                minimize(norm(A'*F_RF*blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4,Q_BB5,Q_BB6,Q_BB7,Q_BB8,Q_BB9,Q_BB10,Q_BB11,Q_BB12,Q_BB13,Q_BB14,Q_BB15)*F_RF'*A-H_obj,'fro'));
            cvx_end
            Q_BB = blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4,Q_BB5,Q_BB6,Q_BB7,Q_BB8,Q_BB9,Q_BB10,Q_BB11,Q_BB12,Q_BB13,Q_BB14,Q_BB15);
            
        case 14
            cvx_begin quiet
                variable Q_BB1(N_RF,N_RF) complex semidefinite
                variable Q_BB2(N_RF,N_RF) complex semidefinite
                variable Q_BB3(N_RF,N_RF) complex semidefinite
                variable Q_BB4(N_RF,N_RF) complex semidefinite
                variable Q_BB5(N_RF,N_RF) complex semidefinite
                variable Q_BB6(N_RF,N_RF) complex semidefinite
                variable Q_BB7(N_RF,N_RF) complex semidefinite
                variable Q_BB8(N_RF,N_RF) complex semidefinite
                variable Q_BB9(N_RF,N_RF) complex semidefinite
                variable Q_BB10(N_RF,N_RF) complex semidefinite
                variable Q_BB11(N_RF,N_RF) complex semidefinite
                variable Q_BB12(N_RF,N_RF) complex semidefinite
                variable Q_BB13(N_RF,N_RF) complex semidefinite
                variable Q_BB14(N_RF,N_RF) complex semidefinite
                minimize(norm(A'*F_RF*blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4,Q_BB5,Q_BB6,Q_BB7,Q_BB8,Q_BB9,Q_BB10,Q_BB11,Q_BB12,Q_BB13,Q_BB14)*F_RF'*A-H_obj,'fro'));
            cvx_end
            Q_BB = blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4,Q_BB5,Q_BB6,Q_BB7,Q_BB8,Q_BB9,Q_BB10,Q_BB11,Q_BB12,Q_BB13,Q_BB14);
            
        case 13
            cvx_begin quiet
                variable Q_BB1(N_RF,N_RF) complex semidefinite
                variable Q_BB2(N_RF,N_RF) complex semidefinite
                variable Q_BB3(N_RF,N_RF) complex semidefinite
                variable Q_BB4(N_RF,N_RF) complex semidefinite
                variable Q_BB5(N_RF,N_RF) complex semidefinite
                variable Q_BB6(N_RF,N_RF) complex semidefinite
                variable Q_BB7(N_RF,N_RF) complex semidefinite
                variable Q_BB8(N_RF,N_RF) complex semidefinite
                variable Q_BB9(N_RF,N_RF) complex semidefinite
                variable Q_BB10(N_RF,N_RF) complex semidefinite
                variable Q_BB11(N_RF,N_RF) complex semidefinite
                variable Q_BB12(N_RF,N_RF) complex semidefinite
                variable Q_BB13(N_RF,N_RF) complex semidefinite
                minimize(norm(A'*F_RF*blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4,Q_BB5,Q_BB6,Q_BB7,Q_BB8,Q_BB9,Q_BB10,Q_BB11,Q_BB12,Q_BB13)*F_RF'*A-H_obj,'fro'));
            cvx_end
            Q_BB = blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4,Q_BB5,Q_BB6,Q_BB7,Q_BB8,Q_BB9,Q_BB10,Q_BB11,Q_BB12,Q_BB13);
            
        case 12
            cvx_begin quiet
                variable Q_BB1(N_RF,N_RF) complex semidefinite
                variable Q_BB2(N_RF,N_RF) complex semidefinite
                variable Q_BB3(N_RF,N_RF) complex semidefinite
                variable Q_BB4(N_RF,N_RF) complex semidefinite
                variable Q_BB5(N_RF,N_RF) complex semidefinite
                variable Q_BB6(N_RF,N_RF) complex semidefinite
                variable Q_BB7(N_RF,N_RF) complex semidefinite
                variable Q_BB8(N_RF,N_RF) complex semidefinite
                variable Q_BB9(N_RF,N_RF) complex semidefinite
                variable Q_BB10(N_RF,N_RF) complex semidefinite
                variable Q_BB11(N_RF,N_RF) complex semidefinite
                variable Q_BB12(N_RF,N_RF) complex semidefinite
                minimize(norm(A'*F_RF*blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4,Q_BB5,Q_BB6,Q_BB7,Q_BB8,Q_BB9,Q_BB10,Q_BB11,Q_BB12)*F_RF'*A-H_obj,'fro'));
            cvx_end
            Q_BB = blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4,Q_BB5,Q_BB6,Q_BB7,Q_BB8,Q_BB9,Q_BB10,Q_BB11,Q_BB12);
            
        case 11
            cvx_begin quiet
                variable Q_BB1(N_RF,N_RF) complex semidefinite
                variable Q_BB2(N_RF,N_RF) complex semidefinite
                variable Q_BB3(N_RF,N_RF) complex semidefinite
                variable Q_BB4(N_RF,N_RF) complex semidefinite
                variable Q_BB5(N_RF,N_RF) complex semidefinite
                variable Q_BB6(N_RF,N_RF) complex semidefinite
                variable Q_BB7(N_RF,N_RF) complex semidefinite
                variable Q_BB8(N_RF,N_RF) complex semidefinite
                variable Q_BB9(N_RF,N_RF) complex semidefinite
                variable Q_BB10(N_RF,N_RF) complex semidefinite
                variable Q_BB11(N_RF,N_RF) complex semidefinite
                minimize(norm(A'*F_RF*blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4,Q_BB5,Q_BB6,Q_BB7,Q_BB8,Q_BB9,Q_BB10,Q_BB11)*F_RF'*A-H_obj,'fro'));
            cvx_end
            Q_BB = blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4,Q_BB5,Q_BB6,Q_BB7,Q_BB8,Q_BB9,Q_BB10,Q_BB11);
            
        case 10
            cvx_begin quiet
                variable Q_BB1(N_RF,N_RF) complex semidefinite
                variable Q_BB2(N_RF,N_RF) complex semidefinite
                variable Q_BB3(N_RF,N_RF) complex semidefinite
                variable Q_BB4(N_RF,N_RF) complex semidefinite
                variable Q_BB5(N_RF,N_RF) complex semidefinite
                variable Q_BB6(N_RF,N_RF) complex semidefinite
                variable Q_BB7(N_RF,N_RF) complex semidefinite
                variable Q_BB8(N_RF,N_RF) complex semidefinite
                variable Q_BB9(N_RF,N_RF) complex semidefinite
                variable Q_BB10(N_RF,N_RF) complex semidefinite
                minimize(norm(A'*F_RF*blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4,Q_BB5,Q_BB6,Q_BB7,Q_BB8,Q_BB9,Q_BB10)*F_RF'*A-H_obj,'fro'));
            cvx_end
            Q_BB = blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4,Q_BB5,Q_BB6,Q_BB7,Q_BB8,Q_BB9,Q_BB10);
            
        case 9
            cvx_begin quiet
                variable Q_BB1(N_RF,N_RF) complex semidefinite
                variable Q_BB2(N_RF,N_RF) complex semidefinite
                variable Q_BB3(N_RF,N_RF) complex semidefinite
                variable Q_BB4(N_RF,N_RF) complex semidefinite
                variable Q_BB5(N_RF,N_RF) complex semidefinite
                variable Q_BB6(N_RF,N_RF) complex semidefinite
                variable Q_BB7(N_RF,N_RF) complex semidefinite
                variable Q_BB8(N_RF,N_RF) complex semidefinite
                variable Q_BB9(N_RF,N_RF) complex semidefinite
                minimize(norm(A'*F_RF*blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4,Q_BB5,Q_BB6,Q_BB7,Q_BB8,Q_BB9)*F_RF'*A-H_obj,'fro'));
            cvx_end
            Q_BB = blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4,Q_BB5,Q_BB6,Q_BB7,Q_BB8,Q_BB9);
            
        case 8
            cvx_begin quiet
                variable Q_BB1(N_RF,N_RF) complex semidefinite
                variable Q_BB2(N_RF,N_RF) complex semidefinite
                variable Q_BB3(N_RF,N_RF) complex semidefinite
                variable Q_BB4(N_RF,N_RF) complex semidefinite
                variable Q_BB5(N_RF,N_RF) complex semidefinite
                variable Q_BB6(N_RF,N_RF) complex semidefinite
                variable Q_BB7(N_RF,N_RF) complex semidefinite
                variable Q_BB8(N_RF,N_RF) complex semidefinite
                minimize(norm(A'*F_RF*blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4,Q_BB5,Q_BB6,Q_BB7,Q_BB8)*F_RF'*A-H_obj,'fro'));
            cvx_end
            Q_BB = blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4,Q_BB5,Q_BB6,Q_BB7,Q_BB8);
            
        case 7
            cvx_begin quiet
                variable Q_BB1(N_RF,N_RF) complex semidefinite
                variable Q_BB2(N_RF,N_RF) complex semidefinite
                variable Q_BB3(N_RF,N_RF) complex semidefinite
                variable Q_BB4(N_RF,N_RF) complex semidefinite
                variable Q_BB5(N_RF,N_RF) complex semidefinite
                variable Q_BB6(N_RF,N_RF) complex semidefinite
                variable Q_BB7(N_RF,N_RF) complex semidefinite
                minimize(norm(A'*F_RF*blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4,Q_BB5,Q_BB6,Q_BB7)*F_RF'*A-H_obj,'fro'));
            cvx_end
            Q_BB = blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4,Q_BB5,Q_BB6,Q_BB7);
            
        case 6
            cvx_begin quiet
                variable Q_BB1(N_RF,N_RF) complex semidefinite
                variable Q_BB2(N_RF,N_RF) complex semidefinite
                variable Q_BB3(N_RF,N_RF) complex semidefinite
                variable Q_BB4(N_RF,N_RF) complex semidefinite
                variable Q_BB5(N_RF,N_RF) complex semidefinite
                variable Q_BB6(N_RF,N_RF) complex semidefinite
                minimize(norm(A'*F_RF*blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4,Q_BB5,Q_BB6)*F_RF'*A-H_obj,'fro'));
            cvx_end
            Q_BB = blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4,Q_BB5,Q_BB6);
            
        case 5
            cvx_begin quiet
                variable Q_BB1(N_RF,N_RF) complex semidefinite
                variable Q_BB2(N_RF,N_RF) complex semidefinite
                variable Q_BB3(N_RF,N_RF) complex semidefinite
                variable Q_BB4(N_RF,N_RF) complex semidefinite
                variable Q_BB5(N_RF,N_RF) complex semidefinite
                minimize(norm(A'*F_RF*blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4,Q_BB5)*F_RF'*A-H_obj,'fro'));
            cvx_end
            Q_BB = blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4,Q_BB5);
            
        case 4
            cvx_begin quiet
                variable Q_BB1(N_RF,N_RF) complex semidefinite
                variable Q_BB2(N_RF,N_RF) complex semidefinite
                variable Q_BB3(N_RF,N_RF) complex semidefinite
                variable Q_BB4(N_RF,N_RF) complex semidefinite
                minimize(norm(A'*F_RF*blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4)*F_RF'*A-H_obj,'fro'));
            cvx_end
            Q_BB = blkdiag(Q_BB1,Q_BB2,Q_BB3,Q_BB4);
            
        case 3
            cvx_begin quiet
                variable Q_BB1(N_RF,N_RF) complex semidefinite
                variable Q_BB2(N_RF,N_RF) complex semidefinite
                variable Q_BB3(N_RF,N_RF) complex semidefinite
                minimize(norm(A'*F_RF*blkdiag(Q_BB1,Q_BB2,Q_BB3)*F_RF'*A-H_obj,'fro'));
            cvx_end
            Q_BB = blkdiag(Q_BB1,Q_BB2,Q_BB3);
            
        case 2
            cvx_begin quiet
                variable Q_BB1(N_RF,N_RF) complex semidefinite
                variable Q_BB2(N_RF,N_RF) complex semidefinite
                minimize(norm(A'*F_RF*blkdiag(Q_BB1,Q_BB2)*F_RF'*A-H_obj,'fro'));
            cvx_end
            Q_BB = blkdiag(Q_BB1,Q_BB2);
            
        case 1
            cvx_begin quiet
                variable Q_BB1(N_RF,N_RF) complex semidefinite
                minimize(norm(A'*F_RF*Q_BB1*F_RF'*A-H_obj,'fro'));
            cvx_end
            Q_BB = Q_BB1;
        
    end
    
    for b = 1:B
        [U, Lambada] = eig(Q_BB((b-1)*N_RF+1:b*N_RF,(b-1)*N_RF+1:b*N_RF));
        F_BB((b-1)*N_RF+1:b*N_RF,(b-1)*N_s+1:b*N_s) = U(:,1:N_s)*Lambada(1:N_s,1:N_s).^(0.5);
    end
    
end
            