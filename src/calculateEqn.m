function [matrix] = calculateEqn(L_eta, alpha, gamma, epsilon, phi, k, time_count, dt, X, NP)
    matrix=X;
    for t=1:time_count
        for phase=1:size(NP)
            third_term = zeros(size(X{1}));
            for j=1:size(NP)
                if phase==j
                    continue
                else
                    third_term=third_term+X{j}.^2
                end
            end
            third_term = 2*gamma*X{phase}.*third_term;
            nablax = nabla_FDM(X{phase});
            matrix{t} =matrix{t} + dt*(-L_eta(t)*(-alpha*X{t}(:,:)+alpha*X{t}(:,:).^3+third_term+2*epsilon*X(phase)*phi^2-k(phase)*nablax));
        end
    end
return

function nabla_X=nabla_FDM(X)
    starti = 2;
    endi = size(X,1)-1
% % % % %     Periodic boundary conditions
    X(1,starti:endi) = X(endi,starti:endi);
    X(endi+1,starti:endi) = X(starti,starti:endi);
    X(starti:endi,1) = X(starti:endi,endi);
    X(starti:endi,endi+1) = X(starti:endi,starti);
% % % % %     Periodic boundary conditions end
    nabla_X=X;
    for i=starti:endi
        for j=starti:endi
            nabla_X(i,j) = (X(i+1,j)+X(i-1,j)+X(i,j+1)+X(i,j-1)-4*X(i,j))/1;
        end
    end
return