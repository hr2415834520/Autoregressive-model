function[x]= AR(X, K, num)
%function return time-series vector that follows to the vector X. Vector X
%has a period K.
%X- forecasting time-series vector
%K- period
%num- number of forecasting periods.
L= floor(length(X)/K);
X1= reshape(X(end-K*L+1:end), K, L);
X2= X1';
Xcomb= zeros(L, K);  %autoregression matrix
for i= 1:L
    for j=1:K
        Xcomb(i, j)= X2(L+1-i, K+1-j);
    end
end
funcPr= inline('[x, sin(x), x.*sin(x)]', 'x');   %add functions
Xpr = zeros(1,K);          
for s=1:num
for i= 1:K
    Xtmp= Xcomb;
    xi= Xtmp(:,i );
    Xtmp(:, i)= [];
    xi(i)=[];
    xj= Xtmp(1,:);
    Xtmp(1, :)= [];
    Gj= funcPr(Xtmp);
    w= lsqr(Gj, xi);
    Xpr(1,i) = funcPr(xj)*w;
end
H= Xcomb;
Xcomb= zeros(L+s, K);
Xcomb= cat(1, Xpr, H);
end
Xans= Xcomb(1:num, :);
for s=1:num
    for i=1:K
        x(i+(s-1)*K)= Xans(num+1-s,K+1-i);  %forecasting time series
    end
end
