function [SumRate,Precoder]=Algorithm1(H,P,MAX_ITER,InitialPrecoder)
N=size(H{1},2);
Nr=[];
K=length(H);
V=cell(K,1);
Precoder=InitialPrecoder;
Hkbar=[];
for k=1:K
    Nr=[Nr size(H{k},1)];
    if k==1
        V{k}=eye(N);
    else
        Hkbar=[Hkbar;H{k-1}];
        V{k}=null(Hkbar);
    end
end
Nr_cum_sum=cumsum([0,Nr]);
Nr_bar=N-Nr_cum_sum;
idx=cumsum([0 Nr_bar]);
SumRate = zeros(MAX_ITER,1);
%% CVX-based Iterative method
cvx_expert true;
for iter=1:MAX_ITER
    cvx_begin quiet
    variable S(sum(Nr_bar),sum(Nr_bar)) complex;
    obj=0;
    power=0;
    for k=1:K
        Phi_k=eye(Nr(k));
        S(idx(k)+1:idx(k+1),idx(k)+1:idx(k+1))==hermitian_semidefinite(Nr_bar(k));
        power=power+real(trace(S(idx(k)+1:idx(k+1),idx(k)+1:idx(k+1))));
        A=eye(Nr(k));
        for m=1:k-1
            Phi_k=Phi_k+H{k}*V{m}*Precoder{m}*V{m}'*H{k}';
            A=A+H{k}*V{m}*S(idx(m)+1:idx(m+1),idx(m)+1:idx(m+1))*V{m}'*H{k}';
        end
        obj=obj+log_det(A+H{k}*V{k}*S(idx(k)+1:idx(k+1),idx(k)+1:idx(k+1))*V{k}'*H{k}');
        for m=1:k-1
            obj=obj-real(trace(inv(Phi_k)*H{k}*V{m}*S(idx(m)+1:idx(m+1),idx(m)+1:idx(m+1))*V{m}'*H{k}'));
        end
    end
    maximize (obj)
    subject to
    power <= P;
    cvx_end
    
    % Update the precoders
    for k=1:K
        Precoder{k}=S(idx(k)+1:idx(k+1),idx(k)+1:idx(k+1));
    end
    for k=1:K
        A=eye(Nr(k));
        for m=1:k-1
            A=A+H{k}*V{m}*Precoder{m}*V{m}'*H{k}';
        end
        SumRate(iter)=SumRate(iter)+real(log2(det(A+H{k}*V{k}*Precoder{k}*V{k}'*H{k}'))-log2(det(A)));
    end
end
