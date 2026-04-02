clear 
addpath("../utils/")
addpath("../code/")
arp_startup

num_ks = 10;
ks = round(logspace(1,3,num_ks));

n = 1e4;
A = diag((1:n).^(-2))*randn(n);

ck_times_dense = zeros(num_ks,1);
rej_times_dense = zeros(num_ks,1);
skarp_times_dense = zeros(num_ks,1);
skpp_times_dense = zeros(num_ks,1);
rpqr_times_dense = zeros(num_ks,1);

for i = 1:length(ks)
    k = ks(i);
    ck_times_dense(i) = timeit(@() arp(A,k,"ck"))
    rej_times_dense(i) = timeit(@() arp(A,k,"rejection"))
    skarp_times_dense(i) = timeit(@() arp(A,k,"rejection","sketchy"))
    skpp_times_dense(i) = timeit(@() sketchy_pivoting(A,k))
    rpqr_times_dense(i) = timeit(@() acc_rpqr(A,k))
end

save("../data/arp_dense.mat")

%% Plot

figure
loglog(ks,rpqr_times_dense,":x","Color",black,"LineWidth",3); hold on
loglog(ks,ck_times_dense,":","Color",purple,"LineWidth",3)
loglog(ks,skpp_times_dense,"--","Color",blue,"LineWidth",3)
loglog(ks,skarp_times_dense,"Color",orange,"LineWidth",3)
loglog(ks,rej_times_dense,"-.","Color",pink,"LineWidth",3)
xlabel("Rank $k$")
ylabel("Runtime (sec)")
title("dense")
legend({"RPQR","ARP [CK]","SkQR","SkARP (ours)","ARP (ours)"},"Location","Northwest")
mysave(gcf,"../figs/arp_dense_time")