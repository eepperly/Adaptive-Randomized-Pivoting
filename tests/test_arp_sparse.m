clear 
addpath("../utils/")
addpath("../code/")
arp_startup

num_ks = 10;
ks = round(logspace(1,3,num_ks));

m = 1e6;
n = 1e4;
A = spdiags(((1:m).^(-2))',0,m,m) * random_sparse(m,n,30);

ck_times_sparse = zeros(num_ks,1);
rej_times_sparse = zeros(num_ks,1);
skarp_times_sparse = zeros(num_ks,1);
skpp_times_sparse = zeros(num_ks,1);

for i = 1:length(ks)
    k = ks(i);
    if max(ck_times_sparse) < 80
        ck_times_sparse(i) = timeit(@() arp(A,k,"ck"))
    end
    rej_times_sparse(i) = timeit(@() arp(A,k,"rejection"))
    skarp_times_sparse(i) = timeit(@() arp(A,k,"rejection","sketchy"))
    if max(skpp_times_sparse) < 80
        skpp_times_sparse(i) = timeit(@() sketchy_pivoting(A,k))
    end
    save("../data/arp_sparse.mat")
end
ck_times_sparse(ck_times_sparse==0) = nan;
skpp_times_sparse(skpp_times_sparse==0) = nan;
save("../data/arp_sparse.mat")

%% Plot
figure
loglog(ks,ck_times_sparse,":","Color",purple,"LineWidth",3); hold on
loglog(ks,skpp_times_sparse,"--","Color",blue,"LineWidth",3)
loglog(ks,skarp_times_sparse,"Color",orange,"LineWidth",3)
loglog(ks,rej_times_sparse,"-.","Color",pink,"LineWidth",3)
xlabel("Rank $k$")
ylabel("Runtime (sec)")
title("sparse")
saveas(gcf,"../figs/arp_sparse_time.fig")
exportgraphics(gcf,"../figs/arp_sparse_time.png")