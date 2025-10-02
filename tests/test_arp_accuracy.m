clear 
addpath("../utils/")
addpath("../code/")
arp_startup

n = 100;
x = linspace(0,2,2*n+1); x(end) = [];
x1 = x(1:n)'; x2 = x(n+1:end)';
y = linspace(0,1,n+1)'; y(end) = [];
X = [kron(x1,ones(n,1)) kron(ones(n,1),y)];
Y = [kron(x2,ones(n,1)) kron(ones(n,1),y)];
A1 = 1 ./ pdist2(X,Y);

load("../data/processed_data.mat")
A2 = data';

trials = 100;
num_ks = 10;

As = {A1,A2};

rpqr_errs = zeros(num_ks,trials);
arp_errs = zeros(num_ks,trials);
skarp_errs = zeros(num_ks,trials);
optarp_errs = zeros(num_ks,trials);
skpp_errs = zeros(num_ks,trials);

for A_idx = 1:2
    A = As{A_idx};
    if A_idx == 1
        ks = 50:50:500;
    else
        ks = 5:5:50;
    end

    Afro = norm(A,"fro");
    for i = 1:length(ks)
        k = ks(i);
        for trial = 1:trials
            [idx,W] = arp(A,k);
            arp_errs(i,trial) = norm(A - W*A(idx,:),"fro") / Afro
            [idx,W] = arp(A,k,[],"sketchy");
            skarp_errs(i,trial) = norm(A - W*A(idx,:),"fro") / Afro
            [idx,W] = arp(A,k,[],"optimal");
            optarp_errs(i,trial) = norm(A - W*A(idx,:),"fro") / Afro
            [idx,W] = sketchy_pivoting(A,k);
            skpp_errs(i,trial) = norm(A - W*A(idx,:),"fro") / Afro
            [idx,W] = acc_rpqr(A,k);
            rpqr_errs(i,trial) = norm(A - W*A(idx,:),"fro") / Afro
        end
    end
    
    if A_idx == 1
        save("../data/arp_kernel.mat","arp_errs","skarp_errs","optarp_errs","skpp_errs","rpqr_errs","A_idx","ks")
    else
        save("../data/arp_genetics.mat","arp_errs","skarp_errs","optarp_errs","skpp_errs","rpqr_errs","A_idx","ks")
    end
    
    figure
    merr = mean(arp_errs,2); upper = max(arp_errs,[],2); lower = min(arp_errs,[],2);
    plot_shaded(ks',merr,lower,upper,purple,"Linestyle",":","LineWidth",3)
    merr = mean(skpp_errs,2); upper = max(skpp_errs,[],2); lower = min(skpp_errs,[],2);
    plot_shaded(ks',merr,lower,upper,blue,"Linestyle","--","LineWidth",3)
    merr = mean(skarp_errs,2); upper = max(skarp_errs,[],2); lower = min(skarp_errs,[],2);
    plot_shaded(ks',merr,lower,upper,orange,"LineWidth",3);
    merr = mean(optarp_errs,2); upper = max(optarp_errs,[],2); lower = min(optarp_errs,[],2);
    plot_shaded(ks',merr,lower,upper,pink,"Linestyle","-.","LineWidth",3);
    merr = mean(rpqr_errs,2); upper = max(rpqr_errs,[],2); lower = min(rpqr_errs,[],2);
    plot_shaded(ks',merr,lower,upper,black,"Marker","x","LineWidth",3,"Linestyle",":");
    set(gca,"YScale","log")
    xlabel("Rank $k$")
    ylabel("Accuracy $\|\mbox{\boldmath $A$} - \mbox{\boldmath $\widehat{A}$}\|_{\mathrm{F}} / \|\mbox{\boldmath $A$}\|_{\mathrm{F}}$")

    if A_idx == 1
        title("kernel")
        legend("","ARP","","SkQR","","SkARP","","OptARP","","RPQR")
        saveas(gcf,"../figs/arp_kernel_acc.fig")
        exportgraphics(gcf,"../figs/arp_kernel_acc.png")
    else
        title("genetics")
        saveas(gcf,"../figs/arp_genetics_acc.fig")
        exportgraphics(gcf,"../figs/arp_genetics_acc.png")
    end
    
end