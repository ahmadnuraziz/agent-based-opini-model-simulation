%initial values
clc
clear
close all
rng(0)
N = 900;
iteration = 120;
mean = 0.497;
degree = mean.*randn(N,1) + 0.35;
attr   = rand(N,3);
C = zeros(N,N);
low_opi = -1; %min
high_opi = 1; %max
low_cpi = -1;
high_cpi = 1;
k_opi = (high_opi-low_opi).*rand(N,1) + low_opi;
k_cpi = (high_cpi-low_cpi).*rand(N,1) + low_cpi;
degree(degree<0) = 0;
degree(degree>1) = 1;
entropy = zeros(5,iteration);

n = 35;
for i=1:n
    j = randi([1,N]);
    while k_opi(j)==0
        j=randi([1,N]);
    end
    k_opi(j) = 0;
    k_cpi(j) = 0;
end
% for i=1:iteration
%     while -0.01<k_opi(i)<0.01
%         (high_opi-low_opi).*rand(N,1) + low_opi;
%     end
%     
%     while -0.01<k_cpi(i)<0.01
%         (high_cpi-low_cpi).*rand(N,1) + low_cpi;
%     end
% end
D = zeros(N,N);
%connection
for j=1:N
    for k=1:N
        if j~=k
            strength = attr(j,1)*attr(k,1) + attr(j,2)*attr(k,2) + attr(j,3)*attr(k,3);
            D(j,k) = strength;
            if (1.5<strength)&&(strength<=3)
                C(j,k) = 2;
            elseif (0<strength)&&(strength)<=1.5
                C(j,k) = 1;
            end
        end
    end
end
alpha_matrix = zeros(N,iteration);
beta_matrix = zeros(N,iteration);
gamma_matrix = zeros(N,iteration);
degree_matrix = zeros(N,iteration);
degree_history = zeros(N,iteration);
clf_history = zeros(N,iteration);
clf = zeros(N,1);

history = zeros(5,iteration);
average = zeros(N,iteration);
%simulation
for i=1:iteration
     
     for j=1:N
        if degree(j)>0.8
            clf(j) = 5;
        elseif (0.6<degree(j))&&(degree(j)<=0.8)
            clf(j) = 4;
        elseif (0.4<degree(j))&&(degree(j)<=0.6)
            clf(j) = 3;
        elseif (0.2<degree(j))&&(degree(j)<=0.4)
            clf(j) = 2;
        elseif degree(j)<0.2
            clf(j) = 1;
        end
     end
    
    history(1,i) = sum(clf==1);
    history(2,i) = sum(clf==2);
    history(3,i) = sum(clf==3);
    history(4,i) = sum(clf==4);
    history(5,i) = sum(clf==5);
    clf_history(:,i) = clf(:);
    degree_history(:,i) = degree(:);
    
    if sum(history(:,i)) ~= N
        disp(i)
    end
    
    if  history(1,i) == 0
        entropy(1,i) = 0;
    else 
        entropy(1,i) = -(history(1,i)/N)*log10(history(1,i)/N);
    end
    
    if  history(2,i) == 0
        entropy(2,i) = 0;
    else 
        entropy(2,i) = -(history(2,i)/N)*log10(history(2,i)/N);
    end
    
    if  history(3,i) == 0
        entropy(3,i) = 0;
    else 
        entropy(3,i) = -(history(3,i)/N)*log10(history(3,i)/N);
    end
    
    if  history(4,i) == 0
        entropy(4,i) = 0;
    else 
        entropy(4,i) = -(history(4,i)/N)*log10(history(4,i)/N);
    end
    
    if  history(5,i) == 0
        entropy(5,i) = 0;
    else 
        entropy(5,i) = -(history(5,i)/N)*log10(history(5,i)/N);
    end
    
    
    alpha = 0;
    beta = 0;
    gamma = 0;
    avg = 0;
    average(:,i) = sum(C,2);
    
    for j=1:N
        for k=1:N
            if j~=k
                alpha = alpha + C(j,k)*attr(j,1);
                beta = beta + C(j,k)*attr(j,2);
                gamma = gamma + C(j,k)*attr(j,3);
                avg = avg + C(j,k)*degree(k);
            end
        end
        alpha = alpha/average(j,i);
        beta  = beta/average(j,i);
        gamma = gamma/average(j,i);
        avg   = avg/average(j,i);
        alpha_matrix(j,i) = alpha;
        beta_matrix(j,i)  = beta;
        gamma_matrix(j,i) = gamma;
        degree_matrix(j,i) = avg;
    end
    
    degree = degree + k_opi.*(degree_matrix(:,i) - degree);
    attr(:,1) = attr(:,1) + k_cpi.*(alpha_matrix(:,i) - attr(:,1));
    attr(:,2) = attr(:,2) + k_cpi.*(beta_matrix(:,i) - attr(:,2));
    attr(:,3) = attr(:,3) + k_cpi.*(gamma_matrix(:,i) - attr(:,3));
    
    for j=1:N
        for k=1:N
            if j~=k
                strength = attr(j,1)*attr(k,1) + attr(j,2)*attr(k,2) + attr(j,3)*attr(k,3);
                D(j,k) = strength;
                if (1.5<strength)&&(strength<=3)
                    C(j,k) = 2;
                elseif (0<strength)&&(strength)<=1.5
                    C(j,k) = 1;
                end
            end
        end
    end
    attr(attr<0) = 0; attr(attr>1) = 1;
    degree(degree<0) = 0; degree(degree>1)=1;
end
x = linspace(1,iteration,iteration);
figure(1)
plot(x,history(1,:),'linewidth',1.2)
hold on
plot(x,history(2,:),'linewidth',1.2)
hold on
plot(x,history(3,:),'linewidth',1.2)
hold on
plot(x,history(4,:),'linewidth',1.2)
hold on
plot(x,history(5,:),'linewidth',1.2)
legend('I','II','III','IV','V')
title('Number of Agents')
xlabel('Iteration')
ylabel('Number of Agents')
ylim([0,N])
savefig('number agents case 1.fig')

figure(2)
plot(x,entropy(1,:),'linewidth',1.2)
hold on
plot(x,entropy(2,:),'linewidth',1.2)
hold on
plot(x,entropy(3,:),'linewidth',1.2)
hold on
plot(x,entropy(4,:),'linewidth',1.2)
hold on
plot(x,entropy(5,:),'linewidth',1.2)
legend('I','II','III','IV','V')
title('entropy')
xlabel('Iteration')
ylabel('Sk')
ylim([0,0.2])
savefig('entropy case 1.fig')

figure(3)

plot(x,sum(entropy),'linewidth',1.2)
title('total entropy')
xlabel('Iteration')
ylabel('S')
ylim([0,1])
savefig('Total entropy case 1.fig')

history(1,iteration)+history(2,iteration)+history(3,iteration)+history(4,iteration)+history(5,iteration)