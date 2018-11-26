tic
load incomp.mat

NodeNum=16;
[len,d]=size(incomp);
total_number=len;
score=zeros(total_number,NodeNum);
a=randperm(total_number);
 
count0=zeros(1,total_number);% online l_2 norm Hodge mismatch--w0
count1=zeros(1,total_number); %online l_2 norm Hodge mismatch--w1
 
count2=zeros(1,total_number); %online l_1 norm Hodge mismatch--w2
count3=zeros(1,total_number); %online l_1 norm hodge mismatch--w3
 
video_better=incomp(:,1);
video_worse=incomp(:,2); 

w0 = zeros(total_number,NodeNum);
w1 = w0;
w2 = w0;
w3 = w0;

for t=1:total_number-1,
    % l2_norm online hodge
    w0(t+1,:) = w0(t,:);
    w0(t+1,video_better(a(t))) = w0(t,video_better(a(t))) - 1/(t+1000) * (w0(t,video_better(a(t))) - w0(t,video_worse(a(t))) - 1); %\gamma_t = 1/(t+1)
    w0(t+1,video_worse(a(t))) = w0(t,video_worse(a(t))) + 1/(t+1000) * (w0(t,video_better(a(t))) - w0(t,video_worse(a(t))) - 1); %\gamma_t = 1/(t+1)
    
    w1(t+1,:) = w1(t,:);
    w1(t+1,video_better(a(t))) = w1(t,video_better(a(t))) - 1/(total_number) * (w1(t,video_better(a(t))) - w1(t,video_worse(a(t))) - 1); %\gamma_t = 1/T
    w1(t+1,video_worse(a(t))) = w1(t,video_worse(a(t))) + 1/(total_number) * (w1(t,video_better(a(t))) - w1(t,video_worse(a(t))) - 1); %\gamma_t = 1/T
    % l_1 norm online hodge
    w2(t+1,:) = w2(t,:);
    w2(t+1,video_better(a(t))) = w2(t,video_better(a(t))) - 1/(t+1000) * sign((w2(t,video_better(a(t))) - w2(t,video_worse(a(t))) - 1)); %\gamma_t = 1/(t+1)
    w2(t+1,video_worse(a(t))) = w2(t,video_worse(a(t))) + 1/(t+1000) * sign((w2(t,video_better(a(t))) - w2(t,video_worse(a(t))) - 1)); %\gamma_t = 1/(t+1)
    
    w3(t+1,:) = w3(t,:);
    w3(t+1,video_better(a(t))) = w3(t,video_better(a(t))) - 1/(total_number) * sign((w3(t,video_better(a(t))) - w3(t,video_worse(a(t))) - 1)); %\gamma_t = 1/T
    w3(t+1,video_worse(a(t))) = w3(t,video_worse(a(t))) + 1/(total_number) * sign((w3(t,video_better(a(t))) - w3(t,video_worse(a(t))) - 1)); %\gamma_t = 1/T
end
 %%%%%Mismatch ratio computation
for t=1:total_number
    for j=1:total_number
        if w0(t,video_better(j))-w0(t,video_worse(j))<0
            count0(t)=count0(t)+1;
        end
    end
    
    for j=1:total_number
        if w1(t,video_better(j))-w1(t,video_worse(j))<0
            count1(t)=count1(t)+1;
        end
    end
  
    for j=1:total_number
        if w2(t,video_better(j))-w2(t,video_worse(j))<0
            count2(t)=count2(t)+1;
        end
    end
    
    for j=1:total_number
        if w3(t,video_better(j))-w3(t,video_worse(j))<0
            count3(t)=count3(t)+1;
        end
    end
    
end

k0=zeros(1,total_number);
k1=zeros(1,total_number);
k2=zeros(1,total_number);
k3=zeros(1,total_number);
 
for i=1:total_number,
    k0(i)=count0(i)./total_number;
    k1(i)=count1(i)./total_number;
    k2(i)=count2(i)./total_number;
    k3(i)=count3(i)./total_number;    
end
 
plot(1:1:total_number,k0,'r-o',1:1:total_number,k1,'-b*',1:1:total_number,k2,'-g+',1:1:total_number,k3,'-y*','LineWidth',2,'MarkerSize',1)
ylim([0 0.5]);  
legend('Online l_2 HodgeRank \gamma_t=1/(t+1000)','Online l_2 HodgeRank \gamma_t=1/T','Online l_1 HodgeRank \gamma_t=1/(t+1000)','Online l_1 HodgeRank \gamma_t=1/T') ;
toc