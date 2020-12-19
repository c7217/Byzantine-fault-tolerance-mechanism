xm=200; 
ym=200; 
n=100; 
r=1500; 
count=0; 

sink.x=0.5*xm; 
sink.y=0.5*ym; 

E_node=0.5; 
E_cluster=5; 
E_sink=10; 
% E_concsume=0.1; 
% E_concsume=50*0.000000001; 
ETX=50*0.000000001; 
ERX=50*0.000000001; 
Efs=10*0.000000000001; 
Emp=0.0013*0.000000000001; 
EDA=5*0.000000001; 

do=sqrt(Efs/Emp); 

D_sink=0.; 
Data=4000; 

dead_node=0 ;
Dead=0; 

TT=0; 
YY=0; 

packets_TO_BS=0; 
packets_TO_CH=0; 

flag_first_dead=0; 
F=zeros(n,2); 

t1=clock; 

for x=1:1:n;
	for y=1:1:2;
        F(x,y)=rand(1,1)*xm;
    end
	S(x).xd=F(x,1);
	S(x).yd=F(x,2);
    S(x).E=E_node;
    S(n+1).E=E_cluster;
    S(n+2).E=E_cluster;
    S(n+3).E=E_sink; 
end

K = 5; 
CENTS=F( ceil(rand(K,1)*size(F,1)),:);	   
old_CENTS=zeros(5,2);
DAL=zeros(size(F,1),K+2);	
CV='okoyobogoc';	

    while old_CENTS ~= CENTS;
        old_CENTS = CENTS;
        for x = 1:size(F,1)
            for y = 1:K  
                DAL(x,y) = norm(F(x,:)-CENTS(y,:)); 
            end
        [Distance, label] = min(DAL(x,1:K)); 
        DAL(x,K+1) = label; 
        DAL(x,K+2) = Distance; 
        DAL(x,K+3) = S(x).E; 
        end

        for x = 1:K
            A = (DAL(:,K+1) == x); 
            CENTS(x,:) = mean(F(A,:));  
            if sum(isnan(CENTS(:))) ~= 0 
                NC = find(isnan(CENTS(:,1)) == 1);       
                for Ind = 1:size(NC,1)
                    CENTS(NC(Ind),:) = F(randi(size(F,1)),:);
                end
            end
        end
    count=count+1;
    end

figure(5);

hold on;
grid on; 

    for x = 1:K
        PT = F(DAL(:,K+1) == x,:);  
        plot(PT(:,1),PT(:,2),CV(2*x-1:2*x),'LineWidth',1);  
        plot(CENTS(:,1),CENTS(:,2),'*r','LineWidth',1);  
    end

S(n+1).xd=CENTS(1,1);     
S(n+1).yd=CENTS(1,2);  
S(n+2).xd=CENTS(2,1);    
S(n+2).yd=CENTS(2,2); 
S(n+3).xd=CENTS(3,1);     
S(n+3).yd=CENTS(3,2); 
S(n+4).xd=CENTS(4,1);     
S(n+4).yd=CENTS(4,2); 
S(n+5).xd=CENTS(5,1);     
S(n+5).yd=CENTS(5,2); 
S(n+6).xd=sink.x; 
S(n+6).yd=sink.y; 
plot(S(n+6).xd,S(n+6).yd,'rx'); 
xlabel('Meter'),ylabel('Meter');

hold on; 

t = timer('TimerFcn', 'stat=false; disp(''Timer!'')',... 
                 'StartDelay',10);
start(t)

stat=true;
while(stat==true)
  disp('.')
  pause(1)
end

for j=1:1:r;
    k=0; 
    
    for  x=1:1:n;
        
        f(x)=rand(1,1); 
            
        if f(x)>=0.5;
            f(x)=1;
        else
            f(x)=0;
        end
        
        DAL(x,K+4) = f(x); 

        distance=DAL(x,K+2); 
        distance1=sqrt( (S(n+6).xd-(S(n+1).xd) )^2 + (S(n+6).yd-(S(n+1).yd) )^2 ); 
        distance2=sqrt( (S(n+6).xd-(S(n+2).xd) )^2 + (S(n+6).yd-(S(n+2).yd) )^2 ); 
        distance3=sqrt( (S(n+6).xd-(S(n+3).xd) )^2 + (S(n+6).yd-(S(n+3).yd) )^2 ); 
        distance4=sqrt( (S(n+6).xd-(S(n+4).xd) )^2 + (S(n+6).yd-(S(n+4).yd) )^2 ); 
        distance5=sqrt( (S(n+6).xd-(S(n+5).xd) )^2 + (S(n+6).yd-(S(n+5).yd) )^2 ); 
        
        if (distance>do) 
			S(x).E=S(x).E-((ETX+EDA)*(4000) + Emp*4000*(distance*distance*distance*distance)); 
            S(n+1).E=S(x).E-((ETX+EDA)*(4000) + Emp*4000*(distance1*distance1*distance1*distance1)); 
            S(n+2).E=S(x).E-((ETX+EDA)*(4000) + Emp*4000*(distance2*distance2*distance2*distance2)); 
            S(n+3).E=S(x).E-((ETX+EDA)*(4000) + Emp*4000*(distance3*distance3*distance3*distance3)); 
            S(n+4).E=S(x).E-((ETX+EDA)*(4000) + Emp*4000*(distance4*distance4*distance4*distance4)); 
            S(n+5).E=S(x).E-((ETX+EDA)*(4000) + Emp*4000*(distance5*distance5*distance5*distance5)); 
        else
            S(x).E=S(x).E-((ETX+EDA)*(4000) + Efs*4000*(distance*distance)); 
            S(n+1).E=S(x).E-((ETX+EDA)*(4000) + Efs*4000*(distance1*distance1)); 
            S(n+2).E=S(x).E-((ETX+EDA)*(4000) + Efs*4000*(distance2*distance2)); 
            S(n+3).E=S(x).E-((ETX+EDA)*(4000) + Efs*4000*(distance3*distance3)); 
            S(n+4).E=S(x).E-((ETX+EDA)*(4000) + Efs*4000*(distance4*distance4)); 
            S(n+5).E=S(x).E-((ETX+EDA)*(4000) + Efs*4000*(distance5*distance5)); 
        end
        
        if S(x).E<=0 ; 
            
            plot(S(x).xd,S(x).yd,'r.','LineWidth',2); 
            
            if dead_node==0 && flag_first_dead==0;
                first_dead=j; 
                S(x).dead_node=1; 
                dead_node=dead_node+1; 
                flag_first_dead=flag_first_dead+1; 
                disp(['the first dead node round',num2str(x),'happen',num2str(etime(clock,t1)),'min',num2str(j),'round']);
                S(x).deadtime=j; 
            end
            if dead_node>=1 && S(x).dead_node==0;
               S(x).dead_node=1; 
                dead_node=dead_node+1;  
                S(x).deadtime=j; 
            end

        else
        S(x).dead_node=0;     
        end

        if cumsum(S(x).dead_node)>n/2;
            break
        else
            k=k+f(x); 
        end
    end
    
    hold on;

STATISTICS(j).ENERGY=0; 
STATISTICS(j).DEAD=dead_node; 

    for x=1:1:n; 
    
        if S(x).E > 0 
            STATISTICS(j).ENERGY = STATISTICS(j).ENERGY+S(x).E; 
        end
    end
    
data{j+1,1}=n;
data{j+1,2}=sum(f(x));
data{j+1,3}=k;
    
G1=f(DAL(:,6)==1); 
G2=f(DAL(:,6)==2); 
G3=f(DAL(:,6)==3); 
G4=f(DAL(:,6)==4); 
G5=f(DAL(:,6)==5); 

G11=sum(G1); 
G22=sum(G2); 
G33=sum(G3); 
G44=sum(G4); 
G55=sum(G5); 

    timer1=tic; 
    timer2=tic; 
        
	if cumsum(k)<=n/2; 
        Ans='NO';
        toc(timer1); 
        data{j+1,4}=Ans;
        data{j+1,5}=toc(timer1);
        data{j+1,6}='NULL';
        data{j+1,7}=0;
        data{j+1,8}=toc(timer1);
    end
    
    if cumsum(k)>n/2; 
        Ans='YES';
        toc(timer1); 
        data{j+1,4}=Ans;
        timer2=tic; 
        
        if G11>size(G1)/2;
            G1_head=1;
        else
            G1_head=0;
        end
        
        if G22>size(G2)/2;
            G2_head=1;
        else
            G2_head=0;
        end
        
        if G33>size(G3)/2;
            G3_head=1;
        else
            G3_head=0;
        end
        
        if G44>size(G4)/2;
            G4_head=1;
            else
            G4_head=0;
        end
        
        if G55>size(G5)/2;
            G5_head=1;
            else
            G5_head=0;
        end
    
    
        SOP = (G1_head & G2_head & G3_head) | (G1_head & G2_head & G4_head) | (G1_head & G2_head & G5_head)| (G1_head & G3_head & G4_head)| (G1_head & G3_head & G5_head)| (G1_head & G4_head & G5_head)| (G2_head & G3_head & G4_head)| (G2_head & G3_head & G5_head)| (G2_head & G4_head & G5_head)| (G3_head & G4_head & G5_head);
        SOP_sum = sum(SOP);    
                
        if SOP_sum == 1;
            Answer = 'YES';
        else
            Answer = 'NO';
        end;
        
        toc(timer2); 
        
        data{j+1,6}=Answer;
        data{j+1,5}=toc(timer1);
        data{j+1,7}=toc(timer2);
        data{j+1,8}=toc(timer1)+toc(timer2);    

    end    
end

for x=1:1:n; 

    Dead=Dead+S(x).dead_node; 
end

disp(['execute',num2str(j),'round',num2str(Dead),'dead nodes']);
    
for j=2:r 
    mylive(j) = n-STATISTICS(j).DEAD; 
	myenergy(j) = STATISTICS(j).ENERGY; 
end

figure(2); 
hold on;

mylive(1)=100;
myenergy(1)=S(1).E+(n-1)*E_node;
hist(myenergy,10);

figure(3); 
hold on; 
plot(mylive,'color','r'); 
grid on;
xlabel('round');
ylabel('surviving of nodes');
title('Number of surviving nodes diagram');

figure(4); 
hold on; 
plot(myenergy,'color','r'); 
grid on;
xlabel('round');
ylabel('energy of nodes');
title('Entire network life cycle diagram');

timer3=tic;  
    
for j=1:1:r;
    for x=1:1:n;
        if(S(x).E>0)
            if( data{j+1,3} <= 50)
                S(x).packets_TO_CH=packets_TO_CH+Data; 
                S(n+1).packets_TO_BS=0; 
                S(n+2).packets_TO_BS=0; 
                S(n+3).packets_TO_BS=0; 
                S(n+4).packets_TO_BS=0; 
                S(n+5).packets_TO_BS=0; 
                S(n+6).D_sink=0; 
            end
            
            if( data{j+1,3} > 50)
                S(x).packets_TO_CH=packets_TO_CH+Data; 
                S(n+1).packets_TO_BS=packets_TO_BS+Data; 
                S(n+2).packets_TO_BS=packets_TO_BS+Data; 
                S(n+3).packets_TO_BS=packets_TO_BS+Data; 
                S(n+4).packets_TO_BS=packets_TO_BS+Data; 
                S(n+5).packets_TO_BS=packets_TO_BS+Data; 
                S(n+6).D_sink=D_sink+S(n+1).packets_TO_BS+S(n+2).packets_TO_BS+S(n+3).packets_TO_BS+S(n+4).packets_TO_BS+S(n+5).packets_TO_BS;
            end
    TT=TT+S(x).packets_TO_CH;        
        end
        
	
    
    toc(timer3);
     
    YY=YY+S(n+6).D_sink; 

	data{j+1,9}=TT;
    data{j+1,10}=YY;   
    data{j+1,11}=toc(timer3);	
    
    end    
end

data{j+1,12}=num2str(etime(clock,t1));

[status, message] = xlswrite(xlsFile, data, sheetName);
dos(['start ' xlsFile]);