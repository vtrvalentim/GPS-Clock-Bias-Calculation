%Defining ER in terms of meters(Earth Radius)
ER=6370000;

%Generating synthetic data - position of satelites and receivier in terms of ER (Eart Radius)
S0=[0.95310;0.8576;0.2775]*ER;
S1=[3.5852;2.07;0]*ER;
S2=[2.9274;2.9274;0]*ER;
S3=[2.6612;0;3.1712]*ER;
S4=[1.4159;0;3.8904]*ER;

%Declaring clock bias and X[S, b] in terms of ER
b=2.354788068*10^(-3)*ER;
X0=[S0;b];

%Defining c as the speed of light
c=299792.458;

%Declaring the initial estimates in terms of ER
Sh0=[0.93310;0.25;0.258819]*ER;
bh0=0*ER;
Xh0=[Sh0;bh0];
%Now generate synthetic data, only for the zero noise
y1=norm(S0-S1)+b;
y2=norm(S0-S2)+b;
y3=norm(S0-S3)+b;
y4=norm(S0-S4)+b;
yb=[y1; y2; y3; y4];


ak=0.3; %declaring step size
Skh=Sh0; %initializing variable Skh with the initial guess value Shatzero
bkh=bh0; %initializing variable bkh with the initial guess value bhatzero
test=1; %initializing test variable
i=1;%initializing counter

while test>0.000000001
%Declaring the current value of Xkhat
Xkh=[Skh;bkh];

%Calculating the jacobian transposed:
r1T=transpose(Skh-S1)/norm(Skh-S1);%calculating unit vector pointing from S1 to current guess
r2T=transpose(Skh-S2)/norm(Skh-S1);%calculating unit vector pointing from S2 to current guess
r3T=transpose(Skh-S3)/norm(Skh-S1);%calculating unit vector pointing from S3 to current guess
r4T=transpose(Skh-S4)/norm(Skh-S1);%calculating unit vector pointing from S4 to current guess
H=[r1T, 1;r2T, 1;r3T, 1;r4T, 1];%calculating the jacobian
HT=transpose(H);

%Calculating the value of h(x)=Rl(S)+b
h1=norm(Skh-S1)+bkh; %for l=1
h2=norm(Skh-S2)+bkh; %for l=2
h3=norm(Skh-S3)+bkh; %for l=3
h4=norm(Skh-S4)+bkh; %for l=4
h=[h1;h2;h3;h4]; %the h vector

%calculating loop test variable = gradient
tt=ak*HT*(yb-h);
test=abs(tt(1,1))+abs(tt(2,1))+abs(tt(3,1))+abs(tt(4,1));

%computing error for plot on step4
eba(i,1)=norm(S0-Skh);
eba(i,2)=norm(b-bkh);

%Iterating Steepest Descent Algorithm: 
Xkh = Xkh+ak*HT*(yb-h);
Skh=Xkh(1:3,1); %initializing variable Skh with the initial guess value Shatzero
bkh=Xkh(4,1); %initializing variable bkh with the initial guess value bhatzero

%run counter (how many times did the loop run?)
i=i+1;

end

%error in meeters for each variable:
errorSa= S0-Skh

%error on clock bias:
errorba= b-bkh

%number of loop iterations
runsa= i


%-----------------------------------------------------------------------------



ak=1; %declaring step size
Skh=Sh0; %initializing variable Skh with the initial guess value Shatzero
bkh=bh0; %initializing variable bkh with the initial guess value bhatzero
testb=1; %initializing test variable
j=1;%initializing counter

while testb>0.000001
%Declaring the current value of Xkhat
Xkh=[Skh;bkh];

%Calculating the inverse of the jacobian:
r1T=transpose(Skh-S1)/norm(Skh-S1);%calculating unit vector pointing from S1 to current guess
r2T=transpose(Skh-S2)/norm(Skh-S1);%calculating unit vector pointing from S2 to current guess
r3T=transpose(Skh-S3)/norm(Skh-S1);%calculating unit vector pointing from S3 to current guess
r4T=transpose(Skh-S4)/norm(Skh-S1);%calculating unit vector pointing from S4 to current guess
H=[r1T, 1;r2T, 1;r3T, 1;r4T, 1];%calculating the jacobian
Hinv=inv(H);

%Calculating the value of h(x)=Rl(S)+b
h1=norm(Skh-S1)+bkh; %for l=1
h2=norm(Skh-S2)+bkh; %for l=2
h3=norm(Skh-S3)+bkh; %for l=3
h4=norm(Skh-S4)+bkh; %for l=4
h=[h1;h2;h3;h4]; %the h vector

%calculating loop test variable = gradient
ttb=ak*Hinv*(yb-h);
testb=abs(ttb(1,1))+abs(ttb(2,1))+abs(ttb(3,1))+abs(ttb(4,1));

%computing error for plot on step4
ebg(j,1)=norm(S0-Skh);
ebg(j,2)=norm(b-bkh);

%Iterating Steepest Descent Algorithm: 
Xkh = Xkh+ak*Hinv*(yb-h);
Skh=Xkh(1:3,1); %initializing variable Skh with the initial guess value Shatzero
bkh=Xkh(4,1); %initializing variable bkh with the initial guess value bhatzero

%run counter (how many times did the loop run?)
j=j+1;

end

%error in meeters for each variable:
errorSb = S0-Skh

%error on clock bias:
errorbb = b-bkh

%number of loop iterations
runsb = j

%---------------
%step 4 plots:

%plots for S error
figure
subplot(2,1,1)       
plot(eba(:,1),'b')
title('S error for Steepest Descent')
xlabel('loop iteration')
ylabel('error in meters')

subplot(2,1,2)     
plot(ebg(:,1),'r') 
title('S error for Gauss-Newton')
xlabel('loop iteration')
ylabel('error in meters')


%plots for b error
figure
subplot(2,1,1)       
plot(eba(:,2),'b')
title('clock bias error for Steepest Descent')
xlabel('loop iteration')
ylabel('error in meters')

subplot(2,1,2)     
plot(ebg(:,2),'r') 
title('clok bias error for Gauss-Newton')
xlabel('loop iteration')
ylabel('error in meters')