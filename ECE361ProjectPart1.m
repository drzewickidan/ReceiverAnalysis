%Group 2, Karanveer Singh, Daniel Drzewicki, Stephen Shetzline
%% Part 1
clc
clear all
close all
data = xlsread('shankar_project_spring#2.xls').';
noTarg = sort(data(1:70));
targ = sort(data(71:100));
z = [zeros(70,1).'  ones(30,1).']; %array of 70's for noTarg and 30 1's for Targ
dataMat=[z;data];
B=(sortrows(dataMat',-2))'; % 2x100 matrix in descending order with each value assigned 1 for Targ and 0 for noTarg

num0=[];  % a cumulative array for number of 0's 
num1=[];  % a cumulative array for number of 1's 
temp = linspace(0, 100000000000000, 101)
zeroMat=B(1,:);
l =[];
for i=1:numel(data)+1    % This loop is used to find the values for ROC
    temp(i);
    num0=[num0 sum(temp<1)];
    num1=[num1 sum(temp>0)];
    num0 = num0./70;    
    num1 = num1./30;
    num0
end

plot(num0,num1)   % This is the plot for the ROC
xlabel('Pf')
ylabel('Pd')
title('Group #2 ROC Plot')
%%
min=999;    % This variable will become the minimum distance
for i=1:numel(num0)  % This for loop is for calculating the minimum distance from ROC curve to (0,1)
    temp=[0,1;num0(i),num1(i)];
    dist = pdist(temp);
    if dist<min
        min=dist;
        minx=num0(i);   % X value of minimum distance
        miny=num1(i);   % X value of minimum distance
    end
    
end
hold on
plot([0, minx], [1, miny])
plot(minx,miny,'*')
plot(0.0571,0.7667,'*')

text(0.5,0.5, ' Min Dist to (0,1) = 0.2021')
text(0.5,0.45, ' Intersection Dist to (0,1) = 0.2402')
plot([0,0.0571],[1,0.76667],'m')
legend('ROC','Minimum Distance','Optimal Point at (0.1143,0.8333)','Intersection Point at (0.0571,0.7667)','Distance Intersection to (0,1)')
areaUnderROC=trapz(num0,num1);   % area under ROC curve
areaUnderROCCum=cumtrapz(num0,num1);
threshold = 3.7813 % index 34
% This section is for the density plots
figure
ksdensity(targ)
hold on
ksdensity(noTarg)
[fnoTarg,xnoTarg] = ksdensity(noTarg);
[fTarg,xTarg] = ksdensity(targ);
legend('PDF of Target Present','PDF of Target Absent')
%threshold=3.6804;
%plot([3.6804,3.6804],[0,0.175])
legend('PDF of Target Present','PDF of Target Absent')
xlabel('Input Data')
ylabel('Estimated Density')
xnoTarg2=xnoTarg>threshold;   %Trying to find values that are False Alarm
xTarg2=xTarg<threshold;       %Trying to find values that are Miss  
a=area(xnoTarg(69:100), fnoTarg(69:100),'FaceColor','flat');
a1=area(xTarg(1:38), fTarg(1:38),'FaceColor','flat');
a.FaceColor=[0 1 0.1];
a1.FaceColor=[1 0.5 0.5];
text(2,0.02, ' P_M')
text(4.2,0.02, ' P_F')
text(6,0.25,'Threshold (optimum)=3.7813')
title('Group #2')
figure
ksdensity(targ)
hold on

xlabel('Input Data')
ylabel('Estimated Density')


ksdensity(noTarg)

notargMat=[xnoTarg;fnoTarg];
targMat=[xTarg;fTarg];
intersection = InterX(notargMat,targMat);
plot(intersection(1,2),intersection(2,2),'*')
legend('PDF of Target Present','PDF of Target Absent','Threshold=3.9858')
title('Group #2')
%Code for shading area of False Alarm
xxx=intersection(1,2):0.01:6;
polyNoTarg=polyfit(xnoTarg,fnoTarg,20);
p11=polyval(polyNoTarg,xxx);
a2=area(xxx,p11);
a2.FaceColor=[0 1 0.1];


%Code for shading area of Miss
xxx2 = -1:0.01:intersection(1,2);
polyTarg=polyfit(xTarg,fTarg,20);
p22=polyval(polyTarg,xxx2);
a3=area(xxx2,p22);
a3.FaceColor=[1 0.5 0.5];

text(2,0.02, ' P_M')
text(4.2,0.02, ' P_F')

thresholdIntersection = intersection(1,2);
%a2=area(xnoTarg(71:100), fnoTarg(71:100),'FaceColor','flat');
%a3=area(xTarg(1:40), fTarg(1:40),'FaceColor','flat');

%calculating PPV for Optimum Threshold
sense = sum(targ>threshold)/30
spec = sum(noTarg<threshold)/70
pFalseAlarm = sum(noTarg>threshold)/70
PTP = sense*(30/100)+((pFalseAlarm)*70/100);
PPV=(sense*(30/100))/(PTP)


%calculating PPV for Intersection Threshold
senseInt = sum(targ>thresholdIntersection)/30
specInt = sum(noTarg<=thresholdIntersection)/70
pFalseAlarmInt = sum(noTarg>thresholdIntersection)/70
PTPInt = senseInt*(30/100)+((pFalseAlarmInt)*70/100);
PPVInt=(senseInt*(30/100))/(PTPInt)



Nf = sum(noTarg>threshold);
Nc = sum(targ>threshold);
NfInt = sum(noTarg>=thresholdIntersection);
NcInt = sum(targ>=thresholdIntersection);




fprintf('\n\n\n-----------------------CONFUSION MATRIX OPTIMUM------------------------\n')
fprintf('Data Collected   Target Detected   Target Not Detected   Total Count\n')
fprintf('-----------------------------------------------------------------------\n')
fprintf('Target Absent        %.f                     %.f                %.f\n',Nf,70-Nf,70)
fprintf('Target Present       %.f                     %.f                %.f\n',Nc,30-Nc,30)
fprintf('Total Count          %.f                    %.f                %.f\n',Nc+Nf,70-Nf+30-Nc,100)
fprintf('-----------------------------------------------------------------------\n')
fprintf('PPV=%.3f',PPV)

fprintf('\n\n\n-----------------------CONFUSION MATRIX INTERSECTION-------------------\n')
fprintf('Data Collected   Target Detected   Target Not Detected   Total Count\n')
fprintf('-----------------------------------------------------------------------\n')
fprintf('Target Absent        %.f                     %.f                %.f\n',NfInt,70-NfInt,70)
fprintf('Target Present       %.f                     %.f                %.f\n',NcInt,30-NcInt,30)
fprintf('Total Count          %.f                    %.f                %.f\n',NcInt+NfInt,70-NfInt+30-NcInt,100)
fprintf('-----------------------------------------------------------------------\n')
fprintf('PPV=%.3f \n',PPVInt)

figure
xAxis = [30:1:70];
Az=areaUnderROC;
A1=Az/(2-Az);
A2=(2*Az^2)/(1+Az);

%yAxis = sqrt((Az.*(1-Az)+(xAxis-1).*(A1-Az^2)+(70-1)*(A2-Az^2))/(xAxis.*70));
yAxis=[];
for i=1:numel(xAxis)
    temp=sqrt((Az*(1-Az)+(xAxis(i)-1).*(A1-Az^2)+(70-1)*(A2-Az^2))/(xAxis(i)*70));
    yAxis=[yAxis, temp];
end
 plot(xAxis,yAxis)


xlabel('Number of samples with (Target Present)')
ylabel('Standard Deviation of ROC area')
hold on
text(45,0.0286, ' Area under ROC=0.9357')
text(45,0.0293, ' Samples Size (Target Absent)=70')
text(45,0.03, ' Samples Size (Target Present)=30')
title('Improvement in Performance-decline in \sigma(A_z)')

performanceIndex=abs(mean(noTarg)-mean(targ))/sqrt(var(noTarg)+var(targ))

function P = InterX(L1,varargin)
%INTERX Intersection of curves
%   P = INTERX(L1,L2) returns the intersection points of two curves L1 
%   and L2. The curves L1,L2 can be either closed or open and are described
%   by two-row-matrices, where each row contains its x- and y- coordinates.
%   The intersection of groups of curves (e.g. contour lines, multiply 
%   connected regions etc) can also be computed by separating them with a
%   column of NaNs as for example
%
%         L  = [x11 x12 x13 ... NaN x21 x22 x23 ...;
%               y11 y12 y13 ... NaN y21 y22 y23 ...]
%
%   P has the same structure as L1 and L2, and its rows correspond to the
%   x- and y- coordinates of the intersection points of L1 and L2. If no
%   intersections are found, the returned P is empty.
%
%   P = INTERX(L1) returns the self-intersection points of L1. To keep
%   the code simple, the points at which the curve is tangent to itself are
%   not included. P = INTERX(L1,L1) returns all the points of the curve 
%   together with any self-intersection points.
%   
%   Example:
%       t = linspace(0,2*pi);
%       r1 = sin(4*t)+2;  x1 = r1.*cos(t); y1 = r1.*sin(t);
%       r2 = sin(8*t)+2;  x2 = r2.*cos(t); y2 = r2.*sin(t);
%       P = InterX([x1;y1],[x2;y2]);
%       plot(x1,y1,x2,y2,P(1,:),P(2,:),'ro')

%   Author : NS
%   Version: 3.0, 21 Sept. 2010

%   Two words about the algorithm: Most of the code is self-explanatory.
%   The only trick lies in the calculation of C1 and C2. To be brief, this
%   is essentially the two-dimensional analog of the condition that needs
%   to be satisfied by a function F(x) that has a zero in the interval
%   [a,b], namely
%           F(a)*F(b) <= 0
%   C1 and C2 exactly do this for each segment of curves 1 and 2
%   respectively. If this condition is satisfied simultaneously for two
%   segments then we know that they will cross at some point. 
%   Each factor of the 'C' arrays is essentially a matrix containing 
%   the numerators of the signed distances between points of one curve
%   and line segments of the other.

    %...Argument checks and assignment of L2
    error(nargchk(1,2,nargin));
    if nargin == 1,
        L2 = L1;    hF = @lt;   %...Avoid the inclusion of common points
    else
        L2 = varargin{1}; hF = @le;
    end
       
    %...Preliminary stuff
    x1  = L1(1,:)';  x2 = L2(1,:);
    y1  = L1(2,:)';  y2 = L2(2,:);
    dx1 = diff(x1); dy1 = diff(y1);
    dx2 = diff(x2); dy2 = diff(y2);
    
    %...Determine 'signed distances'   
    S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
    S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);
    
    C1 = feval(hF,D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
    C2 = feval(hF,D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';

    %...Obtain the segments where an intersection is expected
    [i,j] = find(C1 & C2); 
    if isempty(i),P = zeros(2,0);return; end;
    
    %...Transpose and prepare for output
    i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
    L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
    i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0
    
    %...Solve system of eqs to get the common points
    P = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
                dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows')';
              
    function u = D(x,y)
        u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);
    end
end

