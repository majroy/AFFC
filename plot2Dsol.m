% Z = plot2Dsol(VARS,EDGES) returns z-axis intersections 
% Plots 2D projection of contact condition described by VARS and returns z 
% variables to Z, where Z=[z12 z13 zl zu]. 
%
% If EDGES is true, then edges will be shown on roller outline primaries.
% Eq. and Table numbers refer to identifiers in original publication
%
% Companion subfunction to AFFC.m
% M.J. Roy, 2016
function [retZ]=plot2Dsol(vars,Ed)

Rm=vars(1);
to=vars(2);
tf=vars(3);
Rr=vars(4);
R=vars(5);
a=vars(6);
b=vars(7);
P=vars(8);

%intermediate variables
Ri=Rm+to;
d=Rm+tf+R+Rr;
phi=Rr+R/cosd(a);

xu=R*(1-cosd(b)); %Eq. 1
zu=R*sind(b); %Eq. 2

xl=R*(1-cosd(a)); %Eq. 3
zl=-R*sind(a); %Eq. 4

%set up arrays for drawing roller primaries
EEAngles=[(90-a)*pi/180 (90-b)*pi/180];
RNoseS=-EEAngles(1)-6*pi/12;
RNoseE=EEAngles(2)-18*pi/12;

Roller_theta=linspace(RNoseE,RNoseS)';
Roller_NosePnts=R*[cos(Roller_theta) sin(Roller_theta)];
Roller_NosePnts=[Roller_NosePnts(:,1)+d-Rr...
    Roller_NosePnts(:,2)];



if xl<=xu
    xi=(-R*sind(b)+sind(a)*sind(b)*P-R*sind(a)+R*sind(a)*cosd(b)...
        +sind(b)*R*cosd(a))/(cosd(b)*sind(a)+sind(b)*cosd(a)); % Eq. 5
    zi=(-R*cosd(b)+cosd(b)*sind(a)*P+R*cosd(a))/(cosd(b)*sind(a)...
        +sind(b)*cosd(a)); % Eq. 6
else
    xi=(-R*sind(b)+sind(a)*sind(b)*P-R*sind(a)+R*sind(a)*cosd(b)...
        +sind(b)*R*cosd(a))/(cosd(b)*sind(a)+sind(b)*cosd(a)); % Eq. 5
    zi=(-R*cosd(b)+cosd(b)*sind(a)*P+R*cosd(a))/(cosd(b)*sind(a)...
        +sind(b)*cosd(a)); % Eq. 6
end



fprintf('\n');
%Conditions pertaining to Table 1, Eq. 7-10
if zl+P<=zi && zu>zi && xl>xi && xu>xi
    z13=P/2; %Eq. 7
        fprintf('Condition A');
end

if zl+P>=zi && zu<zi && xl<xi && xu<xi
        fprintf('Condition B');
end

if zl+P>=zi && zu>zi && xl<xi && xu>xi
    fprintf('Condition C');

%momentarily convert a to radians
a=a*pi/180;
%hard coded z13 analytical solution, solved with Maple
z13(1)=(-(2*zl*tan(a)^2*P-2*R*zl*tan(a)-2*R*P*tan(a)+2*xl*zl*tan(a)+2*xl*P*tan(a)+zl^2*tan(a)^2+xl^2+P^2*tan(a)^2-2*R*xl-1/(tan(a)^2+1)*(2*tan(a)^2*zl+2*tan(a)*xl-2*R*tan(a)+2*tan(a)^2*P+2*(-2*xl*P*tan(a)+2*R*P*tan(a)-2*xl*zl*tan(a)+R^2*tan(a)^2-zl^2*tan(a)^2-P^2*tan(a)^2-xl^2+2*R*xl-2*zl*tan(a)^2*P+2*R*zl*tan(a))^(1/2))*tan(a)^2*zl-1/(tan(a)^2+1)*(2*tan(a)^2*zl+2*tan(a)*xl-2*R*tan(a)+2*tan(a)^2*P+2*(-2*xl*P*tan(a)+2*R*P*tan(a)-2*xl*zl*tan(a)+R^2*tan(a)^2-zl^2*tan(a)^2-P^2*tan(a)^2-xl^2+2*R*xl-2*zl*tan(a)^2*P+2*R*zl*tan(a))^(1/2))*tan(a)*xl+R/(tan(a)^2+1)*(2*tan(a)^2*zl+2*tan(a)*xl-2*R*tan(a)+2*tan(a)^2*P+2*(-2*xl*P*tan(a)+2*R*P*tan(a)-2*xl*zl*tan(a)+R^2*tan(a)^2-zl^2*tan(a)^2-P^2*tan(a)^2-xl^2+2*R*xl-2*zl*tan(a)^2*P+2*R*zl*tan(a))^(1/2))*tan(a)-1/(tan(a)^2+1)*(2*tan(a)^2*zl+2*tan(a)*xl-2*R*tan(a)+2*tan(a)^2*P+2*(-2*xl*P*tan(a)+2*R*P*tan(a)-2*xl*zl*tan(a)+R^2*tan(a)^2-zl^2*tan(a)^2-P^2*tan(a)^2-xl^2+2*R*xl-2*zl*tan(a)^2*P+2*R*zl*tan(a))^(1/2))*tan(a)^2*P)/(tan(a)^2+1))^(1/2);
z13(2)=(-(2*zl*tan(a)^2*P-2*R*zl*tan(a)-2*R*P*tan(a)+2*xl*zl*tan(a)+2*xl*P*tan(a)+zl^2*tan(a)^2+xl^2+P^2*tan(a)^2-2*R*xl-1/(tan(a)^2+1)*(2*tan(a)^2*zl+2*tan(a)*xl-2*R*tan(a)+2*tan(a)^2*P-2*(-2*xl*P*tan(a)+2*R*P*tan(a)-2*xl*zl*tan(a)+R^2*tan(a)^2-zl^2*tan(a)^2-P^2*tan(a)^2-xl^2+2*R*xl-2*zl*tan(a)^2*P+2*R*zl*tan(a))^(1/2))*tan(a)^2*zl-1/(tan(a)^2+1)*(2*tan(a)^2*zl+2*tan(a)*xl-2*R*tan(a)+2*tan(a)^2*P-2*(-2*xl*P*tan(a)+2*R*P*tan(a)-2*xl*zl*tan(a)+R^2*tan(a)^2-zl^2*tan(a)^2-P^2*tan(a)^2-xl^2+2*R*xl-2*zl*tan(a)^2*P+2*R*zl*tan(a))^(1/2))*tan(a)*xl+R/(tan(a)^2+1)*(2*tan(a)^2*zl+2*tan(a)*xl-2*R*tan(a)+2*tan(a)^2*P-2*(-2*xl*P*tan(a)+2*R*P*tan(a)-2*xl*zl*tan(a)+R^2*tan(a)^2-zl^2*tan(a)^2-P^2*tan(a)^2-xl^2+2*R*xl-2*zl*tan(a)^2*P+2*R*zl*tan(a))^(1/2))*tan(a)-1/(tan(a)^2+1)*(2*tan(a)^2*zl+2*tan(a)*xl-2*R*tan(a)+2*tan(a)^2*P-2*(-2*xl*P*tan(a)+2*R*P*tan(a)-2*xl*zl*tan(a)+R^2*tan(a)^2-zl^2*tan(a)^2-P^2*tan(a)^2-xl^2+2*R*xl-2*zl*tan(a)^2*P+2*R*zl*tan(a))^(1/2))*tan(a)^2*P)/(tan(a)^2+1))^(1/2);
%and then convert it back to degrees
a=a*180/pi;
z13=min(z13(z13>0));
end

if zl+P<=zi && zu<zi && xl>xi && xu<xi
        fprintf('Condition D');

%momentarily convert b to radians
b=b*pi/180;    
z13(1)=P-(-(zu^2*tan(b)^2+2*R*zu*tan(b)+tan(b)^2*P^2+xu^2-2*tan(b)^2*P*zu+2*tan(b)*P*xu-2*zu*tan(b)*xu-2*R*xu-2*R*tan(b)*P-tan(b)^2*P/(tan(b)^2+1)*(2*tan(b)^2*P+2*tan(b)*xu-2*R*tan(b)-2*tan(b)^2*zu+2*(2*R*tan(b)*P-zu^2*tan(b)^2-tan(b)^2*P^2-xu^2+2*R*xu+2*tan(b)^2*P*zu-2*tan(b)*P*xu+2*zu*tan(b)*xu-2*R*zu*tan(b)+R^2*tan(b)^2)^(1/2))-tan(b)/(tan(b)^2+1)*(2*tan(b)^2*P+2*tan(b)*xu-2*R*tan(b)-2*tan(b)^2*zu+2*(2*R*tan(b)*P-zu^2*tan(b)^2-tan(b)^2*P^2-xu^2+2*R*xu+2*tan(b)^2*P*zu-2*tan(b)*P*xu+2*zu*tan(b)*xu-2*R*zu*tan(b)+R^2*tan(b)^2)^(1/2))*xu+R*tan(b)/(tan(b)^2+1)*(2*tan(b)^2*P+2*tan(b)*xu-2*R*tan(b)-2*tan(b)^2*zu+2*(2*R*tan(b)*P-zu^2*tan(b)^2-tan(b)^2*P^2-xu^2+2*R*xu+2*tan(b)^2*P*zu-2*tan(b)*P*xu+2*zu*tan(b)*xu-2*R*zu*tan(b)+R^2*tan(b)^2)^(1/2))+tan(b)^2/(tan(b)^2+1)*(2*tan(b)^2*P+2*tan(b)*xu-2*R*tan(b)-2*tan(b)^2*zu+2*(2*R*tan(b)*P-zu^2*tan(b)^2-tan(b)^2*P^2-xu^2+2*R*xu+2*tan(b)^2*P*zu-2*tan(b)*P*xu+2*zu*tan(b)*xu-2*R*zu*tan(b)+R^2*tan(b)^2)^(1/2))*zu)/(tan(b)^2+1))^(1/2);
z13(2)=P-(-(zu^2*tan(b)^2+2*R*zu*tan(b)+tan(b)^2*P^2+xu^2-2*tan(b)^2*P*zu+2*tan(b)*P*xu-2*zu*tan(b)*xu-2*R*xu-2*R*tan(b)*P-tan(b)^2*P/(tan(b)^2+1)*(2*tan(b)^2*P+2*tan(b)*xu-2*R*tan(b)-2*tan(b)^2*zu-2*(2*R*tan(b)*P-zu^2*tan(b)^2-tan(b)^2*P^2-xu^2+2*R*xu+2*tan(b)^2*P*zu-2*tan(b)*P*xu+2*zu*tan(b)*xu-2*R*zu*tan(b)+R^2*tan(b)^2)^(1/2))-tan(b)/(tan(b)^2+1)*(2*tan(b)^2*P+2*tan(b)*xu-2*R*tan(b)-2*tan(b)^2*zu-2*(2*R*tan(b)*P-zu^2*tan(b)^2-tan(b)^2*P^2-xu^2+2*R*xu+2*tan(b)^2*P*zu-2*tan(b)*P*xu+2*zu*tan(b)*xu-2*R*zu*tan(b)+R^2*tan(b)^2)^(1/2))*xu+R*tan(b)/(tan(b)^2+1)*(2*tan(b)^2*P+2*tan(b)*xu-2*R*tan(b)-2*tan(b)^2*zu-2*(2*R*tan(b)*P-zu^2*tan(b)^2-tan(b)^2*P^2-xu^2+2*R*xu+2*tan(b)^2*P*zu-2*tan(b)*P*xu+2*zu*tan(b)*xu-2*R*zu*tan(b)+R^2*tan(b)^2)^(1/2))+tan(b)^2/(tan(b)^2+1)*(2*tan(b)^2*P+2*tan(b)*xu-2*R*tan(b)-2*tan(b)^2*zu-2*(2*R*tan(b)*P-zu^2*tan(b)^2-tan(b)^2*P^2-xu^2+2*R*xu+2*tan(b)^2*P*zu-2*tan(b)*P*xu+2*zu*tan(b)*xu-2*R*zu*tan(b)+R^2*tan(b)^2)^(1/2))*zu)/(tan(b)^2+1))^(1/2);
%and then convert it back to degrees
b=b*180/pi;
z13=max(z13(z13>0));
end



if d-Rr+xl-R<=Ri %Eq. 11

    z12=-(to-tf+R*(secd(a)-1))/tand(a); %Eq. 12
    fprintf('I...\n');
    condII=0;
else
    condII=true;
    z12=-sqrt(-to^2+2*to*tf-tf^2+2*R*to-2*R*tf); %Eq. 13
    fprintf('II...\n');

end

%Generate key points for drawing the roller primaries

%if contact condition II contact exists
if condII
    RollerPnts=Roller_NosePnts;
    if d-Rr+xu-R>Ri && d-Rr+xl-R<=Ri
        RollerPnts=[d-Rr+xu-R zu;
            Roller_NosePnts;
            Ri (-Ri+Roller_NosePnts(end,1)+tand(a)*Roller_NosePnts(end,2))/tand(a)];
    end
    if d-Rr+xl-R>Ri && d-Rr+xu-R<=Ri
        RollerPnts=[Ri (Ri-Roller_NosePnts(1,1)+tand(b)*Roller_NosePnts(1,2))/tand(b);
            Roller_NosePnts;
            d-Rr+xl-R zl];
    end
else %condI prevails
    if d-Rr+xu-R<=Ri && d-Rr+xl-R<=Ri
        RollerPnts=[Ri (Ri-Roller_NosePnts(1,1)+tand(b)*Roller_NosePnts(1,2))/tand(b);
        Roller_NosePnts;
        Ri (-Ri+Roller_NosePnts(end,1)+tand(a)*Roller_NosePnts(end,2))/tand(a)];
    end
    if d-Rr+xu-R>Ri && d-Rr+xl-R<=Ri
        RollerPnts=[d-Rr+xu-R zu;
            Roller_NosePnts;
            Ri z12];
    end
    if d-Rr+xl-R>Ri && d-Rr+xu-R<=Ri
        RollerPnts=[(Ri-Roller_NosePnts(1,1)+tand(b)*Roller_NosePnts(1,2))/tand(b);
            Roller_NosePnts;
            d-Rr+xl-R zl];
    end
end

%draw roller primaries
PrevRoll=patch(RollerPnts(:,1),RollerPnts(:,2)+P,[0.6 0.6 0.6]);
hold on
CurrRoll=patch(RollerPnts(:,1),RollerPnts(:,2),[0.7 0.7 0.7]);
if ~Ed
set(CurrRoll,'linewidth',1.2,'Edgecolor','none')
set(PrevRoll,'linewidth',1.2,'Edgecolor','none')
end
axis equal
ax=axis;

plot([Ri Ri],ax(3:4),'k-.','linewidth',1);
CurrRoll=patch(RollerPnts(:,1),RollerPnts(:,2),[0.7 0.7 0.7]);
if ~Ed
set(CurrRoll,'linewidth',1.2,'Edgecolor','none')
end
XLine=linspace(Rm-1,ax(2));
YLine=linspace(ax(3),ax(4));
plot(XLine,z12.*ones(size(XLine)),'k-.')
plot(XLine,z13.*ones(size(XLine)),'k-.')

%plot key lines
for j=1:length(YLine)

QQ(j)=d-Rr-R*cosd(b)+...
    tand(b)*(YLine(j)-R*sind(b)-P);

QQQ(j)=d+tand(a)*(-phi/tand(a)-YLine(j)+P);
end

plot(QQ(QQ>(Rm-1))...
    ,YLine(QQ>(Rm-1))-P,'b--')
plot(QQQ(QQQ>(Rm-1)&QQQ<(Ri+1))...
    ,YLine(QQQ>(Rm-1)&QQQ<(Ri+1)),'b--')

plot(d-Rr+xl-R,zl,'ko','markersize',7,'linewidth',1.1,...
    'markerfacecolor','w')%y=0
plot(d-Rr+xl-R,zl,'k+','markersize',7,'linewidth',1.1) %y=0

plot(d-Rr+xu-R,zu,'ko','markersize',7,'linewidth',1.1,...
    'markerfacecolor','w')%y=0
plot(d-Rr+xu-R,zu,'k+','markersize',7,'linewidth',1.1)%y=0

plot(d-Rr+xi-R,zi,'ko','markersize',7,'linewidth',1.1,...
    'markerfacecolor','w')%y=0
plot(d-Rr+xi-R,zi,'k+','markersize',7,'linewidth',1.1) %y=0

plot(Ri,z13,'ko','markersize',7,'linewidth',1.1,...
    'markerfacecolor','w')%y=0
plot(Ri,z13,'k+','markersize',7,'linewidth',1.1)

plot(Ri,z12,'ko','markersize',7,'linewidth',1.1,...
    'markerfacecolor','w')%y=0
plot(Ri,z12,'k+','markersize',7,'linewidth',1.1)

%Place labels on the plot

%get normalized offset for text labels
xLoff=(ax(2)-ax(1))*0.01;
yLoff=(ax(4)-ax(3))*0.025;
t(1)=text(d-Rr+xl-R+xLoff,zl+yLoff,'x_l,z_l');
t(2)=text(d-Rr+xu-R+xLoff,zu+yLoff,'x_u,z_u');
t(3)=text(d-Rr+xi-R+xLoff,zi+yLoff,'x_i,z_i');
t(4)=text(Ri+xLoff,z13+yLoff,'R_i,z_{13}');
t(5)=text(Ri+xLoff,z12+yLoff,'R_i,z_{12}');

xlabel('x, r','interpreter','none','fontsize',13);
ylabel('z','interpreter','none','fontsize',13)
grid on
box on
plot([Rm Rm],ax(3:4),'k-.','linewidth',1);
retZ=[z12 z13 zl zu];
