function zgy()
%% global variables
close all; clear; clc
global n_link L1 L2 L3 th ang R;
global xo yo x1 y1 x2 y2 x3 y3;
global link1 link2 link3;
global p p_curr p_pre ;
global flag dampflag secondflag;
global nb xc yc;
global deleteind ndelete;
%% plot the control buttons and axis box
fig=gcf;
set(gcf,'position',[100 100 800 600]);
txtbox = uicontrol(fig,'style','text','Position',[200, 500, 500,50],'FontName','Times New Roman','FontSize', 12);
string = {'A mini game to eliminate balloons','Right click on any position in the axis, the linkage will move to the clicked position',...
    'If the clicked position is within the balloon, the balloon and the thread fastens it will be eliminated'};
[outstring,newpos]=textwrap(txtbox,string);
set(txtbox,'String',outstring,'Position',newpos);
drawnow;
r0 = uicontrol(fig,'Style','pushbutton','Units','normalized', ...
    'Position',[0.025 .7 0.17 0.05],'string','EXIT', ...
    'Callback','close','visible','on','BackgroundColor',[0.8 0.8 0.8]);
bg = uibuttongroup(fig,'Visible','on','Position',[0 0.3 0.2 0.3],'SelectionChangeFcn',@bselection);
r1 = uicontrol(bg,'Style','radiobutton','String','PINV','Position',[10 110 120 30], 'HandleVisibility','off');
r2 = uicontrol(bg,'Style','radiobutton','String','Jacobian transpose','Position',[10 70 120 30], 'HandleVisibility','off');
r3 = uicontrol(bg,'Style','radiobutton','String','CCD','Position',[10 30 120 30],'HandleVisibility','off');
set(bg,'SelectedObject',r3);
r4 = uicontrol(fig,'Style','checkbox','String','Damping','Position',[10 120 120 30],'HandleVisibility','off','callback',@dampcallback);
r5 = uicontrol(fig,'Style','checkbox','String','Secondary task','Position',[10 70 120 30],'HandleVisibility','off','callback',@secondcallback);
gca = axes('Position',[0.25,0.1,0.7,0.7]);
axis([-12, 12, -12, 12]);
axis equal;
set(gcf,'WindowButtonDownFcn',@clickcallback);
%% basic parameters
% the length of the three links L1,L2,L3, three angles th(1),th(2),th(3)
L1 = 3; L2 = 4; L3 = 5; n_link = 3; xo = 0; yo = 0; R = 0.5;
error = 0.1; %the position error of the end effector
lamda = 0.2; %the damping parameter
flag = 3;%choose method, 1 means PINV, 2 means Jacobian Transpose, 3 means CCD
% the positions of the four points
th = zeros(3,1);
getposition();
p_pre = p(:,4); %the previous position of the end effector
%% plot the links
scatter(p(1,1),p(2,1),3,'filled');%plot the original point
hold on;
d1 = fill(link1(1,:),link1(2,:),'b','erasemode','xor');
axis equal;hold on
d2 = fill(link2(1,:),link2(2,:),'g','erasemode','xor');
axis equal;hold on
d3 = fill(link3(1,:),link3(2,:),'r','erasemode','xor');
axis equal;hold on
%% draw heart balloons
nb = 9;
ndelete = 0;
deleteind = zeros(1,9);
%xc = [-6, 0, 6, -6, 0, 6, -6, 0, 6];
%yc = [4, 4, 4, 6, 6, 6, 8, 8, 8];
xc = [-2, 2, -4, 0, 4, -6, -2, 2, 6];
yc = [4, 4, 6, 6, 6, 8, 8, 8, 8];
ori = [0,-5];
t = [0:pi/12:2*pi]';
%{
minor = [0.5, 0.5, 0.6, 0.5, 0.6, 0.5, 0.6, 0.5, 0.6];
major = [0.8, 0.7, 0.9, 0.7, 0.9, 0.8, 0.8, 0.9, 0.7];
xb = cos(t)*minor+ones(25,1)*xc;
yb = sin(t)*major+ones(25,1)*yc;
ab = [0.3, 0.3, 0.4, 0.3, 0.4, 0.3, 0.4, 0.3, 0.4];
%}
ab = rand(1,9)*0.1+0.2;
xb=(2*sin(t)-sin(t*2))*ab+ones(25,1)*xc;
yb=(2*cos(t)-cos(t*2))*ab+ones(25,1)*yc;
b1 = fill(xb(:,1),yb(:,1),rand(1,3));hold on
b2 = fill(xb(:,2),yb(:,2),rand(1,3));hold on
b3 = fill(xb(:,3),yb(:,3),rand(1,3));hold on
b4 = fill(xb(:,4),yb(:,4),rand(1,3));hold on
b5 = fill(xb(:,5),yb(:,5),rand(1,3));hold on
b6 = fill(xb(:,6),yb(:,6),rand(1,3));hold on
b7 = fill(xb(:,7),yb(:,7),rand(1,3));hold on
b8 = fill(xb(:,8),yb(:,8),rand(1,3));hold on
b9 = fill(xb(:,9),yb(:,9),rand(1,3));hold on

thread1 = plot([ori(1,1), xc(1,1)],[ori(1,2),yc(1,1)-3*ab(1,1)],'Color',[0.5,0.5,0.5], 'LineWidth',1);hold on
thread2 = plot([ori(1,1), xc(1,2)],[ori(1,2),yc(1,2)-3*ab(1,2)],'Color',[0.5,0.5,0.5], 'LineWidth',1);hold on
thread3 = plot([ori(1,1), xc(1,3)],[ori(1,2),yc(1,3)-3*ab(1,3)],'Color',[0.5,0.5,0.5], 'LineWidth',1);hold on
thread4 = plot([ori(1,1), xc(1,4)],[ori(1,2),yc(1,4)-3*ab(1,4)],'Color',[0.5,0.5,0.5], 'LineWidth',1);hold on
thread5 = plot([ori(1,1), xc(1,5)],[ori(1,2),yc(1,5)-3*ab(1,5)],'Color',[0.5,0.5,0.5], 'LineWidth',1);hold on
thread6 = plot([ori(1,1), xc(1,6)],[ori(1,2),yc(1,6)-3*ab(1,6)],'Color',[0.5,0.5,0.5], 'LineWidth',1);hold on
thread7 = plot([ori(1,1), xc(1,7)],[ori(1,2),yc(1,7)-3*ab(1,7)],'Color',[0.5,0.5,0.5], 'LineWidth',1);hold on
thread8 = plot([ori(1,1), xc(1,8)],[ori(1,2),yc(1,8)-3*ab(1,8)],'Color',[0.5,0.5,0.5], 'LineWidth',1);hold on
thread9 = plot([ori(1,1), xc(1,9)],[ori(1,2),yc(1,9)-3*ab(1,9)],'Color',[0.5,0.5,0.5], 'LineWidth',1);hold on
drawnow;
axis equal
axis([-12, 12, -12, 12]);
    %% call back functions
    % method selection callback
    function bselection(source,eventdata)
        str = get(eventdata.NewValue,'String');
        if str(1) == 'P'
            flag = 1;
        elseif str(1) == 'J'
            flag = 2;
        else
            flag = 3;
        end
    end
    % mouse click call back
    function [xt, yt] = clickcallback(hObject, eventdata)
        crp=get(gca,'CurrentPoint');
        if crp(1,1)^2+crp(1,2)^2>(L1+L2+L3)^2
            xt = (L1+L2+L3)*crp(1,1)/sqrt(crp(1,1)^2+crp(1,2)^2);
            yt = (L1+L2+L3)*crp(1,2)/sqrt(crp(1,1)^2+crp(1,2)^2);
        else
            xt = crp(1,1);
            yt = crp(1,2);
        end
        moving(xt,yt);
    end
    % whether to use damping
    function dampcallback(hObject, eventdata, handles)
        if get(hObject,'Value') == get(hObject,'Max')
            dampflag = 1;
        else
            dampflag = 0;
        end
    end
    % whether to use secondary task
    function secondcallback(hObject, eventdata, handles)
        if get(hObject,'Value') == get(hObject,'Max')
            secondflag = 1;
        else
            secondflag = 0;
        end
    end
    %% calculation in moving
    function moving(xt, yt)
        i=1;
        while (sqrt( (xt-p_curr(1))^2+(yt-p_curr(2))^2 )>error)
            %% inverse kinematics
            if (flag==1)
                %discrete expected positions
                m = 2000;
                dx = [(xt-p_curr(1))/(m-i); (yt-p_curr(2))/(m-i)];
                J = [ -L1*sin(th(1)) - L2*sin(th(1)+th(2)) - L3*sin(th(1)+th(2)+th(3)),   -L2*sin(th(1)+th(2)) - L3*sin(th(1)+th(2)+th(3)),    -L3*sin(th(1)+th(2)+th(3));
                    L1*cos(th(1)) + L2*cos(th(1)+th(2)) + L3*cos(th(1)+th(2)+th(3)),    L2*cos(th(1)+th(2)) + L3*cos(th(1)+th(2)+th(3)),     L3*cos(th(1)+th(2)+th(3))];
                if dampflag == 1
                    INVJ = pinv(J'*J+lamda^2*eye(3))*J';
                else
                    INVJ = pinv(J); %damping, pseudo inverse
                end
                if secondflag == 1
                    dthinput = [0; pi/360; 0]; %secondary task, deta theta input
                    dth = INVJ*dx+(eye(3)-INVJ*J)*dthinput;  %secondary task
                else
                    dth = INVJ*dx;
                end
                th = th + dth; %update joint rotations
                getposition();
                drawmove();
                i=i+1;
                if i>=2000
                    break;
                end
            %% Jacobian Transpose
            elseif (flag==2)
                %discrete expected positions
                dx = [xt-p_curr(1); yt-p_curr(2)];
                J = [ -L1*sin(th(1)) - L2*sin(th(1)+th(2)) - L3*sin(th(1)+th(2)+th(3)),   -L2*sin(th(1)+th(2)) - L3*sin(th(1)+th(2)+th(3)),    -L3*sin(th(1)+th(2)+th(3));
                    L1*cos(th(1)) + L2*cos(th(1)+th(2)) + L3*cos(th(1)+th(2)+th(3)),    L2*cos(th(1)+th(2)) + L3*cos(th(1)+th(2)+th(3)),     L3*cos(th(1)+th(2)+th(3))];
                tt = J*J'*dx;
                alpha = dx'*tt/(tt'*tt);
                if dampflag == 1
                    INVJ = pinv(J'*J+lamda^2*eye(3))*J';
                else
                    INVJ = pinv(J); %damping, pseudo inverse
                end
                if secondflag == 1
                    dthinput = [0; pi/360; 0]; %secondary task, deta theta input
                    dth = alpha*J'*dx+(eye(3)-INVJ*J)*dthinput;  %secondary task
                else
                    dth = alpha*J'*dx;
                end
                th = th + dth;
                getposition();
                drawmove();
            %% CCD
            elseif (flag==3)
                iteration = n_link;
                while (iteration > 0)
                    % e is the point to be rotated, c is the rotate center
                    pe = p(:,n_link+1);
                    pc = p(:,iteration);
                    pt = [xt;yt];
                    %calculate rotation angle
                    a = (pe - pc)/norm(pe-pc);
                    b = (pt - pc)/norm(pt-pc);
                    dth = acos(dot(a, b));
                    %determine whether the angle is positive or negative
                    direction = cross([a(1) a(2) 0],[b(1) b(2) 0]);
                    if direction(3) < 0
                        dth = -dth;
                    end
                    for iter = iteration:n_link
                        th(iter) = th(iter) + dth;
                    end
                    %rotate
                    iteration = iteration-1;
                    getposition();
                    drawmove();
                end
            end
        end
        %% final position of one move
        set (d1,'xdata',link1(1,:),'ydata',link1(2,:));
        set (d2,'xdata',link2(1,:),'ydata',link2(2,:));
        set (d3,'xdata',link3(1,:),'ydata',link3(2,:));
        % determine which balloon is eliminated
        if (ndelete < nb)
            for ib = 1:nb
                dx = p_curr(1)-xc(1,ib);
                dy = p_curr(2)-yc(1,ib);
                if ( deleteind(ib)==0 && (dx^2+dy^2-ab(1,ib)^2)^2-4*ab(1,ib)^2*((dy-ab(1,ib))^2+dx^2)<=0 )
                %if ( deleteind(ib)==0 && (dx)^2/minor(1,ib)^2+(dy)^2/major(1,ib)^2<=1 )
                    bib = eval(strcat('b',num2str(ib)));
                    bthread = eval(strcat('thread',num2str(ib)));
                    delete(bib);
                    delete(bthread);
                    deleteind(ib)=1;
                    ndelete = ndelete+1;
                    break;
                end
            end
        end
        drawnow;
        axis equal
        axis([-12, 12, -12, 12]);
    end
    %% get position, angle and links' coordinates
    function getposition()
        % get new position and new angle
        ang = [th(1);th(1)+th(2);th(1)+th(2)+th(3)];
        x1 = xo+L1*cos(ang(1));
        y1 = yo+L1*sin(ang(1));
        x2 = xo+L1*cos(ang(1))+L2*cos(ang(2));
        y2 = yo+L1*sin(ang(1))+L2*sin(ang(2));
        x3 = xo+L1*cos(ang(1))+L2*cos(ang(2))+L3*cos(ang(3));
        y3 = yo+L1*sin(ang(1))+L2*sin(ang(2))+L3*sin(ang(3));
        p = [xo, x1, x2, x3;
            yo, y1, y2, y3];
        p_curr= p(:,4);
        
        % get the position of links
        leftangle = linspace( (ang(1)+pi/2), (ang(1)+3*pi/2), 12);
        rightangle = fliplr(linspace( (ang(1)+pi/2), (ang(1)+3*pi/2), 12));
        leftarc = [ ones(1,12)*xo+R*cos(leftangle); ones(1,12)*yo+R*sin(leftangle) ];
        rightacr = [ ones(1,12)*(x1)+R*cos(rightangle); ones(1,12)*(y1)+R*sin(rightangle) ];
        upline = [fliplr( linspace( xo+R*cos(ang(1)+pi/2), x1+R*cos(ang(1)+pi/2), 24) );
            fliplr( linspace( yo+R*sin(ang(1)+pi/2), y1+R*sin(ang(1)+pi/2), 24) )];
        bottomline = [linspace( xo+R*cos(ang(1)-pi/2), x1+R*cos(ang(1)-pi/2), 24);
            linspace( yo+R*sin(ang(1)-pi/2), y1+R*sin(ang(1)-pi/2), 24)];
        link1 = [leftarc,bottomline,rightacr,upline];
        
        leftangle = linspace( (ang(2)+pi/2), (ang(2)+3*pi/2), 12);
        rightangle = fliplr(linspace( (ang(1)+pi/2), (ang(1)+3*pi/2), 12));
        leftarc = [ ones(1,12)*x1+R*cos(leftangle); ones(1,12)*y1+R*sin(leftangle) ];
        rightacr = [ ones(1,12)*(x2)+R*cos(rightangle); ones(1,12)*(y2)+R*sin(rightangle) ];
        upline = [fliplr( linspace( x1+R*cos(ang(2)+pi/2), x2+R*cos(ang(2)+pi/2), 32) );
            fliplr( linspace( y1+R*sin(ang(2)+pi/2), y2+R*sin(ang(2)+pi/2), 32) )];
        bottomline = [linspace( x1+R*cos(ang(2)-pi/2), x2+R*cos(ang(2)-pi/2), 32) ;
            linspace( y1+R*sin(ang(2)-pi/2), y2+R*sin(ang(2)-pi/2), 32) ];
        link2 = [leftarc,bottomline,rightacr,upline];
        
        leftangle = linspace( (ang(3)+pi/2), (ang(3)+3*pi/2), 12);
        leftarc = [ ones(1,12)*x2+R*cos(leftangle); ones(1,12)*y2+R*sin(leftangle) ];
        upline = [linspace( x2+R*cos(ang(3)+pi/2), x3, 40);
            linspace( y2+R*sin(ang(3)+pi/2), y3, 40)];
        bottomline = [linspace( x2+R*cos(ang(3)-pi/2), x3, 40);
            linspace( y2+R*sin(ang(3)-pi/2), y3, 40)];
        link3 = [leftarc,bottomline,upline];
    end
    %% move the links
    function drawmove()
        % draw links
        if (p_curr(1)-p_pre(1))^2+(p_curr(2)-p_pre(2))^2>0.5
            set (d1,'xdata',link1(1,:),'ydata',link1(2,:));
            set (d2,'xdata',link2(1,:),'ydata',link2(2,:));
            set (d3,'xdata',link3(1,:),'ydata',link3(2,:));
            drawnow;
            axis equal
            axis([-12, 12, -12, 12]);
            p_pre = p_curr;
        end
    end
end