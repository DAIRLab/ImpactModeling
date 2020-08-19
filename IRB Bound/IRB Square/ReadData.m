

%ReadData
%Junior Team

%This code takes as an input a csv file that contains information about the
%object's position and velocity in the plane. The code will identify the
%regions of the data where an impact has occurred, and will output the
%pre-impact and post-impact set position and velocity of the object as well
%formulate tangential (d) and normal (n) jacobian vectors for the impact. 


function [Pre, Post, d, n] = ReadData(i, impact)
%load csv data
D = readtable(['traj_' num2str(i) '.csv']);
%go from table to matrix
D = D{:, :};

%Column 1: x position 
x = D(:,1);
%Column 2: z position
z = D(:,2);
%Column 3: theta
th = D(:,3);
%Column 4: xDot
xDot = D(:,4);
%Column 5: zDot
zDot = D(:,5);
%Column 6: thetaDot
thDot = D(:,6);

l = length(zDot);

for jj = 1:l
    
     the = th(jj);
     cc1 = -x(1);
     cc2 = z(1);
     ax = 0.03;
     ay = 0.03;
     bx = 0.03;
     by = - 0.03;
     cx = - 0.03;
     cy = - 0.03;
     dx = - 0.03;
     dy = 0.03;

     a = [ax;ay];
     b = [bx;by];
     c = [cx;cy];
     d = [dx;dy];
     
     vx = [a(1),b(1),c(1),d(1),a(1)];
     vy = [a(2),b(2),c(2),d(2),a(2)];

     R = [cos(the),sin(the);-sin(the),cos(the)];

     a2 = R*a;
     b2 = R*b;
     c2 = R*c;
     d2 = R*d;

     vx = [a2(1),b2(1),c2(1),d2(1),a2(1)]+cc1;
     vy = [a2(2),b2(2),c2(2),d2(2),a2(2)]+cc2;
    
end

%To identify the positions at which an impact occurred, we will take a look
%at the changes in zDot over time - in essence, we will try to
%find the parts where zDot changes sign. When an object hits the ground
%from above, it velocity should reverse. Hence, we implement the
%following code:

imp = [];

%This for loop finds a vector imp, which tracks all the times zDot changes
%sign
for k = 2:l
    if sign(zDot(k)) ~= sign (zDot(k-1))
        imp = [k, imp];
    end  
end

l2 = length(imp);
imp = imp - 1;
impR = [];

%This for loop compares the previous imp vector against the z positions, to
%find only the values at which z is at a local minima
for h = 1:l2
    a = imp(h);
    if z(a) < z(a+1) && z(a) < z(a-1) 
        impR = [a,impR];    
    end
end

%output: two matrices that contains all the information for the pre and
%post impacts
impR2 = impR + 1;
miny = [];
minx = [];

for hh = 1:length(impR)
     
     i = impR(hh);
     the = th(i);
     c1 = -x(i);
     c2 = z(i);
     ax = 0.03;
     ay = 0.03;
     bx = 0.03;
     by = - 0.03;
     cx = - 0.03;
     cy = - 0.03;
     dx = - 0.03;
     dy = 0.03;

     a = [ax;ay];
     b = [bx;by];
     c = [cx;cy];
     d = [dx;dy];
     
     vx = [a(1),b(1),c(1),d(1),a(1)];
     vy = [a(2),b(2),c(2),d(2),a(2)];

     R = [cos(the),sin(the);-sin(the),cos(the)];

     a_2 = R*a;
     b_2 = R*b;
     c_2 = R*c;
     d_2 = R*d;

     vx = [a_2(1),b_2(1),c_2(1),d_2(1),a_2(1)]+c1;
     vy = [a_2(2),b_2(2),c_2(2),d_2(2),a_2(2)]+c2;
     
     %find contact position of square
     miny = [miny, min(vy)];
     aa = find(vy == min(vy));
     minx = [minx, vx(min(aa))];
     
end

Pre = [x(impR(impact)), z(impR(impact)), th(impR(impact))];
Post = [x(impR2(impact)), z(impR2(impact)), th(impR2(impact))];

%using line fits and derivatives, find the pre post velocities
frame = impR(impact); %data index of frame
range = frame-6:frame+6;
rate = 1/250;
time = (frame-6)*rate:rate:(frame+6)*rate;    %time vector
for var = 1:3
    position = D(range, var);                    %position vector
    %linearly fit  pre impact data
    if var == 2
        %For y position, we can't just use a simple linear fit. Instead, we
        %subtract out 1/2 g * (t - t_0)^2 to account for the parobolic
        %freelfall and then perform a linear fit
        position = position + 9.8 * 0.5 * (time(:)- ones(length(time),1)*time(1)).^2;  
        
    end 
    
    
    p1 = polyfit(time(1:7)', position(1:7), 1);   
    %add velocity to pre vector
    Pre(3+var)= p1(1);
    %linearly fit post impact data
    p2 = polyfit(time(8:end)', position(8:end), 1);
    %add velocity to post vector
    Post(3+var) = p2(1);   
    
    if var == 2 
        Pre(3 + var) = Pre(3 + var) - 9.8*(time(7)- time(1));
        Post(3 + var) = Post(3 + var) - 9.8*(time(8) - time(1));
    end 
    
end 

%need to transform contact position from global coordinates 
d = [1, 0, z(impR(impact))- miny(impact)]; 
n = [0, 1, minx(impact) + x(impR(impact))];



end
