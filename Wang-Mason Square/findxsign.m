%this program dettermines the position of x1 with respect to the contact
%point. it determines whether it is on the left or on the right. it does
%that by calculating the remainder of the division between the angle and
%pi. if the remainder is between 0 and .25 then the angle is between 0 and
%90 and the the x will be on the right if it is turning clockwise, and so
%on. I used the remainder technique because sometimes the theta is greater
%than pi.

sign =0;
if theta1>0 %turning anticlockwise
    r = rem(theta1,pi)/pi;
    if r >0 && r<0.25 || r >0.5 && r<0.75 %angle between 0 to 90 degrees and 180-270
        sign = -1; %x1 is to the right of contact point
    elseif r >0.25 && r<0.5 || r>0.75 %angle between 90 and 180 and 270-360
        sign = 1; %x1 will be on the left
    else %when it is right on top
        sign = 1;
    end
elseif theta1 <0 %turning clockwise
        r = rem(abs(theta1),pi)/pi;
    if r >0 && r<0.25 || r >0.5 && r<0.75 %angle between 0 to 90 degrees and 180-270
        sign = 1; %x1 is to the left of contact point
    elseif r >0.25 && r<0.5 || r>0.75 %angle between 90 and 180 and 270-360
        sign = -1; %x1 will be on the right
    else %when it is right on top
        sign = 1;
    end
elseif theta1 == 0 %doesn't even rotate
    sign = 1;
end
