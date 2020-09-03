function [xvec,yvec] = ellipse_visual(vec)
x0 = vec(1);
y0 = vec(2);
ang = vec(3);
% Ellipse adds ellipses to the current plot
%
% ellipse_visual(ang,x0,y0) adds an ellipse with semimajor axis of ra,
% a semimajor axis of radius rb, a semimajor axis of ang, centered at
% the point x0,y0.
%
% ra and rb are already pre-established as ra = 0.035 and rb = 0.025
%
% h=ellipse_visual(...) returns the handles to the ellipses.
figure;
hold on 
% Setting up variables
ra = 0.035;
rb = 0.025;
Nb = 1000;
 
% work on the variable sizes
x0=x0(:);
y0=y0(:);
ra=ra(:);
rb=rb(:);
ang=ang(:);
Nb=Nb(:);
if length(ra)~=length(rb),
  error('length(ra)~=length(rb)');
end;
if length(x0)~=length(y0),
  error('length(x0)~=length(y0)');
end;
% how many inscribed elllipses are plotted
if length(ra)~=length(x0)
  maxk=length(ra)*length(x0);
else
  maxk=length(ra);
end;

% Storing values of the ellipse
 xvec = [];
 yvec = [];

% drawing loop
for k=1:maxk
  
  if length(x0)==1
    xpos=x0;
    ypos=y0;
    radm=ra(k);
    radn=rb(k);
    if length(ang)==1
      an=ang;
    else
      an=ang(k);
    end;
  elseif length(ra)==1
    xpos=x0(k);
    ypos=y0(k);
    radm=ra;
    radn=rb;
    an=ang;
  elseif length(x0)==length(ra)
    xpos=x0(k);
    ypos=y0(k);
    radm=ra(k);
    radn=rb(k);
    an=ang(k);
  else
    rada=ra(fix((k-1)/size(x0,1))+1);
    radb=rb(fix((k-1)/size(x0,1))+1);
    an=ang(fix((k-1)/size(x0,1))+1);
    xpos=x0(rem(k-1,size(x0,1))+1);
    ypos=y0(rem(k-1,size(y0,1))+1);
  end;
  
  co=cos(an);
  si=sin(an);
  the=linspace(0,2*pi,Nb(rem(k-1,size(Nb,1))+1,:)+1);
%  x=radm*cos(the)*co-si*radn*sin(the)+xpos;
%  y=radm*cos(the)*si+co*radn*sin(the)+ypos;
  p=line(radm*cos(the)*co-si*radn*sin(the)+xpos,radm*cos(the)*si+co*radn*sin(the)+ypos);
  xvec = [radm*cos(the)*co-si*radn*sin(the)+xpos, xvec];
  yvec = [radm*cos(the)*si+co*radn*sin(the)+ypos, yvec];
  axis equal
  if nargout > 0
    h(k)=p;
  end
  hold on 
  plot(ones(1, 100)*x0, linspace(0, 0.07));
  plot(xvec(yvec == min(yvec)), min(yvec), 'r*');
end;
