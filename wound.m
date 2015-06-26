function wound(Ln)
global CmMax CmMin R Ci Step L Lm Lnl Lc Lo Ln N n0 Dm Do CoTheta Mu CoMin CoNorm
CmMax = 1;  %������������ � ����������� ������ ���� � �����
CmMin = 0.1; %����� ������� ������������ ��������������� ������� �����, ��� �� ��������� ��� ����������� ����������. 
R = 1;
N = 20; %��� ����� �� �������
Ci = 1; 
Step = 1./(R*N);
L = 1; %������������ ������� �� ����������� ������
n0 = 0.1; %�������� ������� � ���������� ������
Lc = 0.7; %0.9 %�������� ���������� ������� � ������ ��������� ����� �������
Lo = 0.02; %�������� ���������� �����
Ln = 1; %������������������ �����
Lm = 10; %�������� ����������� ��������������� ������� ����� �����������,
Lnl = 1; %������� ������ ���������
Do = 0.01;
Dm = 0.1;
CoTheta = 5; %������ ��������, ���� ����� ��������� �'�������� ���������
CoMin = 1; %���� �������� � ���������� ������ ���� �� ������� ����
CoNorm = 5;
Mu = 0.83; %�������� ���������� ����� �������
Nt = 10;
T = 7;
m = 0;
x = 0:Step:R;
t = 0:5*Step:T;

sol = pdepe(m,@pdex4pde,@pdex4ic,@pdex4bc,x,t);
u1 = sol(:,:,1);
u2 = sol(:,:,2);



figure;
surf(x,t,u1);
title('������������ ������� �����');
xlabel(' x');
ylabel(' t');

figure;
surf(x,t,u2);
title('������������ �����');
xlabel('�');
ylabel('t');

Um = (u1(:,1))';
Umt = [Um; t];
Umt1 = prod(Umt);

nt = (L*n0)./(n0+(L-n0)*exp(-Mu*Umt1));

figure;
plot(t, nt, 'b'); hold on;
title('Ncap(x,t)');
title('ٳ������ ������� � ����� ����');
xlabel('t');
ylabel('n');

end



% --------------------------------------------------------------------------

function [c,f,s] = pdex4pde(x,t,u,DuDx)
global Dm Do Lc Lm n0 L Mu Ln CoTheta Lnl Lo CmMin
c = [1; 1];
f = [Dm; Do] .* DuDx;

F1 = -Lc*u(1)./CmMin +Lm*(1-u(2)./CoTheta)-Lnl*u(1);
F2 = (Ln*L*n0)./(n0+(L-n0)*exp(-Mu*u(1).*t)) - Lo*u(2)./CoTheta;
s = [F1; F2];
end
% --------------------------------------------------------------------------
function phiM = t0M(x) %������������ �������������� ������� ����� � ���������� ������ ����
global CmMax CmMin R
phiM = 0.1;
end

function phiOx = t0Ox(x)
global CoMin CoTheta R
phiOx = CoMin + ((CoTheta-CoMin)./R)*x.^2;
end 



function u0 = pdex4ic(x)
u0 = [t0M(x); t0Ox(x)];
end 


function [pl,ql,pr,qr] = pdex4bc(xl,ul,xr,ur,t)
global CoTheta
pl = [1; 0];
ql = [1; 1];
pr = [0; ur(2)-CoTheta];
qr = [1 ; 0];
end

