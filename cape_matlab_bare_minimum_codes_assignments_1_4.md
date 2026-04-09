# CAPE Laboratory — Bare Minimum MATLAB Codes (Assignments 1–4)

This file contains **bare minimum working MATLAB codes** for all 4 assignments.

---

# Assignment 1 — Nonlinear Equation

## MATLAB Built‑in (fzero)

```matlab
clc; clear

T = 250 + 273;
P = 10;

Tc = 407.5;
Pc = 111.3;
R = 0.08206;

 a = 27*R^2*Tc^2/(64*Pc);
 b = R*Tc/(8*Pc);

f = @(v) (P + a./v.^2).*(v-b) - R*T;

v = fzero(f,1);

fprintf('Volume = %f\n',v)
```

---

# Assignment 2 — CSTR (fsolve)

```matlab
clc; clear

x0 = [1 300 300];

sol = fsolve(@cstr,x0);

disp(sol)

function F = cstr(x)

CA = x(1);
T = x(2);
Tj = x(3);

k0 = 36e6;
E = 12000;
R = 1.987;

r = k0*exp(-E/(R*T))*CA;

F = zeros(3,1);

F(1) = 1*(10-CA) - r;
F(2) = 500*(298-T) + 6500*r -150*(T-Tj);
F(3) = 600*(298-Tj) +150*(T-Tj);

end
```

---

# Assignment 3 — ode45

```matlab
clc; clear

[t,y] = ode45(@model,[0 50],[10 300 300]);

plot(t,y)
legend('CA','T','Tj')

function dydt = model(t,y)

CA = y(1);
T = y(2);
Tj = y(3);

k0 = 36e6;
E = 12000;
R = 1.987;

r = k0*exp(-E/(R*T))*CA;

dydt = zeros(3,1);

dydt(1) = 1*(10-CA) - r;
dydt(2) = 500*(298-T) + 6500*r -150*(T-Tj);
dydt(3) = 600*(298-Tj) +150*(T-Tj);

end
```

---

# Assignment 3 — BVP (bvp4c)

```matlab
clc; clear

x = linspace(0,10,20);
solinit = bvpinit(x,[300 0]);

sol = bvp4c(@ode,@bc,solinit);

plot(sol.x,sol.y(1,:))

function dydx = ode(x,y)

alpha = 0.05;
beta = 2.7e-9;
Tinf = 200;

dydx = [y(2);
alpha*(y(1)-Tinf)+beta*(y(1)^4-Tinf^4)];

end

function res = bc(ya,yb)

res = [ya(1)-300
       yb(1)-400];

end
```

---

# Assignment 4 — Explicit

```matlab
clc; clear

x = linspace(0,1,50);
dx = x(2)-x(1);

alpha = 1;

dt = 0.0001;
r = alpha*dt/dx^2;

T = 350*ones(size(x));
T(1)=300;
T(end)=400;

for t=0:dt:100

Tnew = T;

for i=2:length(x)-1
Tnew(i)=T(i)+r*(T(i+1)-2*T(i)+T(i-1));
end

T=Tnew;

end

plot(x,T)
```

---

# Assignment 4 — PDEPE

```matlab
alpha = [1 10 100];

x = linspace(0,1,50);
t = [1 5 10 50 100];

for a = alpha

sol = pdepe(0,@(x,t,u,DuDx)pde(x,t,u,DuDx,a),@ic,@bc,x,t);

T = sol(:,:,1);

figure
plot(x,T)

end

function [c,f,s] = pde(x,t,u,DuDx,a)

c=1;
f=a*DuDx;
s=0;

end

function u0 = ic(x)

u0 = 350;

end

function [pl,ql,pr,qr] = bc(xl,ul,xr,ur,t)

pl = ul-300;
ql = 0;

pr = ur-400;
qr = 0;

end
```

---

# END

