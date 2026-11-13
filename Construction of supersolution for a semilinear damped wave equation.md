# Construction of supersolution for a semilinear damped wave equation

One-dimensional damped wave equation, $2\\lt p\\le 3$, negative total mass

## 1. The problem sd

We study mild solutions of $\\partial_t^2u+\\partial_tu-\\partial_x^2u=|u|^p$ on $\\R$, with data $\\ep(u_0,u_1)$:

\\\[u(t)=\\ep L(t,u_0,u_1)+\\int_0^t S(t-\\sigma)|u(\\sigma)|^p\\,d\\sigma,\\qquad L:=(\\pt_t+1)S(t)u_0+S(t)u_1.\\\]

The propagator is

\\\[S(t)f(x)=\\frac{e^{-t/2}}{2}\\int\_{x-t}^{x+t}I_0\\!\\Big(\\tfrac12\\sqrt{t^2-(x-y)^2}\\Big)f(y)\\,dy,\\\]

with a nonnegative kernel $K$. Hence $\\ep L\\le u$ for all $0\\le t\\lt T\_{\\max}$.

Assumptions: $2\\lt p\\le3$, $u_0,u_1\\in L^\\infty$ supported in $\[-R_0,R_0\]$, $R_0\\ge1$.

**Theorem (global supersolution and global existence).** Assume in addition $H_0:=\\int(u_0+u_1)\\,dx\\lt0$, and let $\\Hneg:=\\int\\big((u_0)\_-+(u_1)\_-\\big)dx$. Then there is $\\ep_0\\gt0$ such that for every $0\\lt\\ep\\le\\ep_0$ there are $T_U\\gt0$ and $a\\in C^2(\[0,\\infty))$ with \\\[a(0)=0,\\quad \\dot a\\ge0,\\quad \\ddot a\\ge-2^{-3}\\dot a,\\quad 0\\le a(t)\\lt\\frac{-\\ep H_0}{2R_0}.\\\] With $\\Psi(t,x):=S(t+2^2)\\mathbf 1\_{\[-R_0,R_0\]}(x)$ the function \\\[U:=\\ep L(t,u_0,u_1)+a(t)\\Psi\\\] is a global supersolution and \\\[\\ep L(t,u_0,u_1)\\le u\\le U,\\qquad T\_{\\max}=\\infty.\\\] Moreover $\\|u(t)\\|\_{L^\\infty}\\le C\\ep(1+t)^{-1/2}$.

## 2. Heuristic construction

Fix a time shift $t_0\\ge2^2$ (the theorem uses $t_0=2^2$) and define \\\[\\Psi(t,x):=S(t+t_0)\\mathbf 1\_{\[-R_0,R_0\]}(x)=\\int\_{-R_0}^{R_0}K(t+t_0,x-y)\\,dy.\\\] Since the kernel $K$ is nonnegative, $\\Psi(t,x)\\ge0$ for $t\\ge0$, $x\\in\\R$. By construction $\\Psi$ is a time-shifted solution of the homogeneous damped wave equation: \\\[\\partial_t^2\\Psi+\\partial_t\\Psi-\\partial_x^2\\Psi=0.\\\] Write $\\mathcal L:=\\partial_t^2+\\partial_t-\\partial_x^2$.

### Positivity of the shifted profile

The purpose of the time shift is to guarantee positivity of the combination $\\Psi+2\\partial_t\\Psi$.

**Lemma (positive kernel combination).** Let $K(s,z)=\\frac12e^{-s/2}I_0(\\rho/2)$ with $\\rho=\\rho(s,z)=\\sqrt{s^2-z^2}$, $s\\gt0$, $|z|\\lt s$. Then \\\[K(s,z)+2\\partial_sK(s,z)\\ge\\frac{s}{s+2^2}K(s,z).\\\] *Proof.* Since $\\partial_s\\rho=s/\\rho$ and $I_0'=I_1$, direct differentiation gives $\\partial_sK=-\\frac12K+\\frac12e^{-s/2}I_1(\\rho/2)\\frac{s}{2\\rho}$, hence \\\[K+2\\partial_sK=K\\,\\frac s\\rho\\,\\frac{I_1(\\rho/2)}{I_0(\\rho/2)}.\\\] By Lemma A (Bessel quotient, in the Appendix), \\\[K+2\\partial_sK\\ge K\\,\\frac{s}{2\\big(1+\\sqrt{1+\\rho^2/4}\\big)}.\\\] Since $\\rho\\le s$ and $\\sqrt{1+y^2}\\le1+y$, we have $2\\big(1+\\sqrt{1+\\rho^2/4}\\big)\\le2(2+s/2)=s+2^2$. $\\square$

Integrate over $y\\in\[-R_0,R_0\]$ and take $s=t+t_0$: \\\[\\Psi+2\\partial_t\\Psi\\ge\\frac{t+t_0}{t+t_0+2^2}\\Psi.\\\] More generally this gives $\\Psi+2\\partial_t\\Psi\\ge c\_{t_0}\\Psi$ with $c\_{t_0}:=\\frac{t_0}{t_0+2^2}\\gt0$, for every $t_0\\gt0$. With $t_0\\ge2^2$, $\\frac{t+t_0}{t+t_0+2^2}\\ge\\frac{t_0}{t_0+2^2}\\ge\\frac12$, and therefore \\\[\\Psi+2\\partial_t\\Psi\\ge\\tfrac12\\Psi.\\\]

### Ansatz

We seek a supersolution of the form $U:=\\ep L(t,u_0,u_1)+a(t)\\Psi$, where $a\\in C^2(\[0,\\infty))$, $a\\ge0$, $a(0)=0$. The condition $a(0)=0$ ensures $U(0,x)=\\ep u_0(x)$. Moreover $\\partial_tU(0,x)=\\ep u_1(x)+\\dot a(0)\\Psi(0,x)$, so $\\partial_tU(0,x)\\ge\\ep u_1(x)$ whenever $\\dot a(0)\\ge0$.

Both $L$ and $\\Psi$ solve the homogeneous damped wave equation, so $\\mathcal LU=\\mathcal L(a\\Psi)$. By the product rule, $\\mathcal L(a\\Psi)=\\ddot a\\Psi+2\\dot a\\partial_t\\Psi+\\dot a\\Psi+a\\,\\mathcal L\\Psi=\\ddot a\\,\\Psi+\\dot a(\\Psi+2\\partial_t\\Psi)$. We impose \\\[\\dot a(t)\\ge0,\\qquad\\ddot a(t)\\ge-\\tfrac18\\dot a(t).\\\] Combining with the positivity property $\\Psi+2\\partial_t\\Psi\\ge\\frac12\\Psi$, \\\[(\\pt_t^2+\\pt_t-\\pt_x^2)(a\\Psi)=\\ddot a\\,\\Psi+\\dot a(\\Psi+2\\pt_t\\Psi)\\ge\\big(-\\tfrac18+\\tfrac12\\big)\\dot a\\Psi\\ge\\tfrac14\\dot a\\Psi.\\\] Consequently $\\mathcal LU\\ge\\frac14\\dot a\\Psi$, and the construction reduces to finding a nondecreasing $a$ with $\\ddot a\\ge-\\frac18\\dot a$ such that the pointwise condition \\\[|U|^p\\le\\tfrac14\\,\\dot a\\,\\Psi.\\tag{$\\ast$}\\\] Then $\\mathcal LU\\ge|U|^p$, i.e. $U$ is a genuine supersolution. Since $U(0,x)=\\ep u_0(x)$ and $\\pt_tU(0,x)\\ge\\ep u_1(x)$, the comparison principle (Lemma D), together with the verification lemma (Lemma D$'$), yields $u\\le U$ for all $t\\ge0$, provided $\\ep\\gt0$ is sufficiently small. On the other side, the kernel of $S(t)$ is nonnegative and $|u|^p\\ge0$, so $u-\\ep L=\\int_0^tS(t-\\sigma)|u(\\sigma)|^p\\,d\\sigma\\ge0$, i.e. $\\ep L\\le u$ for $0\\le t\\lt T\_{\\max}$.

### Two regions

Fix a sufficiently large time $T_U\\gt0$. The verification of $(\\ast)$ is divided into two complementary spacetime regions.

- **Step I: short times and the exterior region** ($t\\le T_U$ or $|x|\\ge \\frac{t}{2(R_0+1)}$): we prove $|U|^p\\le2^{-2}A(t)\\Psi$.
- **Step II: large times in the interior region** ($t\\ge T_U$, $|x|\\lt\\frac{t}{2(R_0+1)}$): we prove $|U|^p\\le2^{-2}B(t,d(t))\\Psi$. Here \\\[d(t):=-\\ep H_0-2R_0a(t),\\qquad H_0=\\int\_{\\R}(u_0+u_1)\\,dx\\lt0.\\\] The quantity $d(t)$ measures the remaining net negative mass in the construction. Since $\\dot d(t)=-2R_0\\dot a(t)$, it decreases as $a(t)$ increases.

We then define $a$ by the ODE $\\dot a=A(t)+B(t,d(t))$, $a(0)=0$. Since $A$ and $B$ are nonnegative, $\\dot a\\ge0$, $A\\le\\dot a$ and $B(t,d(t))\\le\\dot a$. Therefore the Step I and Step II estimates imply $(\\ast)$ in their respective regions. Concretely, with \\\[A(t)=C_c\\ep^p(3+t)^p\\exp\\!\\Big(-\\tfrac{(p-1)t}{2^7(R_0+1)^2}\\Big),\\qquad B(t,d)=2^{\\frac{p+3}{2}}\\tfrac{e}{R_0}(t+T_U)^{-\\frac{p-1}{2}}\[d\]\_+^p.\\\] The estimates for $d(t)$ (Section 5) show that the ODE has a global solution with $d(t)\\gt0$ and $0\\le a(t)\\lt\\frac{-\\ep H_0}{2R_0}$ for every $t\\ge0$. Consequently $U$ is globally defined, and $(\\ast)$ holds everywhere. For every $T\\lt T\_{\\max}$ the comparison principle gives \\\[\\ep L(t,u_0,u_1)(x)\\le u(t,x)\\le U(t,x)\\quad\\text{on }\[0,T\]\\times\\R.\\\] The resulting $L^\\infty$ bound is independent of $T\\lt T\_{\\max}$, so the continuation criterion yields $T\_{\\max}=\\infty$.

## 3. Short times and exterior region

Standing assumptions: $T_U\\ge2^2(R_0+1)^2$ and $0\\le a(t)\\lt\\frac{-\\ep H_0}{2R_0}$ for $t\\ge0$. Write $M_0:=\\|u_0\\|\_{L^\\infty}+\\|u_1\\|\_{L^\\infty}$ and \\\[C_c:=2^2\\big\[e^2M_0\\big\]^pR_0^{p-1}\\exp\\!\\Big(\\tfrac{(p-1)T_U}{2^7(R_0+1)^2}\\Big).\\\]

**Lemma.** If $t\\le T_U$ or $|x|\\ge\\frac{t}{2(R_0+1)}$, then $|U|^p\\le 2^{-2}A(t)\\Psi$.

*Proof.* **Step 1 (linear part).** By positivity of the kernels, $|L(t,u_0,u_1)|\\le L(t,|u_0|,|u_1|)$. The linear profile estimate (Appendix) then gives \\\[|L(t,u_0,u_1)(x)|\\le e^2\\Big(\\tfrac52+\\tfrac t8\\Big)M_0\\,\\Psi(t,x).\\\]

**Step 2 (size of $a$).** From $H_0\\lt0$ and the definition of $\\Hneg$, \\\[-H_0=-\\int(u_0+u_1)\\le\\Hneg\\le\\int\_{-R_0}^{R_0}\\big(|u_0|+|u_1|\\big)\\le2R_0M_0.\\\] Hence $0\\le a(t)\\lt\\frac{\\ep\\Hneg}{2R_0}\\le\\ep M_0$, and \\\[|U|\\le\\ep|L|+a\\Psi\\le\\Big\[e^2\\Big(\\tfrac52+\\tfrac t8\\Big)+1\\Big\]\\ep M_0\\Psi\\le e^2(3+t)\\,\\ep M_0\\,\\Psi.\\\]

**Step 3 (the factor $\\Psi^{p-1}$).**

- If $t\\le T_U$: $\\Psi\\le R_0$, and since $T_U-t\\ge0$, \\\[\\Psi^{p-1}\\le R_0^{p-1}\\exp\\!\\Big(\\tfrac{(p-1)(T_U-t)}{2^7(R_0+1)^2}\\Big).\\\]
- If $t\\ge T_U$ and $|x|\\ge\\frac{t}{2(R_0+1)}$: the exterior decay lemma (valid since $t\\ge2^2(R_0+1)^2$) gives $\\Psi\\le R_0\\exp\\!\\big(-\\frac{t}{2^7(R_0+1)^2}\\big)$, so \\\[\\Psi^{p-1}\\le R_0^{p-1}e^{-\\frac{(p-1)t}{2^7(R_0+1)^2}}\\le R_0^{p-1}\\exp\\!\\Big(\\tfrac{(p-1)(T_U-t)}{2^7(R_0+1)^2}\\Big).\\\]

**Step 4 (conclusion).** In both regions, \\begin{align\*}|U|^p&\\le\\big(e^2(3+t)\\ep M_0\\big)^p\\,\\Psi^{p-1}\\,\\Psi\\\\&\\le(e^2M_0)^p\\ep^p(3+t)^pR_0^{p-1}\\,e^{\\frac{(p-1)T_U}{2^7(R_0+1)^2}}\\,e^{-\\frac{(p-1)t}{2^7(R_0+1)^2}}\\,\\Psi\\;=\\;2^{-2}A(t)\\Psi,\\end{align\*} by the definitions of $C_c$ and $A$. $\\square$

## 4. Large-time interior estimate

Here $0\\lt-H_0\\le\\Hneg$, so $\\Hneg\\gt0$. We strengthen the assumption on $T_U$ to $T_U\\ge2^3(R_0+1)^2$ and assume the front condition for all $t\\ge T_U$: \\\[\\frac{t\\,d(t)}{\\ep\\Hneg}\\ge2^5(R_0+1)^2\\Big(\\frac{2^2pt\\log t}{p-1}\\Big)^{1/2}.\\tag{F}\\\]

Write $U=P-N$ with $P=\\ep L(t,(u_0)\_+,(u_1)\_+)+a\\Psi$ and $N=\\ep L(t,(u_0)\_-,(u_1)\_-)$. Their masses are $\\ep\\Hneg-d$ and $\\ep\\Hneg$. With $\\theta=\\frac{(R_0+1)(|x|+2^2(R_0+1))}{t}\\in\[0,1\]$ the product-form bounds give \\\[U\\le K\\big\[e^{\\theta}(\\ep\\Hneg-d)-e^{-\\theta}\\ep\\Hneg\\big\],\\qquad -U\\le K\\big\[e^{\\theta}\\ep\\Hneg-e^{-\\theta}(\\ep\\Hneg-d)\\big\].\\\]

**Lemma.** Under the front condition \\\[\\frac{t\\,d(t)}{\\ep\\Hneg}\\ge2^5(R_0+1)^2\\Big(\\frac{2^2pt\\log t}{p-1}\\Big)^{1/2}\\quad(t\\ge T_U),\\\] one has $|U|^p\\le\\frac14B(t,d(t))\\Psi$ for $|x|\\lt\\frac{t}{2(R_0+1)}$.

*Proof.* Put $\\rho(t):=\\big(\\frac{2^2pt\\log t}{p-1}\\big)^{1/2}$. Since $t\\ge T_U\\ge2^3(R_0+1)^2\\ge32$, we have $\\rho(t)\\ge1$.

**Preliminaries.**

- *$\\theta\\le1$.* For $|x|\\lt\\frac{t}{2(R_0+1)}$, $\\theta(t,x)\\lt\\frac12+\\frac{2^2(R_0+1)^2}{t}\\le1$, because $t\\ge2^3(R_0+1)^2$.
- *Masses.* With $H\_{\\mathrm{pos}}:=\\int\\big((u_0)\_++(u_1)\_+\\big)$ we have $H_0=H\_{\\mathrm{pos}}-\\Hneg$, so $\\ep H\_{\\mathrm{pos}}+2R_0a=\\ep\\Hneg-d$. This is the mass of $P$, and it is $\\ge0$. The mass of $N$ is $\\ep\\Hneg$.
- *Rough bound.* Since $0\\le\\ep\\Hneg-d\\le\\ep\\Hneg$ and $\\theta\\le1$, both right-hand sides in the product-form bounds are at most $e\\,\\ep\\Hneg K$. Hence $|U|\\le e\\,\\ep\\Hneg\\,K(t,x)$.
- *Comparison of $\\Psi$ and $K$.* For these $(t,x)$: $\\Psi\\ge2e^{-1}R_0K$ and $K\\le\\frac12t^{-1/2}e^{-x^2/4t}$.

Two cases:

- **Case 1: $U(t,x)\\gt0$.** The upper product bound forces $e^{-\\theta}\\ep\\Hneg\\lt e^{\\theta}(\\ep\\Hneg-d)$, i.e. \\\[\\frac{d}{\\ep\\Hneg}\\lt1-e^{-2\\theta}\\le2\\theta,\\qquad\\text{so}\\qquad|x|\\gt\\frac{t\\,d}{2(R_0+1)\\ep\\Hneg}-2^2(R_0+1).\\\] By (F) this gives $|x|\\gt2^4(R_0+1)\\rho-2^2(R_0+1)\\gt\\rho(t)$, using $R_0+1\\ge2$ and $\\rho\\ge1$. So $U\\gt0\\Rightarrow|x|\\gt\\rho(t)$. \
  For $|x|\\ge\\rho(t)$ the definition of $\\rho$ gives $\\exp\\!\\big(-\\frac{(p-1)x^2}{4t}\\big)\\le t^{-p}$. Moreover (F) and $\\rho\\ge1$ give $\\frac{t\\,d}{\\ep\\Hneg}\\ge2^5(R_0+1)^2\\gt\\frac e2$, that is $\\frac{e\\ep\\Hneg}{2t}\\le d$. Therefore \\\[\\frac{|U|^p}{\\Psi}\\le\\frac{e}{2^pR_0}t^{-\\frac{p-1}2}e^{-\\frac{(p-1)x^2}{4t}}(e\\ep\\Hneg)^p\\le\\frac{e}{R_0}t^{-\\frac{p-1}2}\\Big(\\frac{e\\ep\\Hneg}{2t}\\Big)^p\\le\\frac e{R_0}t^{-\\frac{p-1}2}d^p.\\\]
- **Case 2: $|x|\\lt\\rho(t)$.** By Case 1, $U\\le0$, so $|U|=U\_-$. The lower product bound gives \\\[U\_-\\le K\\big\[e^{-\\theta}d+(e^{\\theta}-e^{-\\theta})\\ep\\Hneg\\big\]\\le K\\Big\[d+\\frac{2e(R_0+1)(|x|+2^2(R_0+1))}{t}\\ep\\Hneg\\Big\],\\\] using $e^{-\\theta}\\le1$ and $e^\\theta-e^{-\\theta}\\le2e\\theta$ for $0\\le\\theta\\le1$. Since $R_0+1\\ge2$, $\\rho\\ge1$ and $|x|\\lt\\rho$, we get $|x|+2^2(R_0+1)\\le\\frac92(R_0+1)\\rho$, hence \\\[U\_-\\le K\\Big\[d+9e(R_0+1)^2\\frac{\\rho\\,\\ep\\Hneg}{t}\\Big\]\\le2d\\,K,\\\] where the last step uses $9e\\lt2^5$ and (F). Then \\\[\\frac{|U|^p}{\\Psi}\\le\\frac{(2d)^pK^p}{2e^{-1}R_0K}=\\frac{e\\,2^{p-1}}{R_0}d^pK^{p-1}\\le\\frac{e}{R_0}t^{-\\frac{p-1}2}d^p,\\\] because $2^{p-1}K^{p-1}\\le t^{-(p-1)/2}$ by the Gaussian bound.

**Conclusion.** In both cases $|U|^p/\\Psi\\le\\frac{e}{R_0}t^{-(p-1)/2}d(t)^p$. Since $t\\ge T_U$ we have $t+T_U\\le2t$, so $t^{-(p-1)/2}\\le2^{(p-1)/2}(t+T_U)^{-(p-1)/2}$. Therefore \\\[\\frac{|U|^p}{\\Psi}\\le2^{\\frac{p-1}2}\\frac e{R_0}(t+T_U)^{-\\frac{p-1}2}d^p=\\tfrac14B(t,d(t)).\\ \\square\\\]

## 5. Existence of $a$ and lower bound for $d$

We construct $a$ as the solution of $\\dot a=A(t)+B(t,d(t))$, $a(0)=0$, with $d=-\\ep H_0-2R_0a$.

**Lemma.** Assume (T), (FA) and (E1)–(E3) below. Then the ODE has a unique solution $a\\in C^2(\[0,\\infty))$ with $\\dot a\\ge0$, $0\\le a\\lt\\frac{-\\ep H_0}{2R_0}$, $\\ddot a\\ge-\\frac18\\dot a$, the front condition (F) holds for $t\\ge T_U$, and for all $t\\ge0$ \\\[d(t)\\gt d\_\*(t):=\\tfrac12\\Big((-2^{-1}\\ep H_0)^{1-p}+\\Lambda(t+T_U)\\Big)^{-1/(p-1)}\\gt0,\\\] where $\\Lambda(\\sigma)=e\\,2^{(5+p)/2}(p-1)\\int_1^\\sigma\\tau^{-(p-1)/2}d\\tau$.

**Preliminaries on $\\Lambda$ and $d\_\*$.** Explicitly, \\\[\\Lambda(\\sigma)=\\frac{e\\,2^{(7+p)/2}(p-1)}{3-p}\\big(\\sigma^{\\frac{3-p}2}-1\\big)\\ \\ (2\\lt p\\lt3),\\qquad \\Lambda(\\sigma)=e\\,2^5\\log\\sigma\\ \\ (p=3).\\\] $\\Lambda$ is increasing, so $\\Lambda(t+T_U)\\le\\Lambda(2t)$ for $t\\ge T_U$. Since $X+Y\\le2\\max\\{X,Y\\}$ and $2^{-1/(p-1)}\\ge2^{-1}$ for $p\\ge2$, \\\[d\_\*(t)\\ge2^{-2}\\min\\Big(-\\tfrac{\\ep H_0}{2},\\ \\Lambda(t+T_U)^{-1/(p-1)}\\Big).\\tag{D}\\\]

**Conditions on $T_U$.** We impose \\\[T_U\\ge2^3(R_0+1)^2\\tag{T}\\\] \\\[\\frac1{\\Hneg}\\,2^{-2}\\min\\Big(-\\tfrac{H_0}{2},\\ \\Lambda(2t)^{-1/(p-1)}\\Big)\\ge2^5(R_0+1)^2\\Big(\\frac{2^2p\\log t}{(p-1)t}\\Big)^{1/2},\\quad t\\ge T_U.\\tag{FA}\\\] Such a $T_U$ exists. Since $p\\ge2$, $\\tau^{-(p-1)/2}\\le\\tau^{-1/2}$ for $\\tau\\ge1$, so $\\int_1^\\sigma\\tau^{-(p-1)/2}d\\tau\\le2\\sigma^{1/2}$ and \\\[\\Lambda(\\sigma)\\le e\\,2^{(7+p)/2}(p-1)\\sigma^{1/2}.\\\] Hence $\\Lambda(2t)^{-1/(p-1)}\\ge C\\,t^{-1/(2(p-1))}$. After multiplying (FA) by $t$, the left side is bounded below by a positive multiple of $t^{\\frac{2p-3}{2(p-1)}}$ and the right side is a positive multiple of $(t\\log t)^{1/2}$. The exponent gap is positive because \\\[\\tfrac{2p-3}{2(p-1)}-\\tfrac12=\\tfrac{p-2}{2(p-1)}\\gt0.\\\] This is the only place where $p\\gt2$ is used.

**Smallness of $\\ep$.** After fixing $T_U$, we impose \\\[\\ep\\le1,\\tag{E1}\\\] \\\[2R_0\\ep^{p-1}\\int_0^\\infty\\eta(t)\\,dt\\lt\\log2,\\tag{E2}\\\] \\\[p\\Big(2R_0\\ep^{p-1}\\sup\_{t\\ge0}\\eta(t)+e\\,2^{(5+p)/2}(\\ep\\Hneg)^{p-1}\\Big)\\le2^{-4},\\tag{E3}\\\] \\\[p\\Big(2^3R_0\\ep M_0\\Big)^{p-1}\\le2^{-2},\\tag{E4}\\\] where \\\[\\eta(t):=2^3C_c(3+t)^pe^{-\\frac{(p-1)t}{2^7(R_0+1)^2}}\\max\\Big((-H_0)^{-1},\\ \\Lambda(t+T_U)^{\\frac1{p-1}}\\Big).\\\] $\\eta$ does not depend on $\\ep$. It is bounded and integrable, since $\\Lambda(t+T_U)^{1/(p-1)}\\le C(1+t)^{1/(2(p-1))}$ and the exponential factor dominates. Since $p-1\\gt0$, (E1)–(E4) hold for all $0\\lt\\ep\\le\\ep_0$ with $\\ep_0$ small.

*Proof of the lemma.*

**Step 1 (global existence).** Since $p\\gt1$, $z\\mapsto\[z\]\_+^p$ is $C^1$, so the right-hand side of the ODE is continuous in $t$ and locally Lipschitz in $a$. This gives a unique local solution. The right-hand side is $\\ge0$, so $\\dot a\\ge0$ and $a\\ge0$. Then $d=-\\ep H_0-2R_0a\\le-\\ep H_0$, so $0\\le\[d\]\_+\\le-\\ep H_0$. Hence $B(t,d(t))$, like $A(t)$, is bounded on every finite interval, uniformly in $a$. So $\\dot a$ stays bounded on finite intervals, there is no finite-time blow-up, and the solution is global. The right-hand side is $C^1$ in $(t,a)$, so $a\\in C^2$.

**Step 2 (lower bound $d\\gt d\_\*$).** First, $d\_\*(0)\\le\\frac12\\big((-2^{-1}\\ep H_0)^{1-p}\\big)^{-1/(p-1)}=-\\frac{\\ep H_0}{4}\\lt-\\ep H_0=d(0)$.

Suppose $d\\ge d\_\*$ on $\[0,T\]$. Then $d\\gt0$, $\[d\]\_+=d$, and \\\[2R_0B(t,d)=e\\,2^{(5+p)/2}(t+T_U)^{-\\frac{p-1}2}d^p.\\tag{B}\\\] By (D) and $\\ep\\le1$, \\\[\\frac1{d\_\*}\\le2^2\\max\\Big\\{\\big(-\\tfrac{\\ep H_0}2\\big)^{-1},\\Lambda(t+T_U)^{\\frac1{p-1}}\\Big\\}\\le2^3\\ep^{-1}\\max\\Big\\{(-H_0)^{-1},\\Lambda(t+T_U)^{\\frac1{p-1}}\\Big\\}.\\\] By the definitions of $A$ and $\\eta$ this gives $\\frac{A}{d}\\le\\frac{A}{d\_\*}\\le\\ep^{p-1}\\eta$, i.e. \\\[A(t)\\le\\ep^{p-1}\\eta(t)\\,d(t).\\tag{A}\\\] Since $\\dot d=-2R_0(A+B)$, (A) and (B) give \\\[\\dot d\\ge-2R_0\\ep^{p-1}\\eta\\,d-e\\,2^{(5+p)/2}(t+T_U)^{-\\frac{p-1}2}d^p.\\\] Let $y:=d^{1-p}$. Then $\\dot y=-(p-1)d^{-p}\\dot d$, so \\\[\\dot y\\le2(p-1)R_0\\ep^{p-1}\\eta\\,y+(p-1)e\\,2^{(5+p)/2}(t+T_U)^{-\\frac{p-1}2}.\\\] With the integrating factor $\\mu(t)=\\exp\\!\\big(-2(p-1)R_0\\ep^{p-1}\\int_0^t\\eta\\big)\\in(0,1\]$ we get \\\[\\tfrac{d}{dt}(\\mu y)\\le(p-1)e\\,2^{(5+p)/2}\\mu\\,(t+T_U)^{-\\frac{p-1}2}.\\\] Integrate on $\[0,t\]$ and use $\\mu\\le1$, $T_U\\ge1$ (so $\\int_0^t(s+T_U)^{-\\frac{p-1}2}ds\\le\\int_1^{t+T_U}\\sigma^{-\\frac{p-1}2}d\\sigma$) and $y(0)=(-\\ep H_0)^{1-p}\\le(-2^{-1}\\ep H_0)^{1-p}$: \\\[y(t)\\le\\exp\\!\\Big(2(p-1)R_0\\ep^{p-1}\\!\\int_0^t\\eta\\Big)\\Big\[(-2^{-1}\\ep H_0)^{1-p}+\\Lambda(t+T_U)\\Big\].\\\] By (E2) the exponential is $\\lt2^{p-1}$. Since $x\\mapsto x^{-1/(p-1)}$ is strictly decreasing, \\\[d(t)=y^{-\\frac1{p-1}}\\gt\\tfrac12\\Big\[(-2^{-1}\\ep H_0)^{1-p}+\\Lambda(t+T_U)\\Big\]^{-\\frac1{p-1}}=d\_\*(t).\\\] Continuity argument: let $T^\*:=\\sup\\{T\\gt0:\\ d\\ge d\_\*\\text{ on }\[0,T\]\\}$. Then $T^\*\\gt0$ by the strict inequality at $0$. If $T^\*\\lt\\infty$, the estimate above upgrades $d\\ge d\_\*$ on $\[0,T^\*\]$ to $d\\gt d\_\*$ there. By compactness, $d-d\_\*$ has a positive minimum on $\[0,T^\*\]$, so $d\\ge d\_\*$ persists past $T^\*$, a contradiction. Hence $d\\gt d\_\*\\gt0$ for all $t\\ge0$. Finally $d\\gt0$ means $a\\lt\\frac{-\\ep H_0}{2R_0}$.

**Step 3 (estimate for $\\ddot a$).** As $A,B\\gt0$ and $\\ddot a=\\dot A+\\frac{d}{dt}B(t,d(t))$, it suffices to show $\\dot A/A\\ge-\\frac18$ and $\\frac{d}{dt}\\log B\\ge-\\frac18$.

First, $\\frac{\\dot A}{A}=\\frac{p}{t+3}-\\frac{p-1}{2^7(R_0+1)^2}\\ge-\\frac18$, using $p-1\\le2$ and $R_0\\ge1$.

Next, $\\frac{\\dot d}{d}=-\\frac{2R_0(A+B)}{d}$. By (A) and (B), and since $0\\lt d\\le-\\ep H_0\\le\\ep\\Hneg$ and $t+T_U\\ge1$, \\\[\\frac{2R_0(A+B)}{d}\\le2R_0\\ep^{p-1}\\sup_s\\eta(s)+e\\,2^{(5+p)/2}(\\ep\\Hneg)^{p-1}.\\\] Now $\\frac{d}{dt}\\log B=-\\frac{p-1}{2(t+T_U)}+p\\frac{\\dot d}{d}$, and by (T), $\\frac{p-1}{2(t+T_U)}\\le\\frac{p-1}{2T_U}\\le\\frac1{16}$. With (E3), \\\[\\frac{d}{dt}\\log B\\ge-\\tfrac1{16}-p\\Big(2R_0\\ep^{p-1}\\sup\\eta+e\\,2^{\\frac{5+p}2}(\\ep\\Hneg)^{p-1}\\Big)\\ge-\\tfrac18.\\\] Therefore $\\ddot a\\ge-\\frac18A-\\frac18B=-\\frac18\\dot a$.

**Step 4 (front condition).** Let $t\\ge T_U$. By Step 2 and (D), using $t+T_U\\le2t$ and monotonicity of $\\Lambda$, \\\[d(t)\\ge2^{-2}\\min\\Big\\{-\\tfrac{\\ep H_0}2,\\Lambda(2t)^{-\\frac1{p-1}}\\Big\\}\\ge\\ep\\,2^{-2}\\min\\Big\\{-\\tfrac{H_0}2,\\Lambda(2t)^{-\\frac1{p-1}}\\Big\\},\\\] where the last step uses $\\min\\{\\ep X,Y\\}\\ge\\ep\\min\\{X,Y\\}$ for $0\\lt\\ep\\le1$. Multiplying by $t/(\\ep\\Hneg)$ and using (FA) gives exactly (F). $\\square$

## 6. Proof of the main theorem

1. **Step 1 (choice of $T_U$ and $\\ep$).** Choose $T_U$ so that (T) and (FA) hold. Then $C_c$ and $\\eta$ are fixed, independent of $\\ep$. Choose $\\ep_0$ so that (E1)–(E4) hold for $0\\lt\\ep\\le\\ep_0$, and fix $0\\lt\\ep\\le\\ep_0$. By the ODE lemma there is a unique global $a\\in C^2$ with $\\dot a\\ge0$, $\\ddot a\\ge-\\frac18\\dot a$, $d\\gt d\_\*\\gt0$ and $0\\le a\\lt\\frac{-\\ep H_0}{2R_0}$; and (F) holds for $t\\ge T_U$.
2. **Step 2 (supersolution inequality).** Both $L$ and $\\Psi$ solve the homogeneous equation, so \\\[(\\pt_t^2+\\pt_t-\\pt_x^2)U=\\ddot a\\Psi+\\dot a(\\Psi+2\\pt_t\\Psi)\\ge\\big(-\\tfrac18+\\tfrac12\\big)\\dot a\\Psi=\\tfrac38\\dot a\\Psi\\ge\\tfrac14\\dot a\\Psi,\\\] using $\\Psi+2\\pt_t\\Psi\\ge\\frac12\\Psi$ and $\\ddot a\\ge-\\frac18\\dot a$. It remains to show $|U|^p\\le\\frac14\\dot a\\Psi$ everywhere. \
   If $t\\le T_U$ or $|x|\\ge\\frac{t}{2(R_0+1)}$: the standing assumptions of the short-time lemma hold by Step 1, so $|U|^p\\le\\frac14A\\Psi\\le\\frac14\\dot a\\Psi$, because $B\\ge0$. \
   If $t\\ge T_U$ and $|x|\\lt\\frac{t}{2(R_0+1)}$: $d\\gt0$ by Step 1, and $\\ep\\Hneg-d=\\ep\\Hneg+\\ep H_0+2R_0a=\\ep H\_{\\mathrm{pos}}+2R_0a\\ge0$. With (F), all hypotheses of the interior lemma hold, so $|U|^p\\le\\frac14B(t,d)\\Psi\\le\\frac14\\dot a\\Psi$, because $A\\ge0$. \
   The two regions cover $\[0,\\infty)\\times\\R$, hence $\\mathcal LU\\ge|U|^p$ and $U$ is a global supersolution.
3. **Step 3 (initial ordering and comparison).** Since $a(0)=0$ and $L(0)=u_0$, $\\pt_tL(0)=u_1$, \\\[U(0)=\\ep u_0,\\qquad \\pt_tU(0)=\\ep u_1+\\dot a(0)\\Psi(0,\\cdot)\\ge\\ep u_1.\\\] Let $u$ be the maximal mild solution on $\[0,T\_{\\max})$ and fix $0\\lt T\\lt T\_{\\max}$. Positivity of $S(t)$ in the mild formulation gives $\\ep L\\le u$. The verification lemma, together with (E4), then yields $q\\ge-\\frac14$ for the quotient $q$ of the comparison principle. The comparison principle, together with Step 2 and the initial ordering above, gives $u\\le U$ on $\[0,T\]$. As $T\\lt T\_{\\max}$ was arbitrary, \\\[\\ep L(t,u_0,u_1)(x)\\le u(t,x)\\le U(t,x),\\qquad0\\le t\\lt T\_{\\max}.\\\]
4. **Step 4 (global continuation).** Suppose $T\_{\\max}\\lt\\infty$. $L(t,u_0,u_1)$ is bounded in $L^\\infty$ on finite time intervals, $0\\le a\\lt\\frac{-\\ep H_0}{2R_0}$, and $\\Psi$ is bounded on bounded time intervals. Hence $\\sup\_{t\\lt T\_{\\max}}\\|U(t)\\|\_{L^\\infty}\\lt\\infty$. From $\\ep L\\le u\\le U$ we get $|u|\\le\\max\\{|\\ep L|,|U|\\}$, so $\\sup\_{t\\lt T\_{\\max}}\\|u(t)\\|\_{L^\\infty}\\lt\\infty$. This contradicts the blow-up alternative for the maximal mild solution. Therefore $T\_{\\max}=\\infty$, the two-sided bound holds for every $t\\ge0$, and the theorem is proved. $\\square$

## Appendix: tools

- **Lemma A (Bessel quotient).** For every $\\rho\\gt0$, \\\[\\frac{I_1(\\rho/2)}{I_0(\\rho/2)}\\ge\\frac{\\rho}{2\\big(1+\\sqrt{1+\\rho^2/4}\\big)},\\qquad 0\\le1-\\frac{I_1(\\rho/2)}{I_0(\\rho/2)}\\le\\frac2\\rho.\\\]

  *Proof.* Put $x=\\rho/2$ and $r(x)=I_1(x)/I_0(x)$. From $I_0'=I_1$ and $I_1'=I_0-\\frac1xI_1$, \\\[r'=1-\\frac rx-r^2.\\tag{R}\\\] Let $g(x):=\\frac{x}{1+\\sqrt{1+x^2}}$ and $s=\\sqrt{1+x^2}$. A direct computation gives \\\[1-\\frac gx-g^2-g'=\\frac{s-1}{s(s+1)}\\gt0.\\\] Near the origin $r=\\frac x2-\\frac{x^3}{16}+O(x^5)$ and $g=\\frac x2-\\frac{x^3}8+O(x^5)$, so $r\\gt g$ there. If a first contact $r(x_0)=g(x_0)$ occurred, (R) would give $(r-g)'(x_0)=1-\\frac gx-g^2-g'\\gt0$, which is impossible for a first contact from above. Hence $r\\ge g$, which is the first inequality.

  For the second, the integral representations of $I_0,I_1$ give $0\\le r\\le1$. Moreover $g=\\frac{s-1}{x}$ (since $(s-1)(s+1)=x^2$), and $s\\ge x$, so $g\\ge1-\\frac1x$. Thus $0\\le1-r\\le1-g\\le\\frac1x=\\frac2\\rho$. $\\square$
- **Lemma B (product form in the interior region).** Let $R_0\\ge1$, $t\\ge2^2(R_0+1)^2$, $|x|\\le\\frac{t}{2(R_0+1)}$, and let $f\\ge0$ be bounded with support in $\[-R_0,R_0\]$. If $F$ is any of $S(t)f(x)$, $(\\pt_t+1)S(t)f(x)$, $S(t+2^2)f(x)$, then \\\[e^{-\\theta}K(t,x)\\int f\\le F(t,x)\\le e^{\\theta}K(t,x)\\int f,\\qquad\\theta=\\frac{(R_0+1)(|x|+2^2(R_0+1))}{t}.\\\]

  *Proof.* Write $\\rho(t,z)=\\sqrt{t^2-z^2}$, so $K(t,z)=\\frac12e^{-t/2}I_0(\\rho(t,z)/2)$.

  **Step 1 (geometry).** The assumptions give $R_0\\le2^{-4}t$, $R_0^2\\le2^{-2}t$, $|x|\\le2^{-2}t$. So for $|y|\\le R_0$, $|x-y|\\lt\\frac t2$ and $\\rho(t,x),\\rho(t,x-y)\\ge\\frac{\\sqrt3}{2}t$. Also $\[-R_0,R_0\]\\subset(x-t,x+t)$, so the boundary terms in the formula for $(\\pt_t+1)S(t)$ vanish; its kernel is $K+\\partial_tK$.

  **Step 2 (increments of $\\log I_0$).** Put $\\ell(\\rho)=\\log I_0(\\rho/2)$, so $\\ell'=\\frac12\\frac{I_1}{I_0}(\\rho/2)$. For $0\\lt\\rho\\le\\rho'$, Lemma A gives \\\[0\\le\\frac{\\rho'-\\rho}{2}-\\big\[\\ell(\\rho')-\\ell(\\rho)\\big\]=\\frac12\\int\_\\rho^{\\rho'}\\Big(1-\\frac{I_1}{I_0}\\Big)ds\\le\\log\\frac{\\rho'}\\rho\\le\\frac{\\rho'-\\rho}{\\rho},\\\] and, since $0\\le I_1/I_0\\le1$, $|\\ell(\\rho')-\\ell(\\rho)|\\le\\frac{|\\rho'-\\rho|}{2}$.

  **Step 3 (spatial translation, case $F=S(t)f$).** For $|y|\\le R_0$, $\\rho(t,x-y)^2-\\rho(t,x)^2=2xy-y^2$. Hence \\\[\\Big|\\log\\frac{K(t,x-y)}{K(t,x)}\\Big|\\le\\frac{|2xy-y^2|}{2(\\rho(t,x-y)+\\rho(t,x))}\\le\\frac{R_0(|x|+R_0)}{t},\\\] so $e^{-R_0(|x|+R_0)/t}K(t,x)\\le K(t,x-y)\\le e^{R_0(|x|+R_0)/t}K(t,x)$. Integrating against $f$ gives the claim, since $R_0(|x|+R_0)\\le(R_0+1)(|x|+4(R_0+1))$.

  **Step 4 (the operator $(\\pt_t+1)S(t)$).** With $z=x-y$, \\\[\\frac{K+\\partial_tK}{K}=\\frac12+\\frac{t}{2\\rho}\\frac{I_1}{I_0}(\\rho/2).\\\] Lower bound: by Lemma A and $t\\ge\\rho$, this is $\\ge\\frac12+\\frac t{2\\rho}\\big(1-\\frac2\\rho\\big)\\ge1-\\frac t{\\rho^2}\\ge1-\\frac{4}{3t}\\ge\\exp\\!\\big(-\\frac{8}{3t}\\big)$, using $1-s\\ge e^{-2s}$ for $0\\le s\\le\\frac14$. \
  Upper bound: it is $\\le\\frac12+\\frac t{2\\rho}=1+\\frac{z^2}{2\\rho(t+\\rho)}\\le1+\\frac{z^2}{3t^2}\\le1+\\frac{|z|}{6t}\\le\\exp\\!\\big(\\frac{|x|+R_0}{6t}\\big)$. \
  Combine with Step 3 and use $R_0^2+\\frac83\\le2^2(R_0+1)^2$.

  **Step 5 (the time shift $S(t+2^2)$).** For $z=x-y$ let $\\rho=\\rho(t,z)$, $\\rho'=\\rho(t+2^2,z)$. Then $\\log\\frac{K(t+2^2,z)}{K(t,z)}=-2+\\ell(\\rho')-\\ell(\\rho)$ and $(\\rho')^2-\\rho^2=2^3t+2^4$. By Step 1, $0\\le\\frac{\\rho'-\\rho}{\\rho}\\le\\frac8t$, and an elementary computation using $\\rho\\ge\\frac{\\sqrt3}2t$, $|z|\\le\\frac t2$ gives $0\\le\\frac{\\rho'-\\rho}2-2\\le\\frac{2^2z^2}{3t^2}\\le\\frac{2|z|}{3t}\\le\\frac{2(|x|+R_0)}{3t}$. By Step 2, \\\[-\\frac8t\\le\\log\\frac{K(t+2^2,z)}{K(t,z)}\\le\\frac{2(|x|+R_0)}{3t}.\\\] Combine with Step 3 and use $R_0^2+8\\le2^2(R_0+1)^2$. $\\square$
- **Lemma C (Gaussian bound for the kernel).** For $t\\ge1$ and $|x|\\lt t/2$, \\\[K(t,x)\\le\\frac12t^{-1/2}\\exp\\!\\Big(-\\frac{x^2}{4t}\\Big).\\\]

  *Proof.* Set $\\rho=\\sqrt{t^2-x^2}$ and $\\delta=t-\\rho$. The standard bound $I_0(r)\\le e^r/\\sqrt{1+2r}$ gives \\\[K(t,x)\\le\\frac{1}{2\\sqrt{1+\\rho}}e^{-\\delta/2}.\\\] Since $x^2=\\delta(2t-\\delta)$, we have $-\\frac\\delta2=-\\frac{x^2}{4t}-\\frac{\\delta^2}{4t}$. So it suffices to prove \\\[\\sqrt{\\frac{t}{1+\\rho}}\\,e^{-\\delta^2/(4t)}\\le1.\\\] If $0\\le\\delta\\le1$, then $1+\\rho=t+1-\\delta\\ge t$ and we are done. If $\\delta\\gt1$, let $z=(\\delta-1)/t$. Since $|x|\\lt t/2$, $\\rho\\gt\\frac{\\sqrt3}2t$, so $0\\lt z\\lt1-\\frac{\\sqrt3}{2}$, and \\\[\\log\\frac t{1+\\rho}=-\\log(1-z)\\le\\frac{2z}{\\sqrt3}=\\frac{2(\\delta-1)}{\\sqrt3\\,t}\\le\\frac{\\delta^2}{2t},\\\] the last inequality because $\\delta^2-\\frac4{\\sqrt3}\\delta+\\frac4{\\sqrt3}\\gt0$ (negative discriminant). Hence $\\sqrt{t/(1+\\rho)}\\le e^{\\delta^2/(4t)}$. $\\square$
- **Lemma D (comparison principle).** Let $p\\gt1$, $T\\gt0$, $\\mathcal L=\\pt_t^2+\\pt_t-\\pt_x^2$. Let $u,U\\in C(\[0,T\];L^\\infty)$ be mild/distributional solutions of $\\mathcal Lu=|u|^p$ and $\\mathcal LU=|U|^p+h$ with $0\\le h\\in L^1((0,T);L^\\infty)$, and $U(0)\\ge u(0)$, $\\pt_tU(0)\\ge\\pt_tu(0)$. Define \\\[q:=\\frac{|U|^p-|u|^p}{U-u}\\ \\ (U\\ne u),\\qquad q:=0\\ \\ (U=u).\\\] If $q\\ge-\\frac14$, then $u\\le U$ on $\[0,T\]\\times\\R$.

  *Proof.* Put $V=U-u$. Then $\\mathcal LV=h+qV$, equivalently \\\[\\Big(\\pt_t^2+\\pt_t-\\pt_x^2+\\tfrac14\\Big)V=h+\\Big(q+\\tfrac14\\Big)V.\\\] The sine propagator of the shifted operator is $\\widetilde S(t)f(x)=\\frac12e^{-t/2}\\int\_{x-t}^{x+t}f(y)\\,dy$, and both $\\widetilde S(t)$ and $(\\pt_t+1)\\widetilde S(t)$ preserve nonnegativity. By Duhamel, \\\[V=V^{(0)}+\\mathcal KV,\\qquad(\\mathcal Kg)(t):=\\int_0^t\\widetilde S(t-\\sigma)\\Big\[\\Big(q(\\sigma)+\\tfrac14\\Big)g(\\sigma)\\Big\]d\\sigma,\\\] where $V^{(0)}\\ge0$ (free evolution of ordered data, plus the contribution of $h\\ge0$). The operator $\\mathcal K$ is positive, so $V\\ge\\mathcal KV$ and, by induction, $V\\ge\\mathcal K^NV$ for all $N\\ge1$. Since $u,U$ are bounded, so is $q$; let $M_q:=\\|q+\\frac14\\|\_{L^\\infty}$. As $\\|\\widetilde S(\\tau)\\|\_{L^\\infty\\to L^\\infty}\\le T$ on $\[0,T\]$, \\\[\\|\\mathcal K^Ng(t)\\|\_{L^\\infty}\\le\\frac{(M_qTt)^N}{N!}\\sup\_{\[0,T\]}\\|g\\|\_{L^\\infty}.\\\] Applying this to $g=V$ and letting $N\\to\\infty$ gives $V\\ge0$. $\\square$

  **Lemma D$'$ (verification of $q\\ge-\\frac14$).** Let $u\\ge\\ep L(t,u_0,u_1)$ and $U=\\ep L+a\\Psi$ with $a\\ge0$. If $p\\big\[2^3R_0\\ep(\\|u_0\\|\_\\infty+\\|u_1\\|\_\\infty)\\big\]^{p-1}\\le2^{-2}$ (this is (E4)), then $q\\ge-2^{-2}$.

  *Proof.* By positivity, $u\_-,U\_-\\le\\ep L(t,(u_0)\_-,(u_1)\_-)\\le\\ep L(t,|u_0|,|u_1|)$. If $U\\ne u$, the mean value theorem gives $q=p|\\zeta|^{p-2}\\zeta$ for some $\\zeta$ between $u$ and $U$. If $\\zeta\\ge0$ then $q\\ge0$. If $\\zeta\\lt0$, then $|\\zeta|\\le\\ep L(t,|u_0|,|u_1|)\\le2^3R_0\\ep(\\|u_0\\|\_\\infty+\\|u_1\\|\_\\infty)$, by the uniform decay bound $\\|L(t,|u_0|,|u_1|)\\|\_{L^\\infty}\\le2^3R_0(\\|u_0\\|\_\\infty+\\|u_1\\|\_\\infty)(1+t)^{-1/2}$. So $q\\ge-p|\\zeta|^{p-1}\\ge-2^{-2}$. Where $U=u$, $q=0$. $\\square$

Interactive rendering of the LaTeX source. Proofs are condensed to their key steps; the full derivations remain in the source.
