# Projekt 2, del 1 – Finita differensmetoden (facit)

Koden ligger i samma mapp: `U1.m`, `U2.m` (+ `bvp_fdm.m`), `U3.m` (+ `heat_euler.m`), `U4.m` (+ `heat_cn.m`).
Kör respektive `U*.m`-skript; varje skript är uppdelat i sektioner (`%%`) per deluppgift.

---

## U1

### a) Härledning av (6)

Taylorutveckla kring $x_j$:

$$
u_{j+1} = u_j + h u'_j + \tfrac{h^2}{2} u''_j + \tfrac{h^3}{6} u'''_j + \tfrac{h^4}{24} u^{(4)}_j + \dots
$$

$$
u_{j+2} = u_j + 2h u'_j + 2h^2 u''_j + \tfrac{4h^3}{3} u'''_j + \tfrac{2h^4}{3} u^{(4)}_j + \dots
$$

Sök $a,b,c$ så att $a u_j + b u_{j+1} + c u_{j+2} \approx u'_j$. Koefficienterna framför $u_j$, $u'_j$ och $u''_j$ ger

$$
a + b + c = 0,\qquad h(b + 2c) = 1,\qquad \tfrac{h^2}{2}(b + 4c) = 0
\;\Rightarrow\; b = \tfrac{4}{2h},\; c = -\tfrac{1}{2h},\; a = -\tfrac{3}{2h}.
$$

Insättning:

$$
-u_{j+2} + 4u_{j+1} - 3u_j = 2h\,u'_j + \left(4\cdot\tfrac{1}{6} - \tfrac{4}{3}\right)h^3 u'''_j + \left(\tfrac{4}{24}-\tfrac{2}{3}\right)h^4 u^{(4)}_j + \dots
= 2h\,u'_j - \tfrac{2}{3}h^3 u'''_j - \tfrac{1}{2}h^4 u^{(4)}_j + \dots
$$

Dividera med $2h$:

$$
u'(x_j) = \frac{-u_{j+2} + 4u_{j+1} - 3u_j}{2h} + \underbrace{\frac{h^2}{3}u'''(x_j) + \frac{h^3}{4}u^{(4)}(x_j) + \dots}_{O(h^2)}
$$

Den ledande feltermen är alltså $\dfrac{h^2}{3}u'''(x_j)$.

### b) Fel för $u(x)=\sin(e^x)$ vid $x=0.75$

Exakt derivata: $u'(x) = e^x\cos(e^x)$. Vi beräknar $|D_h u - u'(0.75)|$ för $h = 2^{-k}$, $k=1,\dots,8$, och plottar i loglog. Referenskurvan är $h^2$: i en loglog-plot blir $C h^p$ en rät linje med lutning $p$, så om felkurvorna är parallella med $h^2$-linjen avtar felet som $O(h^2)$.

| $k$ | Central (2) | Framåt (6) | Bakåt (7) |
|---|---|---|---|
| 1 | $2.01\cdot10^{-1}$ | $2.33$ | $5.29\cdot10^{-1}$ |
| 4 | $4.96\cdot10^{-3}$ | $8.40\cdot10^{-3}$ | $1.08\cdot10^{-2}$ |
| 8 | $1.95\cdot10^{-5}$ | $3.87\cdot10^{-5}$ | $3.92\cdot10^{-5}$ |

Felet minskar med faktor $\approx 4$ när $h$ halveras, dvs. ordning 2 för alla tre. Den centrala formeln har ungefär halva felet, eftersom dess ledande felterm är $\tfrac{h^2}{6}u'''$ jämfört med $\tfrac{h^2}{3}u'''$ för de ensidiga. För stora $h$ (k = 1, 2) är vi inte i det asymptotiska området än.

---

## U2

### a) Analytisk lösning

$u'' + 16u = d\cos(2x)$. Homogen lösning: $u_h = C_1\cos(4x) + C_2\sin(4x)$.
Partikulärlösning $u_p = A\cos(2x)$: $-4A + 16A = d \Rightarrow A = d/12$.

$$
u(x) = C_1\cos(4x) + C_2\sin(4x) + \tfrac{d}{12}\cos(2x)
$$

Randvillkor:
- $u'(0) = 4C_2 = \gamma \Rightarrow C_2 = \gamma/4$
- $u(2\pi) = C_1 + d/12 = \delta \Rightarrow C_1 = \delta - d/12$

$$
\boxed{u(x) = \left(\delta - \tfrac{d}{12}\right)\cos(4x) + \tfrac{\gamma}{4}\sin(4x) + \tfrac{d}{12}\cos(2x)}
$$

Med $d=5,\ \gamma=\delta=1$: $u(x) = \tfrac{7}{12}\cos 4x + \tfrac14\sin 4x + \tfrac{5}{12}\cos 2x$.

### b) Diskretisering

$h = 2\pi/N$, $x_j = jh$, $f_j = f(x_j)$. I de inre punkterna $j = 1,\dots,N-1$:

$$
\frac{w_{j+1} - 2w_j + w_{j-1}}{h^2} + 16 w_j = f_j .
$$

Vänster rand med (6): $\dfrac{-w_2 + 4w_1 - 3w_0}{2h} = \gamma \Rightarrow w_0 = \dfrac{4w_1 - w_2 - 2h\gamma}{3}$.

Insatt i ekvationen för $j=1$:

$$
\left(-\frac{2}{3h^2} + 16\right)w_1 + \frac{2}{3h^2}w_2 = f_1 + \frac{2\gamma}{3h}.
$$

Höger rand: $w_N = \delta$ flyttas till högerledet i ekvationen för $j = N-1$. Systemet $A\mathbf{w} = \mathbf{q}$ blir

$$
A = \begin{pmatrix}
-\frac{2}{3h^2}+16 & \frac{2}{3h^2} & & & \\
\frac{1}{h^2} & -\frac{2}{h^2}+16 & \frac{1}{h^2} & & \\
& \ddots & \ddots & \ddots & \\
& & \frac{1}{h^2} & -\frac{2}{h^2}+16 & \frac{1}{h^2}\\
& & & \frac{1}{h^2} & -\frac{2}{h^2}+16
\end{pmatrix},
\qquad
\mathbf{q} = \begin{pmatrix} f_1 + \frac{2\gamma}{3h} \\ f_2 \\ \vdots \\ f_{N-2} \\ f_{N-1} - \frac{\delta}{h^2} \end{pmatrix}.
$$

Efter lösning fås $w_0$ ur formeln ovan och $w_N = \delta$.

### c), d)
Se `bvp_fdm.m` och `U2.m`. För $N = 400$ är maxfelet $\approx 2.6\cdot10^{-3}$.

### e) Konvergensstudie

Felet mäts med den normaliserade 2-normen (11) i de inre punkterna. Noggrannhetsordningen skattas ur två på varandra följande $N$ (där $h$ halveras):

$$
p \approx \log_2\frac{\|r\|_{N}}{\|r\|_{2N}} .
$$

| $N$ | $\|r\|_2$ | $p$ |
|---|---|---|
| 50 | $6.72\cdot10^{-2}$ | – |
| 100 | $1.88\cdot10^{-2}$ | 1.84 |
| 200 | $4.94\cdot10^{-3}$ | 1.93 |
| 400 | $1.26\cdot10^{-3}$ | 1.97 |
| 800 | $3.19\cdot10^{-4}$ | 1.99 |

Metoden är av ordning 2, vilket stämmer med att både (3) och (6) är $O(h^2)$.

### f) $f(x) = e^{-(x-\pi)^2}$

Problemet är linjärt, så skillnaden mellan fall (ii) och (i) löser det homogena problemet $v'' + 16v = 0$, $v'(0) = 1$, $v(2\pi) = 0$, dvs.

$$
u_{(ii)} - u_{(i)} = \tfrac14\sin(4x).
$$

Lösning (i) är en ren respons på den lokaliserade källan kring $x=\pi$ (med horisontell tangent i $x=0$), medan (ii) dessutom har en överlagrad svängning $\tfrac14\sin 4x$ som ger lutning 1 i $x=0$. Numeriskt avviker skillnaden från $\tfrac14\sin 4x$ med $\approx 10^{-3}$ för $N=400$ (diskretiseringsfelet).

---

## U3

### a) Diskretisering (centrala differenser + Euler framåt)

$\Delta x = 1/N$, $x_j = j\Delta x$, $t_n = n\Delta t$. För $j = 1,\dots,N-1$:

$$
w_{j,n+1} = w_{j,n} + \Delta t\left( D\,\frac{w_{j+1,n} - 2w_{j,n} + w_{j-1,n}}{\Delta x^2} + f(x_j, t_n)\right),
$$

med $w_{0,n} = T_L(t_n)$, $w_{N,n} = T_R(t_n)$ och $w_{j,0} = u(x_j, 0)$. På matrisform:

$$
\mathbf{w}_{n+1} = (I + \Delta t\,A)\,\mathbf{w}_n + \Delta t\,\mathbf{b}_n,
\qquad
A = \frac{D}{\Delta x^2}\,\mathrm{tridiag}(1,-2,1),
\qquad
\mathbf{b}_n = \mathbf{f}_n + \frac{D}{\Delta x^2}\begin{pmatrix} T_L(t_n)\\0\\\vdots\\0\\T_R(t_n)\end{pmatrix}.
$$

### b) Se `heat_euler.m`.

### c) Analytisk lösning

Variabelseparation $u = X(x)\,G(t)$ ger $G'/G = D X''/X = -\lambda$. Med $X(0) = X(1) = 0$ blir $X_k = \sin(k\pi x)$, $\lambda_k = D k^2\pi^2$. Begynnelsedata $\sin(5\pi x)$ innehåller bara $k=5$:

$$
\boxed{u(x,t) = e^{-25\pi^2 D t}\sin(5\pi x)}
$$

För $N = 100$, $\Delta t = 10^{-5}$, $T = 0.01$ blir maxfelet $1.7\cdot10^{-4}$ ✓.

### d) Stabilitet

Egenvärdena till $A$ är $\lambda_k = -\dfrac{4D}{\Delta x^2}\sin^2\!\left(\dfrac{k\pi}{2N}\right)$, $k = 1,\dots,N-1$. Euler framåt är stabil om $|1 + \Delta t\,\lambda_k| \le 1$ för alla $k$, dvs. $\Delta t\,|\lambda_k| \le 2$. Det största $|\lambda_k|$ är $\approx 4D/\Delta x^2$, vilket ger

$$
\boxed{\Delta t \le \Delta t_{\max} = \frac{\Delta x^2}{2D}}
$$

(exakt: $\Delta x^2 / (2D\cos^2(\pi/2N))$, marginellt större).

Över gränsen får den högfrekventa moden ($k \approx N$) förstärkningsfaktor $< -1$; den växer och byter tecken varje steg, vilket ger en sicksackformad lösning som exploderar. Den startas av avrundningsfel ($\sim10^{-16}$), så tillväxten syns först efter tillräckligt många steg. Test med $M = 1000$ steg ($D=1$, exakt lösning $\approx 0$ vid sluttiden):

| $N$ | $0.99\,\Delta t_{\max}$ | $1.01\,\Delta t_{\max}$ | $1.05\,\Delta t_{\max}$ |
|---|---|---|---|
| 20 | $1.3\cdot10^{-21}$ | $1.5\cdot10^{-13}$ | $6\cdot10^{19}$ |
| 50 | $4.8\cdot10^{-18}$ | $2.2\cdot10^{-9}$ | $6\cdot10^{23}$ |
| 100 | $4.7\cdot10^{-6}$ | $3.7\cdot10^{-6}$ | $3\cdot10^{24}$ |

Teorin stämmer: under gränsen dämpas allt, strax över växer avrundningsfelet (långsamt vid 1 %, explosivt vid 5 %). (För $N=100$ vid $1.01\,\Delta t_{\max}$ har störningen ännu inte hunnit bli större än själva lösningen.) Eftersom $\Delta t_{\max} \propto \Delta x^2$ blir villkoret mycket snävt för fina nät.

---

## U4

### a) Crank–Nicolson

Med samma $A$ och $\mathbf{b}_n$ som i U3a blir (16) på matrisform

$$
\left(I - \frac{\Delta t}{2}A\right)\mathbf{w}_{n+1} = \left(I + \frac{\Delta t}{2}A\right)\mathbf{w}_n + \frac{\Delta t}{2}\left(\mathbf{b}_n + \mathbf{b}_{n+1}\right),
$$

där $\mathbf{b}_n = \mathbf{f}_n + \dfrac{D}{\Delta x^2}\,(T_L(t_n), 0, \dots, 0, T_R(t_n))^\top$. I varje tidssteg löses alltså ett tridiagonalt linjärt system. Randvärdena sätts till $w_{0,n} = T_L(t_n)$, $w_{N,n} = T_R(t_n)$.

### b) Se `heat_cn.m`.

### c) Data (18), $N = 100$, $\Delta t = 10^{-3}$, $T = 0.1$

Maxfel $4.3\cdot10^{-7}$ ✓.

### d) Stabilitet

Förstärkningsfaktorn för mod $k$ är

$$
g_k = \frac{1 + \tfrac{\Delta t}{2}\lambda_k}{1 - \tfrac{\Delta t}{2}\lambda_k},\qquad \lambda_k < 0 \;\Rightarrow\; |g_k| < 1 \text{ för alla } \Delta t > 0,
$$

så metoden är ovillkorligt stabil. Experiment med $N=100$, $T=1$ (Euler skulle kräva $\Delta t \le 5\cdot10^{-5}$):

| $\Delta t$ | 1 | 0.2 | 0.1 | 0.01 |
|---|---|---|---|---|
| maxfel | $4.2\cdot10^{-2}$ | $1.9\cdot10^{-3}$ | $4.5\cdot10^{-4}$ | $2.9\cdot10^{-5}$ |

Ingen instabilitet ens med ett enda tidssteg ($\Delta t = 2\cdot10^4\,\Delta t_{\max}$), och felet avtar ungefär som $\Delta t^2$ så länge tidsfelet dominerar (från $\Delta t = 0.2$ till $0.1$ minskar det med faktor $\approx 4$). Vid $\Delta t = 0.01$ börjar rumsfelet ($\sim\Delta x^2$) dominera: för $N=100$ planar felet ut på $\approx 3\cdot10^{-5}$ även om $\Delta t$ minskas ytterligare. Men för stora $\Delta t$ går $g_k \to -1$ för höga moder, så de dämpas dåligt och kan ge oscillationer (t.ex. för diskontinuerliga begynnelsedata) – stabilitet garanterar inte noggrannhet.

### e) Konvergensstudie med $\Delta t = \Delta x$, $T = 1$

| $N$ | $\|r\|_2$ | $p$ |
|---|---|---|
| 10 | $2.20\cdot10^{-3}$ | – |
| 20 | $5.34\cdot10^{-4}$ | 2.04 |
| 40 | $1.32\cdot10^{-4}$ | 2.02 |
| 80 | $3.27\cdot10^{-5}$ | 2.01 |
| 160 | $8.16\cdot10^{-6}$ | 2.00 |
| 320 | $2.04\cdot10^{-6}$ | 2.00 |

Crank–Nicolson har felet $O(\Delta t^2 + \Delta x^2)$. Med $\Delta t \propto \Delta x$ är båda bidragen $O(\Delta x^2)$, så vi förväntar ordning 2 – vilket stämmer.

### f) Den konstruerade lösningen

$u(x,t) = \sin(e^t x)$ har "vågtal" $e^t$: vid $t = 1, 2, 4$ blir det $2.7,\ 7.4,\ 54.6$, dvs. ungefär $0.4,\ 1.2,\ 8.7$ perioder på staven, och antalet svängningar växer exponentiellt. Amplituden är hela tiden 1.

Värmeledning jämnar normalt ut variationer, så för att upprätthålla en allt mer oscillerande profil krävs en källterm $f \approx e^{2t}\sin(e^t x)$ som växer som $e^{2t}$ och som omväxlande värmer och kyler med allt finare rumslig struktur (samt en högerrand vars temperatur oscillerar allt snabbare, $\sin(e^t)$). Det är inte realistiskt för en riktig stav. Numeriskt kräver det dessutom $\Delta x \ll e^{-t}$ för att lösningen ska kunna representeras.
