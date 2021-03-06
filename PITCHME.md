## scipy.integrate.odeint

Hibiki Mizuno 

---

## scipy.integrate.odeint

```python
from scipy.integrate import odeint

def f(uv, t, a, b, c, I): # u, v:v[0], y'=z:v[1]
    dudt = c*(-uv[1] + uv[0] - uv[0]**3 /3 + I)
    dvdt = uv[0] - b*uv[1] + a
    return [dudt, dvdt]

if __name__ == "__main__":
    t_list = np.linspace(0.0, 3000.0, 20000)
    a = 0.7;b = 0.8;c = 10;
    var_init = [1, 2]
    for i in range(3):
        I = 1.32 + i/50 + 0.07
        var_list = odeint(f, var_init, t_list, args=(a, b, c, I))
        print(var_list)
        plt.plot(var_list[:, 0], var_list[:, 1], label="I={}".format(round(I,3)))

```




---

## scipy.integrate.odeint

https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html

>Solve a system of ordinary differential equations using lsoda from the FORTRAN library odepack.

- FORTRANで書かれたODEPACKのLSODAを使って計算
だそうです。　

---
## odepack/LSODA
https://people.sc.fsu.edu/~jburkardt/f77_src/odepack/odepack.html
>LSODA automatically switches between nonstiff and stiff solvers depending on the behavior of the problem;

---

## odepack/LSODA
https://people.sc.fsu.edu/~jburkardt/f77_src/odepack/odepack.html
>LSODA, written jointly with L. R. Petzold, solves systems dy/dt = f with a dense or banded Jacobian when the problem is stiff, but it automatically selects between nonstiff (Adams) and stiff (BDF) methods. It uses the nonstiff method initially, and dynamically monitors data in order to decide which method to use.

---

## odepack/LSODA
https://people.sc.fsu.edu/~jburkardt/f77_src/odepack/odepack.html

- Adams法とBDFを動的に切り替える


---
## odepack/LSODA

- Adams法とは<<<
	- 線形多段法

- BDFとは
	- 後退微分法

---


## Adams法
- オイラー法、ルンゲクッタ法などは一段法
	- 一つ手前の分点での$ x(n), f(x(n))$の値から次の分点の値を求める
- 線形多段法
	- アダムス・バッシュフォース法 adams-bashforth <<<<
	- アダムス・ムルトン法 adams-moulton

---

## Adams-Bashforth法

$$
\mbox{常微分方程式}\ \ \frac{dy}{dx} = f(x, y)
$$

は

`\[y(x_{k+1}) = y(x_k) + \int_{x_k}^{x_{k+1}} f(x, y(x)) dx\]`


と同値である。したがって、知っている情報を用いて、右側の積分を近似的に表現することにより、次のステップの$y$に対する近似値 $y_{k+1}$を求めることができる.

---

### Lagrangeの補間多項式

予備知識。

$(n+1)$ 個の異なる $x$ 座標 $x_0, x_1, ..., x_n$に対して $y$座標$y_0, y_1, ..., y_n$が与えられた時, これらすべての点を通過する$n$次の多項式$g_n(x)$, すなわち


`\[g_n(x) = a_0 + a_1 x + a_1 x^2 + \cdots + a_n x^n\ \ \ \ (a_n \not =0)\\
g_n(x_i) = y_i\ \ \ \ (i = 0,1,...,n)
\]`

を満たすような$g_n(x)$がただ一つ存在する.

-----

`\[
\begin{align}
L_i(x) &= \prod^n_{j=0\\j\not=i} \frac{(x - x_j)}{(x_i - x_j)}\\
&= \frac{(x - x_0)\cdots(x - x_{i-1})(x - x_{i+1})\cdots (x- x_{n})}{(x_i - x_0)\cdots (x_i - x_{i-1})(x_i - x_{i+1})\cdots (x_i - x_n)}\\
\end{align}
\]`
`\[(i = 0,1,...,n)\]`

を定義する. 

---

$L_i (x)$は明らかに$n$次の多項式で

`\[
\begin{equation*}
L_i(x_k) = \delta_{ik}=
      \left\{
      \begin{aligned}
             1 \hspace{5mm}\mbox{if}\ i = k\\
             0 \hspace{5mm} \mbox{if}\ i\not= k\\
      \end{aligned}
      \right.
  \end{equation*}
\]`

が成り立つ.

---

このように定義された $L_i (x)\ \ (i = 0,1,...,n)$を用いると, 求める補間多項式$ f_n(x)$は

`\[
g_n(x)= \sum_{i=0}^n y_i L_i(x)\\
\]`

のように表すことができる. この$g_n(x)$はLagrangeの補間多項式(Lagrange interpolation polyno-mial)と呼ばれる.

---

### 2次のAdams-Bashforth法
$x\_k$と $x\_{k+1}$の間の$f(x, y(x))$を現在いる第$k$ステップの情報$x\_k,\ f\_k( = f(x\_k, y \_k))$と, １つ前の第k-1ステップの情報$x\_{k-1}, f\_{k-1}( = f(x\_{k-1}, y\_{k-1}))$を使った１次のLagrange補間多項式 $g\_1(x)$で近似することを考える.

---

すなわち
`\[
\begin{align}
y_{k+1} &= y_k + \int_{x_k}^{x_{k+1}} g_1(x) dx\\
g_1(x) &= f_{k-1} L_{k-1}(x) + f_k L_k (x),\\
L_{k-1}(x) &= \frac{x - x_k}{x_{k-1} - x_k}\hspace{10mm} L_{k}(x) = \frac{x - x_{k-1}}{x_{k} - x_{k-1}}\\
\end{align}
\]`

---

ここで簡単のため刻み$h$は一定とし, また積分の便宜のために $x = x\_k + h\xi$により $\xi$を導入すると, $x_{k\pm1} = x_k \pm h$より

`\[
\int_{x_k}^{x_{k+1}} L_{k-1}(x) dx = \int_{x_k}^{x_{k+1}} \frac{x - x_k}{x_{k-1} - x_k} dx \\
= \int_0^1 \frac{h\xi}{-h} hd\xi = -h \left[\frac{1}{2} \xi^2\right]^1_0 = -\frac{1}{2}h,\\
\]`

---


`\[
\int_{x_k}^{x_{k+1}} L_{k}(x) dx = \int_{x_k}^{x_{k+1}} \frac{x - x_{k-1}}{x_{k} - x_{k-1}} dx\\ 
= \int_0^1 \frac{h(\xi+1)}{h} hd\xi = h \left[\frac{1}{2} \xi^2 + \xi\right]^1_0 = \frac{3}{2}h,\\
\]`

となって１つの積分スキーム
`\[
y_{k+1} = y_k + \frac{h}{2}(-f_{k-1} + 3f_k)\\
\]`
を得る.

---

終わり

---

刻み幅 $x\_{k+1} - x{k} = h$で, $x = x\_k + h\xi$ より
`\[
\frac{x - x_k}{x_{k-1} - x_k} = \frac{h\xi}{x_{k-1} - x_k}
\]`

$x\_{k} - x\_{k-1} = h$より $x\_{k-1} - x_k = -h$

`\[
\frac{h\xi}{x_{k-1} - x_k} = \frac{h\xi}{-h}
\]`

$x = x_k + h\xi$より, 
`\[
\frac{dx}{d\xi} = h,\ \ \ \ x: x_k \to x_{k+1}\\
\xi: 0 \to 1\ \ \ \because x_{k+1} - x_k = h
\]`

---

このスキームの局所打ち切り誤差は
`\[
\begin{align}
\Delta_{k+1} &= y(x_{k+1}) - \left\{y(x_k) + \frac{h}{2}[-f(x_{k-1}, y(x_{k-1})) + 3f(x_k, y(x_k))]\right\}\\
&= y(x_{k+1}) -  y(x_{k}) - \frac{h}{2}[-y^\prime(x_{k-1})) + 3y^\prime(x_k)]\\
&= \left\{ y + hy^\prime + \frac{1}{2}h^2 y^{\prime\prime} + \frac{1}{6}h^3y^{\prime\prime\prime} + \cdots \right\} - \frac{h}{2}\left[2y^\prime + hy^{\prime\prime} - \frac{1}{2}h^2hy^{\prime\prime\prime} + \cdots \right]\\
&= \frac{5}{12}h^3hy^{\prime\prime\prime} + O(h^4)\\
\end{align}
\]`