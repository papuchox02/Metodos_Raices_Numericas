"""Método de Newton-Raphson - Laboratorio No. 4.

Casos:
  1) f(x) = x^3 - x - 2  (converge, cuadrática)
  2) f(x) = x^(1/3)      (diverge: cada paso da x_{n+1} = -2 x_n)
  3) Newton amortiguado (backtracking) sobre el caso 2
"""
import math
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def newton_raphson(f, df, x0, tol=1e-5, max_iter=50):
    """x_{n+1} = x_n - f(x_n)/f'(x_n).

    Retorna: (x, iteraciones, historial, convergio)
        historial: lista de tuplas (i, x_i, f(x_i), f'(x_i), x_{i+1}, Er%).
    """
    if tol <= 0 or max_iter < 1:
        raise ValueError("tol debe ser positiva y max_iter debe ser mayor que cero.")
    xi = x0
    historial = []
    for i in range(1, max_iter + 1):
        dfx = df(xi)
        if abs(dfx) < 1e-12:
            raise ZeroDivisionError(f"Derivada casi nula en x = {xi}")
        x_sig = xi - f(xi) / dfx
        er = abs((x_sig - xi) / x_sig) * 100 if x_sig != 0 else float("inf")
        historial.append((i, xi, f(xi), dfx, x_sig, er))
        if abs(x_sig - xi) < tol or abs(f(x_sig)) < 1e-12:
            return x_sig, i, historial, True
        xi = x_sig
    return xi, max_iter, historial, False


def newton_amortiguado(f, df, x0, tol=1e-5, max_iter=100):
    """Aplica Newton con backtracking hasta reducir el residuo.

    Retorna el mismo formato que :func:`newton_raphson`. Si no encuentra un
    paso que reduzca el residuo, termina con ``convergio=False``.
    """
    if tol <= 0 or max_iter < 1:
        raise ValueError("tol debe ser positiva y max_iter debe ser mayor que cero.")
    xi = x0
    historial = []
    for i in range(1, max_iter + 1):
        fx, dfx = f(xi), df(xi)
        if fx == 0:
            return xi, i, historial, True
        if abs(dfx) < 1e-12:
            raise ZeroDivisionError(f"Derivada casi nula en x = {xi}")
        paso = fx / dfx
        lam = 1.0
        while abs(f(xi - lam * paso)) >= abs(fx) and lam > 1e-6:
            lam /= 2.0
        if abs(f(xi - lam * paso)) >= abs(fx):
            return xi, i, historial, False
        x_sig = xi - lam * paso
        historial.append((i, xi, fx, dfx, x_sig, lam))
        if abs(x_sig - xi) < tol:
            return x_sig, i, historial, True
        xi = x_sig
    return xi, max_iter, historial, False


def orden_convergencia(xs, raiz):
    """Estima $p$ con $p≈ln(e[n+1]/e[n])/ln(e[n]/e[n-1])$."""
    e = [abs(x - raiz) for x in xs if abs(x - raiz) > 1e-15]
    return [math.log(e[k + 1] / e[k]) / math.log(e[k] / e[k - 1])
            for k in range(1, len(e) - 1)
            if e[k] != e[k - 1] and e[k + 1] > 0]


# Funciones de los casos
f1 = lambda x: x**3 - x - 2
df1 = lambda x: 3 * x**2 - 1
RAIZ1 = 1.5213797068045676

f2 = lambda x: np.cbrt(x)                       # raíz real de x^(1/3)
df2 = lambda x: 1.0 / (3.0 * np.cbrt(x) ** 2)


if __name__ == "__main__":
    # ---- Caso 1: convergencia ----
    print("=== CASO 1: f(x) = x^3 - x - 2, x0 = 1.5 ===")
    x, n, hist, ok = newton_raphson(f1, df1, 1.5, tol=1e-10)
    for i, xi, fx, dfx, xs, er in hist:
        print(f"i={i}  x={xi:.10f}  f={fx:.3e}  f'={dfx:.4f}  Er={er:.3e}%")
    print(f"Raíz = {x:.10f} | iteraciones = {n}")
    p = orden_convergencia([1.5] + [h[4] for h in hist], RAIZ1)
    print("Orden de convergencia estimado (debería acercarse a 2):", [round(v, 2) for v in p])

    # ---- Caso 2: divergencia ----
    print("\n=== CASO 2: f(x) = x^(1/3), x0 = 0.5 (DIVERGE) ===")
    x, n, hist2, ok = newton_raphson(f2, df2, 0.5, max_iter=12)
    for i, xi, fx, dfx, xs, er in hist2:
        print(f"i={i:>2}  x={xi:>12.4f}  f={fx:>8.4f}  f'={dfx:>10.5f}  x_sig={xs:>12.4f}")
    print(f"¿Convergió? {ok}  (cada paso hace x_(n+1) = -2 x_n)")

    # ---- Caso 3: solución (amortiguado) ----
    print("\n=== CASO 3: Newton amortiguado sobre x^(1/3), x0 = 0.5 ===")
    x, n, hist3, ok = newton_amortiguado(f2, df2, 0.5)
    print(f"¿Convergió? {ok} | x = {x:.3e} | iteraciones = {n} | primer lambda usado = {hist3[0][5]}")

    # ---- Gráfica ----
    xx = np.linspace(-3, 3, 500)
    plt.figure(figsize=(8, 5))
    plt.plot(xx, np.cbrt(xx), label="f(x) = x^(1/3)")
    plt.axhline(0, color="black", linewidth=0.8)
    trayectoria = [h[1] for h in hist2[:6]]
    plt.plot(trayectoria, [np.cbrt(t) for t in trayectoria], "ro-", label="Newton (diverge)")
    plt.xlim(-10, 10); plt.legend(); plt.grid(True)
    plt.title("Newton-Raphson sobre x^(1/3)")
    Path("graficas").mkdir(parents=True, exist_ok=True)
    plt.savefig("graficas/newton.png", dpi=150, bbox_inches="tight")
    plt.show()