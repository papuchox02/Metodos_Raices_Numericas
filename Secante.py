"""Método de la Secante y comparación con Newton-Raphson - Laboratorio No. 4.

Función: f(x) = cos(x) - x   (raíz ≈ 0.7390851332)
"""
import math
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from newton import newton_raphson, orden_convergencia

RAIZ = 0.7390851332151607


def f(x):
    """Evalúa la función de prueba ``f(x) = cos(x) - x``."""
    return math.cos(x) - x


def df(x):
    """Evalúa la derivada de la función de prueba."""
    return -math.sin(x) - 1


def secante(f, x0, x1, tol=1e-5, max_iter=50):
    """x_{n+1} = x_n - f(x_n)(x_n - x_{n-1}) / (f(x_n) - f(x_{n-1})).

    Retorna: (x, iteraciones, historial, convergio)
        historial: lista de tuplas (i, x_{i-1}, x_i, x_{i+1}, Er%).
    """
    if tol <= 0 or max_iter < 1:
        raise ValueError("tol debe ser positiva y max_iter debe ser mayor que cero.")
    if x0 == x1:
        raise ValueError("x0 y x1 deben ser distintos para iniciar la secante.")
    historial = []
    for i in range(1, max_iter + 1):
        f0, f1 = f(x0), f(x1)
        if abs(f1 - f0) < 1e-14:
            raise ZeroDivisionError("Secante horizontal: f(x_n) ≈ f(x_{n-1}).")
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        er = abs((x2 - x1) / x2) * 100 if x2 != 0 else float("inf")
        historial.append((i, x0, x1, x2, er))
        if abs(x2 - x1) < tol or abs(f(x2)) < 1e-14:
            return x2, i, historial, True
        x0, x1 = x1, x2
    return x1, max_iter, historial, False


if __name__ == "__main__":
    tol = 1e-10

    xs, ns, hs, oks = secante(f, 0.0, 1.0, tol)
    xn, nn, hn, okn = newton_raphson(f, df, 1.0, tol)

    print("=== SECANTE (x0 = 0, x1 = 1) ===")
    for i, a, b, c, er in hs:
        print(f"i={i}  x_(i-1)={a:.8f}  x_i={b:.8f}  x_(i+1)={c:.10f}  Er={er:.3e}%")
    print(f"Raíz = {xs:.10f} | iteraciones = {ns}")

    print("\n=== NEWTON-RAPHSON (x0 = 1) ===")
    for i, xi, fx, dfx, x_sig, er in hn:
        print(f"i={i}  x_i={xi:.8f}  x_(i+1)={x_sig:.10f}  Er={er:.3e}%")
    print(f"Raíz = {xn:.10f} | iteraciones = {nn}")

    ps = orden_convergencia([0.0, 1.0] + [h[3] for h in hs], RAIZ)
    pn = orden_convergencia([1.0] + [h[4] for h in hn], RAIZ)
    print(f"\nOrden estimado - Secante (≈1.618): {[round(v, 2) for v in ps]}")
    print(f"Orden estimado - Newton  (≈2):     {[round(v, 2) for v in pn]}")

    # Gráfica de error vs iteración
    es = [abs(h[3] - RAIZ) for h in hs]
    en = [abs(h[4] - RAIZ) for h in hn]
    plt.figure(figsize=(7, 5))
    plt.semilogy(range(1, len(es) + 1), [max(e, 1e-17) for e in es], "o-", label="Secante")
    plt.semilogy(range(1, len(en) + 1), [max(e, 1e-17) for e in en], "s-", label="Newton-Raphson")
    plt.xlabel("Iteración"); plt.ylabel("|x - r|"); plt.grid(True, which="both"); plt.legend()
    plt.title("cos(x) - x = 0: Secante vs Newton")
    Path("graficas").mkdir(parents=True, exist_ok=True)
    plt.savefig("graficas/secante_vs_newton.png", dpi=150, bbox_inches="tight")
    plt.show()