"""Método del Punto Fijo - Laboratorio No. 4.

Función: f(x) = e^(-x) - x   (raíz ≈ 0.5671432904)
Se prueban varias formas x = g(x) y se evalúa |g'(r)| < 1 en la raíz.
"""
import math
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

RAIZ = 0.5671432904097838   # solución de e^(-x) = x, usada solo para evaluar |g'(r)|


def f(x):
    """Evalúa $f(x)=e^{-x}-x$, cuya raíz es aproximadamente ``RAIZ``."""
    return math.exp(-x) - x


# Formas de g(x)
def g1(x):
    """Primera transformación de punto fijo: $g_1(x)=e^{-x}$."""
    return math.exp(-x)


def g2(x):
    """Segunda transformación de punto fijo: ``g2(x) = -log(x)``."""
    return -math.log(x)


def g3(x):
    """Tercera transformación de punto fijo: $g_3(x)=(x+e^{-x})/2$."""
    return (x + math.exp(-x)) / 2

dg1 = lambda x: -math.exp(-x)
dg2 = lambda x: -1.0 / x
dg3 = lambda x: (1 - math.exp(-x)) / 2


def punto_fijo(g, x0, tol=1e-5, max_iter=100):
    """Itera x_{i+1} = g(x_i).

    Retorna: (x, iteraciones, historial, estado)
        estado: "convergio", "no convergio (max_iter)" o "fallo: <motivo>".
        historial: lista de tuplas (i, x_i, g(x_i), Er%).
    """
    if tol <= 0 or max_iter < 1:
        raise ValueError("tol debe ser positiva y max_iter debe ser mayor que cero.")
    xi = x0
    historial = []
    for i in range(1, max_iter + 1):
        try:
            x_sig = g(xi)
        except (ValueError, OverflowError, ZeroDivisionError) as e:
            return xi, i, historial, f"fallo: {e} (xi = {xi:.6g})"
        if not math.isfinite(x_sig):
            return x_sig, i, historial, "fallo: el valor generado no es finito"
        er = abs((x_sig - xi) / x_sig) * 100 if x_sig != 0 else float("inf")
        historial.append((i, xi, x_sig, er))
        if abs(x_sig) > 1e10:
            return x_sig, i, historial, "fallo: divergió a infinito"
        if abs(x_sig - xi) < tol:
            return x_sig, i, historial, "convergio"
        xi = x_sig
    return xi, max_iter, historial, "no convergio (max_iter)"


def telarana(g, x0, n=10, nombre="g(x)", archivo="graficas/punto_fijo.png", xlim=(0.05, 1.2)):
    """Genera y guarda el diagrama de telaraña de la iteración ``x[n+1]=g(x[n])``."""
    Path(archivo).parent.mkdir(parents=True, exist_ok=True)
    x = np.linspace(*xlim, 400)
    gv = np.array([g(v) if v > 0 else np.nan for v in x])
    plt.figure(figsize=(6, 6))
    plt.plot(x, gv, label=nombre)
    plt.plot(x, x, "k--", label="y = x")
    xi = x0
    try:
        for _ in range(n):
            xn = g(xi)
            plt.plot([xi, xi, xn], [xi, xn, xn], "r-", linewidth=1)
            xi = xn
    except (ValueError, OverflowError, ZeroDivisionError):
        pass
    plt.xlim(xlim); plt.ylim(xlim)
    plt.legend(); plt.grid(True); plt.title(f"Telaraña - {nombre}")
    plt.savefig(archivo, dpi=150, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    x0 = 0.5
    casos = [("g1(x) = e^(-x)", g1, dg1, "graficas/pf_g1.png"),
             ("g2(x) = -ln(x)", g2, dg2, "graficas/pf_g2.png"),
             ("g3(x) = (x + e^(-x))/2", g3, dg3, "graficas/pf_g3.png")]

    for nombre, g, dg, arch in casos:
        x, n, hist, estado = punto_fijo(g, x0)
        print(f"\n=== {nombre} ===")
        print(f"|g'(r)| en la raíz = {abs(dg(RAIZ)):.4f}  ->  criterio |g'|<1: {abs(dg(RAIZ)) < 1}")
        for i, xi, gx, er in hist[:8]:
            print(f"  i={i:>2}  x={xi:>10.6f}  g(x)={gx:>10.6f}  Er={er:>9.4f}%")
        print(f"Resultado: {estado} | x = {x:.8f} | iteraciones = {n}")
        telarana(g, x0, nombre=nombre, archivo=arch)