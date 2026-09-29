"""Falsa Posición (Regula Falsi) y variante Illinois - Laboratorio No. 4.

Función de prueba: f(x) = x^3 + x^2 - 1  (raíz real ≈ 0.7549)
Es convexa en [0, 2], por lo que la versión clásica deja un extremo "congelado".
"""
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt


def f(x):
    """Evalúa la función de prueba ``f(x) = x**3 + x**2 - 1``."""
    return x**3 + x**2 - 1


def falsa_posicion(f, a, b, tol=1e-5, max_iter=1000, metodo="clasica"):
    """Regula Falsi en [a, b].

    metodo: "clasica" o "illinois" (divide entre 2 el f del extremo que se
            repite dos veces seguidas, rompiendo el estancamiento).

    Retorna:
        (raiz, iteraciones, historial, convergio, extremo_fijo)
        historial: lista de tuplas (i, a, b, xr, f(xr), Ea).
        extremo_fijo: cuántas iteraciones seguidas terminó igual el extremo más estancado.
    """
    if tol <= 0 or max_iter < 1:
        raise ValueError("tol debe ser positiva y max_iter debe ser mayor que cero.")
    if metodo not in {"clasica", "illinois"}:
        raise ValueError("metodo debe ser 'clasica' o 'illinois'.")
    fa, fb = f(a), f(b)
    if fa == 0:
        return a, 0, [], True, 0
    if fb == 0:
        return b, 0, [], True, 0
    if fa * fb > 0:
        raise ValueError("No hay cambio de signo en [a, b]: Bolzano no se cumple.")

    historial = []
    xr_ant = None
    lado = 0           # -1: se movió b la última vez, +1: se movió a
    racha = maxracha = 0

    for i in range(1, max_iter + 1):
        denominador = fa - fb
        if abs(denominador) < 1e-15:
            raise ZeroDivisionError("Regula Falsi: denominador casi nulo.")
        xr = b - fb * (a - b) / denominador
        fr = f(xr)
        ea = abs(xr - xr_ant) if xr_ant is not None else float("inf")
        historial.append((i, a, b, xr, fr, ea))

        if fr == 0 or abs(fr) < tol or ea < tol:
            return xr, i, historial, True, maxracha

        if fa * fr < 0:                 # raíz en [a, xr] -> se mueve b
            b, fb = xr, fr
            if metodo == "illinois" and lado == -1:
                fa /= 2.0
            racha = racha + 1 if lado == -1 else 1
            lado = -1
        else:                           # raíz en [xr, b] -> se mueve a
            a, fa = xr, fr
            if metodo == "illinois" and lado == +1:
                fb /= 2.0
            racha = racha + 1 if lado == +1 else 1
            lado = +1
        maxracha = max(maxracha, racha)
        xr_ant = xr
    return xr, max_iter, historial, False, maxracha


def graficar(hist_clasica, hist_illinois, archivo="graficas/falsa_posicion.png"):
    """Genera y guarda la comparación gráfica entre ambas variantes."""
    Path(archivo).parent.mkdir(parents=True, exist_ok=True)
    x = np.linspace(-0.2, 2.1, 400)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    ax[0].plot(x, f(x), label="f(x)")
    ax[0].axhline(0, color="black", linewidth=0.8)
    ax[0].plot([h[3] for h in hist_clasica[:6]], [h[4] for h in hist_clasica[:6]], "ro-", label="xr (clásica, 6 primeras)")
    ax[0].set_title("f(x) = x³ + x² - 1")
    ax[0].legend(); ax[0].grid(True)

    ax[1].semilogy([h[0] for h in hist_clasica], [h[5] if h[5] != float("inf") else np.nan for h in hist_clasica], "r.-", label="Clásica")
    ax[1].semilogy([h[0] for h in hist_illinois], [h[5] if h[5] != float("inf") else np.nan for h in hist_illinois], "g.-", label="Illinois")
    ax[1].set_title("Error absoluto estimado por iteración")
    ax[1].set_xlabel("Iteración"); ax[1].set_ylabel("Ea"); ax[1].legend(); ax[1].grid(True, which="both")
    plt.savefig(archivo, dpi=150, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    a, b, tol = 0.0, 2.0, 1e-5

    r1, n1, h1, ok1, fijo1 = falsa_posicion(f, a, b, tol, metodo="clasica")
    r2, n2, h2, ok2, fijo2 = falsa_posicion(f, a, b, tol, metodo="illinois")

    print("--- Falsa posición CLÁSICA ---")
    print(f"{'i':>3} {'a':>10} {'b':>10} {'xr':>12} {'f(xr)':>12} {'Ea':>10}")
    for i, ai, bi, xr, fr, ea in h1[:10]:
        print(f"{i:>3} {ai:>10.6f} {bi:>10.6f} {xr:>12.8f} {fr:>12.3e} {ea:>10.2e}")
    print(f"Raíz: {r1:.8f} | Iteraciones: {n1} | Iteraciones seguidas con un extremo fijo: {fijo1}")

    print("\n--- Falsa posición con ILLINOIS ---")
    print(f"Raíz: {r2:.8f} | Iteraciones: {n2} | Iteraciones seguidas con un extremo fijo: {fijo2}")

    graficar(h1, h2)