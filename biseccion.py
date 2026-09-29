"""Método de Bisección - Laboratorio No. 4 (Métodos Numéricos, UTP).

Función de prueba: f(x) = x^3 - 7x + 6  (raíces exactas: -3, 1 y 2)
"""
import math
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def f(x):
    """Evalúa la función de prueba $f(x)=x^3-7x+6$."""
    return x**3 - 7 * x + 6


def df(x):
    """Evalúa la derivada de la función de prueba."""
    return 3 * x**2 - 7


def biseccion(f, a, b, tol=1e-5, max_iter=100):
    """Encuentra una raíz de f en [a, b] por bisección.

    Parámetros:
        f: función continua.
        a, b: extremos del intervalo, con f(a)*f(b) < 0 (Teorema de Bolzano).
        tol: tolerancia sobre el error absoluto estimado (b - a)/2.
        max_iter: límite de iteraciones.

    Retorna:
        (raiz, iteraciones, historial, convergio)
        historial: lista de tuplas (i, a, b, m, f(m), Ea, Er%).
    """
    if tol <= 0 or max_iter < 1:
        raise ValueError("tol debe ser positiva y max_iter debe ser mayor que cero.")
    fa, fb = f(a), f(b)
    if not math.isfinite(fa) or not math.isfinite(fb):
        raise ValueError("f(a) y f(b) deben ser valores finitos.")
    if fa == 0:
        return a, 0, [], True
    if fb == 0:
        return b, 0, [], True
    if fa * fb > 0:
        raise ValueError("No hay cambio de signo en [a, b]: Bolzano no se cumple.")

    historial = []
    m_ant = None
    for i in range(1, max_iter + 1):
        m = (a + b) / 2.0
        fm = f(m)
        ea = (b - a) / 2.0
        er = abs((m - m_ant) / m) * 100 if (m_ant is not None and m != 0) else float("nan")
        historial.append((i, a, b, m, fm, ea, er))

        if fm == 0 or ea < tol:
            return m, i, historial, True

        if fa * fm < 0:
            b = m
            fb = fm
        else:
            a = m
            fa = fm
        m_ant = m
    return m, max_iter, historial, False


def biseccion_hibrida(f, df, a, b, tol=1e-5, max_iter=100):
    """Combina Newton con bisección para conservar el intervalo acotado.

    Si el paso de Newton sale del intervalo o la derivada es nula, usa el
    punto medio y conserva el cambio de signo.
    """
    if tol <= 0 or max_iter < 1:
        raise ValueError("tol debe ser positiva y max_iter debe ser mayor que cero.")
    fa = f(a)
    fb = f(b)
    if fa == 0:
        return a, 0, [], True
    if fb == 0:
        return b, 0, [], True
    if fa * fb > 0:
        raise ValueError("No hay cambio de signo en [a, b]: Bolzano no se cumple.")

    historial = []
    m = (a + b) / 2.0
    m_ant = None
    for i in range(1, max_iter + 1):
        fm = f(m)
        ea = (b - a) / 2.0
        er = abs((m - m_ant) / m) * 100 if (m_ant is not None and m != 0) else float("nan")
        historial.append((i, a, b, m, fm, ea, er))

        if fm == 0 or abs(fm) < tol or ea < tol:
            return m, i, historial, True

        derivada = df(m)
        candidato = m - fm / derivada if derivada != 0 else None
        if candidato is None or not (a < candidato < b):
            candidato = (a + b) / 2.0

        fc = f(candidato)
        if fc == 0:
            return candidato, i, historial, True
        if fa * fc < 0:
            b, fb = candidato, fc
        else:
            a, fa = candidato, fc
        m_ant = m
        m = candidato

    return m, max_iter, historial, False


def graficar(f, a, b, raiz, xlim=(-4, 3), archivo="graficas/biseccion.png"):
    """Genera y guarda la gráfica de la función y del intervalo final."""
    Path(archivo).parent.mkdir(parents=True, exist_ok=True)
    x = np.linspace(xlim[0], xlim[1], 500)
    plt.figure(figsize=(8, 5))
    plt.plot(x, f(x), label="f(x) = x³ - 7x + 6")
    plt.axhline(0, color="black", linewidth=0.8)
    plt.axvspan(a, b, alpha=0.15, color="orange", label=f"Intervalo [{a}, {b}]")
    plt.plot(raiz, f(raiz), "ro", label=f"Raíz ≈ {raiz:.6f}")
    plt.grid(True)
    plt.legend()
    plt.title("Método de Bisección")
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.savefig(archivo, dpi=150, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    a, b, tol = 1.5, 3.0, 1e-5   # cambia el intervalo para hallar las otras raíces (-3 y 1)
    print(f"f(a) = {f(a):.4f}, f(b) = {f(b):.4f}  ->  f(a)*f(b) < 0: {f(a)*f(b) < 0}")

    raiz, n, hist, ok = biseccion_hibrida(f, df, a, b, tol)

    print(f"\n{'i':>3} {'a':>10} {'b':>10} {'m':>12} {'f(m)':>12} {'Ea':>10} {'Er(%)':>10}")
    for i, ai, bi, m, fm, ea, er in hist:
        print(f"{i:>3} {ai:>10.6f} {bi:>10.6f} {m:>12.8f} {fm:>12.3e} {ea:>10.2e} {er:>10.4f}")

    n_teorico = math.ceil(math.log2((b - a) / tol))
    print(f"\nRaíz aproximada: {raiz:.8f}")
    print(f"Iteraciones realizadas: {n} | Iteraciones teóricas n >= log2((b-a)/tol): {n_teorico}")

    graficar(f, a, b, raiz)