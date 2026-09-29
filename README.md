Tabla de resultados 
### 📊 Tabla Comparativa de Resultados Reales

|      Método                   | Función Evaluada         | Parámetros Iniciales     | Raíz Obtenida | Iteraciones           | Comportamiento / Solución con Copilot |
|                               |                          |                          |                |            | :--- |
| **Secante**                   | \(f(x) = \cos(x) - x\)   | \(x_0 = 0\), \(x_1 = 1\) | `0.7390851332` | **6** | Convergencia súper-lineal (orden estimado \(\approx 1.6\)). |
| **Newton-Raphson**            | \(f(x) = \cos(x) - x\)   | \(x_0 = 1\)              | `0.7390851332` | **4** | Más rápido que la secante. Convergencia cuadrática (orden \(\approx 2.0\)). |
| **Newton Clásico**            | \(f(x) = x^3 - x - 2\)   | \(x_0 = 1.5\)            | `1.5213797068` | **3** | Alta velocidad por cercanía al valor real (orden \(\approx 2.01\)). |
| **Newton (Caso Falla)**       | \(f(x) = x^{1/3}\)       | \(x_0 = 0.5\)             | *Ninguna* | *Diverge* | **Falla.** Entra en ciclo infinito donde cada paso hace \(x_{n+1} = -2x_n\). |
| **Newton Amortiguado**        | \(f(x) = x^{1/3}\)       | \(x_0 = 0.5\)            | `1.907e-06` | **18** | **Solución:** Copilot aplicó un factor lambda \(= 0.5\) que logró la convergencia. |
| **Punto Fijo**                | \(g_1(x) = e^{-x}\)      | \(x_0 = 0.5\)            | `0.56714076` | **18** | Converge de manera lineal debido a que \(\Vert{}g'(r)\Vert{} < 1\) (True). |
| **Falsa Posición Clásica**    | \(f(x) = x^3 + x^2 - 1\) | \([a, b] = [0, 2]\)      | `0.75486526` | **26** | **Estancamiento.** El extremo \(b=2.0\) se congeló durante 25 iteraciones. |
| **Falsa Posición (Illinois)** | \(f(x) = x^3 + x^2 - 1\) | \([a, b] = [0, 2]\)       | `0.75487767` | **9** | **Solución:** Modificación propuesta por Copilot. Redujo las iteraciones fijas a 3. |
