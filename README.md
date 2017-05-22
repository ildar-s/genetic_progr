# genetic_progr
Genetic programming algorithm

String representation:
```
(pow2(tanh((-9.00)+x[1])*(sin(x[0])*((-4.00)+(-8.00))))-((pow2(tanh(x[1]))+(pow2(x[0])+(x[2]*x[1])))*(pow3(tanh(x[2]))/pow3(1.00))))
```

Tree representation:
```
  ┌─────pow2
  │       │       ┌─────tanh
  │       │       │       │       ┌─────(-9.00)
  │       │       │       └─────+
  │       │       │               └─────x[1]
  │       └─────*
  │               │       ┌─────sin
  │               │       │       └─────x[0]
  │               └─────*
  │                       │       ┌─────(-4.00)
  │                       └─────+
  │                               └─────(-8.00)
-
  │               ┌─────pow2
  │               │       └─────tanh
  │               │               └─────x[1]
  │       ┌─────+
  │       │       │       ┌─────pow2
  │       │       │       │       └─────x[0]
  │       │       └─────+
  │       │               │       ┌─────x[2]
  │       │               └─────*
  │       │                       └─────x[1]
  └─────*
          │       ┌─────pow3
          │       │       └─────tanh
          │       │               └─────x[2]
          └─────/
                  └─────pow3
                          └─────1.00
```
