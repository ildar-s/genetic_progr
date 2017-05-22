# genetic_progr
Genetic programming algorithm

Implementation of the genetic programming algorithm. An application to the classification problem is added too for demo purposes.

There is an option for vizualizing trees, converting the string expression to the tree representation and vice versa. An example is shown below.

String representation:
```
(pow2(tanh((-9.00)+x[1])*(sin(x[0])*((-4.00)+(-8.00))))-((pow2(tanh(x[1]))+(pow2(x[0])+(x[2]*x[1])))*(pow3(tanh(x[2]))/pow3(1.00))))
```

Tree representation (vizualized by the code):
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
