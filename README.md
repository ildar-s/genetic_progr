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


```
All options:
  --help                          help message

Parameters of GP algorithm:
  --ntrees arg (=1000)            number of trees
  --max_depth arg (=2)            max depth of the tree, in the initial forest
  --max_depth_ex arg (=5)         max depth of the tree, after transformations
  --ndim arg (=0)                 dimensionality of the parameter space. Only 
                                  need to set if no data is provided, e.g. for 
                                  random tree vizualization
  --ngen_max arg (=30000)         max generations
  --err_thr arg (=150)            error threshold used in stop criteria
  --lambda arg (=0.1)             is a parameter controlling depth of 
                                  cross-over points: lower lambdas mean higher 
                                  chance of picking lower-depth points thus 
                                  compensating for the natural abundance of 
                                  high-depth nodes in trees (e.g. leafs)
  --ratio_ss arg (=0.01)          ratio of trees considered for cross-over to 
                                  the total number of trees: the selection 
                                  among only a subset of the forest gives 
                                  better chance to the not-so-good trees to be 
                                  selected too.
  --func_prob arg (=0.85)         probability of creating a function node (as 
                                  opposed to creating a terminal node) in 
                                  forest initialization
  --term_const_prob arg (=0.2)    probability of creating a constant terminal 
                                  node (as opposed to creating a variable 
                                  terminal node) in forest initialization
  --consts_min arg (=-10.0)       minimum value of terminal constants
  --consts_max arg (=10.0)        maximum value of terminal constants
  --consts_n arg (=20)            number of different terminal constants

I/O control:
  --cfg arg                       GP configuration file, if not present, will 
                                  attempt to read "default.cfg", if not present
                                  either, will use default values
  -d [ --data ] arg               CSV file with data, first column always 
                                  ignored, if it's training data then last 
                                  column - labels (target)
  -f [ --final_regressor ] arg    final regressor file
  -i [ --initial_regressor ] arg  initial regressor file
  -p [ --predictions ] arg        file to output predictions to
  --csv_has_header arg (=yes)     yes or no

Workflow control:
  --logging arg (=info)           logging level: "silent", "info", or "debug"
  --n_iter arg (=1)               number of fitting procedures and resulting 
                                  decision trees
  --bs                            use bootstrapped subsample for each fitting 
                                  procedure
  -t [ --testonly ]               ignore label information and just test
  --random_seed arg               random seed
  --random_tree_demo              vizualize a random tree and quit
```
