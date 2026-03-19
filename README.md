# Physics-Informed smoothing 
This repo implements a class of physics-informed smoothing models to integrate empirical data with physical apriori knowledge about the data generating mechanism. In this repo, the physics will be expressed through a wide range of possibly nonlinear ordinary differential equations (ODEs). The general estimation problem looks as follows:

$$
J(u,\theta) = \sum_{j}(y_{j}-y(t_{j}))^{2}+ \lambda \int_{[0,T]} ||u||^{2}dt 
\\
\\
\mathrm{s.t.}\quad 
\begin{cases}
    \frac{dy}{dt}(t)=f(t,y(t),\theta)+u(t) \\
    y(0) = y_{0}
\end{cases}
$$
where $\theta$ is a vector parameter used to parametrise the regularising ODE model. The $\lambda$ parameter allows the model to span from a pure data-fitting framework ($\lambda \to 0$) to pure parameter estimation ($\lambda \to \infty$). 


## Project structure
- `src/solvers` contains the actual implementations of the solver used to solve the estimation problem. Multiple estimation strategies are implemented, more details are available in the dedicated section.
- `tests` provides some tests to inspect if the behaviour of the solvers complies with the expected one. These are used to test if new modifications break the current implementation.
- `examples`: contains a set of examples used to inspect the functioning of the implemented solvers
- `simulations` will contain a series of simulation studies used to extensively compare the various implemented methods

## Estimation algorithms
The repo implements two major approaches to solve the estimation problems, which are the Parameter Cascading approach and the Tracking Estimation approach. Other approaches will be implemented as well for comparison purposes (e.g. nonlinear least squares).

### Parameter cascading
The Parameter Cascading approach [(Ramsay,2007)](https://doi.org/10.1111/j.1467-9868.2007.00610.x) implemented in `src/solvers/parameter_cascading.R` consists of a multi-criteria approach. The optimisation problem is then slightly different from the estimation problem presented before. In particular, we distinguish between an outer optimisation criterion with respect to $\theta$
$$
H(\theta) = \sum_{j} (y_{j}-y(t_{j}))^{2},
$$
and an inner one to be optimised with respect to $u$ (with $\theta$ fixed):
$$
J(u,\theta) = \sum_{j}(y_{j}-y(t_{j}))^{2}+ \lambda \int_{[0,T]} ||u||^{2}dt 
\\
\\
\mathrm{s.t.}\quad 
\begin{cases}
    \frac{dy}{dt}(t)=f(t,y(t),\theta)+u(t) \\
    y(0) = y_{0}
\end{cases}
$$

The inner problem is solved using an iterative scheme that employs the adjoint method to compute the gradient of the optimal control problem governed by the ODE. The gradient computation is based on the followigng optimality system:
$$
\begin{cases}
    \frac{dy}{dt} = f(t,y(t), \theta) + u(t) \\
    -\frac{dp}{dt} = \frac{df}{dy}p(t) + 2\sum_{j}(y_{j}-y(t_{j})) \\
    u + \frac{1}{2\lambda}p \geq0
\end{cases}
$$
The inner gradient is then fed into an appropriate optimisation algorithm to converge to the optimum. This optimisation scheme is implemented within `src/solvers/general_ode_system_solver.R`, more details on this solver will be provided later.

On the other hand, the outer gradient is computed relying on the implicit function theorem. To see why, consider gradient of the outer error criterion expanded thorugh the chain rule
$$
\frac{dH}{d\theta} = \frac{\partial H}{\partial y} \frac{d y}{d \theta} 
=  \frac{\partial H}{\partial y} (\frac{\partial y}{\partial \theta} +\frac{\partial y}{\partial u}\frac{d u}{d \theta}).
$$

Note that such computation involves the sensitivity of $y$ and $u$ to $\theta$. Such quantities can be computed by differentiating the inner optimality system with respect to $\theta$

$$
\begin{cases}
    \frac{dy_\theta}{dt} = f_{\theta}(t,y(t), \theta) + f(t,y(t), \theta)y_{\theta}(t) + u_{\theta}(t) \\
    -\frac{dp_{\theta}}{dt} = \frac{df}{dy}p_{\theta}(t) + 2\sum_{j}(y_{j}-y(t_{j})) \\
    u + \frac{1}{2\lambda}p \geq0
\end{cases}
$$


### Tracking estimator
The tracking estimator is based on the works of [(Brunel and Clairon, 2015)](https://doi.org/10.48550/arXiv.1410.7554), [(Clairon and Brunel,2018)](https://doi.org/10.1080/01621459.2017.1319841) and here implemented in a modified version in `src/solvers/tracking_ode_solver.R`. The approach is quite similar to parameter cascading. Actually, it can be considered as a parameter cascading having the same optimisation criterion $J(u,\theta)$ both for the inner and the outer problem, namely $H(\theta)=\min_{u}J(u,\theta)$. For this reason, the solution of the inner problem is delegated to the same solver (`src/solvers/general_ode_system_solver.R`) used by parameter cascading. 

The slight modification made to the outer optimisation cost allows for simplified gradient computation and a nice interpretation of the optimisation procedure. Let's start by looking at the gradient: 

$$
\frac{dH}{d\theta} = \frac{\partial H}{\partial y} \frac{d y}{d \theta} 
=  \frac{\partial H}{\partial y} (\frac{\partial y}{\partial \theta} +\frac{\partial y}{\partial u}\frac{\partial u}{\partial \theta}) 
= \frac{\partial H}{\partial y} \frac{\partial y}{\partial \theta} +\frac{\partial H}{\partial u}\frac{\partial u}{\partial \theta}
$$

However note that, since $H(\theta)=\min_{u} J(u,\theta)$, as the inner optimisation converges, $\frac{\partial H}{\partial u} = 0$, implying
$\frac{dH}{d\theta} = \frac{\partial H}{\partial y} \frac{\partial y}{\partial \theta}$. Thus, we can avoid the computation of $\frac{dy}{d\theta}$ by solving the sensitivity system as done in the parameter cascading approach, but we can compute $\frac{\partial y}{\partial \theta}$ by a single forward pass. Additional speedups can obtained by relying on the adjoint method, exploiting the already computed adjoint variable $p$.

[TO BE COMPLETED]

Another intersting aspect is the interpretation of the optimisation procedure as a joint optimisation approach. Indeed, since the outer cost function is the same as the inner one, $\min_{u,\theta} J(u,\theta) = \min_{\theta}\min_{u}J(u,\theta) = \min_{\theta}H(\theta)$.

### The inner solver

[TO BE COMPLETED]

## Work In Progress
- One of the goal of the project is to quantify the uncertainty related to the estimate. This is particularly difficult due to the presence of non-linearity. The idea is to consider a linearization around the current estimate, using techiniques as the Laplace approximation and delta method
- Some difficulites that may be encountered in this context: relating the uncertainty in the estimation of the nonparametric term $y$ and the uncertainty related to parameter $\theta$ estimation, which role plays $\lambda$ in UQ, ...
- 