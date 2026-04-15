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
The Parameter Cascading approach [(Ramsay,2007)](https://doi.org/10.1111/j.1467-9868.2007.00610.x) implemented in `src/solvers/inverse-solvers/parameter_cascading.R` consists of a multi-criteria approach. The optimisation problem is then slightly different from the estimation problem presented before. In particular, we distinguish between an outer optimisation criterion with respect to $\theta$

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

The inner problem is solved using an iterative scheme that employs the adjoint method to compute the gradient of the optimal control problem governed by the ODE. The gradient computation is based on the following optimality system:

$$
\begin{cases}
    \frac{dy}{dt} = f(t,y(t), \theta) + u(t) \\
    -\frac{dp}{dt} = \frac{df}{dy}p(t) + 2\sum_{j}(y_{j}-y(t_{j})) \\
    u + \frac{1}{2\lambda}p = 0
\end{cases}
$$

The inner gradient is then fed into an appropriate optimisation algorithm to converge to the optimum. This optimisation scheme is implemented within `src/solvers/forward-solvers/dto_fwd_solver.R`, more details on this solver are provided below.

On the other hand, the outer gradient is computed relying on the implicit function theorem. To see why, consider gradient of the outer error criterion expanded through the chain rule

$$
\frac{dH}{d\theta} = \frac{\partial H}{\partial y} \frac{d y}{d \theta} 
=  \frac{\partial H}{\partial y} (\frac{\partial y}{\partial \theta} +\frac{\partial y}{\partial u}\frac{d u}{d \theta}).
$$

Note that such computation involves the sensitivity of $y$ and $u$ to $\theta$. Such quantities can be computed by differentiating the inner optimality system with respect to $\theta$

$$
\begin{cases}
    \frac{dy_\theta}{dt} = f_{\theta}(t,y(t), \theta) + f(t,y(t), \theta)y_{\theta}(t) + u_{\theta}(t) \\
    -\frac{dp_{\theta}}{dt} = \frac{df}{dy}p_{\theta}(t) + 2\sum_{j}(y_{j}-y(t_{j})) \\
    u_{\theta} + \frac{1}{2\lambda}p_{\theta} = 0
\end{cases}
$$


### Tracking estimator
The tracking estimator is based on the works of [(Brunel and Clairon, 2015)](https://doi.org/10.48550/arXiv.1410.7554), [(Clairon and Brunel,2018)](https://doi.org/10.1080/01621459.2017.1319841) and is implemented in `src/solvers/inverse-solvers/tracking_ode_solver.R`. The approach is quite similar to parameter cascading. Actually, it can be considered as a parameter cascading having the same optimisation criterion $J(u,\theta)$ both for the inner and the outer problem, namely $H(\theta)=\min_{u}J(u,\theta)$. For this reason, the solution of the inner problem is delegated to the same solver (`src/solvers/forward-solvers/dto_fwd_solver.R`) used by parameter cascading. 

The slight modification made to the outer optimisation cost allows for simplified gradient computation and a nice interpretation of the optimisation procedure. Let's start by looking at the gradient: 

$$
\frac{dH}{d\theta} = \frac{\partial H}{\partial y} \frac{d y}{d \theta} 
=  \frac{\partial H}{\partial y} (\frac{\partial y}{\partial \theta} +\frac{\partial y}{\partial u}\frac{\partial u}{\partial \theta}) 
= \frac{\partial H}{\partial y} \frac{\partial y}{\partial \theta} +\frac{\partial H}{\partial u}\frac{\partial u}{\partial \theta}
$$

However note that, since $H(\theta)=\min_{u} J(u,\theta)$, as the inner optimisation converges, $\frac{\partial H}{\partial u} = 0$, implying
$\frac{dH}{d\theta} = \frac{\partial H}{\partial y} \frac{\partial y}{\partial \theta}$. Thus, we can avoid the computation of $\frac{dy}{d\theta}$ by solving the sensitivity system as done in the parameter cascading approach, but we can compute $\frac{\partial y}{\partial \theta}$ by a single forward pass. Additional speedups can be obtained by relying on the adjoint method, exploiting the already computed adjoint variable $p$.

Another interesting aspect is the interpretation of the optimisation procedure as a joint optimisation approach. Indeed, since the outer cost function is the same as the inner one, $\min_{u,\theta} J(u,\theta) = \min_{\theta}\min_{u}J(u,\theta) = \min_{\theta}H(\theta)$.

### The inner solver

The inner problem is solved by the discretize-then-optimize (DtO) solver in `src/solvers/forward-solvers/dto_fwd_solver.R`, with scheme implementations in `src/solvers/numerical-solvers/ode_solvers.R`.

#### Discrete control defect and its linearization (core fix)

Let $t_1 < \dots < t_{n_s}$, $\Delta t_t = t_{t+1}-t_t$. For each interval $t=1,\dots,n_s-1$, the solver defines a discrete control defect

$$
u_t = \mathcal{D}_t(y_t, y_{t+1}, \theta),
$$

which depends on the selected time-discretization scheme.

The key abstraction introduced in the numerical layer is the exact local linearization

$$
\delta u_t
= D^{(t)}_{y_t}\,\delta y_t
+ D^{(t)}_{y_{t+1}}\,\delta y_{t+1}
+ D^{(t)}_{\theta}\,\delta\theta,
$$

with matrices returned by the API:

- $D^{(t)}_{y_t}$ as `Du_dyc`
- $D^{(t)}_{y_{t+1}}$ as `Du_dyn`
- $D^{(t)}_{\theta}$ as `Du_dtheta`

This is exposed through:

- scheme-level method `differentiate_discrete_control_map(...)`
- forward-solver wrapper `get_discrete_control_jacobian(...)`

The API also provides split accessors `dcontrol_dy_curr`, `dcontrol_dy_next`, `dcontrol_dtheta` (and forward-solver wrappers `get_dcontrol_dy_curr`, `get_dcontrol_dy_next`, `get_dcontrol_dtheta`) for cases where only one derivative block is needed.

and is consumed by both tracking and cascading outer gradients.

#### Scheme-specific defects

Let $f_t = f(y_t,t_t,\theta)$.

1. Euler

$$
\mathcal{D}_t = \frac{y_{t+1}-y_t}{\Delta t_t} - f(y_t,t_t,\theta).
$$

2. Crank-Nicolson

$$
\mathcal{D}_t = \frac{y_{t+1}-y_t}{\Delta t_t}
- \frac{1}{2}\left[f(y_t,t_t,\theta)+f(y_{t+1},t_{t+1},\theta)\right].
$$

3. GL1 (midpoint)

$$
\mathcal{D}_t = \frac{y_{t+1}-y_t}{\Delta t_t}
- f\!\left(\frac{y_t+y_{t+1}}{2},\frac{t_t+t_{t+1}}{2},\theta\right).
$$

4. GL2 (2-stage Gauss-Legendre)

Stage states satisfy

$$
Y_i = y_t + \Delta t_t\sum_{j=1}^2 A_{ij} f(Y_j, t_t+c_j\Delta t_t,\theta),\quad i=1,2,
$$

and

$$
\mathcal{D}_t = \frac{y_{t+1}-y_t}{\Delta t_t}
- \frac{1}{2}\left[f(Y_1,t_1,\theta)+f(Y_2,t_2,\theta)\right].
$$

For GL2, derivatives are obtained by differentiating the stage system and solving the induced linear systems, so that `Du_dyc`, `Du_dyn`, `Du_dtheta` are consistent with the actual implicit stage equations.

#### Tracking gradient path (now scheme-consistent)

State sensitivities satisfy

$$
S_{t+1} = -\left(D^{(t)}_{y_{t+1}}\right)^{-1}
\left(D^{(t)}_{y_t}S_t + D^{(t)}_{\theta}\right),
\quad
S_1 = \frac{\partial y_0}{\partial\theta}.
$$

Then

$$
\nabla_\theta H
= \sum_{t,v}\frac{2}{n_s}\,r_{t,v}\,S_{t,v,:},
\quad r_{t,v}=y_{t,v}-y^{obs}_{t,v}.
$$

This recursion is implemented using the linearization API, not Euler-specific formulas.

#### Cascading IFT path (now scheme-consistent)

The inner objective in discrete form is

$$
J_{in}(Y,\theta)
= \frac{1}{n_s}\sum_{t,v} m_{t,v}(y_{t,v}-y^{obs}_{t,v})^2
+ \lambda\sum_{t=1}^{n_s-1} w_t\lVert u_t\rVert^2,
$$

with observation mask $m_{t,v}$ and trapezoidal weights

$$
w_t = \frac{\Delta t_{t-1}+\Delta t_t}{2},\quad \Delta t_0=0.
$$

Differentiating the inner optimality condition gives

$$
A S = -B,
$$

where $A=\partial^2 J_{in}/\partial Y^2$ and $B=\partial^2 J_{in}/(\partial Y\partial\theta)$.

For each interval $t$, with $D_y=[D^{(t)}_{y_t},D^{(t)}_{y_{t+1}}]$:

- $A$ receives the weighted Gauss-Newton block contribution
    $2\lambda w_t D_y^\top D_y$,
- $B$ receives
    $2\lambda w_t (D_y^\top D^{(t)}_{\theta})$
    plus the missing cross term
    $2\lambda w_t (\partial D_y/\partial\theta)^\top u_t$.

Initial-condition sensitivity is imposed as a hard boundary condition,

$$
S_1 = \frac{\partial y_0}{\partial\theta},
$$

by solving the reduced system on free time blocks and injecting the constrained first block afterwards.

The outer gradient is then

$$
\nabla_\theta H
= \sum_{t,v} 2\,r_{t,v}\,S_{t,v,:},
$$

with $H(\theta)=\sum_{t,v}m_{t,v}(y_{t,v}-y^{obs}_{t,v})^2$.

#### Practical implication

The outer sensitivity/gradient code no longer relies on Euler-only approximations or finite-difference fallback shortcuts in the inverse layer. All methods (`euler`, `cn`, `gl1`, `gl2`) use the same low-level linearization interface, ensuring method-consistent derivatives.

## Work In Progress
- One of the goal of the project is to quantify the uncertainty related to the estimate. This is particularly difficult due to the presence of non-linearity. The idea is to consider a linearization around the current estimate, using techniques such as the Laplace approximation and delta method
- Some difficulties that may be encountered in this context: relating the uncertainty in the estimation of the nonparametric term $y$ and the uncertainty related to parameter $\theta$ estimation, which role plays $\lambda$ in UQ, ...
- 