# utl-symbolic-solution-for-the-gradient-of-the-cumulative-bivariate-normal-using-erf-and-sympy
Symbolic solution for the gradient of the cumulative bivariate normal using erf and sympy
    %let pgm=utl-symbolic-solution-for-the-gradient-of-the-cumulative-bivariate-normal-using-erf-and-sympy;

    Symbolic solution for the gradient of the cumulative bivariate normal using erf and sympy

    For a better explanation of the numerical solution see Rick's blog
    https://blogs.sas.com/content/iml/2013/09/20/gradient-of-the-bivariate-normal-cumulative-distribution.html

    github
    https://tinyurl.com/y8574sps
    https://github.com/rogerjdeangelis/utl-symbolic-solution-for-the-gradient-of-the-cumulative-bivariate-normal-using-erf-and-sympy

    The symbolic solution uses the error fumction, erf(x)

    Note the error function erf

    erf(t) = 2/sqrt(pi)*integral(exp(-t**2), t=0 z)
    looks like a CDF

    So the gradient solution is not exacly a close form solution.
    Howere SAS, R and Python provide the error function erf(x) so why
    not use it.

    Sympy can provide general formulation for the gradient of
    the cumulative bivariate normal.


       Four Solutiuons
             1. Symbolic formulation (substitute (x-muX)/sigmaX and (y-muY)/sigmaY
                   for general souion. May be useful for other mutivariate CDFs)
             2. scipy numerical solution
             3. mvtnorm r
             4 approximation r
             5 ehy the gradient is usefull

    This is the final symbolic solution in SAS.
    This matches Rick's numerical solutions.
    The equations are from sympy.

    /*----                                                                   ----*/
    /*---- for a general solution us the translation below                   ----*/
    /*---- use x = (x-muX)/sigmaX and y = (y-muY)/sigmaY                     ----*/
    /*----                                                                   ----*/

    data final;
      length gradx grady $200;
      input x y rho;

      pi=constant('pi');

      dF_dy= sqrt(2)*(1 - erf(sqrt(2)*(rho*y - x)/(2*sqrt(1 - rho**2))))*exp(-y**2/2)/(4*sqrt(pi));
      df_dx =sqrt(2)*(1 - erf(sqrt(2)*(rho*x - y)/(2*sqrt(1 - rho**2))))*exp(-x**2/2)/(4*sqrt(pi));

      gradx = catx(' ','Gradient partial df/dx at (x=',dF_dx,'y=',dF_dy,') with rho=0.5');
      grady = catx(' ','Gradient partial df/dx at (x=',dF_dx,'y=',dF_dy,') with rho=0.5');

      put // gradx / grady;
    cards4;
    1 0 .5
    -1 0 .5
    0 1 .5
    ;;;;
    run;quit;


    /*                       _           _ _                  _       _   _
    / |  ___ _   _ _ __ ___ | |__   ___ | (_) ___   ___  ___ | |_   _| |_(_) ___  _ __
    | | / __| | | | `_ ` _ \| `_ \ / _ \| | |/ __| / __|/ _ \| | | | | __| |/ _ \| `_ \
    | | \__ \ |_| | | | | | | |_) | (_) | | | (__  \__ \ (_) | | |_| | |_| | (_) | | | |
    |_| |___/\__, |_| |_| |_|_.__/ \___/|_|_|\___| |___/\___/|_|\__,_|\__|_|\___/|_| |_|
             |___/

    */

    %utl_pybegin;
    parmcards4;
    import sympy as sp
    from sympy import symbols, exp, pi, erf, diff, sqrt, pprint
    x, y, rho = symbols('x y rho')
    def phi(t):
        return exp(-t**2/2) / sqrt(2*pi)

    def Phi(t):
        return (1 + erf(t/sqrt(2))) / 2

    def F(x, y, rho):
        return Phi(x) * Phi((y - rho*x) / sqrt(1 - rho**2))

    dF_dx = diff(F(x, y, rho), x)
    dF_dy = diff(F(x, y, rho), y)

    dF_dx = sp.simplify(dF_dx)
    dF_dy = sp.simplify(dF_dy)

    expected_dF_dx = phi(x) * Phi((y - rho*x) / sqrt(1 - rho**2))
    expected_dF_dy = phi(y) * Phi((x - rho*y) / sqrt(1 - rho**2))

    print("dF/dy:", sp.simplify(expected_dF_dy))
    print("dF/dx:", sp.simplify(expected_dF_dx))

    pprint(sp.simplify(expected_dF_dy))
    pprint(sp.simplify(expected_dF_dx))
    ;;;;
    %utl_pyend;

    data final;
      length gradx grady $200;
      input x y rho;

      pi=constant('pi');

      dF_dy= sqrt(2)*(1 - erf(sqrt(2)*(rho*y - x)/(2*sqrt(1 - rho**2))))*exp(-y**2/2)/(4*sqrt(pi));
      df_dx =sqrt(2)*(1 - erf(sqrt(2)*(rho*x - y)/(2*sqrt(1 - rho**2))))*exp(-x**2/2)/(4*sqrt(pi));

      grad = catx(' ','df/dx=',dF_dx,'df/dy=',dF_dy,') at x=',x,'y=',y,'rho=',rho);

      put / grad ;
    cards4;
    1 0 .5
    -1 0 .5
    0 1 .5
    ;;;;
    run;quit;


    /*           _               _
      ___  _   _| |_ _ __  _   _| |_
     / _ \| | | | __| `_ \| | | | __|
    | (_) | |_| | |_| |_) | |_| | |_
     \___/ \__,_|\__| .__/ \__,_|\__|
                    |_|
    */

    /**************************************************************************************************************************/
    /*                                                                                                                        */
    /*  PYTHON OUTPUT                                                                                                         */
    /*  -------------                                                                                                         */
    /*                                                                                                                        */
    /*  dF/dy: sqrt(2)*(1 - erf(sqrt(2)*(rho*y - x)/(2*sqrt(1 - rho**2))))*exp(-y**2/2)/(4*sqrt(pi))                          */
    /*  dF/dx: sqrt(2)*(1 - erf(sqrt(2)*(rho*x - y)/(2*sqrt(1 - rho**2))))*exp(-x**2/2)/(4*sqrt(pi))                          */
    /*                                                                                                                        */
    /*                                                                                                                        */
    /*  LETS APPLY THE EQUATIONS IN SAS. RESULT AGREES WITH RICK                                                              */
    /*  --------------------------------------------------------                                                              */
    /*                                                                                                                        */
    /*  df/dx= 0.0681997949 df/dy= 0.3494309345 ) at x= 1 y= 0 rho= 0.5                                                       */
    /*                                                                                                                        */
    /*  df/dx= 0.1737709296 df/dy= 0.0495113459 ) at x= -1 y= 0 rho= 0.5                                                      */
    /*                                                                                                                        */
    /*  df/dx= 0.3494309345 df/dy= 0.0681997949 ) at x= 0 y= 1 rho= 0.5                                                       */
    /*                                                                                                                        */
    /*                                                                                                                        */
    /*   PYTHON OUTPUT                                                                                                        */
    /*                                                  2                                                                     */
    /*                                                -y                                                                      */
    /*                  /       /  ___            \\  ----                                                                    */
    /*              ___ |       |\/ 2 *(rho*y - x)||   2                                                                      */
    /*            \/ 2 *|1 - erf|-----------------||*e                                                                        */
    /*                  |       |      __________ ||                                                                          */
    /*                  |       |     /        2  ||                                                                          */
    /*   dF             \       \ 2*\/  1 - rho   //                                                                          */
    /*   --   =   ----------------------------------------                                                                    */
    /*   dx                           ____                                                                                    */
    /*                            4*\/ pi                                                                                     */
    /*                                                  2                                                                     */
    /*                                                -x                                                                      */
    /*                  /       /  ___            \\  ----                                                                    */
    /*              ___ |       |\/ 2 *(rho*x - y)||   2                                                                      */
    /*            \/ 2 *|1 - erf|-----------------||*e                                                                        */
    /*                  |       |      __________ ||                                                                          */
    /*   dF             |       |     /        2  ||                                                                          */
    /*   --   =         \       \ 2*\/  1 - rho   //                                                                          */
    /*   dy          ----------------------------------------                                                                 */
    /*                             ____                                                                                       */
    /*                            4*\/ pi                                                                                     */
    /*                                                                                                                        */
    /**************************************************************************************************************************/

    /*___             _                                              _           _
    |___ \   ___  ___(_)_ __  _   _  _ __  _   _ _ __ ___   ___ _ __(_) ___ __ _| |
      __) | / __|/ __| | `_ \| | | || `_ \| | | | `_ ` _ \ / _ \ `__| |/ __/ _` | |
     / __/  \__ \ (__| | |_) | |_| || | | | |_| | | | | | |  __/ |  | | (_| (_| | |
    |_____| |___/\___|_| .__/ \__, ||_| |_|\__,_|_| |_| |_|\___|_|  |_|\___\__,_|_|
                       |_|    |___/
    */

    %utl_pybegin;
    parmcards4;
    import numpy as np
    from scipy.stats import norm
    def bivariate_normal_cdf_gradient(x, y, rho):
        pdf_x = norm.pdf(x)
        pdf_y = norm.pdf(y)
        cdf_y_given_x = norm.cdf((y - rho * x) / np.sqrt(1 - rho**2))
        cdf_x_given_y = norm.cdf((x - rho * y) / np.sqrt(1 - rho**2))
        dF_dx = pdf_x * cdf_y_given_x
        dF_dy = pdf_y * cdf_x_given_y
        return dF_dx, dF_dy
    rho = .5
    x, y = 0, 0
    gradient = bivariate_normal_cdf_gradient(x, y, rho)
    print(f"\n Gradient at (x={x}, y={y}) with rho={rho}:")
    print(f" dF/dx = {gradient[0]:.6f}")
    print(f" dF/dy = {gradient[1]:.6f}")

    x, y = -1, 0
    gradient = bivariate_normal_cdf_gradient(x, y, rho)
    print(f"\n Gradient at (x={x}, y={y}) with rho={rho}:")
    print(f" dF/dx = {gradient[0]:.6f}")
    print(f" dF/dy = {gradient[1]:.6f}")

    x, y = 0, 1
    gradient = bivariate_normal_cdf_gradient(x, y, rho)
    print(f"\n Gradient at (x={x}, y={y}) with rho={rho}:")
    print(f" dF/dx = {gradient[0]:.6f}")
    print(f" dF/dy = {gradient[1]:.6f}" )
    ;;;;
    %utl_pyend;

    /**************************************************************************************************************************/
    /*                                                                                                                        */
    /*  Gradient at (x=0, y=0) with rho=0.5:                                                                                  */
    /*  dF/dx = 0.199471                                                                                                      */
    /*  dF/dy = 0.199471                                                                                                      */
    /*                                                                                                                        */
    /*  Gradient at (x=-1, y=0) with rho=0.5:                                                                                 */
    /*  dF/dx = 0.173771                                                                                                      */
    /*  dF/dy = 0.049511                                                                                                      */
    /*                                                                                                                        */
    /*  Gradient at (x=0, y=1) with rho=0.5:                                                                                  */
    /*  dF/dx = 0.349431                                                                                                      */
    /*  dF/dy = 0.068200                                                                                                      */
    /*                                                                                                                        */
    /**************************************************************************************************************************/

    /*____                  _
    |___ /   _ __ _____   _| |_ _ __   ___  _ __ _ __ ___    _ __
      |_ \  | `_ ` _ \ \ / / __| `_ \ / _ \| `__| `_ ` _ \  | `__|
     ___) | | | | | | \ V /| |_| | | | (_) | |  | | | | | | | |
    |____/  |_| |_| |_|\_/  \__|_| |_|\___/|_|  |_| |_| |_| |_|

    */

    %utl_rbeginx;
    parmcards4;
    library(mvtnorm)
    # Function to calculate the gradient of bivariate normal CDF
    grad_bvn_cdf <- function(x, y, rho) {
      # Standard normal density function (PDF)
      phi <- function(z) dnorm(z)

      # Standard normal cumulative distribution function (CDF)
      Phi <- function(z, mean, sd) pnorm(z, mean = mean, sd = sd)

      # Calculate partial derivatives
      dF_dx <- phi(x) * Phi(y, mean = rho * x, sd = sqrt(1 - rho^2))
      dF_dy <- phi(y) * Phi(x, mean = rho * y, sd = sqrt(1 - rho^2))
      df<-as.data.frame(c(dF_dx = dF_dx, dF_dy = dF_dy))
      return(df)
    }
    # Set parameters
    grad1 <- grad_bvn_cdf(0, 0, .5)
    grad2 <- grad_bvn_cdf(-1, 0, .5)
    grad3 <- grad_bvn_cdf(0, 1, .5)
    res<-rbind(grad1,grad2,grad3)
    res;
    ;;;;
    %utl_rendx;


    /**************************************************************************************************************************/
    /*                                                                                                                        */
    /* dF_dx      0.19947114                                                                                                  */
    /* dF_dy      0.19947114                                                                                                  */
    /*                                                                                                                        */
    /* dF_dx1     0.17377093                                                                                                  */
    /* dF_dy1     0.04951135                                                                                                  */
    /*                                                                                                                        */
    /* dF_dx2     0.34943093                                                                                                  */
    /* dF_dy2     0.06819979                                                                                                  */
    /*                                                                                                                        */
    /**************************************************************************************************************************/

    /*  _                                     _                 _   _
    | || |     __ _ _ __  _ __  _ __ _____  _(_)_ __ ___   __ _| |_(_) ___  _ __   _ __
    | || |_   / _` | `_ \| `_ \| `__/ _ \ \/ / | `_ ` _ \ / _` | __| |/ _ \| `_ \ | `__|
    |__   _| | (_| | |_) | |_) | | | (_) >  <| | | | | | | (_| | |_| | (_) | | | || |
       |_|    \__,_| .__/| .__/|_|  \___/_/\_\_|_| |_| |_|\__,_|\__|_|\___/|_| |_||_|
                   |_|   |_|
    */
    %utl_rbeginx;
    parmcards4;
    library(mvtnorm)
    eps <- 1e-6
    rho = .5
    x=1;
    y=0;
    F <- function(x, y) pmvnorm(lower = c(-Inf, -Inf), upper = c(x, y),
                                mean = c(0, 0), sigma = matrix(c(1, rho, rho, 1), 2, 2))
    numerical_grad_x <- (F(x + eps, y) - F(x - eps, y)) / (2 * eps)
    numerical_grad_y <- (F(x, y + eps) - F(x, y - eps)) / (2 * eps)
    print(c(x=x,y=y,rho=rho,numerical_dF_dx = numerical_grad_x, numerical_dF_dy = numerical_grad_y))
    ;;;;
    %utl_rendx;

    /**************************************************************************************************************************/
    /*                                                                                                                        */
    /* Matches Rick                                                                                                           */
    /*                                                                                                                        */
    /* x   y   rho  numerical_dF_dx numerical_dF_dy                                                                           */
    /* 1   0    .5  0.06819979      0.34943093                                                                                */
    /*                                                                                                                        */
    /**************************************************************************************************************************/

    /*___             _             _   _
    | ___|  __      _| |__  _   _  | |_| |__   ___
    |___ \  \ \ /\ / / `_ \| | | | | __| `_ \ / _ \
     ___) |  \ V  V /| | | | |_| | | |_| | | |  __/
    |____/    \_/\_/ |_| |_|\__, |  \__|_| |_|\___|
                            |___/
                         _ _            _
      __ _ _ __ __ _  __| (_) ___ _ __ | |_
     / _` | `__/ _` |/ _` | |/ _ \ `_ \| __|
    | (_| | | | (_| | (_| | |  __/ | | | |_
     \__, |_|  \__,_|\__,_|_|\___|_| |_|\__|
     |___/
    */
    The gradient of the bivariate normal cumulative distribution function (CDF) has several practical
    applications in statistics, finance, and other fields:

    1. Sensitivity analysis: The gradient provides information about how sensitive the cumulative
    probability is to changes in the input variables. This is useful for understanding which
    factors have the most significant impact on the probability[1].

    2. Optimization problems: In financial modeling and portfolio optimization, the gradient can
    be used to find optimal solutions that maximize or minimize certain objectives
    related to bivariate normal distributions[1].

    3. Confidence regions: The gradient can be used to calculate confidence regions for
    bivariate normal parameters, which is important in statistical inference and hypothesis testing[1].

    4. Monte Carlo simulations: The gradient can be used to improve the efficiency of
    Monte Carlo simulations involving bivariate normal distributions,
    particularly in importance sampling techniques[1].

    5. Risk assessment: In finance and insurance, the gradient can help quantify how changes
    in correlated risk factors affect the overall risk profile of a portfolio or insurance policy[1].

    6. Statistical modeling: The gradient is useful in fitting bivariate normal models to data,
    especially in maximum likelihood estimation procedures[1].

    7. Machine learning: In some machine learning algorithms, particularly those involving Gaussian
    processes or probabilistic models, the gradient of the bivariate normal
     CDF can be used to optimize model parameters[1].

    8. Quality control: In manufacturing processes where two correlated variables need
    to be monitored simultaneously, the gradient can help in designing
    efficient control charts and tolerance regions[1].

    9. Copula theory: The gradient is useful in copula-based modeling, which is widely
    used in finance for modeling dependence structures between random variables[1].

    10. Numerical methods: The gradient can be used to implement numerical methods for
    bivariate normal probability calculations, improving the accuracy and efficiency of
    computational algorithms[1].

    These applications demonstrate the importance of the gradient of the bivariate normal
    CDF in various fields, particularly where understanding the relationship between
    correlated variables is crucial for decision-making or analysis.

    Citations:
    [1] https://blogs.sas.com/content/iml/2013/09/20/gradient-of-the-bivariate-normal-cumulative-distribution.html
    [2] https://blogs.sas.com/content/iml/2013/09/25/compute-contours-of-the-bivariate-normal-cdf.html
    [3] https://github.com/tensorflow/probability/issues/422
    [4] https://stats.stackexchange.com/questions/71976/partial-derivative-of-bivariate-normal-cdf-and-pdf
    [5] https://math.stackexchange.com/questions/4373566/derivatives-of-multivariate-normal-cdf-in-terms-of-lower-dimension-pdfs-and-cdfs


     /*              _
      ___ _ __   __| |
     / _ \ `_ \ / _` |
    |  __/ | | | (_| |
     \___|_| |_|\__,_|

    */
