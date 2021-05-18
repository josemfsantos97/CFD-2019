# Finite Differences and Finite Volume Methods
The first assignment is divided into two parts. The first part consisted of an implementation of a finite difference approximation and to get insights about sources of error. In the second part the linearized potential flow equations and the small perturbations approach was solved with the Finite Volume Methods.

### Finite Differences - Implementation and error decay analysis
The first derivative of a function is approximated using a fourth-order central difference scheme and a second-order forward difference scheme.

![Analytical Error and Truncation Error](https://github.com/josemfsantos97/CFD-2019/blob/main/Homework01/images/finite_diff.png)

### Finite Volume Methods
The potential equation is obtained from the assumption of isentropic flow. In transonic and supersonic, shock waves are present hence, realistically, flow is not isentropic. Moreover, the small perturbation theory is applied. This renders a linear subsonic and supersonic equation and a non-linear transonic
equation, respectively shown below.

![Linear Subsonic and Supersonic](https://github.com/josemfsantos97/CFD-2019/blob/main/Homework01/images/form1.jpg)


![Non-Linear Transonic](https://github.com/josemfsantos97/CFD-2019/blob/main/Homework01/images/form2.jpg)


#### Subsonic regime - visualizing results

![Mach Subsonic](https://github.com/josemfsantos97/CFD-2019/blob/main/Homework01/images/machsub.jpg)

#### Supersonic regime - visualizing results

![Mach Subsonic](https://github.com/josemfsantos97/CFD-2019/blob/main/Homework01/images/machsup.jpg)

