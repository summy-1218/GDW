## 2 Wake Modeling
The user has two different options for calculating the effect of the wake on the turbine rotor aerodynamics: either the classic blade element momentum theory or the more recently developed generalized dynamic wake model. Each model is explained in more detail in the following sections.

### 2.1 Blade Element Momentum
Blade element momentum (BEM) theory is one of the oldest and most commonly used methods for calculating induced velocities on wind turbine blades. This theory is an extension of actuator disk theory, first proposed by the pioneering propeller work of Rankine and Froude in the late 19th century. The BEM theory, generally attributed to Betz and Glauert (1935), actually originates from two different theories: blade element theory and momentum theory (see Leishman 2000). Blade element theory assumes that blades can be divided into small elements that act independently of surrounding elements and operate aerodynamically as two-dimensional airfoils whose aerodynamic forces can be calculated based on the local flow conditions. These elemental forces are summed along the span of the blade to calculate the total forces and moments exerted on the turbine. The other half of BEM, the momentum theory, assumes that the loss of pressure or momentum in the rotor plane is caused by the work done by the airflow passing through the rotor plane on the blade elements. Using the momentum theory, one can calculate the induced velocities from the momentum lost in the flow in the axial and tangential directions. These induced velocities affect the inflow in the rotor plane and therefore also affect the forces calculated by blade element theory. This coupling of two theories ties together blade element momentum theory and sets up an iterative process to determine the aerodynamic forces and also the induced velocities near the rotor.

In practice, BEM theory is implemented by breaking the blades of a wind turbine into many elements along the span. As these elements rotate in the rotor plane, they trace out annular regions, shown in Figure 1, across which the momentum balance takes place. These annular regions are also where the induced velocities from the wake change the local flow velocity at the rotor plane. BEM can also be used to analyze stream tubes through the rotor disk, which can be smaller than the annular regions and provide more computational fidelity. However, as currently written, AeroDyn only allows analysis using annular regions.

Because of its simplicity, BEM theory does have its limitations. One primary assumption is that the calculations are static; it is assumed that the airflow field around the airfoil is always in equilibrium and that the passing flow accelerates instantaneously to adjust to the changes in the vorticity in the wake. In practice, it has been shown that the airfoil response takes time to adjust to a changing wake resulting from new inflow or turbine operating conditions (Snel and Schepers 1995). In order to model this time lag effect correctly, we recommended that the user utilize the generalized dynamic wake model described below. One other limitation is that BEM theory breaks down when the blades experience large deflections out of the rotor plane. Because the theory assumes that momentum is balanced in a plane parallel to the rotor, any deflections of the rotor will lead to errors in the aerodynamic modeling. Another limitation of BEM theory comes from blade element theory. This theory is based on the assumption that the forces acting on the blade element are essentially two-dimensional, meaning that spanwise flow is neglected. This assumption also implies that there is very little spanwise pressure variation (which would create spanwise flow), and the theory is therefore less accurate for heavily loaded rotors with large pressure gradients across the span. Some other limitations of the original theory include no modeling of tip or hub vortex influence on the induced velocities and an inability to account for skewed inflow. However, corrections to the original theory have provided some methods to model these aerodynamic effects and will be explained in more detail below. In spite of the limitations listed above, BEM theory has been used widely as a reliable model for calculating the induced velocity and elemental forces on wind turbine blades, and it has been retained as a useful model in AeroDyn.

The advantage of the BEM theory is that each blade element is modeled as a two-dimensional airfoil. Figure 2 is an example of an airfoil with the velocities and angles that determine the forces on the element and also the induced velocities from the wake influence. Figure 3 shows the resultant aerodynamic forces on the element and their components perpendicular and parallel to the rotor plane. These are the forces that dictate the thrust (perpendicular) and torque (parallel) of the rotor, which are the dominant forces for turbine design. In Figure 3, the angle relating the lift and drag of the airfoil element to the thrust and torque forces is the local inflow angle, $\varphi$ (or $\phi$ in the figures). As shown in Figure 2, this inflow angle is the sum of the local pitch angle of the blade, $\beta$ , and the angle of attack, $\alpha$ . The local pitch angle is dependent on the static blade geometry, elastic deflections, and the active or passive blade pitch control system. The angle of attack is a function of the local velocity vector, which is in turn constrained by the incoming local wind speed, rotor speed, blade element velocities and induced velocities. Note in Figure 2 that the velocities of the element from blade deflections ( $v_{e-o p}$ and $v_{e-i p}$ ) affect the inflow angle and angle of attack, but are not directly affected by the induced velocities from the wake. This assumption is consistent with momentum theory, but it might not be the appropriate physical model for the element-wake coupling.

Because we are required to obtain the angle of attack to determine the aerodynamic forces on an element, we must first determine the inflow angle based on the two components of the local velocity vector. Assuming that the blade motion is very small, the resulting equation is dependent on the induced velocities in both the axial and tangential directions as well as the local tip speed ratio:
$$
tan \varphi=\frac{U_{\infty}(1-a)}{\Omega r\left(1+a'\right)}=\frac{1-a}{\left(1+a'\right) \lambda_{r}} . \tag{1}
$$

However, if the blade motion is significant we must include the local velocities in the calculation of the inflow angle, as follows:
$$
tan \varphi=\frac{U_{\infty}(1-a)+v_{e-op}}{\Omega r\left(1+a'\right)+v_{e-ip}} . \tag{2}
$$

This equation holds for all elements of the blade along the span, although typically the inflow angle changes with element location.

The induced velocity components in Equations 1 and 2 are a function of the forces on the blades and we use BEM theory to calculate them. A thorough derivation of these equations can be found in most wind turbine design handbooks (Manwell et al. 2002; Burton et al. 2001), and so it will only be summarized here. From blade element theory and Figure 3, the thrust distributed around an annulus of width $dr$ (see Figure 1) is equivalent to
$$
dT=B \frac{1}{2} \rho V_{total }^{2}\left(C_{l} cos \varphi+C_{d} sin \varphi\right) c d r, \tag{3}
$$
and the torque produced by the blade elements in the annulus is equivalent to
$$
d Q=B \frac{1}{2} \rho V_{total }^{2}\left(C_{l} sin \varphi-C_{d} cos \varphi\right) c r d r . \tag{4}
$$

Now, to relate the induced velocities in the rotor plane to the elemental forces of Equations 3 and 4 we must incorporate the momentum part of the theory, which states that the thrust extracted by each rotor annulus is equivalent to
$$
dT=4 \pi r \rho U_{\infty}^{2}(1-a) a d r, \quad \tag{5}
$$
and the torque extracted from each annular section is equivalent to
$$
dQ=4 \pi r^{3} \rho U_{\infty} \Omega(1-a) a' d r . \tag{6}
$$

Thus, when we include two-dimensional airfoil tables of lift and drag coefficient as a function of the angle of attack, $\alpha$ , we have a set of equations that can be iteratively solved for the induced velocities and the forces on each blade element. However, before we solve our system of equations, we would like to take into account several corrections to the BEM theory. These corrections include tip- and hub-loss models to account for vortices shed at these locations, the Glauert correction to account for large induced velocities ($a >0.4$), and the skewed wake correction to model the effects of incoming flow that is not perpendicular to the rotor plane. Each of these will be described in a section below.

Note that Equations 1-6 do not include terms for coning angle or teeter angle of the rotor plane. AeroDyn assumes that the built-in coning, effective coning of the rotor blades from large aeroelastic deflections, and teeter do not significantly change the aerodynamics of the rotor in operation. This assumption is tenuous, particularly for large deflection angles that will change the shape of the wake by introducing an effective skew angle. Because of this concern, the effective coning and teeter may be introduced in a future version of the code.

#### 2.1.1 Tip-Loss Model
One of the major limitations of the original blade element momentum theory is that there is no influence of vortices shed from the blade tips into the wake on the induced velocity field. These tip vortices create multiple helical structures in the wake, as seen in Figure 4, and they play a major role in the induced velocity distribution at the rotor. The effect on induced velocity in the rotor plane is most pronounced near the tips of the blades, an area that also has the greatest influence on the power produced by the turbine. To compensate for this deficiency in BEM theory, AeroDyn uses a theory originally developed by Prandtl (see Glauert 1935). Prandtl simplified the wake of the turbine by modeling the helical vortex wake pattern as vortex sheets that are convected by the mean flow and have no direct effect on the wake itself. This theory is summarized by a correction factor to the induced velocity field, $F$ , and can be expressed simply by the following:
$$
F=\frac{2}{\pi} cos ^{-1} e^{-f}, \tag{7}
$$
where,
$$
f=\frac{B}{2} \frac{R-r}{r sin \varphi} . \quad \tag{8}
$$
<sup>1</sup> Prandtl’s original tip-loss model was based on the inflow angle at the blade tip, $\varphi_{tip }$ . He later revised this (Glauert 1935) to the local inflow angle, making the calculation easier to implement with minimal loss in accuracy.

This correction factor is used to modify the momentum part of the blade element momentum equations, replacing Equations 5 and 6 with the following:
$$
d T=4 \pi r \rho U_{\infty}^{2}(1-a) a F d r \quad \tag{9}
$$
$$
d Q=4 \pi r^{3} \rho U_{\infty} \Omega(1-a) a' F d r . \quad \tag{10}
$$

Because of its reasonable accuracy for most operating conditions and easy formulaic implementation, the Prandtl model is often used in engineering codes such as AeroDyn. However, like most engineering models it has limitations that affect its accuracy. One limitation of this model is that it assumes the wake does not expand, limiting its validity to lightly loaded rotors. Also, Glauert (1935) showed that the accuracy, relative to the more accurate and computationally expensive Goldstein solution (1929), of this model decreases with lower numbers of blades (less than three) and higher tip speed ratios.

Figure 5 is an example of the radial distribution of the tip-loss correction for a blade that is operating such that the inflow angle, $\varphi$ , is constant along the span at 10°. When the tip-loss model is employed, the tip-loss factor sharply decreases as the radial position along the blade approaches the blade tip. This corresponds to a dramatic increase in the induction factor near the tip. As the induction factor increases, the resultant relative wind speed for a given blade segment decreases along with the angle of attack (see Figure 2). As a result, the loading (lift and drag forces) decreases near the tip.

In addition to the Prandtl model, users of AeroDyn also have the option of using an empirical relationship for the tip loss based on the Navier-Stokes solutions of Xu and Sankar (2002) as described in the following equations:
$$
F_{new }=\frac{F_{Prandl }^{0.85}+0.5}{2} \quad for \ 0.7 \leq r / R \leq 1,
$$
or
$$
F_{new }=1-\left(\frac{r}{R}\right) \frac{1-F_{Prandl (r / R=0.7)}}{0.7} \quad for \ r / R<0.7 .
$$

These relationships are a correction for the Prandtl model and must be used in conjunction with Equations 7 and 8. Note, however, that this correction was based on a specific turbine design (UAE Phase 6, Hand et al. 2001) at one wind speed and may not be applicable to all turbine configurations. It also results in a tip-loss factor greater than zero at the tip, which is physically unrealistic at the tip blade station.

#### 2.1.2 Hub-Loss Model
Much like the tip-loss model, the hub-loss model serves to correct the induced velocity resulting from a vortex being shed near the hub of the rotor. The hub-loss model uses a nearly identical implementation of the Prandtl tip-loss model to describe the effect of this vortex, replacing Equation 8 with the following:
$$
f=\frac{B}{2} \frac{r-R_{hub}}{r sin \varphi} . \tag{13}
$$

For a given element, the local aerodynamics may be affected by both the tip loss and hub loss, in which case the tip-loss and hub-loss correction factors are multiplied to create the total loss factor used in Equations 9 and 10.

#### 2.1.3 Glauert Correction
Another limitation of the BEM theory is that when the induction factor is greater than about 0.4, the basic theory becomes invalid. This occurs with turbines operating at high tip speed ratios (e.g. constant speed turbine at low wind speeds), as the rotor enters what is known as the turbulent wake state ($a>0.5$). According to momentum theory, this operating state results when some of the flow in the far wake starts to propagate upstream, which is a violation of the basic assumptions of BEM theory. Physically, this flow reversal cannot occur, and what actually happens is more flow entrains from outside the wake and the turbulence increases. The flow behind the rotor slows down, but the thrust on the rotor disk continues to increase. To compensate for this effect, Glauert (1926) developed a correction to the rotor thrust coefficient based on experimental measurements of helicopter rotors with large induced velocities. While this model was originally developed as a correction to the thrust coefficient of an entire rotor, it has also been used to correct the local coefficient of the individual blade elements when used with BEM theory. Because of this, it is important to understand the Glauert correction's relationship to the tip-loss model. When the losses near the tip are high, the induced velocities are large; therefore, the possibility of a turbulent wake near the tips increases. Thus, for each element the total induced velocity calculation must use a combination of the tip-loss and Glauert corrections. Buhl (2004) derived a modification to the Glauert empirical relation that included the tip-loss correction as follows:
$$
C_T=\frac{8}{9}+\left(4 F-\frac{40}{9}\right) a+\left(\frac{50}{9}-4 F\right) a^{2}, \quad \tag{14}
$$
or, solving for the induction factor,
$$
a=\frac{18 F-20-3 \sqrt{C_{T}(50-36 F)+12 F(3 F-4)}}{36 F-50} . \tag{15}
$$

This empirical relationship is different from those in the models of other authors (Manwell 2002; Burton 2001). But, this relationship is necessary to eliminate a numerical instability when using the Glauert correction to calculate the elemental thrust in conjunction with the tip-loss correction model.

Figure 6 shows an example of the Glauert correction when the tip-loss factor is equal to one. When the induction factor, $a$ , is 0.4 the BEM theory and Glauert correction produce the same value for thrust coefficient of 0.96. The slopes are also equivalent at this induction factor. When the tip-loss factor is less than one (e.g. 0.75, as in Figure 7), the BEM theory predicts a much lower thrust coefficient for most induction factors. Thus, to prevent numerical instability in AeroDyn, the Glauert correction must also adjust so that the value and slopes again match at the induction factor of 0.4. This figure demonstrates the sensitivity of the induction factor to the tiploss factor seen in Equation 15.

Again, note that the Glauert correction was developed as a correction to an entire rotor disk; the original researchers did not intend it to be applied to a rotor annulus. However, because of a limited amount of experimental data, an alternative model with BEM theory does not currently exist.

#### 2.1.4 Skewed Wake Correction
Another disadvantage of blade element momentum theory is that it was originally designed for axisymmetric flow. Often, however, wind turbines operate at yaw angles relative to the incoming wind, which produces a skewed wake behind the rotor. The BEM model needs to be corrected to account for this skewed wake effect. The formulation used in AeroDyn is based on an equation originally developed by Glauert (1926) who was primarily interested in the autogyro. The basic formula of the skewed wake correction he derived is
$$
a_{skew }=a\left[1+K \frac{r}{R} cos \psi\right], \tag{16}
$$
where the constant K is a function of the skew angle.

Many skewed wake correction models are derived from this formulation. The one implemented in AeroDyn is based on a method developed by Pitt and Peters (1981) (see also Snel and Schepers 1995). Assuming steady inflow conditions, the skewed wake formulation is
$$
a_{skew }=a\left[1+\frac{15 \pi}{32} \frac{r}{R} tan \frac{\chi}{2} cos \psi\right], \tag{17}
$$
where $\psi$ is defined as the azimuth angle that is zero at the most downwind position of the rotor plane, after accounting for both tilt and yaw (see Figure 8). This position has the greatest amount of induced velocity, whereas the most upwind position ( $cos \psi=-1$ ) has the least induced velocity.

Although Glauert's model originally assumed a was the induced velocity for the entire rotor, AeroDyn uses this correction for each element. Therefore, a and $a_{skew }$ apply to the local elemental induced velocities.

Notice that the constant K in Equation 16 is a function of the wake skew angle, $\chi$ , rather than the rotor yaw angle, $\gamma.$ The wake angle is the actual flow angle leaving the turbine and is slightly larger than the skew angle, which is defined as the difference between the incoming flow and rotor plane, as seen in Figure 8. Using the analysis of Coleman et al. (1945), we can relate the wake skew angle to the yaw angle in the following formula:
$$
tan \chi=\frac{U_{\infty}\left(sin \gamma-a tan \frac{\chi}{2}\right)}{U_{\infty}(cos \gamma-a)}, \tag{18}
$$
which can be approximated by the following relationship according to Burton et al. (2001):
$$
\chi =(0.6a+1)\gamma . \tag{19}
$$

As with previous models, the skewed wake correction has limitations. The major limitation of this model is that it assumes a cylindrical wake, which is valid only for lightly loaded rotors. Also, there is no firm theoretical basis for applying this correction to BEM theory, even though it has improved results for predicting yaw moments and motions when compared with those of the standard theory (Hansen 1992; Snel and Schepers 1995). In addition, recent research has found that this correction may be too large in some situations (see Eggers et al. 2000). We believe that, for wind turbines consistently operating in yaw, the generalized dynamic wake model described below is a better alternative for making more accurate predictions of the turbine aerodynamics.

#### 2.1.5 Other Corrections
Previous researchers (Wilson and Patton 1978) have suggested various other corrections to the BEM theory. These corrections include accounting for the blade thickness effect on local angle of attack, cascade width for high solidity turbines, and spanwise gaps for partial span pitch control. Blade thickness and cascade effects can be aerodynamically significant near the rotor hub and may affect the in-plane yaw forces on the rotor. At this time, AeroDyn does not model these effects, but future research may necessitate adding these corrections at some later time. Spanwise gaps are not modeled in AeroDyn because partial span pitch control is not used in most modern turbine designs.

#### 2.1.6 Final Iteration Procedure for Blade Element Momentum Theory
Now that all of the equations for BEM theory have been established, we will identify the iteration procedure used in AeroDyn to calculate the induced velocities, angles of attack, and thrust coefficients for each blade element along the span of a blade. To begin the calculation we must first estimate the axial induction factor. One efficient way to do this is to use Equations 21 and 27 below, assume that the inflow angle $\varphi$ is small ($sin \varphi ≈\varphi$), the tangential induction $a'$ is zero, the tip and hub-loss corrections F are one, the drag coefficient $C_{d}$ is zero, the lift coefficient, $C_{l}=2 \pi \alpha$ , and finally, $\alpha=\varphi-\beta$ . After some rearranging, we arrive at the initial estimate of the axial induction factor:
$$
a=\frac{1}{4}\left[2+\pi \lambda_{r} \sigma'-\sqrt{4-4 \pi \lambda_{r} \sigma'+\pi \lambda_{r}^{2} \sigma'\left(8 \beta+\pi \sigma'\right)}\right] . \tag{20}
$$

From here we can estimate the inflow angle using an initial assumption of zero for the tangential induction, $a'$ and
$$
tan \varphi=\frac{U_{\infty}(1-a)+v_{e-op}}{\Omega r\left(1+a'\right)+v_{e-ip}} . \tag{21}
$$

Next, AeroDyn determines the thrust coefficient for the element using the following:
$$
C_{T}=\left[1+\frac{\sigma'(1-a)^{2}\left(C_{l} cos \varphi+C_{d} sin \varphi\right)}{sin ^{2} \varphi}\right] . \tag{22}
$$

Then, the tip- and hub-loss corrections are calculated as follows:
$$
F_{tip }=\frac{2}{\pi} cos ^{-1} e^{-\frac{B}{2} \frac{R-r}{r sin \varphi}} \tag{23}
$$
$$
F_{hub }=\frac{2}{\pi} cos ^{-1} e^{-\frac{B}{2}\frac{r-R_{hub}}{r sin \varphi}} \tag{24}
$$
$$
F=F_{hub} F_{tip} . \tag{25}
$$

Now, if $C_{T}>0.96 ~F$, the element is highly loaded and the modified Glauert correction will be used to determine the new axial induction factor:
$$
a=\frac{18 F-20-3 \sqrt{C_{T}(50-36 F)+12 F(3 F-4)}}{36 F-50} . \tag{26}
$$

If $C_{T} ≤0.96 ~F$, the standard BEM theory is used to calculate the axial induction:
$$
a=\left[1+\frac{4 F sin ^{2} \varphi}{\sigma'\left(C_{l} cos \varphi+C_{d} sin \varphi\right)}\right]^{-1} . \tag{27}
$$

The tangential induction factor is calculated using
$$
a'=\left[-1+\frac{4 F sin \varphi cos \varphi}{\sigma'\left(C_{l} sin \varphi-C_{d} cos \varphi\right)}\right]^{-1} . \tag{28}
$$

And finally, the effect of skew is included using the skewed wake correction factor:
$$
a_{skew }=a\left[1+\frac{15 \pi}{32} \frac{r}{R} tan \frac{\chi}{2} cos \psi\right] . \tag{29}
$$

This process is then repeated for each element, starting again with Equation 21 and iterated until the values of induction factors and inflow angle have converged to their final values.

In AeroDyn, the user has some control over how the induced velocities are calculated. Four of the options in the calculation routine are (1) whether to include drag in the induction calculations (axial or tangential), (2) whether to include tip losses, (3) whether to include hub losses, and (4) whether to calculate rotational induction. If the user chooses not to include drag in the induction calculation (as recommended by Wilson and Lissaman 1974) the $C_{d}$ term in the above equations is set to zero. Similarly, if the tangential induction is neglected, AeroDyn will ignore Equation 28 and assume all induction is in the axial direction. Finally, if both tip and hub losses are ignored, the parameter F will be set to one for all of the above equations. If tip losses are desired but not hub losses, the parameter F will be calculated only near the tip, and likewise when only hub losses are modeled.

Currently in AeroDyn, these equations are not directly coupled with the dynamic stall routines explained below. In this iteration process, only static coefficients of lift and drag are used to calculate the properties of the wake. Once all of the induced velocities and angles of attack are calculated, the dynamic stall routines are called and the forces returned to the aeroelastic code are calculated. This decoupling was assumed for two reasons. First, the wake cannot fully respond to short-term dynamic stall events, so they should not always affect the wake. Second, the code is faster and simpler without this coupling. However, some dynamic stall events such as those due to persistent skewed flow can affect the entire wake, so this assumption is not always valid. The decoupling should be refined after future research. This is another reason that the blade element momentum method is not preferred for unsteady or highly skewed flows.

### 2.2 Generalized Dynamic Wake
The generalized dynamic wake (GDW) model of AeroDyn is based on the work of Peters and He (1989) and was implemented in the code by Suzuki (2000) for his Ph.D. thesis at the University of Utah. This model was originally developed for the helicopter industry, and it is also known as the acceleration potential method. An advantage of this method is that it allows for a more general distribution of pressure across a rotor plane than BEM theory. It is an extension of the often-used Pitt and Peters (1981) model, with more flow states and a fully nonlinear implementation to account for turbulence and spatial variation of the inflow.

The GDW method is based on a potential flow solution to Laplace's equation. Kinner (1937) used this solution to develop the equations for the pressure distributions in the rotor plane, which consist of an infinite series of Legendre functions in the radial direction and trigonometric functions in the azimuthal direction. In his derivation, Kinner started from the Euler equations (inviscid and incompressible flow), assumed that the induced velocities were small in comparison to the mean wind speed and regarded the rotor as an infinite number of slender blades, to keep the solidity low.

The main advantages of the generalized dynamic wake method over blade element momentum theory include inherent modeling of the dynamic wake effect, tip losses, and skewed wake aerodynamics. The dynamic wake effect is the time lag in the induced velocities created by vorticity being shed from the blades and being convected downstream. Figure 9 contains an example of the time lag effect on the power output of a turbine. It shows the measured power output of the Tjæreborg turbine (Suzuki 2000) operating in a 10.6 m/s mean wind. During a 60 second period, the turbine blades were pitched from 0.2° to 3.9° over a period of 1.0 seconds and then back to 0.2° again. After each blade pitch change, the power changes quickly resulting in an overshoot and then gradually returns to a new equilibrium value over the next 25+ seconds. Also in the figure are predictions of this event using both BEM and GDW theory. Notice that the BEM theory has no time lag, but the GDW does. However, the time constant is shorter than that exhibited by the data in the figure. The small oscillations in the BEM prediction are due to structural vibrations and not the aerodynamic model itself. We estimate that similar amplitude structural oscillations also appear in the GDW model. Therefore, most of the overshoot seen in Figure 9 is due to the aerodynamic time lag predicted by the model. This time lag is longer than that experienced by helicopters and was calibrated to typical wind turbine response times using experimental data in Suzuki's thesis (2000).

Another advantage of this method is that the induced velocities in the rotor plane are determined from a set of first-order differential equations, which can be solved using a non-iterative technique. The technique used in AeroDyn is the fourth-order Adams-Bashford-Moulton (Press et al. 1982) predictor-corrector method. Because iteration is not required, the model can also be directly incorporated with a dynamic stall model for determining the aerodynamic coefficients of each blade element. Although, as currently written, AeroDyn determines the dynamic stall effects after the GDW equations have been solved.

Like the BEM theory, the GDW method has its limitations. As with most wake models, the generalized dynamic wake was developed for lightly loaded rotors and assumes that the induced velocities are small relative to the mean flow. This basic assumption leads to instability of the method at low wind speeds when the turbulent wake state is approached (Laino and Hansen, 2004). To avoid this computational instability, AeroDyn currently switches to the BEM method when the mean wind speed is below 8 m/s. Another disadvantage of the model is that it does not account for wake rotation. To correct for this, AeroDyn uses the BEM equation to calculate the tangential induction factor, as in Equation 28. Finally, the GDW method assumes that the rotor plane is a flat disk. Therefore, the effect of large aeroelastic deflections or significant coning of the rotor blades on the wake aerodynamics will not be accurately modeled.

#### 2.2.1 Basic Derivation
The basic governing equations of the generalized dynamic wake are derived from the Euler equations. Assuming that the induced velocities are small perturbations relative to the freestream inflow, conservation of momentum simplifies to
$$
\frac{\partial u_{i}}{\partial t}+U_{\infty j} \frac{\partial u_{i}}{\partial x_{j}}=-\frac{1}{\rho} \frac{\partial p}{\partial x_{i}}, \tag{30}
$$
and conservation of mass resulting in
$$
\frac{\partial u_{i}}{\partial x_{i}}=0, \tag{31}
$$
finally leading to Laplace's equation for the pressure distribution:
$$
\nabla^2 p=0 . \quad \tag{32}
$$

It is convenient to non-dimensionalize these equations with the rotor tip speed, which is a widely used convention in rotorcraft aerodynamics, and also the hub-height wind speed, which is common in wind turbine aerodynamics. This results in the following nondimensional quantities:
$$
time: \hat{t}=\Omega t \quad \tag{33}
$$
$$
displacements: \hat{x}_{i}=\frac{x_{i}}{R} \quad \tag{34}
$$
$$
velocities: \hat{u}_{i}=\frac{u_{i}}{\Omega R} \ and \ \hat{U}_{\infty}=\frac{U_{\infty}}{\Omega R} \quad \tag{35}
$$
$$
pressure: \Phi=\frac{p}{\rho \cdot(\Omega R)^{2}} \quad \tag{36}
$$

Note that most of these dimensionless variables are dependent on the rotor speed, $\Omega$ . Because the rotor speed may change over the course of a simulation (e.g. a variable speed turbine), these quantities must be calculated at the beginning of each time step.

The two primary equations for the generalized dynamic wake are then made dimensionless. Laplace’s equation (Equation 32) is also true for the dimensionless pressure:
$$
\nabla^2 \Phi=0, \quad \tag{37}
$$
and the momentum equation becomes
$$
\frac{\partial \hat{u}_{i}}{\partial \hat{t}}+\hat{U}_{\infty j} \frac{\partial \hat{u}_{i}}{\partial \hat{x}_{j}}=-\frac{\partial \Phi}{\partial \hat{x}_{i}} . \tag{38}
$$

The boundary conditions for these differential equations are given by the aerodynamic loading on the rotor blades and the requirement that the pressure return to ambient pressure far from the rotor. Also, the pressure discontinuity across the rotor plane must apply a force equal to the rotor thrust.

Using linear superposition, the pressure field can be divided into two components: one modeling the spatial variation of the pressure distribution, $\Phi^{V}$ , and one modeling the unsteadiness, $\Phi^{A}$ where:
$$
\Phi=\Phi^{V}+\Phi^{A} . \tag{39}
$$

By dividing the pressure field into two components, Equation 38 can also be divided into two separate equation sets as follows:
$$
\frac{\partial \hat{u}_{i}}{\partial \hat{t}}=-\frac{\partial \Phi^{A}}{\partial \hat{x}_{i}}, \tag{40}
$$
and
$$
\hat{U}_{\infty j} \frac{\partial \hat{u}_{i}}{\partial \hat{x}_{j}}=-\frac{\partial \Phi^{V}}{\partial \hat{x}_{i}} . \tag{41}
$$

Assuming that the differential equations (Equations 40 and 41) are linear and can be represented by a set of operators L and E the equations become:
$$
\frac{\partial \hat{u}_{i}}{\partial \hat {t}}=\hat {u}_{i}^{*}=L\left[\Phi^{A}\right], \tag{42}
$$
and
$$
\Phi_{i}=L[\Phi ^{V}] . \tag{43}
$$

As long as the operators L and E are invertible the solution for the dimensionless pressure field is
$$
\Phi=\Phi^{V}+\Phi^{A}=L^{-1}[\hat{u}]+E^{-1}[\hat{u}]^{*} . \tag{44}
$$

If an operator M is defined as the inverse of E , Equation 44 becomes,
$$
M[\hat{u}]^{*}+L^{-1}[\hat{u}]=\Phi . \tag{45}
$$

This is the general form of the governing equation of the generalized dynamic wake model that relates the induced velocity to the pressure field on the rotor disk. The methods used to solve these equations are described below.

#### 2.2.2 Pressure Distribution
Kinner (1937) developed the pressure distribution that satisfies Laplace’s equation (Equation 37) and that gives pressure discontinuity across a circular disk (the rotor). This solution was originally developed for the problem of a circular wing, but with different boundary conditions, also applies to a yawed actuator disk. The pressure distribution is given in an ellipsoidal coordinate system.
$$
\Phi(v, \eta, \psi, \hat{t})=\sum_{m=0}^{\infty} \sum_{n=m+1, m+3, \cdots}^{\infty} P_{n}^{m}(v) Q_{n}^{m}(i \eta)\left[C_{n}^{m}(\hat{t}) cos (m \psi)+D_{n}^{m}(\hat{t}) sin (m \psi)\right], \tag{46}
$$
where and $v , \eta$ and $\psi$ are ellipsoidal coordinates defined by the following relationships:
$$
x=\sqrt{1-v^{2}} \sqrt{1+\eta^{2}} cos \psi \tag{47}
$$
$$
y=\sqrt{1-v^{2}} \sqrt{1+\eta^{2}} sin \psi \quad \tag{48}
$$
$$
z=v \eta . \tag{49}
$$

Note that this xyz coordinate system follows the convention of the helicopter industry, with the $x - y$ plane parallel with the rotor plane, and the positive-z-axis perpendicular to the rotor plane in the upwind direction (see Figure 12). Also, note that the relationship between the dimensionless radius, $\hat{r}$ , and v in the ellipsoidal coordinate is given as $\hat{r}=\sqrt{1-v^{2}}$ .

The $v-\eta-\psi$ coordinate system covers the entire three-dimensional space once and only once, if $v, \eta$, and $\psi$ are restricted to the ranges
$$
-1 \leq v \leq 1 \tag{50}
$$
$$
0 \leq \eta \leq \infty \tag{51}
$$
$$
0 \leq \psi \leq 2 \pi . \tag{52}
$$

Figures 10 and 11 show contours of constant v and $\eta$ , respectively, in the $x-z$ plane (perpendicular to the rotor plane). The constant v surfaces are hyperboloids and the constant $\eta$ surfaces are ellipsoids. Both families of surfaces are azimuthally symmetric about the z-axis. The coordinate $\psi$ is the azimuthal angle measured from the positive x-axis and is positive in the clockwise direction. The $\eta=0$ surface represents both sides of the disk surface.

The pressure field of Equation 46 is discontinuous only within the unit circle (the rotor), where $\eta=0$ . And, because the pressure is perfectly continuous outside of the rotor, this distribution satisfies one of the boundary conditions that the rotor thrust force is zero outside the rotor boundary.

This pressure discontinuity provides thrust force on the rotor that simulates the aerodynamic forces on the blades. Although the actual aerodynamic forces act only on the blades and are discretely distributed, the distribution in Equation 46 gives a continuous distribution. However, the distribution starts to have peaks at the blades and to show the characteristics of discontinuity, as the number of terms of the series solution (flow states) increases.

The rotor disk pressure loading can be obtained as the pressure difference between the upwind and downwind surfaces of the rotor plane (He 1989),
$$
P(\hat{r}, \psi, \hat{t})=-2 \sum_{m=0}^{\infty} \sum_{n=m+1, m+3, \cdots}^{\infty} P_{n}^{m}(v) Q_{n}^{m}(i 0)\left[C_{n}^{m}(\hat{t}) cos (m \psi)+D_{n}^{m}(\hat{t}) sin (m \psi)\right],
$$
or,
$$
P(\hat{r}, \psi, \hat{t})=\sum_{m=0}^{\infty} \sum_{n=m+1, m+3, \cdots}^{\infty} \hat{P}_{n}^{m}(v)\left[\tau_{n}^{m c}(\hat{t}) cos (m \psi)+\tau_{n}^{m s}(\hat{t}) sin (m \psi)\right], \tag{54}
$$
where
$$
\hat{P}_{n}^{m}(v)=(-1)^{m} \frac{P_{n}^{m}(v)}{\rho_{n}^{m}} \tag{55}
$$
$$
\left(\rho_{n}^{m}\right)^{2}=\frac{1}{2 n+1} \frac{(n+m) !}{(n-m) !} \quad \tag{56}
$$
$$
\tau_{n}^{m c}=(-1)^{m+1}2Q_{n}^{m}(i0)\rho_{n}^{m}C_{n}^{m} \tag{57}
$$
$$
\tau_{n}^{m s}=(-1)^{m+1} 2 Q_{n}^{m}(i 0) \rho_{n}^{m} D_{n}^{m} . \tag{58}
$$

The dimensionless pressure quantities, $\tau_{n}^{m c}$ and $\tau_{n}^{m s}$ , couple the pressure distribution to the forces on the blades, as explained below.

The term $\hat{P}_{n}^{m}(v)$ is called the “normalized” associated Legendre function of the first kind, since it satisfies
$$
\int_{0}^{1}\left[\hat{P}_{n}^{m}(v)\right]^{2} d v=1 \quad \tag{59}
$$

#### 2.2.3 Induced Velocity Distribution
Similar to the expansion of the pressure distribution, the induced velocity distribution (the component normal to the rotor plane) can be expressed as an infinite series as shown in Equation 60 (He 1989):
$$
\hat{u}(\hat{r}, \psi, \hat{t})=\sum_{r=0}^{\infty} \sum_{j=r+1, r+3, \cdots}^{\infty} \phi_{j}^{r}(v)\left[\alpha_{j}^{r}(\hat{t}) cos (r \psi)+\beta_{j}^{r}(\hat{t}) sin (r \psi)\right], \tag{60}
$$
where the radial shape functions, $\phi_{j}^{r}(v)$ , are linearly independent and complete for a given harmonic, r . The coefficients $\alpha_{j}^{r}$ and $\beta_{j}^{r}$ can be regarded as the time-dependent states of the induced-velocity field. The shape functions are defined as
$$
\phi_{j}^{r}(\hat{r})=\sqrt{(2 j+1) H_{j}^{r}} \sum_{q=r, r+2, ...}^{j-1} \hat{r}^{q} \frac{(-1)^{\frac{q-r}{2}}(j+q) ! !}{(q-r) ! !(q+r) ! !(j-q-1) ! !} \tag{61}
$$
$$
H_{j}^{r}=\frac{(j+r-1) ! !(j-r-1) ! !}{(j+r) ! !(j-r) ! !} \tag{62}
$$
$n!!$ is a double factorial, defined as
$$
n!! \equiv
\begin{cases}
n \cdot(n-2) \cdot(n-4) \cdots 3 \cdot 1 & (n= odd ) \\
n \cdot(n-2) \cdot(n-4) \cdots 4 \cdot 2 & (n= even )
\end{cases} . \tag{63}
$$

For example, the first four shape functions of $r=2$ are
$$
\phi_{1}^{0}(\hat{r})=\sqrt{3} \tag{64}
$$
$$
\phi_{3}^{0}(\hat{r})=\sqrt{7}-\frac{5 \sqrt{7}}{2} \hat{r}^{2} \tag{65}
$$
$$
\phi_{1}^{2}(\hat{r})=\sqrt{\frac{15}{2}} \hat{r} \quad \tag{66}
$$
$$
\phi_{3}^{2}(\hat{r})=\frac{\sqrt{210}}{4} \hat{r}^{2} \tag{67}
$$

#### 2.2.4 Expansion of Governing Equation
With both the dimensionless pressure at the rotor disk and the induced velocity distribution expressed as infinite series of sines and cosines, they can be combined into the governing equation (45) (He 1989). The two operators are represented by square matrices.

The pressure coefficients in Equation 54, $\tau_{n}^{m c}$ and $\tau_{n}^{m s}$ , and the velocity coefficients in Equation 60, $\alpha_{j}^{r}$ and $\beta_{j}^{r}$ , have a relationship shown in Equations 68 and 69, which separate the cosine terms and sine terms. Equation 68 is the governing equation for the cosine terms and Equation 69 is the governing equation for the sine terms.
$$
\left[M^{c}\right]\left\{\begin{array}{c}\vdots \\ \left\{\alpha_{j}^{r}\right\} \\ \vdots \end{array}\right\}+\left[L^{c}\right]^{-1}\left\{\begin{array}{c}\vdots \\ \left\{\alpha_{j}^{r}\right\} \\ \vdots \end{array}\right\}=\frac{1}{2}\left\{\begin{array}{c}\vdots \\ \left\{\tau_{n}^{m c}\right\} \\ \vdots \end{array}\right\} \tag{68}
$$
$$
\left[M^{s}\right]\left\{\begin{array}{c}\vdots \\ \left\{\beta_{j}^{r}\right\} \\ \vdots\end{array}\right\}+\left[L^{s}\right]^{-1}\left\{\begin{array}{c}\vdots \\ \left\{\beta_{j}^{r}\right\} \\ \vdots\end{array}\right\}=\frac{1}{2}\left\{\begin{array}{c}\vdots \\ \left\{\tau_{n}^{m s}\right\} \\ \vdots\end{array}\right\} \tag{69}
$$
where $[M^{c}]$ and $[M^{s}]$ are the cosine and sine terms of the M-operator (see Equation 45) and similarly with the $[L^{c}]$ and $[L^{s}]$ matrices. The M-operator (apparent mass matrix) is nearly identical for both of the cosine and sine equations, with the exception that $[M^{c}]$ has elements for $m = 0$ and $[M^{s}]$ does not. The M-operators for both sine and cosine terms are given as
$$
[M]=\left[\begin{array}{lll}\ddots & & \\ & K_{n}^{m} & \\ & & \ddots\end{array}\right] \tag{70}
$$
$$
\zeta_{n}^{m}=\frac{2}{\pi} H_{n}^{m}=\frac{2}{\pi} \frac{(n+m-1) ! !(n-m-1) ! !}{(n+m) ! !(n-m) ! !} \tag{71}
$$

This mass matrix is purely diagonal, indicating that there is no coupling in either the harmonic or radial direction. This diagonal structure also simplifies the computation in a time marching scheme.

The L-operator (inflow gain matrix) is different for the sine and cosine equations and can be divided into matrices dependent on the flow parameters, $[\hat{V}]$ , and the wake skew angle, $[\tilde{L}]$ ,as follows:
$$
[L^{c}]^{-1}=[\tilde{L}^{c}]^{-1}[\hat {V}^{c}] \tag{72}
$$
$$
[L^{s}]^{-1}=[\tilde{L}^{s}]^{-1}[\hat {V}^{s}] \, . \tag{73}
$$

The matrices dependent on wake skew angle can be expressed as
$$
\left[\tilde{\sigma}_{j n}^{0 m}\right]^{c}=X^{m}\left[\Gamma_{j n}^{0 m}\right] \tag{74}
$$
$$
\left[\tilde{T}_{j n}^{r m}\right]^{c}=\left[X^{|m-r|}+(-1)^{l} X^{|m+r|}\right]\left[\Gamma_{j n}^{r m}\right] \tag{75}
$$
$$
\left[\tilde{\pi}_{j n}^{r m}\right]^{s}=\left[X^{|m-r|}-(-1)^{l} X^{|m+r|}\right]\left[\Gamma_{j n}^{r m}\right] \tag{76}
$$
with
$$
l=min (r,m) \tag{77}
$$
$$
X=tan \left|\frac{\chi}{2}\right|=\frac{\mu}{V_{T}+|\lambda|} \quad(0 \leq \chi \leq \pi / 2) \tag{78}
$$
$$
\Gamma_{j n}^{r m}=\frac{(-1)^{\frac{n+j-2 r}{2}}}{\sqrt{H_{n}^{m} H_{j}^{r}}} \frac{2 \sqrt{(2 n+1)(2 j+1)}}{(j+n)(j+n+2)\left[(j+n)^{2}-1\right]}, \text{ for } r+m= even \tag{79}
$$
$$
\Gamma_{j n}^{r m}=\frac{\pi}{2 \sqrt{H_{n}^{m} H_{j}^{r}}} \frac{sign(r-m)}{\sqrt{(2 n+1)(2 j+1)}}, \quad \text{ for } r+m= \text{odd and } |j-n|=1 \tag{80}
$$
$$
\Gamma_{j n}=0, \quad \text{ for } r+m= \text{odd and } |j-n| \neq 1. \tag{81}
$$

Note that the wake angle function of Equation 78 is determined from the average wake angle using the flow parameters described below. Also, the terms $X^{m}$ and $X^{|m-r|}$ , in Equations 74 - 76, can become $0^{0}$ . This zero raised to the zeroth power is considered one.

As with the M-operator the matrices dependent on flow parameters, $[\hat{V}]$ , are nearly identical for both of the cosine and sine equations, with the exception of m = 0 as follows:
$$
\left[\hat{V}^{c}\right]=\left[\begin{array}{lll}\ddots & & \\ & \hat{V}_{n}^{m} & \\ & & \ddots\end{array}\right] \text{ for } m=0,1,2,3, ... \tag{82}
$$
$$
\left[\hat{V}^{s}\right]=\left[\begin{array}{lll}\ddots & \hat{V}_{n}^{m} & \\ & & \ddots\end{array}\right] \text{ for } m=1,2,3, ... \tag{83}
$$
where
$$
\hat{V}_{1}^{0}=V_{T} \text{ for } (m, n)=(0,1) \quad \tag{84}
$$
$$
\hat{V}_{n}^{m}=V \text{ for } (m, n) \neq(0,1) . \quad \tag{85}
$$

The flow parameters of Equations 84 and 85 are based on the inflow (including blade motion) and induced velocities described in the next section.

#### 2.2.5 Flow Parameters
The inflow parameter, $V$, accounts for the energy that the rotor subtracts from the flow. It is calculated as follows (He 1989):
$$
V=\frac{\mu^{2}+\left(\lambda+\lambda_{m}\right) \lambda}{\sqrt{\mu^{2}+\lambda^{2}}} \tag{86}
$$
and the total flow parameter
$$
V_{T}=\sqrt{\mu^{2}+\lambda^{2}} . \tag{87}
$$
where μ is the advance ratio (in-plane velocity divided by tip speed), and the induced velocities, λ , are found by
$$
\lambda=\lambda_{m}+\lambda_{f} \quad \tag{88}
$$
$$
\lambda_{f}=\frac{U_{\infty} cos \chi}{\Omega R} \quad \tag{89}
$$

And $\lambda_{m}$ is calculated using Equation 90:
$$
\lambda_{m}=\sqrt{3} \alpha_{1}^{0} . \tag{90}
$$

The directions of each of these quantities are shown in Figure 12. Again, the wake skew angle, $\chi$ , of Equation 89 is the average over the entire rotor. The coupling between $V_{T}$ and $\alpha_{1}^{0}$ makes the theory nonlinear.

#### 2.2.6 Pressure Coefficients
The pressure coefficients, $\tau_{n}^{m c}$ and $\tau_{n}^{m s}$ , need to be coupled with the blade loading, which gives the boundary conditions of the model. Let $L_{i}^{q}$ be equal to the aerodynamic force normal to the rotor plane acting on blade element i of blade q (the element thrust force). If this element thrust force is normalized by the thrust force from the dynamic pressure of the flow:
$$
\text{normalized element thrust force} =\frac{L_{i}^{q}}{\rho A(\Omega R)^{2}} \quad \tag{91}
$$

The pressure coefficients are the normalized total thrust force multiplied by the radial expansion shape function and the azimuthal mode shape (modified from He 1989).
$$
\tau_{n}^{m c}=\frac{1}{\pi \rho \Omega^{2} R^{4}} \sum_{q=1}^{B}\left[\sum_{i=1}^{N_{E}} L_{i}^{q} \phi_{n}^{m}\left(\hat{r}_{i}\right)\right] cos (m \psi) \tag{92}
$$
$$
\tau_{n}^{m s}=\frac{1}{\pi \rho \Omega^{2} R^{4}} \sum_{q=1}^{B}\left[\sum_{i=1}^{N_{E}} L_{i}^{q} \phi_{n}^{m}\left(\hat{r}_{i}\right)\right] sin (m \psi) \tag{93}
$$

The sine terms (Equation 93) do not need to be defined for $m=0$ , because they are always multiplied by zero. The cosine terms for $m=0$ , however, need a slightly different definition because of the way $[\tilde{L}^{c}]$ is defined for $m=0$ in Equation 74. $[\tilde{L}^{c}]$ in Equation 74 is only one half of the $[\tilde{L}^{c}]$ in Equation 75 (modified from He, 1989). Therefore,
$$
\tau_{n}^{0 c}=\frac{1}{2 \pi \rho \Omega^{2} R^{4}} \sum_{q=1}^{B}\left[\sum_{i=1}^{N_{E}} L_{i}^{q} \phi_{n}^{0}\left(\hat{r}_{i}\right)\right] . \tag{94}
$$

#### 2.2.7 Final Governing Equations
The final set of equations for the GDW model using the parameters defined above are:
$$
\left[M^{c}\right]\left\{\begin{array}{c} \vdots \\ \left\{\alpha_{j}^{r}\right\} \\ \vdots \end{array}\right\}^{*}+\left[\tilde{L}^{c}\right]^{-1}\left[\hat{V}^{c}\right]\left\{\begin{array}{c}\vdots \\ \left\{\alpha_{j}^{r}\right\} \\ \vdots\end{array}\right\}=\frac{1}{2}\left\{\begin{array}{c}\vdots \\ \tau_{n}^{m c} \\ \vdots \end{array}\right\} \tag{95}
$$
$$
\left[M^{s}\right]\left\{\begin{array}{c} \vdots \\ \left\{\beta_{j}^{r}\right\} \\ \vdots \end{array}\right\}^{*}+\left[\tilde{L}^{s}\right]^{-1}\left[\hat{V}^{s}\right]\left\{\begin{array}{c} \vdots \\ \left\{\beta_{j}^{r}\right\} \\ \vdots \end{array}\right\}=\frac{1}{2}\left\{\begin{array}{c} \vdots \\ \left\{\tau_{n}^{m s}\right\} \\ \vdots \end{array}\right\} . \tag{96}
$$

The cosine term and the sine terms are not coupled. This indicates that the wake rotation is not considered in the generalized dynamic wake model itself. However, as mentioned above, AeroDyn uses the BEM method (Equation 28) to calculate the tangential components of the induced velocity. This completes the equations necessary for the generalized dynamic wake as implemented in AeroDyn.

#### 2.2.8 Procedure for Generalized Dynamic Wake Calculations
The above equations are written for an infinite number of azimuthal harmonics and radial shape functions. When implemented into a computer algorithm, the number of functions used in the modeling of the pressure distribution and induced velocity field must be truncated. For the induced velocity field,
$$
\hat{u}(\hat{r}, \psi, \hat{t})=\sum_{r=0}^{N} \sum_{j=r+1, r+3, \cdots}^{2 S_{r}+r-1} \phi_{j}^{r}(v)\left[\alpha_{j}^{r}(\hat{t}) cos (r \psi)+\beta_{j}^{r}(\hat{t}) sin (r \psi)\right], \tag{97}
$$
where N is the highest harmonic in the azimuthal direction and $S_{r}$ is the number of radial shape functions for the $r^{th }$ harmonic. The AeroDyn user must choose the values of the number of harmonics and radial shape functions to be modeled based on the structural dynamics and the desired resolution of the pressure or induced velocity distribution. The number of harmonics is often related to the number of blades. For example, He (1989) states that for a time-averaged solution of a four-bladed rotor, the induced velocity distribution can be truncated at the fourth harmonic, with little loss in accuracy. More harmonics may be required for an unsteady calculation.

Based on the number of harmonics, the number of radial shape functions can be determined. Table 1 shows the proper choice of the number of radial shape functions based on the mathematical consistency of the highest polynomial power of $\hat{r}$ for the radial shape function at each harmonic value, m . For example, in order to truncate the induced velocity distribution at the fourth harmonic ( $N=4$ ) with a radial variation up to $\hat{r}^{8}$ , the number of shape functions for each harmonic is then $S_{0}=5$ , $S_{1}=4$ , $S_{2}=4$ , $S_{3}=3$ , and $S_{4}=3$ . Now, remember that for $m=0$ only one inflow state is modeled, while all other values of m model two inflow states (i.e. the sine and cosine terms of Equation 97). We can then calculate the total number of inflow states for this example, $5+2(4+4+3+3)=33$ . Whereas, if we calculate the total number of inflow states using all harmonic values for that given power of $\hat{r}$ , we arrive at 45 inflow states, which is the last column in the table.

**Table 1. Choice for the Number of Inflow Radial Shape Functions**

| Power of $\hat{r}$ | Highest m (harmonic value) |  |  |  |  |  |  |  |  |  |  |  |  | Total Inflow States |
|--------------------|-----------------------------|----|----|----|----|----|----|----|----|----|----|----|----|---------------------|
|                    | 0                           | 1  | 2  | 3  | 4  | 5  | 6  | 7  | 8  | 9  | 10 | 11 | 12 |                     |
| 0                  | 1                           |    |    |    |    |    |    |    |    |    |    |    |    | 1                   |
| 1                  | 1                           | 1  |    |    |    |    |    |    |    |    |    |    |    | 3                   |
| 2                  | 2                           | 1  | 1  |    |    |    |    |    |    |    |    |    |    | 6                   |
| 3                  | 2                           | 2  | 1  | 1  |    |    |    |    |    |    |    |    |    | 10                  |
| 4                  | 3                           | 2  | 2  | 1  | 1  |    |    |    |    |    |    |    |    | 15                  |
| 5                  | 3                           | 3  | 2  | 2  | 1  | 1  |    |    |    |    |    |    |    | 21                  |
| 6                  | 4                           | 3  | 3  | 2  | 2  | 1  | 1  |    |    |    |    |    |    | 28                  |
| 7                  | 4                           | 4  | 3  | 3  | 2  | 2  | 1  | 1  |    |    |    |    |    | 36                  |
| 8                  | 5                           | 4  | 4  | 3  | 3  | 2  | 2  | 1  | 1  |    |    |    |    | 45                  |
| 9                  | 5                           | 5  | 4  | 4  | 3  | 3  | 2  | 2  | 1  | 1  |    |    |    | 55                  |
| 10                 | 6                           | 5  | 5  | 4  | 4  | 3  | 3  | 2  | 2  | 1  | 1  |    |    | 66                  |
| 11                 | 6                           | 6  | 5  | 5  | 4  | 4  | 3  | 3  | 2  | 2  | 1  | 1  |    | 78                  |
| 12                 | 7                           | 6  | 6  | 5  | 5  | 4  | 4  | 3  | 3  | 2  | 2  | 1  | 1  | 91                  |

Once the user determines number of inflow states, the calculation of the induced velocity proceeds as follows. Because the GDW method in AeroDyn is based on a solution of ordinary differential equations in time, it must rely on initial values of various parameters to accurately calculate the effect of the wake. The initial values are based on BEM calculations of the operating turbine over the first second of the time simulation. After one second has passed in the simulation, AeroDyn switches from the BEM method to the GDW using the BEM solution as the initial condition for the GDW method. The blade forces are used to calculate the pressure coefficients, $\tau_{n}^{m c}$ and $\tau_{n}^{m s}$ , in Equations 91 through 94, which form the right hand side of Equations 95 and 96. Equations 92-94 transform the loading on the blades to a pressure distribution around the entire actuator disk. The apparent mass matrix, $[M]$ , is calculated for each harmonic and radial shape function based on Equations 70 and 71. The inflow gain matrix, $[\tilde{L}]$ , is calculated based on the wake skew angle in the formulae given in Equations 74-81. The flow parameter matrix, $[\hat{V}]$ , is assembled based on formulae 82-90. Once all of the matrices are assembled, Equations 95 and 96 are solved using a fourth-order Adams-Bashford-Moulton (Press et al. 1982) predictor-corrector method, for each of the azimuthal harmonics and radial shape functions. The solutions of these equations are the coefficients of Equation 60, $\alpha_{j}^{r}$ and $\beta_{j}^{r}$ . These coefficients are then fed back into Equation 60, along with the radial shape functions of Equations 61-63, in order to calculate the induced velocity field at any point in the rotor plane. These induced velocities are then used to determine the angle of attack for each element. This angle of attack is passed to the airfoil aerodynamics routines that return the elemental force based on either static or dynamic stall airfoil conditions as explained below. These forces are the output of the AeroDyn routines and are passed back to the aeroelastic code for further analysis. Starting at the next time step, the process is repeated using the most recent estimated forces to calculate the pressure coefficients in Equations 92 through 94, which serve as the initial conditions for Equations 95 and 96.