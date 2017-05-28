Theory of Boyle
===============

The Boyle-Computation engine is built on top of what is known as the BioModel (*Reference goes here*). The program definitions and details are written here for clarity.

Variables and Constants
-----------------------

The variables that form part of the system are:

- :math:`T` : Reactor Temperature
- :math:`T_{inc}` : Increase rate of reactor temperature (if it is not constant)
- :math:`T_{base}` : Starting base temperature, if temperature increases
- :math:`Qin_{daily}` : Daily volume of feed flowing in the reactor
- :math:`Qout_{daily}`: Daily volume of effluent flowing out of the reactor
- :math:`P_{int}`: Pumping interval (frequence) is either continuous or intermittent.
- :math:`P_{idx}`: Pumping index referring to the efficiency of the pump
- :math:`runtime`: Time period framing the current simulation cycle.

The constants that form the part of the computation system are:

- :math:`K_{d0}`: Bacterial death rate
- :math:`K_{s}`: Half-saturation constants for bacterial growth on the substrate
- :math:`Ks_{NH3}`: Half-saturation constants for bacterial growth on Ammonia
- :math:`Ki_{carb_enz}`: VFA inhibition of carbohydrate enzymatic hydrolysis
- :math:`Ki_{prot_enz}`: VFA inhibition of protein enzymatic hydrolysis
- :math:`Ki_{HAc_HPr}`: Acetic acid inhibition of propionic acid degradation
- :math:`Ki_{HAc_HBut}`: Acetic acid inhibition of butyric acid degradation
- :math:`Ki_{HAc_HVal}`: Acetic acid inhibition of valeric acid degradation
- :math:`Ki_{NH3_HAc}`: Ammonia inhibition of acetic acid degration
- :math:`Ki_{LCFA}`: LCFA Inhibition involving all levels of anaerobic digestion
- :math:`pK_{low}`: Lower pH dropoff values for the Michaelis pH function
- :math:`pk_{high}`: Higher pH dropoff values for the Michaelis pH function

Calculation of Maximum Enzyme Reaction
--------------------------------------

The maximum enzyme reaction and the bacterial growth rates are a function of the reactor temperature. This compuatation is given by the following equations based on the temperatures.

If the :math:`T < T_{opt}`, then the equation for computation of maximum enzyme reaction is:

.. math:: \mu_{max}[T] = \mu_{max}[T_{opt}] - \alpha \cdot (T - T_{opt})
    :label: minimumenzyme

and if :math:`T > T_{opt}` then the equation becomes:

.. math:: \mu_{max}[T] = (\mu_{max}[T_{opt}] + \alpha \cdot (T_{opt} - T_{0})) \cdot \frac{(T_{max} - T)}{(T_{max} - T_{opt})}
    :label:  maximumenzyme

where :math:`\mu_{max_T0}` is the maximum growth rate of bacteria at :math:`T_{0}` and :math:`\mu_{max}` is the maximum growth rate of bacteriat at :math:`T`.

Calculation of Dissociation Constant
------------------------------------

.. math:: X_{T} = X_{T0} + \Delta T \cdot (a + \Delta T \cdot (b + \Delta T \cdot c))

Newton-Raphson Method for pH Computation
----------------------------------------

.. math:: RHS(1) = \left( \frac{CO_2}{44} \cdot \frac{Ka1_{CO_2} \cdot H}{H^2 + H \cdot Ka1_{CO_2} \cdot Ka2_{CO_2}} \right)
.. math:: RHS(2) = \left( 2 * \cdot \frac{CO_2}{44} \cdot \frac{Ka1_{CO_2} \cdot (4 * H * Ka2_{CO_2} - 2 \cdot Ka1_{CO_2} \cdot Ka2_{CO_2} ) }{ (H^2 + H \cdot Ka1_{CO_2} + Ka1_{CO_2} \cdot Ka2_{CO_2}) } \right)
.. math:: RHS(3) = \left( \frac{HAc}{60} \cdot \frac{Ka_{HAc} }{H + Ka_{HAc} } \right)
.. math:: RHS(4) = \left( \frac{HPr}{74} \cdot \frac{Ka_{HPr} }{H + Ka_{Hpr} } \right)
.. math:: RHS(5) = \left( \frac{HBut}{88} \cdot \frac{Ka_{HBut} }{H + Ka_{HBut} } \right)
.. math:: RHS(6) = \left( \frac{HVal}{102} \cdot \frac{Ka_{HVal} }{H + Ka_{HVal} } \right)
.. math:: RHS(7) = \left( \frac{H_2PO_4}{31} \cdot \frac{H}{H + Ka_{H_2PO_4} } \right)
.. math:: RHS(8) = \left( 2 \cdot \frac{H_2PO_4}{31} \cdot \frac{Ka_{H_2PO_4} }{H + Ka_{H_2PO_4} } \right)
.. math:: RHS(9) = \left( \frac{H_2S}{34} \cdot \frac{Ka_{H_2S}}{H + Ka_{H_2S}} \right)
.. math:: RHS(10) = \left( \frac{LCFA}{282.46} \cdot \frac{Ka_{LCFA}}{H + Ka_{LCFA}} \right)
.. math:: RHS(11) = \left( \frac{A}{35.5} \right)
.. math:: RHS(12) = \left( \frac{NH_4}{14} \cdot \frac{H}{H + Ka_{NH_4}} \right)
.. math:: RHS(13) = \frac{Z}{39}
.. math:: RHS(14) = \frac{K_w}{H}

.. math::
    H_{func} = RHS(1) + RHS(2) + RHS(3) + RHS(4) + RHS(5) \
    + RHS(6) + RHS(7) + RHS(8) + RHS(9) + RHS(10) \
      + RHS(11) - RHS(12) - RHS(13) + RHS(14)


.. math:: RHS(1) = - \left( \frac{CO_2}{44} \cdot \frac{ Ka1_{CO_2} \cdot (H^2 - Ka1_{CO_2} \cdot Ka2_{CO_2}) }{ H^2 + H \cdot Ka1_{CO_2} + Ka1_{CO_2} \cdot Ka2_{CO_2} } \right)
.. math:: RHS(2) = \left( \frac{CO_2}{44} \cdot \frac{ Ka1_{CO_2} \cdot (4 \cdot H \cdot Ka2_{CO_2} - 2 \cdot Ka1_{CO_2} \cdot Ka2_{CO_2}) }{ (H^2 + H \cdot Ka1_{CO_2} + Ka1_{CO_2} \cdot Ka2_{CO_2}) } \right)
.. math:: RHS(3) = \left( \frac{HAc}{60} \cdot \frac{Ka_{HAc}}{H + Ka_{HAc}} \right)
.. math:: RHS(4) = \left( \frac{HPr}{74} \cdot \frac{Ka_{HPr}}{H + Ka_{HPr}} \right)
.. math:: RHS(5) = \left( \frac{HBut}{88} \cdot \frac{Ka_{HBut}}{H + Ka_{HBut}} \right)
.. math:: RHS(6) = \left( \frac{HVal}{102} \cdot \frac{Ka_{HVal}}{H + Ka_{HVal}} \right) 
.. math:: RHS(7) = \left( \frac{H_2PO_4}{31} \cdot \frac{Ka_{H_2PO_4}}{H + Ka_{H_2PO_4}} \right)
.. math:: RHS(8) = \left( \frac{H_2PO_4}{31} \cdot \frac{2 \cdot Ka_{H_2PO_4} }{ (H + Ka_{H_2PO_4})^2 } \right)
.. math:: RHS(9) = \left( \frac{H_2S}{34} \cdot \frac{Ka_{H_2S}}{H + Ka_{H_2S}} \right)
.. math:: RHS(10) = \left( \frac{LCFA}{14123} \cdot \frac{50 \cdot Ka_{LCFA}}{H + Ka_{LCFA}} \right) 
.. math:: RHS(11) = 0
.. math:: RHS(12) = \left( \frac{NH_4}{14} \cdot \frac{Ka_{NH_4}}{ H + Ka_{NH_4} } \right)
.. math:: RHS(13) = 0
.. math:: RHS(14) = \frac{K_w}{H^2}

.. math::
    dHfunc_{dH} = RHS(1) - RHS(2) - RHS(3) - RHS(4) - RHS(5) - RHS(6) + RHS(7) - RHS(8) - RHS(9) - RHS(10) + RHS(11) - RHS(12) + RHS(13) - RHS(14)

Calculation of Ion Concentration based on pH
--------------------------------------------

.. math:: HCO_{3_neg} = \frac{CO_2}{44} \cdot \frac{ Ka1_{CO_2} \cdot H }{ H^2 + H \cdot Ka1_{CO_2} + Ka1_{CO_2} \cdot Ka2_{CO_2} } \cdot 61

.. math:: CO_{3_neg} = \frac{CO_2}{44} \cdot \frac{ Ka1_{CO_2} \cdot H }{ H^2 + H \cdot Ka1_{CO_2} + Ka1_{CO_2} \cdot Ka2_{CO_2} } \cdot 62

.. math:: Ac_{neg} = \frac{HAc}{60} \cdot \frac{Ka_{HAc}}{H + Ka_{HAc}} \cdot 59

.. math:: Pr_{neg} = \frac{HPr}{74} \cdot \frac{Ka_{HPr}}{H + Ka_{HPr}} \cdot 73

.. math:: But_{neg} = \frac{HBut}{88} \cdot \frac{Ka_{HBut}}{H + Ka_{HBut}} \cdot 87

.. math:: Val_{neg} = \frac{HVal}{102} \cdot \frac{Ka_{HVal}}{H + Ka_{HVal}} \cdot 101

.. math:: HPO_{4_neg} = \frac{H_2PO_4}{31} \cdot \frac{Ka_{H_2PO_4}}{H + Ka_{H_2PO_4}} \cdot 30

.. math:: HS_{neg} = \frac{H_2S}{34} \cdot \frac{Ka_{H_2S}}{H + Ka_{H_2S}} \cdot 33

.. math:: LCFA_{neg} = \frac{LCFA}{282.46} \cdot \frac{Ka_{LCFA}}{H + Ka_{LCFA}} \cdot 281.46

.. math:: NH3_{liq} = \frac{NH_4}{14} \cdot \frac{Ka_{NH_4}}{H + Ka_{NH_4}} \cdot 14

Calculation of Bacterial Growth Rates per Unit Mass
---------------------------------------------------

FpH is the normalized Michaelis pH function describing the effect of pH on the bacterial growth rates which is given by:

.. math:: FpH = \frac{ 1 + 2 \cdot 10^{0.5 \cdot (pK_{low} - pK_{high}) }}{ 1 + 10^{(pH - pk_{h})} + 10^{(pK_{low} - pH)} }

Using the FpH the growth rates are computed by :math:`\mu(N)` where N=1:8 which are the growth rates of the acidogenic glucose degrader, acidogenic amino acid degrader, acidogenci lipid degrader, acetogenic LCFA degrader, acetogenic propionate degrader, acetogenic butyrate degrader, acetogenic valerate degrader, and aceticlastic methanogens.

.. math:: \mu[1] = \mu_{max}[T] = \cdot \frac{ CH_4 }{ K_s + CH_4 } \cdot \frac{ NH_4 }{ Ks_{NH_4} + NH_4 } \cdot \frac{ Ki_{LCFA} }{ LCFA + Ki_{LCFA} } \cdot FpH

.. math:: \mu[2] = \mu_{max}[T] \cdot \frac{ Amino }{ K_s + Amino } \cdot \frac{ Ki_{LCFA} }{ LCFA + Ki_{LCFA} } \cdot FpH

.. math:: \mu[3] = \mu_{max}[T] \cdot \frac{ Lipids }{ K_s + Lipids } \cdot \frac{ NH_4 }{ Ks_{NH_3} + NH_4 } \cdot \frac{ Ki_{LCFA} }{ LCFA + Ki_ {LCFA} } \cdot FpH

.. math:: \mu[4] = \mu_{max}[T] \cdot \frac{ LCFA }{ LCFA + K_s + \frac{ LCFA^2 }{ Ki_{LCFA} } } \cdot \frac{ NH_4 }{ Ks_{NH_3} + NH_4 } \cdot FpH

.. math:: \mu[5] = \mu_{max}[T] \cdot \frac{ HPr }{ K_s + HPr } \cdot \frac{ NH_4 }{ Ks{NH_3} + NH_4 } \cdot \frac{ Ki_{HAc, HPr} }{ HAc + Ki_{HAc, HPr} } \cdot \frac{ Ki_{LCFA} }{ LCFA + Ki_{LCFA} } \cdot FpH

.. math:: \mu[6] = \mu_{max}[T] \cdot \frac{ HBut }{ K_s + HBut } \cdot \frac{ NH_4 }{ Ks_{NH_3} + NH_4 } \cdot \frac{ Ki_{HAc, HBut} }{ HAc + Ki_{HAc, HBut} } \cdot \frac{ Ki_{LCFA} }{ LCFA + Ki_{LCFA} } \cdot FpH

.. math:: \mu[7] = \mu_{max}[T] \cdot \frac{ HVal }{ K_s + HVal } \cdot \frac{ NH_4 }{ Ks_{NH_3} + NH_4 } \cdot \frac{ Ki_{HAc, Hval} }{ HAc + Ki_{HAc, HVal} } \cdot \frac{ Ki_{LCFA} }{ LCFA + Ki_{LCFA} } \cdot FpH

.. math:: \mu[8] = \mu_{max}[T] \cdot \frac{ HAc }{ K_s + HAc } \cdot \frac{ NH_4 }{ Ks_{NH_3} + NH_4 } \cdot \frac{ Ki_{NH_3, HAc} }{ \frac{ NH_4 \cdot Ks_{NH_4} }{ H + Ka_{NH_4} } + Ki_{NH_3, HAc} } \cdot \frac{ Ki_{LCFA} }{ LCFA + Ki_{LCFA} } \cdot FpH

Calculation of Gasflow
----------------------

.. math:: \alpha = \frac{ \left[ \frac{ Ka_{NH_4} }{ H + Ka_{NH_4} }; 1; \frac{ H^2 }{ H^2 + H \cdot Ka1_{CO_2} + Ka1_{CO_2} \cdot Ka2_{CO_2} }; \frac{ H }{ H + Ka_{H_2S} } \right] }{ K_H }

.. math:: dadH = \frac{ \left[ - \frac{ Ka_{NH_4} }{ H + Ka_{NH_4} }; 0; \frac{ Ka1_{CO_2} \cdot H \cdot (H + 2 \cdot Ka2_{CO_2}) }{ (H^2 + H \cdot Ka1_{CO_2} + Ka1_{CO_2} \cdot Ka2_{CO_2})^2 }; \frac{ Ka_{H_2S} }{ (H + Ka_{H_2S})^2 } \right] }{ K_H }

