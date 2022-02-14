
# Analytical solution to the $P_{max}$ dependent model



$$ \frac{dB}{dt} = (P_{max}(I,T) - \psi B)B -R(T)B-MB$$

with

$$P_{max}(I,T) = P_{Tmax}\left(\frac{T_{max}-T}{T_{max} - T_{opt}} \right)\left(\frac{T}{T_{opt}}\right)^{(T_{opt}/(T_{max} - T_{opt})}\tanh\left(\frac{I}{I_k}\right),$$

$$ R(T)=R_{max}\left(\frac{R_{Tmax}-T}{R_{Tmax} - R_{Topt}} \right)\left(\frac{T}{R_{Topt}}\right)^{(R_{Topt}/(R_{Tmax} - R_{Topt})} .$$

Let $\alpha = (P_{max} - R -M)/\psi$. Now, this model is a separable one-dimensional ODE, where

$$
\begin{aligned}

 \frac{dB}{dt} &= (P_{max} - \psi B -R-M)B \\
 \implies & \frac{dB}{dt} = (\psi\alpha - \psi B )B \\
 \implies &\int \frac{dB}{\psi(\alpha - B )B} = \int dt \\
 \implies &\int \frac{dB}{(B - \alpha )B} = -\psi t + C \\
 

\end{aligned}
$$

Now, we solve the LHS using integration by partial fractions, where

$$ LHS = \int\left( \frac{E}{(B-\alpha )} + \frac{F}{B}\right)dB $$

with $EB + F(B-\alpha) = 1$. Let $B = 0$, then we have $F = -1/\alpha$. Now, let $B = \alpha$ then we have $E = 1/\alpha$. So,

$$
\begin{aligned}
    LHS &= \int\left( \frac{1}{\alpha(B-\alpha )} - \frac{1}{\alpha B}\right)dB\\
    &= \frac{1}{\alpha} \int\left( \frac{1}{(B-\alpha )} - \frac{1}{ B}\right)dB\\
    &= \frac{1}{\alpha} \left( \log(B-\alpha ) - \log(B)\right)\\
    &= \frac{1}{\alpha} \left( \log\left(\frac{B-\alpha}{B} \right) \right).
\end{aligned}
$$

Thus, we have

$$
\begin{aligned}

 
 \frac{1}{\alpha} \left( \log\left(\frac{B-\alpha}{B} \right) \right)& = -\psi t + C \\
  \implies & \frac{B-\alpha}{B} = G\exp(-\psi \alpha t) \text{ : with } G\in \mathbb{R}\\
  \implies & 1-G\exp(-\psi \alpha t) =  \alpha/B\\
  \implies & B = \frac{\alpha}{1-G\exp(-\psi \alpha t)}.\\

\end{aligned}
$$

Now, we use the initial condition $B(0) = B_0$ to solve for $G$, where

$$ G = \frac{B_0-\alpha}{B_0} = 1 - \alpha/B_0$$

Thus,substituting in $G$ and $\alpha$, the analytical solution for the model is

$$ B(t) =  \frac{(P_{max} - R -M)/\psi}{1-(1 - (P_{max} - R -M)/\psi B_0)\exp(- (P_{max} - R -M) t)}$$



