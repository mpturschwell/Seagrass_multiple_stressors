
# Analytical solution to the Gompertz model

We have the following model:

$$ \frac{dB}{dt} = P_{max}(I,T) \log \left(\frac{B_{max}}{B}\right)B -R(T)B-MB$$

with

$$P_{max}(I,T) = P_{Tmax}\left(\frac{T_{max}-T}{T_{max} - T_{opt}} \right)\left(\frac{T}{T_{opt}}\right)^{(T_{opt}/(T_{max} - T_{opt})}\tanh\left(\frac{I}{I_k}\right)$$

$$ R(T)=R_{max}\left(\frac{R_{Tmax}-T}{R_{Tmax} - R_{Topt}} \right)\left(\frac{T}{R_{Topt}}\right)^{(R_{Topt}/(R_{Tmax} - R_{Topt})}$$

Now rearranging the ODE we get 

$$
\begin{aligned}

    \frac{1}{B}\frac{dB}{dt} &= P_{max} \left( \log (B_{max}) - \log(B) \right) -R-M \\

    \implies & P_{max} \log(B) +\frac{1}{B}\frac{dB}{dt} = P_{max} \log (B_{max}) -R-M \\

\end{aligned}
$$

Now we use a change of variables, where $y = \log(B) \implies dy/dt = 1/B$. Then we have

$$
\begin{aligned}

 P_{max} y +\frac{dy}{dB}\frac{dB}{dt} &= P_{max} \log (B_{max}) -R-M\\

    \implies  P_{max} y +\frac{dy}{dt} &= P_{max} \log (B_{max}) -R-M \text{: by the chain rule}

\end{aligned}
$$

Which is the form of a linear differential equation, so we can solve by multiplying through by an integrating factor $z$, where $z = \exp(\int P_{max}dt)=\exp(P_{max}t)$. This gives

$$
\begin{aligned}

     P_{max} \exp(P_{max}t)y +\frac{dy}{dt}\exp(P_{max}t) & = \left(P_{max} \log (B_{max}) -R-M\right)\exp(P_{max}t)\\

     \implies  \frac{d}{dt}\left(\exp(P_{max}t)y\right)  &= \left(P_{max} \log (B_{max}) -R-M\right)\exp(P_{max}t)\\

     \implies \int \frac{d}{dt}\left(\exp(P_{max}t)y\right) dt &= \int \left(P_{max} \log (B_{max}) -R-M\right)\exp(P_{max}t) dt\\

     \implies \exp(P_{max}t)y &= \frac{\left(P_{max} \log (B_{max}) -R-M\right)}{P_{max}}\exp(P_{max}t) + C \text{ : where } C \in \mathbb R \\

     \implies y &= \log (B_{max}) -  \frac{R+M}{P_{max}} + C\exp(-P_{max}t)\\

     \implies \log(B) &= \log (B_{max}) - \frac{R+M}{P_{max}} + C\exp(-P_{max}t)\\

     B &= \exp\left(\log (B_{max}) - \frac{R+M}{P_{max}} + C\exp(-P_{max}t)\right )\\

     B &= B_{max}\exp\left(-\frac{R+M}{P_{max}} + C\exp(-P_{max}t)\right)

\end{aligned}
$$

Now, we use the initial condition $B(0) = B_0$ to solve for $C$, where

$$ \log(B_0) = \frac{\left(P_{max} \log (B_{max}) -R-M\right)}{P_{max}}+ C $$

$$\implies C = \log(B_0) - \log (B_{max}) + \frac{\left(R+M\right)}{P_{max}}$$

$$\implies C = \log\left(\frac{B_0}{B_{max}}\right)  + \frac{\left(R+M\right)}{P_{max}}$$

Combining these expressions we have that the analytical solution to the Gompertz model is 

$$B(t) = B_{max}\exp\left(-\frac{R+M}{P_{max}} + \left(\log\left(\frac{B_0}{B_{max}}\right)  + \frac{\left(R+M\right)}{P_{max}}\right)\exp(-P_{max}t)\right)$$

