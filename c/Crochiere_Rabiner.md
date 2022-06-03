# Clockwise Commutator

To deconstruct the legs in a clockwise commutator model that starts at 0 and goes down, with $M$ being the number 
of polyphase legs and downsample factor and overbar indicating the clockwise model.

$$ \bar{p}_{\rho} \left( n \right) = h \left( n M - \rho \right) $$

$$ \bar{x}_{\rho} \left( n \right) = x \left( n M - \rho \right) $$

So the matrices would look like the below (if the top row is $\rho = 0$ and it proceeds down from there to 1, 2, 3, ...)

$$
H_{poly} = \left[
\begin{matrix}
0 & M & 2M & \cdots \\
-1 & M-1 & 2M-1 & \cdots \\
-2 & M-2 & 2M-2 & \cdots \\
\vdots & \vdots & \vdots & \cdots \\
-M+1 & 1 & M+1 & \cdots \\
\end{matrix}
\right]
$$

