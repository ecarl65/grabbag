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

The data polyphase decomposition has advances and doesn't "seem" realizable in a block diagram sense. Luckily we're
doing this on a computer with buffers so it should work just fine.

$$
X_{poly} = \left[
\begin{matrix}
0 & M & 2M & \cdots \\
1 & M+1 & 2M+1 & \cdots \\
2 & M+2 & 2M+2 & \cdots \\
\vdots & \vdots & \vdots & \cdots \\
M-1 & 2M-1 & 3M-1 & \cdots \\
\end{matrix}
\right]
$$

# Procedure

This would seem to indicate that the amount of data should be some integer multiple of $M$ plus 1. When chunking up
the input this is important to keep in mind.

More important is getting the right filter length. The procedure might be:

1. Estimate the filter size for the requirements. Call it $N_H$
2. Find integer $L$ for nearest (next highest? Lowest?) result of $2 M L + 1$ to $N_H$
3. $$ N_H = round \left( \frac{N_H - 1}{2 M} \right) 2 M + 1 $$
4. Find the poly filter matrix decomposition length.
$$ \left( \lceil \frac{N_H - 1}{M} \rceil + 1 \right) M $$


