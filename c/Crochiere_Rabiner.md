# Clockwise Commutator

To deconstruct the legs in a clockwise commutator model that starts at the first row being 0 and goes down, with $M$ being the number 
of polyphase legs and downsample factor and overbar indicating the clockwise model.

$$ \bar{p}_{\rho} \left( n \right) = h \left( n M - \rho \right) \quad \rho = 0, 1, \ldots, M-1 $$

$$ \bar{x}_{\rho} \left( n \right) = x \left( n M + \rho \right) \quad \rho = 0, 1, \ldots, M-1 $$

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
$$ \left( \left\lceil \frac{N_H - 1}{M} \right\rceil + 1 \right) M $$

# Channelizer

The Crochiere and Rabiner critically sampled channelizer is a little different from the Sharpin version. This will detail some of the differences.

1. It's a clockwise commutator model.
2. It uses the DFT down the channels instead of the IDFT.
3. The filters are defined as $ \bar{p}_{\rho} ( m ) = h ( m M - \rho ) $
4. The branch input signals are defined as $ x_\rho (m ) = x ( m M + \rho ), \quad \rho = 0, 1, \ldots, M-1 $

# Oversampled Channelizer

For the oversampled version it's even more different. Let's use the case in which $ K = M I $, where $M$ is the downsampling amount now and $K$ is the number of channels. We're interested in $I = 2$, so a 2x oversampled channelizer. This is about right for making the stopband of the filter hit by the time we get to aliasing on the output. 

Now the filters are defined with a _subtle_ difference:

$$ \bar{p}_{\rho} ( m ) = h ( m M - \rho ), \quad \rho = 0, 1, \ldots, K-1 $$

Do you see it? The difference is that the $\rho$ extends to $K-1$ but the downsampling is based on half of that, $M$. This leads to a similar structure as the Sharpin one in which half the channels are repeats of the filter taps of the others but with an extra delay. But it _doesn't_ have the upsampling in the filter structure.

The data into each channel is defined in a more traditional commutator sense, without the repitition the filter has, as

$$ x_{\rho} ( r ) = x ( r K + \rho), \quad \rho = 0, 1, \ldots, K-1 $$

The other main difference is that *after* the commutator but *before* the polyphase filters there is an upsample by $I = 2$ step. I *think* this can maybe be handled by modifying the `ostride` value in FFTW. Which would make it super convenient to not have to do any extra data copies. 

# Input Signal

The python input can just import `chirp` from `scipy.signal`. But we have to make our own for the C version. It's the ideal signal. Pulling from [this](https://en.wikipedia.org/wiki/Chirp) wikipedia page, the instantaneous frequency of a chirp is:

$$f(t)=ct+f_{0}$$

where $f_0$ is the starting frequency at time 0 and $c$ is the chirp rate, assumed to be:

$$ c = \frac{f_1 - f_0}{T}$$

Where $f_1$ is the final frequency and $T$ is the time it takes to sweep from $f_0$ to $f_1$.

Integrating that frequency and putting it as the phase of a sinusoid gives us

$$x(t)=\sin \left[\phi_0 + 2 \pi \left( \frac{c}{2} t^2 + f_0 t \right)\right]$$

We'll probably just assuming the starting phase is zero. If we also discretize the equation with $t = n T_s$ or $t = n / F_s$ then we end up with:

$$x(t)=\sin \left[2 \pi \left( \frac{c T_s^2}{2} n^2 + f_0 n T_s \right)\right]$$

