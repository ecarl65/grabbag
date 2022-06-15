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

# Persistent Memory

Trying to keep some things around from call-to-call when doing this in buffered mode is going to be a challenge. Don't want to re-design
the filter every time and do the FFT. Or the Overlap/Save, perhaps. Although previous processing just treated each buffer as 
unique and maybe that's the way to go. But if we can make it an object that lives in the top-level processor like currently that 
would be ideal. But that can't be done with C. Will have to move to pybind11 or CFFI.

It _is_ worth noting that to design and prepare the channelizer all it requires is:
* The sample rate
* The downsample factor
* The filter overlap amount

That seems simple enough that looking at the pybind11 code it can be done.

# Dimensions

This shows the dimensions of the oversampled filterbank. In that case $K$ is the number of channels, $M$ is the downsample rate, and $I$ is an integer such that $K = MI$. In practice we'll pretty much always make $I = 2$. Note that a lot of these dimensions change in or out of the FFTs because we're assuming that we're starting with real input data. When you do an FFT of real input data that's length $N$ that makes the output be length $N // 2 + 1$, using the python notation for integer/floor division. In that case the output has the first and last samples real and the in-between samples are complex. Doing the inverse gets you back to length $N$. 

## Input

| 0 | 1 | 2 | ... | N<sub>full</sub> - 1 |
| - | - | - | --- | -------------------- |

## Buffer Size

N<sub>buffer</sub> &Lt; N<sub>full</sub>

| 0 | 1 | 2 | ... | N<sub>buffer</sub> - 1 |
| - | - | - | --- | ---------------------- |

## Filter

N<sub>filter</sub> &Lt; N<sub>buffer</sub>


| 0 | 1 | 2 | ... | N<sub>filter</sub> - 1 |
| - | - | - | --- | ---------------------- |

## Polyphase Filter

| 0     | 1 | 2 | 3 | ... | N<sub>buffer</sub> / M - 1 |
|:-----:|---|---|---|-----|----------------------------|
| 1     |   |   |   |     |                            |
| ...   |   |   |   |     |                            |
| K - 1 |   |   |   |     |                            |

## FFT Filter

| 0     | 1 | 2 | 3 | ... | N<sub>buffer</sub> / (2M) |
|:-----:|---|---|---|-----|---------------------------|
| 1     |   |   |   |     |                           |
| ...   |   |   |   |     |                           |
| K - 1 |   |   |   |     |                           |


## Polyphase Data

| 0     | 1 | ... | N<sub>buffer</sub> / K - 1 |
|:-----:|---|-----|----------------------------|
| 1     |   |     |                            |
| ...   |   |     |                            |
| K - 1 |   |     |                            |

## FFT Data

| 0     | 1 | ... | N<sub>buffer</sub> / (2K) |
|:-----:|---|-----|---------------------------|
| 1     |   |     |                           |
| ...   |   |     |                           |
| K - 1 |   |     |                           |

## Circular Convolution

In reality we're supposed to upsample each input polyphase channel of data by $I$, where $K = MI$. But we can cleverly get around that by doing it in the frequency domain. Upsampling in the time domain is equivalent to $I$ copies in the frequency domain. But rather than copying the data over and over we can just do a modulo on the indexing when we multiply and wrap around the values again. I'm a little shocked I didn't know about this trick, but when doing an FFT-based convolution it sure seems like a simple hack. It means we can do the polyphase decomposition and FFTs without rearranging any input data and passing parameters to FFTW to handle strides. The parameters for the FFT of an input array of size N<sub>buffer</sub> is:

```c++
int n_size[] = {N_buffer / M};    // Logical size of FFT
int rank = 1;                     // Multiple 1D FFTs
int howmany = K;                  // How many FFTs to do
int idist = 1;                    // Distance in input between starting each FFT
int odist = N_buffer / (2*K) + 1; // Distance in output between start of each FFT output
int istride = K;                  // Distance between successive inputs in each FFT, this is the polyphase jump
int ostride = 1;                  // Distance between successive outputs in each FFT, this is how the polyphase alters output
int *inembed = NULL;
int *onembed = NULL;
```
    
So the equation would be something like the following:

```c++
for (size_t row = 0; row < K; row++) {
    for (size_t col = 0; col < N_buffer / (2 * M) + 1; col++) {
        fft_mult[row][col] = fft_filt[row][col] * fft_data[row, col % (N_buffer / (2 * K) + 1)];
    }
}
```
## Inverse FFT

The inverse FFT gets us back to real samples, so it alters the dimensions to be:

| 0     | 1 | ... | N<sub>buffer</sub> / M - 1 |
|:-----:|---|-----|----------------------------|
| 1     |   |     |                            |
| ...   |   |     |                            |
| K - 1 |   |     |                            |

## Channel FFT

The next FFT is across the channels and basically performs the modulation to the center frequency of each channel. Because after this step we'll be doing things with discarding invalid samples and handling delays it actually makes those operations *much* easier if we transpose the output. Because C is row-major order that means that each row will correspond to an output time sample and each column is a channel. Again, since we have a real input that means that we go from $K$ channels to $K // 2 + 1$ output channels.

| 0                          | 1 | ... | K / 2 |
|:--------------------------:|---|-----|-------|
| 1                          |   |     |       |
| 2                          |   |     |       |
| 3                          |   |     |       |
| 4                          |   |     |       |
| ...                        |   |     |       |
| N<sub>buffer</sub> / M - 1 |   |     |       |
