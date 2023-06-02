package jp.f_matano44.ja_world;

/**
 * Matlab functions implemented for WORLD
 * Since these functions are implemented as the same function of Matlab,
 * the source code does not follow the style guide (Names of variables
 * and functions).
 * Please see the reference of Matlab to show the usage of functions.
*/
final class MatlabFunctions {
    private MatlabFunctions() {
        throw new IllegalStateException("MatlabFunction isn't allowed to create instance.");
    }


    //-----------------------------------------------------------------------------
    // FilterForDecimate() calculates the coefficients of low-pass filter and
    // carries out the filtering. This function is only used for decimate().
    //-----------------------------------------------------------------------------
    static void FilterForDecimate(final double[] x, int x_length, int r, double[] y) {
        // double a[3], b[2];  // filter Coefficients
        double[] a = new double[3];
        double[] b = new double[2];

        switch (r) {
            case 11:  // fs : 44100 (default)
                a[0] = 2.450743295230728;
                a[1] = -2.06794904601978;
                a[2] = 0.59574774438332101;
                b[0] = 0.0026822508007163792;
                b[1] = 0.0080467524021491377;
                break;
            case 12:  // fs : 48000
                a[0] = 2.4981398605924205;
                a[1] = -2.1368928194784025;
                a[2] = 0.62187513816221485;
                b[0] = 0.0021097275904709001;
                b[1] = 0.0063291827714127002;
                break;
            case 10:
                a[0] = 2.3936475118069387;
                a[1] = -1.9873904075111861;
                a[2] = 0.5658879979027055;
                b[0] = 0.0034818622251927556;
                b[1] = 0.010445586675578267;
                break;
            case 9:
                a[0] = 2.3236003491759578;
                a[1] = -1.8921545617463598;
                a[2] = 0.53148928133729068;
                b[0] = 0.0046331164041389372;
                b[1] = 0.013899349212416812;
                break;
            case 8:  // fs : 32000
                a[0] = 2.2357462340187593;
                a[1] = -1.7780899984041358;
                a[2] = 0.49152555365968692;
                b[0] = 0.0063522763407111993;
                b[1] = 0.019056829022133598;
                break;
            case 7:
                a[0] = 2.1225239019534703;
                a[1] = -1.6395144861046302;
                a[2] = 0.44469707800587366;
                b[0] = 0.0090366882681608418;
                b[1] = 0.027110064804482525;
                break;
            case 6:  // fs : 24000 and 22050
                a[0] = 1.9715352749512141;
                a[1] = -1.4686795689225347;
                a[2] = 0.3893908434965701;
                b[0] = 0.013469181309343825;
                b[1] = 0.040407543928031475;
                break;
            case 5:
                a[0] = 1.7610939654280557;
                a[1] = -1.2554914843859768;
                a[2] = 0.3237186507788215;
                b[0] = 0.021334858522387423;
                b[1] = 0.06400457556716227;
                break;
            case 4:  // fs : 16000
                a[0] = 1.4499664446880227;
                a[1] = -0.98943497080950582;
                a[2] = 0.24578252340690215;
                b[0] = 0.036710750339322612;
                b[1] = 0.11013225101796784;
                break;
            case 3:
                a[0] = 0.95039378983237421;
                a[1] = -0.67429146741526791;
                a[2] = 0.15412211621346475;
                b[0] = 0.071221945171178636;
                b[1] = 0.21366583551353591;
                break;
            case 2:  // fs : 8000
                a[0] = 0.041156734567757189;
                a[1] = -0.42599112459189636;
                a[2] = 0.041037215479961225;
                b[0] = 0.16797464681802227;
                b[1] = 0.50392394045406674;
                break;
            default:
                a[0] = 0.0;
                a[1] = 0.0;
                a[2] = 0.0;
                b[0] = 0.0;
                b[1] = 0.0;
        }
    
        // Filtering on time domain.
        double[] w = new double[3]; //[3] = {0.0, 0.0, 0.0};
        double wt;
        for (int i = 0; i < x_length; ++i) {
            wt = x[i] + a[0] * w[0] + a[1] * w[1] + a[2] * w[2];
            y[i] = b[0] * wt + b[1] * w[0] + b[1] * w[1] + b[0] * w[2];
            w[2] = w[1];
            w[1] = w[0];
            w[0] = wt;
        }
    }


    //-----------------------------------------------------------------------------
    // fftshift() swaps the left and right halves of input vector.
    // http://www.mathworks.com/help/matlab/ref/fftshift.html
    //
    // Input:
    //   x              : Input vector
    //   x_length       : Length of x
    //
    // Output:
    //   y              : Swapped vector x
    //
    // Caution:
    //   Lengths of index and edges must be the same.
    //-----------------------------------------------------------------------------
    static void fftshift(final double[] x, int x_length, double[] y) {
        for (int i = 0; i < x_length / 2; ++i) {
            y[i] = x[i + x_length / 2];
            y[i + x_length / 2] = x[i];
        }
    }


    //-----------------------------------------------------------------------------
    // histc() counts the number of values in vector x that fall between the
    // elements in the edges vector (which must contain monotonically
    // nondecreasing values). n is a length(edges) vector containing these counts.
    // No elements of x can be complex.
    // http://www.mathworks.co.jp/help/techdoc/ref/histc.html
    //
    // Input:
    //   x              : Input vector
    //   x_length       : Length of x
    //   edges          : Input matrix (1-dimension)
    //   edges_length   : Length of edges
    //
    // Output:
    //   index          : Result counted in vector x
    // Caution:
    //   Lengths of index and edges must be the same.
    //-----------------------------------------------------------------------------
    static void histc(final double[] x, int x_length, final double[] edges,
        int edges_length, int[] index) {
        int count = 1;
        int i = 0;

        for (; i < edges_length; ++i) {
            index[i] = 1;
            if (edges[i] >= x[0]) break;
        }

        for (; i < edges_length; ++i) {
            if (edges[i] < x[count]) {
                index[i] = count;
            } else {
                index[i--] = count++;
            }
            if (count == x_length) break;
        }
        count--;

        for (i++; i < edges_length; ++i)
            index[i] = count;
    }


    //-----------------------------------------------------------------------------
    // interp1() interpolates to find yi, the values of the underlying function Y
    // at the points in the vector or array xi. x must be a vector.
    // http://www.mathworks.co.jp/help/techdoc/ref/interp1.html
    //
    // Input:
    //   x          : Input vector (Time axis)
    //   y          : Values at x[n]
    //   x_length   : Length of x (Length of y must be the same)
    //   xi         : Required vector
    //   xi_length  : Length of xi (Length of yi must be the same)
    //
    // Output:
    //   yi         : Interpolated vector
    //-----------------------------------------------------------------------------
    static void interp1(final double[] x, final double[] y, int x_length, final double[] xi,
        int xi_length, double[] yi) {
        double[] h = new double[x_length - 1];
        int[] k = new int[xi_length];
        
        for (int i = 0; i < x_length - 1; ++i) h[i] = x[i + 1] - x[i];
        for (int i = 0; i < xi_length; ++i) {
            k[i] = 0;
        }
        
        histc(x, x_length, xi, xi_length, k);
        
        for (int i = 0; i < xi_length; ++i) {
            double s = (xi[i] - x[k[i] - 1]) / h[k[i] - 1];
            yi[i] = y[k[i] - 1] + s * (y[k[i]] - y[k[i] - 1]);
        }
    }


    //-----------------------------------------------------------------------------
    // decimate() carries out down sampling by both IIR and FIR filters.
    // Filter coeffiencts are based on FilterForDecimate().
    //
    // Input:
    //   x          : Input signal
    //   x_length   : Length of x
    //   r          : Coefficient used for down sampling
    //                (fs after down sampling is fs/r)
    // Output:
    //   y          : Output signal
    //-----------------------------------------------------------------------------
    static void decimate(final double[] x, int x_length, int r, double[] y) {
        final int kNFact = 9;
        double[] tmp1 = new double[x_length + kNFact * 2];
        double[] tmp2 = new double[x_length + kNFact * 2];

        for (int i = 0; i < kNFact; ++i)
            tmp1[i] = 2 * x[0] - x[kNFact - i];
        for (int i = kNFact; i < kNFact + x_length; ++i)
            tmp1[i] = x[i - kNFact];
        for (int i = kNFact + x_length; i < 2 * kNFact + x_length; ++i)
            tmp1[i] = 2 * x[x_length - 1] - x[x_length - 2 - (i - (kNFact + x_length))];

        FilterForDecimate(tmp1, 2 * kNFact + x_length, r, tmp2);
        for (int i = 0; i < 2 * kNFact + x_length; ++i)
            tmp1[i] = tmp2[2 * kNFact + x_length - i - 1];
        FilterForDecimate(tmp1, 2 * kNFact + x_length, r, tmp2);
        for (int i = 0; i < 2 * kNFact + x_length; ++i)
            tmp1[i] = tmp2[2 * kNFact + x_length - i - 1];

        int nout = (x_length - 1) / r + 1;
        int nbeg = r - r * nout + x_length;

        int count = 0;
        for (int i = nbeg; i < x_length + kNFact; i += r)
            y[count++] = tmp1[i + kNFact - 1];
    }


    //-----------------------------------------------------------------------------
    // matlab_round() calculates rounding.
    //
    // Input:
    //   x    : Input value
    //
    // Output:
    //   y    : Rounded value
    //-----------------------------------------------------------------------------
    // static int matlab_round(double x) {
    //     // return x > 0 ? static_cast<int>(x + 0.5) : static_cast<int>(x - 0.5);
    //     return x > 0 ? (int) (x + 0.5) : (int) (x - 0.5);
    // }
    // -> y = (int) Math.round(x)


    //-----------------------------------------------------------------------------
    // diff() calculates differences and approximate derivatives
    // http://www.mathworks.co.jp/help/techdoc/ref/diff.html
    //
    // Input:
    //   x          : Input signal
    //   x_length   : Length of x
    //
    // Output:
    //   y          : Output signal
    //-----------------------------------------------------------------------------
    static void diff(final double[] x, int x_length, double[] y) {
        for (int i = 0; i < x_length - 1; ++i) y[i] = x[i + 1] - x[i];
    }


    //-----------------------------------------------------------------------------
    // interp1Q() is the special case of interp1().
    // We can use this function, provided that All periods of x-axis is the same.
    //
    // Input:
    //   x          : Origin of the x-axis
    //   shift      : Period of the x-axis
    //   y          : Values at x[n]
    //   x_length   : Length of x (Length of y must be the same)
    //   xi         : Required vector
    //   xi_length  : Length of xi (Length of yi must be the same)
    //
    // Output:
    //   yi         : Interpolated vector
    //
    // Caution:
    //   Length of xi and yi must be the same.
    //-----------------------------------------------------------------------------
    static void interp1Q(double x, double shift, final double[] y, int x_length,
        final double[] xi, int xi_length, double[] yi) {
        double[] xi_fraction = new double[xi_length];
        double[] delta_y = new double[x_length];
        int[] xi_base = new int[xi_length];
        
        double delta_x = shift;
        for (int i = 0; i < xi_length; ++i) {
            // xi_base[i] = static_cast<int>((xi[i] - x) / delta_x);
            xi_base[i] = (int) ((xi[i] - x) / delta_x);
            xi_fraction[i] = (xi[i] - x) / delta_x - xi_base[i];
        }
        diff(y, x_length, delta_y);
        delta_y[x_length - 1] = 0.0;
        
        for (int i = 0; i < xi_length; ++i) {
            yi[i] = y[xi_base[i]] + delta_y[xi_base[i]] * xi_fraction[i];
        }
    }


    static final class Random {
        long g_randn_x = 123456789;
        long g_randn_y = 362436069;
        long g_randn_z = 521288629;
        long g_randn_w = 88675123;
        final long uintMax = (long) Math.pow(2, 32) - 1;

        public Random() {
            randn_reseed();
        }

        /**
         * randn() generates pseudorandom numbers based on xorShift method.
         *
         * @return A generated pseudorandom number
         */
        public double randn() {
            long t;
            t = g_randn_x ^ (g_randn_x << 11) & uintMax;
            g_randn_x = g_randn_y & uintMax;
            g_randn_y = g_randn_z & uintMax;
            g_randn_z = g_randn_w & uintMax;
            g_randn_w = ((g_randn_w ^ (g_randn_w >> 19)) ^ ((t ^ (t >> 8)))) & uintMax;

            long tmp = (g_randn_w >> 4) & uintMax;
            for (int i = 0; i < 11; ++i) {
                t = (g_randn_x ^ (g_randn_x << 11)) & uintMax;
                g_randn_x = g_randn_y & uintMax;
                g_randn_y = g_randn_z & uintMax;
                g_randn_z = g_randn_w & uintMax;
                g_randn_w = ((g_randn_w ^ (g_randn_w >> 19)) ^ (t ^ (t >> 8)))  & uintMax;
                tmp += g_randn_w >> 4;            
                tmp &= uintMax;
            }
            return (tmp / 268435456.0 - 6.0);
        }

        //-----------------------------------------------------------------------------
        // randn_reseed() forces to seed the pseudorandom generator using initial
        // values.
        //-----------------------------------------------------------------------------
        void randn_reseed() {
            g_randn_x = 123456789;
            g_randn_y = 362436069;
            g_randn_z = 521288629;
            g_randn_w = 88675123;
        }
    }


    //-----------------------------------------------------------------------------
    // fast_fftfilt() carries out the convolution on the frequency domain.
    //
    // Input:
    //   x                : Input signal
    //   x_length         : Length of x
    //   h                : Impulse response
    //   h_length         : Length of h
    //   fft_size         : Length of FFT
    //   forward_real_fft : Struct to speed up the forward FFT
    //   inverse_real_fft : Struct to speed up the inverse FFT
    //
    // Output:
    //   y                : Calculated result.
    //-----------------------------------------------------------------------------
    void fast_fftfilt(final double[] x, int x_length, final double[] h, int h_length,
        int fft_size, final Common.ForwardRealFFT forward_real_fft,
        final Common.InverseRealFFT inverse_real_fft, double[] y) {
        // fft_complex *x_spectrum = new fft_complex[fft_size];
        double[][] x_spectrum = new double[fft_size][2];

        for (int i = 0; i < x_length; ++i) {
            forward_real_fft.waveform[i] = x[i] / fft_size;
        }
        for (int i = x_length; i < fft_size; ++i) {
            forward_real_fft.waveform[i] = 0.0;
        }
        Fft.fft_execute(forward_real_fft.forward_fft);
        for (int i = 0; i <= fft_size / 2; ++i) {
            x_spectrum[i][0] = forward_real_fft.spectrum[i][0];
            x_spectrum[i][1] = forward_real_fft.spectrum[i][1];
        }

        for (int i = 0; i < h_length; ++i) {
            forward_real_fft.waveform[i] = h[i] / fft_size;
        }
        for (int i = h_length; i < fft_size; ++i) {
            forward_real_fft.waveform[i] = 0.0;
        }
        Fft.fft_execute(forward_real_fft.forward_fft);

        for (int i = 0; i <= fft_size / 2; ++i) {
            inverse_real_fft.spectrum[i][0] =
                x_spectrum[i][0] * forward_real_fft.spectrum[i][0]
                - x_spectrum[i][1] * forward_real_fft.spectrum[i][1];
            inverse_real_fft.spectrum[i][1] =
                x_spectrum[i][0] * forward_real_fft.spectrum[i][1]
                + x_spectrum[i][1] * forward_real_fft.spectrum[i][0];
        }
        Fft.fft_execute(inverse_real_fft.inverse_fft);

        for (int i = 0; i < fft_size; ++i) {
            y[i] = inverse_real_fft.waveform[i];
        }
        
        x_spectrum = null; //delete[] x_spectrum;
    }

    //-----------------------------------------------------------------------------
    // matlab_std() calculates the standard deviation of the input vector.
    //
    // Input:
    //   x          : Input vector
    //   x_length   : Length of x
    //
    // Output:
    //   Calculated standard deviation
    //-----------------------------------------------------------------------------
    double matlab_std(final double[] x, int x_length) {
        double average = 0.0;
        for (int i = 0; i < x_length; ++i) {
            average += x[i];
        }
        average /= x_length;

        double s = 0.0;
        for (int i = 0; i < x_length; ++i) {
            s += Math.pow(x[i] - average, 2.0);
        }
        s /= (x_length - 1);

        return Math.sqrt(s);
    }

}
