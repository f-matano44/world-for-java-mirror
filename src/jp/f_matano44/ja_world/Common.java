package jp.f_matano44.ja_world;

/** 
 * common.cpp includes functions used in at least two files.
 * (1) Common functions
 * (2) FFT, IFFT and minimum phase analysis.

 * In FFT analysis and minimum phase analysis,
 * Functions "Initialize*()" allocate the mamory.
 * Functions "Destroy*()" free the accolated memory.
 * FFT size is used for initialization, and structs are used to keep the memory.
 * Functions "GetMinimumPhaseSpectrum()" calculate minimum phase spectrum.
 * Forward and inverse FFT do not have the function "Get*()",
 * because forward FFT and inverse FFT can run in one step.
 */
final class Common {
    private Common() {
        throw new IllegalStateException("fft isn't allowed to create instance.");
    }


    private static void setParametersForLinearSmoothing(int boundary, int fft_size, int fs,
        double width, final double[] power_spectrum, double[] mirroring_spectrum,
        double[] mirroring_segment, double[] frequency_axis
    ) {
        for (int i = 0; i < boundary; ++i) {
            mirroring_spectrum[i] = power_spectrum[boundary - i];
        }
        for (int i = boundary; i < fft_size / 2 + boundary; ++i) {
            mirroring_spectrum[i] = power_spectrum[i - boundary];
        }
        for (int i = fft_size / 2 + boundary; i <= fft_size / 2 + boundary * 2; ++i) {
            mirroring_spectrum[i] =
                power_spectrum[fft_size / 2 - (i - (fft_size / 2 + boundary))];
        }

        mirroring_segment[0] = mirroring_spectrum[0] * fs / fft_size;
        for (int i = 1; i < fft_size / 2 + boundary * 2 + 1; ++i) {
            mirroring_segment[i] = 
                mirroring_spectrum[i] * fs / fft_size + mirroring_segment[i - 1];
        }

        for (int i = 0; i <= fft_size / 2; ++i) {
            frequency_axis[i] = (double) (i) / fft_size * fs - width / 2.0;
        }
    }


    //-----------------------------------------------------------------------------
    // Structs on FFT
    //-----------------------------------------------------------------------------
    // Forward FFT in the real sequence
    static final class ForwardRealFFT {
        int fft_size;
        double[] waveform;
        double[][] spectrum; //fft_complex *spectrum;
        Fft.Plan forward_fft;

        /**
         * These functions are used to speed up the processing.
         * Forward FFT
         */
        public ForwardRealFFT(int fft_size) {
            this.fft_size = fft_size;
            this.waveform = new double[fft_size];
            this.spectrum = new double[fft_size][2]; //fft_complex[fft_size];
            this.forward_fft = Fft.fft_plan_dft_r2c_1d(fft_size,
            this.waveform, this.spectrum, Fft.FFT_ESTIMATE);
        }
    }

    // Inverse FFT in the real sequence
    static final class InverseRealFFT {
        int fft_size;
        double[] waveform;
        double[][] spectrum; // fft_complex *spectrum;
        Fft.Plan inverse_fft;

        // Inverse FFT
        public InverseRealFFT(int fft_size) {
            this.fft_size = fft_size;
            this.waveform = new double[fft_size];
            this.spectrum = new double[fft_size][2]; // fft_complex[fft_size];
            this.inverse_fft = Fft.fft_plan_dft_c2r_1d(fft_size,
            this.spectrum, this.waveform, Fft.FFT_ESTIMATE);
        }
    }

    // Inverse FFT in the complex sequence
    static final class InverseComplexFFT {
        int fft_size;
        double[][] input;   // fft_complex *
        double[][] output;  // fft_complex *
        Fft.Plan inverse_fft;

        // Inverse FFT (Complex)
        public InverseComplexFFT(int fft_size) {
            this.fft_size = fft_size;
            this.input  = new double[fft_size][2]; /* fft_complex * */ 
            this.output = new double[fft_size][2]; /* fft_complex * */ 
            this.inverse_fft = Fft.fft_plan_dft_1d(fft_size,
            this.input, this.output, Fft.FFT_BACKWARD, Fft.FFT_ESTIMATE);
        }
    }

    // Minimum phase analysis from logarithmic power spectrum
    static final class MinimumPhaseAnalysis {
        int fft_size;
        double[] log_spectrum;
        double[][] minimum_phase_spectrum;  // fft_complex *minimum_phase_spectrum;
        double[][] cepstrum;                // fft_complex *cepstrum;
        Fft.Plan inverse_fft;
        Fft.Plan forward_fft;

        // Minimum phase analysis (This analysis uses FFT)
        public MinimumPhaseAnalysis(int fft_size) {
            this.fft_size = fft_size;
            this.log_spectrum = new double[fft_size];
            this.minimum_phase_spectrum = new double[fft_size][2];  // fft_complex * 
            this.cepstrum = new double[fft_size][2];                // fft_complex *
            this.inverse_fft = Fft.fft_plan_dft_r2c_1d(fft_size,
            this.log_spectrum, this.cepstrum, Fft.FFT_ESTIMATE);
            this.forward_fft = Fft.fft_plan_dft_1d(fft_size,
            this.cepstrum, this.minimum_phase_spectrum,
                Fft.FFT_FORWARD, Fft.FFT_ESTIMATE);
        }
    }


    //-----------------------------------------------------------------------------
    // GetSuitableFFTSize() calculates the suitable FFT size.
    // The size is defined as the minimum length whose length is longer than
    // the input sample.
    //
    // Input:
    //   sample : Length of the input signal
    //
    // Output:
    //   Suitable FFT size
    //-----------------------------------------------------------------------------
    static int getSuitableFFTSize(int sample) {
        return (int) (Math.pow(2.0,
        (int) (Math.log((double) (sample)) / ConstantNumbers.kLog2) + 1.0));
        // return static_cast<int>(pow(2.0,
        //     static_cast<int>(log(static_cast<double>(sample)) / world::kLog2) + 1.0));
    }


    //-----------------------------------------------------------------------------
    // These four functions are simple max() and min() function
    // for "int" and "double" type.
    //-----------------------------------------------------------------------------
    static int myMaxInt(int x, int y) {
        return x > y ? x : y;
    }
    
    static double myMaxDouble(double x, double y) {
        return x > y ? x : y;
    }
    
    static int myMinInt(int x, int y) {
        return x < y ? x : y;
    }
    
    static double myMinDouble(double x, double y) {
        return x < y ? x : y;
    }
    //-----------------------------------------------------------------------------
    // These functions are used in at least two different .cpp files


    //-----------------------------------------------------------------------------
    // DCCorrection interpolates the power under f0 Hz
    // and is used in CheapTrick() and D4C().
    //-----------------------------------------------------------------------------
    static void dcCorrection(
        final double[] input, double f0, int fs, int fft_size, double[] output
    ) {
        int upper_limit = 2 + /* static_cast<int> */ (int) (f0 * fft_size / fs);
        double[] low_frequency_replica = new double[upper_limit];
        double[] low_frequency_axis = new double[upper_limit];
        
        for (int i = 0; i < upper_limit; ++i) {
            low_frequency_axis[i] = (double) (i) * fs / fft_size;
        }
        
        int upper_limit_replica = upper_limit - 1;

        MatlabFunctions.interp1Q(f0 - low_frequency_axis[0],
            (double) (-1) * (fs) / fft_size, input, upper_limit + 1,
            low_frequency_axis, upper_limit_replica, low_frequency_replica);
        
        for (int i = 0; i < upper_limit_replica; ++i) {
            output[i] = input[i] + low_frequency_replica[i];
        }
    }


    //-----------------------------------------------------------------------------
    // LinearSmoothing() carries out the spectral smoothing by rectangular window
    // whose length is width Hz and is used in CheapTrick() and D4C().
    //-----------------------------------------------------------------------------
    static void linearSmoothing(final double[] input, double width, int fs, int fft_size,
        double[] output) {
        int boundary = /* static_cast<int> */ (int) (width * fft_size / fs) + 1;

        // These parameters are set by the other function.
        double[] mirroring_spectrum = new double[fft_size / 2 + boundary * 2 + 1];
        double[] mirroring_segment = new double[fft_size / 2 + boundary * 2 + 1];
        double[] frequency_axis = new double[fft_size / 2 + 1];
        setParametersForLinearSmoothing(boundary, fft_size, fs, width,
            input, mirroring_spectrum, mirroring_segment, frequency_axis);
        
        double[] low_levels = new double[fft_size / 2 + 1];
        double[] high_levels = new double[fft_size / 2 + 1];
        double origin_of_mirroring_axis = -(boundary - 0.5) * fs / fft_size;
        double discrete_frequency_interval
            = /* static_cast<double> */ (double) (fs) / fft_size;
        
        MatlabFunctions.interp1Q(origin_of_mirroring_axis, discrete_frequency_interval,
            mirroring_segment, fft_size / 2 + boundary * 2 + 1, frequency_axis,
            fft_size / 2 + 1, low_levels);
        
        for (int i = 0; i <= fft_size / 2; ++i) {
            frequency_axis[i] += width;
        }
        
        MatlabFunctions.interp1Q(origin_of_mirroring_axis, discrete_frequency_interval,
            mirroring_segment, fft_size / 2 + boundary * 2 + 1, frequency_axis,
            fft_size / 2 + 1, high_levels);
        
        for (int i = 0; i <= fft_size / 2; ++i) {
            output[i] = (high_levels[i] - low_levels[i]) / width;
        }
    }


    //-----------------------------------------------------------------------------
    // NuttallWindow() calculates the coefficients of Nuttall window whose length
    // is y_length and is used in Dio(), Harvest() and D4C().
    //-----------------------------------------------------------------------------
    static void nuttallWindow(int y_length, double[] y) {
        double tmp;
        for (int i = 0; i < y_length; ++i) {
            tmp  = i / (y_length - 1.0);
            y[i] = 0.355768 - 0.487396 * Math.cos(2.0 * Math.PI * tmp)
                + 0.144232 * Math.cos(4.0 * Math.PI * tmp)
                - 0.012604 * Math.cos(6.0 * Math.PI * tmp);
        }
    }


    //-----------------------------------------------------------------------------
    // GetSafeAperiodicity() limit the range of aperiodicity from 0.001 to
    // 0.999999999999 (1 - world::kMySafeGuardMinimum).
    //-----------------------------------------------------------------------------
    static double getSafeAperiodicity(double x) {
        return myMaxDouble(0.001, myMinDouble(0.999999999999, x));
    }


    static void getMinimumPhaseSpectrum(final MinimumPhaseAnalysis minimum_phase) {
        // Mirroring
        for (int i = minimum_phase.fft_size / 2 + 1;
            i < minimum_phase.fft_size; ++i
        ) {
            minimum_phase.log_spectrum[i] =
            minimum_phase.log_spectrum[minimum_phase.fft_size - i];
        }

        // This fft_plan carries out "forward" FFT.
        // To carriy out the Inverse FFT, the sign of imaginary part
        // is inverted after FFT.
        Fft.fft_execute(minimum_phase.inverse_fft);
        minimum_phase.cepstrum[0][1] *= -1.0;
        for (int i = 1; i < minimum_phase.fft_size / 2; ++i) {
            minimum_phase.cepstrum[i][0] *= 2.0;
            minimum_phase.cepstrum[i][1] *= -2.0;
        }
        minimum_phase.cepstrum[minimum_phase.fft_size / 2][1] *= -1.0;
        for (int i = minimum_phase.fft_size / 2 + 1;
            i < minimum_phase.fft_size; ++i) {
            minimum_phase.cepstrum[i][0] = 0.0;
            minimum_phase.cepstrum[i][1] = 0.0;
        }

        Fft.fft_execute(minimum_phase.forward_fft);

        // Since x is complex number, calculation of exp(x) is as following.
        // Note: This FFT library does not keep the aliasing.
        double tmp;
        for (int i = 0; i <= minimum_phase.fft_size / 2; ++i) {
            tmp = Math.exp(minimum_phase.minimum_phase_spectrum[i][0]
                / minimum_phase.fft_size);
            minimum_phase.minimum_phase_spectrum[i][0] =
                tmp * Math.cos(minimum_phase.minimum_phase_spectrum[i][1]
                / minimum_phase.fft_size);
            minimum_phase.minimum_phase_spectrum[i][1] = 
                tmp * Math.sin(minimum_phase.minimum_phase_spectrum[i][1]
                / minimum_phase.fft_size);
        }
    }
}
