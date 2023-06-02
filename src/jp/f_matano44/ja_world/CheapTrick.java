package jp.f_matano44.ja_world;

/**
 * Spectral envelope estimation on the basis of the idea of CheapTrick.
 */
public final class CheapTrick {
    private CheapTrick() {
        throw new IllegalStateException("CheapTrick isn't allowed to create instance.");
    }

    /** Struct for CheapTrick. */
    public static final class Option {
        public double q1;
        public double f0_floor;
        public int fft_size;

        /**
         * CheapTrick's Option.
         *
         * @param fs sampling rate
         */
        public Option(int fs) {
            // q1 is the parameter used for the spectral recovery.
            // Since The parameter is optimized, you don't need to change the parameter.
            q1 = -0.15;
            // f0_floor and fs are used to determine fft_size;
            // We strongly recommend not to change this value unless you have enough
            // knowledge of the signal processing in CheapTrick.
            f0_floor = ConstantNumbers.kFloorF0;
            fft_size = getFFTSizeForCheapTrick(fs, f0_floor);
        }
    }


    /**
     * CheapTrick() calculates the spectrogram that consists of spectral envelopes
     * estimated by CheapTrick.
     *
     *  @param x                    : Input signal
     *  @param fs                   : Sampling frequency
     *  @param temporal_positions   : Time axis
     *  @param f0                   : F0 contour

     *  @return double[][]
     *      spectrogram             : Spectrogram estimated by CheapTrick.
     */
    public static double[][] estimateSp(
        double[] x, double[] f0, double[] temporal_positions, int fs
    ) {
        return estimateSp(x, f0, temporal_positions, fs, new CheapTrick.Option(fs));
    }


    /**
     * CheapTrick() calculates the spectrogram that consists of spectral envelopes
     * estimated by CheapTrick.

     *  @param x                    : Input signal
     *  @param fs                   : Sampling frequency
     *  @param temporal_positions   : Time axis
     *  @param f0                   : F0 contour
     *  @param option               : Struct to order the parameter for CheapTrick

     *  @return double[][]
     *      spectrogram             : Spectrogram estimated by CheapTrick.
     */
    public static double[][] estimateSp(
        double[] x, double[] f0, double[] temporal_positions,
        int fs, CheapTrick.Option option
    ) {
        int fft_size = getFFTSizeForCheapTrick(fs, option.f0_floor);
        double[][] sp = new double[f0.length][fft_size / 2 + 1];

        CheapTrickMain(x.clone(), x.length, fs, temporal_positions.clone(),
            f0.clone(), f0.length, option, sp
        );

        return sp;
    }


    /**
     * GetFFTSizeForCheapTrick() calculates the FFT size based on the sampling
     * frequency and the lower limit of f0 (kFloorF0 defined in constantnumbers.h).
     *
     * @param fs Sampling frequency
     * @param f0_floor Lower f0 limit
     *
     * @return FFT size
     */
    public static int getFFTSizeForCheapTrick(final int fs, final double f0_floor) {
        return (int) (Math.pow(2.0, 1.0 +
            (int) (Math.log(3.0 * fs / f0_floor + 1) / ConstantNumbers.kLog2)));
    }


    /**
     * GetF0FloorForCheapTrick() calculates actual lower f0 limit for CheapTrick
     * based on the sampling frequency and FFT size used. Whenever f0 is below
     * this threshold the spectrum will be analyzed as if the frame is unvoiced
     * (using kDefaultF0 defined in constantnumbers.h).
     *
     * @param fs Sampling frequency
     * @param fft_size FFT size
     *
     * @return Lower f0 limit (Hz)
     */
    public static double getF0FloorForCheapTrick(int fs, int fft_size) {
        return 3.0 * fs / (fft_size - 3.0);
    }


    //-----------------------------------------------------------------------------
    // SmoothingWithRecovery() carries out the spectral smoothing and spectral
    // recovery on the Cepstrum domain.
    //-----------------------------------------------------------------------------
    private static void SmoothingWithRecovery(
        double f0, int fs, int fft_size, double q1,
        final Common.ForwardRealFFT forward_real_fft,
        final Common.InverseRealFFT inverse_real_fft, double[] spectral_envelope
    ) {
        double[] smoothing_lifter = new double[fft_size];
        double[] compensation_lifter = new double[fft_size];

        smoothing_lifter[0] = 1.0;
        compensation_lifter[0] = (1.0 - 2.0 * q1) + 2.0 * q1;
        double quefrency;
        for (int i = 1; i <= forward_real_fft.fft_size / 2; ++i) {
            quefrency = (double) (i) / fs;
            smoothing_lifter[i] = Math.sin(Math.PI * f0 * quefrency) /
            (Math.PI * f0 * quefrency);
            compensation_lifter[i] = (1.0 - 2.0 * q1) + 2.0 * q1 *
            Math.cos(2.0 * Math.PI * quefrency * f0);
        }

        for (int i = 0; i <= fft_size / 2; ++i)
            forward_real_fft.waveform[i] = Math.log(forward_real_fft.waveform[i]);
        for (int i = 1; i < fft_size / 2; ++i)
            forward_real_fft.waveform[fft_size - i] = forward_real_fft.waveform[i];
        Fft.fft_execute(forward_real_fft.forward_fft);

        for (int i = 0; i <= fft_size / 2; ++i) {
            inverse_real_fft.spectrum[i][0] = forward_real_fft.spectrum[i][0] *
            smoothing_lifter[i] * compensation_lifter[i] / fft_size;
            inverse_real_fft.spectrum[i][1] = 0.0;
        }
        Fft.fft_execute(inverse_real_fft.inverse_fft);

        for (int i = 0; i <= fft_size / 2; ++i)
            spectral_envelope[i] = Math.exp(inverse_real_fft.waveform[i]);
    }


    //-----------------------------------------------------------------------------
    // GetPowerSpectrum() calculates the power_spectrum with DC correction.
    // DC stands for Direct Current. In this case, the component from 0 to F0 Hz
    // is corrected.
    //-----------------------------------------------------------------------------
    private static void GetPowerSpectrum(
        int fs, double f0, int fft_size,
        final Common.ForwardRealFFT forward_real_fft
    ) {
        int half_window_length = (int) Math.round(1.5 * fs / f0);

        // FFT
        for (int i = half_window_length * 2 + 1; i < fft_size; ++i)
            forward_real_fft.waveform[i] = 0.0;
        Fft.fft_execute(forward_real_fft.forward_fft);

        // Calculation of the power spectrum.
        double[] power_spectrum = forward_real_fft.waveform;
        for (int i = 0; i <= fft_size / 2; ++i)
            power_spectrum[i] =
                forward_real_fft.spectrum[i][0] * forward_real_fft.spectrum[i][0]
                + forward_real_fft.spectrum[i][1] * forward_real_fft.spectrum[i][1];

        // DC correction
        Common.DCCorrection(power_spectrum, f0, fs, fft_size, power_spectrum);
    }


    //-----------------------------------------------------------------------------
    // SetParametersForGetWindowedWaveform()
    //-----------------------------------------------------------------------------
    private static void SetParametersForGetWindowedWaveform(int half_window_length,
        int x_length, double currnet_position, int fs, double current_f0,
        int[] base_index, int[] safe_index, double[] window
    ) {
        for (int i = -half_window_length; i <= half_window_length; ++i)
            base_index[i + half_window_length] = i;
        int origin = (int) Math.round(currnet_position * fs + 0.001);
        for (int i = 0; i <= half_window_length * 2; ++i)
            safe_index[i] =
                Common.MyMinInt(x_length - 1, Common.MyMaxInt(0, origin + base_index[i]));

        // Designing of the window function
        double average = 0.0;
        double position;
        for (int i = 0; i <= half_window_length * 2; ++i) {
            position = base_index[i] / 1.5 / fs;
            window[i] = 0.5 * Math.cos(Math.PI * position * current_f0) + 0.5;
            average += window[i] * window[i];
        }
        average = Math.sqrt(average);
        for (int i = 0; i <= half_window_length * 2; ++i) window[i] /= average;
    }


    //-----------------------------------------------------------------------------
    // GetWindowedWaveform() windows the waveform by F0-adaptive window
    //-----------------------------------------------------------------------------
    private static void GetWindowedWaveform(final double[] x, int x_length, int fs,
        double current_f0, double currnet_position,
        final Common.ForwardRealFFT forward_real_fft, MatlabFunctions.Random random
    ) {
        int half_window_length = (int) Math.round(1.5 * fs / current_f0);

        int[] base_index = new int[half_window_length * 2 + 1];
        int[] safe_index = new int[half_window_length * 2 + 1];
        double[] window  = new double[half_window_length * 2 + 1];

        SetParametersForGetWindowedWaveform(half_window_length, x_length,
            currnet_position, fs, current_f0, base_index, safe_index, window);

        // F0-adaptive windowing
        double[] waveform = forward_real_fft.waveform;
        for (int i = 0; i <= half_window_length * 2; ++i)
            waveform[i] = x[safe_index[i]] * window[i] +
            random.randn() * ConstantNumbers.kMySafeGuardMinimum;
        double tmp_weight1 = 0;
        double tmp_weight2 = 0;
        for (int i = 0; i <= half_window_length * 2; ++i) {
            tmp_weight1 += waveform[i];
            tmp_weight2 += window[i];
        }
        double weighting_coefficient = tmp_weight1 / tmp_weight2;
        for (int i = 0; i <= half_window_length * 2; ++i)
            waveform[i] -= window[i] * weighting_coefficient;
    }


    //-----------------------------------------------------------------------------
    // AddInfinitesimalNoise()
    //-----------------------------------------------------------------------------
    private static void AddInfinitesimalNoise(
        final double[] input_spectrum, int fft_size, double[] output_spectrum, 
        MatlabFunctions.Random random
    ) {
        for (int i = 0; i <= fft_size / 2; ++i)
            output_spectrum[i] = 
                input_spectrum[i] + Math.abs(random.randn()) * ConstantNumbers.kEps;
    }


    //-----------------------------------------------------------------------------
    // CheapTrickGeneralBody() calculates a spectral envelope at a temporal
    // position. This function is only used in CheapTrick().
    // Caution:
    //   forward_fft is allocated in advance to speed up the processing.
    //-----------------------------------------------------------------------------
    private static void CheapTrickGeneralBody(final double[] x, int x_length, int fs,
        double current_f0, int fft_size, double current_position, double q1,
        final Common.ForwardRealFFT forward_real_fft,
        final Common.InverseRealFFT inverse_real_fft, double[] spectral_envelope
        , MatlabFunctions.Random random
    ) {
        // F0-adaptive windowing
        GetWindowedWaveform(x, x_length, fs, current_f0, current_position,
            forward_real_fft, random);

        // Calculate power spectrum with DC correction
        // Note: The calculated power spectrum is stored in an array for waveform.
        // In this imprementation, power spectrum is transformed by FFT (NOT IFFT).
        // However, the same result is obtained.
        // This is tricky but important for simple implementation.
        GetPowerSpectrum(fs, current_f0, fft_size, forward_real_fft);

        // Smoothing of the power (linear axis)
        // forward_real_fft.waveform is the power spectrum.
        Common.LinearSmoothing(forward_real_fft.waveform, current_f0 * 2.0 / 3.0,
            fs, fft_size, forward_real_fft.waveform);

        // Add infinitesimal noise
        // This is a safeguard to avoid including zero in the spectrum.
        AddInfinitesimalNoise(forward_real_fft.waveform, fft_size,
            forward_real_fft.waveform, random);

        // Smoothing (log axis) and spectral recovery on the cepstrum domain.
        SmoothingWithRecovery(current_f0, fs, fft_size, q1, forward_real_fft,
            inverse_real_fft, spectral_envelope);
    }


    private static void CheapTrickMain(final double[] x, int x_length, int fs,
        final double[] temporal_positions, final double[] f0, int f0_length,
        final CheapTrick.Option option, double[][] spectrogram
    ) {
        MatlabFunctions.Random random = new MatlabFunctions.Random();
        int fft_size = option.fft_size;

        double f0_floor = getF0FloorForCheapTrick(fs, fft_size);
        double[] spectral_envelope = new double[fft_size];

        Common.ForwardRealFFT forward_real_fft = new Common.ForwardRealFFT(fft_size);
        Common.InverseRealFFT inverse_real_fft = new Common.InverseRealFFT(fft_size);

        double current_f0;
        for (int i = 0; i < f0_length; ++i) {
            current_f0 = f0[i] <= f0_floor ? ConstantNumbers.kDefaultF0 : f0[i];
            CheapTrickGeneralBody(x, x_length, fs, current_f0, fft_size,
                temporal_positions[i], option.q1, forward_real_fft,
                inverse_real_fft, spectral_envelope, random);
            for (int j = 0; j <= fft_size / 2; ++j)
            spectrogram[i][j] = spectral_envelope[j];
        }
    }

    /* void InitializeCheapTrickOption(int fs, CheapTrick.Option option) {
        // q1 is the parameter used for the spectral recovery.
        // Since The parameter is optimized, you don't need to change the parameter.
        option.q1 = -0.15;
        // f0_floor and fs are used to determine fft_size;
        // We strongly recommend not to change this value unless you have enough
        // knowledge of the signal processing in CheapTrick.
        option.f0_floor = ConstantNumbers.kFloorF0;
        option.fft_size = GetFFTSizeForCheapTrick(fs, option);
    } */
}