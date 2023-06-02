package jp.f_matano44.ja_world;

import java.util.Arrays;

/** Band-aperiodicity estimation on the basis of the idea of D4C. */
public final class D4C {
    private D4C() {
        throw new IllegalStateException("D4C isn't allowed to create instance.");
    }

    /** Struct for D4C. */
    public static class Option {
        public double threshold;


        /** Option's constructor. */
        public Option() {
            this.threshold = ConstantNumbers.kThreshold;
        }
    }


    /**
     * D4C: Band-aperiodicity estimation.
     *
     * @param x Input signal
     * @param f0 F0 contour
     * @param temporal_positions Temporal positions
     * @param fs Sampling rate
     *
     * @return double[][] Band-aperiodicity
     */
    public static double[][] estimateAp(
        final double[] x, double[] f0, 
        final double[] temporal_positions, final int fs
    ) {
        int fft_size = CheapTrick.getFFTSizeForCheapTrick(fs, ConstantNumbers.kFloorF0);
        D4C.Option option = new D4C.Option();
        double[][] ap = new double[f0.length][fft_size / 2 + 1];
        d4cMain(
            x.clone(), x.length, fs, temporal_positions.clone(),
            f0.clone(), f0.length, fft_size, option, ap
        );
        return ap;
    }


    /**
     * D4C: Band-aperiodicity estimation.
     *
     * @param x Input signal
     * @param f0 F0 contour
     * @param temporal_positions Temporal positions
     * @param fs Sampling rate
     * @param fft_size CheapTrick.getFFTSizeForCheapTrick(fs, f0_floor)
     * @param option D4C.Option
     *
     * @return double[][] Band-aperiodicity
     */
    public static double[][] estimateAp(
        final double[] x, double[] f0, final double[] temporal_positions,
        final int fs, int fft_size, D4C.Option option
    ) {
        double[][] ap = new double[f0.length][fft_size / 2 + 1];
        d4cMain(
            x.clone(), x.length, fs, temporal_positions.clone(),
            f0.clone(), f0.length, fft_size, option, ap
        );
        return ap;
    }


    //-----------------------------------------------------------------------------
    // SetParametersForGetWindowedWaveform()
    //-----------------------------------------------------------------------------
    static void setParametersForGetWindowedWaveform(int half_window_length,
        int x_length, double current_position, int fs, double current_f0,
        int window_type, double window_length_ratio, int[] base_index,
        int[] safe_index, double[] window
    ) {
        for (int i = -half_window_length; i <= half_window_length; ++i) {
            base_index[i + half_window_length] = i;
        }
        final int origin = (int) Math.round(current_position * fs + 0.001);
        for (int i = 0; i <= half_window_length * 2; ++i) {
            safe_index[i] =
                Common.myMinInt(x_length - 1, Common.myMaxInt(0, origin + base_index[i]));
        }

        // Designing of the window function
        double position;
        if (window_type == ConstantNumbers.kHanning) {  // Hanning window
            for (int i = 0; i <= half_window_length * 2; ++i) {
                position = (2.0 * base_index[i] / window_length_ratio) / fs;
                window[i] = 0.5 * Math.cos(Math.PI * position * current_f0) + 0.5;
            }
        } else {  // Blackman window
            for (int i = 0; i <= half_window_length * 2; ++i) {
                position = (2.0 * base_index[i] / window_length_ratio) / fs;
                window[i] = 0.42 + 0.5 * Math.cos(Math.PI * position * current_f0)
                    + 0.08 * Math.cos(Math.PI * position * current_f0 * 2);
            }
        }
    }


    //-----------------------------------------------------------------------------
    // GetWindowedWaveform() windows the waveform by F0-adaptive window
    // In the variable window_type, 1: hanning, 2: blackman
    //-----------------------------------------------------------------------------
    static void getWindowedWaveform(final double[] x, int x_length, int fs,
        double current_f0, double current_position, int window_type,
        double window_length_ratio, double[] waveform, MatlabFunctions.Random random
    ) {
        int half_window_length = (int) Math.round(window_length_ratio * fs / current_f0 / 2.0);

        int[] base_index = new int[half_window_length * 2 + 1];
        int[] safe_index = new int[half_window_length * 2 + 1];
        double[] window  = new double[half_window_length * 2 + 1];

        setParametersForGetWindowedWaveform(half_window_length, x_length,
            current_position, fs, current_f0, window_type, window_length_ratio,
            base_index, safe_index, window);

        // F0-adaptive windowing
        for (int i = 0; i <= half_window_length * 2; ++i) {
            waveform[i] = x[safe_index[i]] * window[i]
                + random.randn() * ConstantNumbers.kMySafeGuardMinimum;
        }

        double tmp_weight1 = 0;
        double tmp_weight2 = 0;
        for (int i = 0; i <= half_window_length * 2; ++i) {
            tmp_weight1 += waveform[i];
            tmp_weight2 += window[i];
        }
        double weighting_coefficient = tmp_weight1 / tmp_weight2;
        for (int i = 0; i <= half_window_length * 2; ++i) {
            waveform[i] -= window[i] * weighting_coefficient;
        }
    }


    //-----------------------------------------------------------------------------
    // GetCentroid() calculates the energy centroid (see the book, time-frequency
    // analysis written by L. Cohen).
    //-----------------------------------------------------------------------------
    static void getCentroid(final double[] x, int x_length, int fs,
        double current_f0, int fft_size, double current_position,
        final Common.ForwardRealFFT forward_real_fft, double[] centroid,
        MatlabFunctions.Random random
    ) {
        for (int i = 0; i < fft_size; ++i) {
            forward_real_fft.waveform[i] = 0.0;
        }
        getWindowedWaveform(
            x, x_length, fs, current_f0,
            current_position, ConstantNumbers.kBlackman, 4.0,
            forward_real_fft.waveform, random
        );
        double power = 0.0;
        for (int i = 0; i <= (int) Math.round(2.0 * fs / current_f0) * 2; ++i) {
            power += forward_real_fft.waveform[i] * forward_real_fft.waveform[i];
        }
        for (int i = 0; i <= (int) Math.round(2.0 * fs / current_f0) * 2; ++i) {
            forward_real_fft.waveform[i] /= Math.sqrt(power);
        }

        Fft.fft_execute(forward_real_fft.forward_fft);
        double[] tmp_real = new double[fft_size / 2 + 1];
        double[] tmp_imag = new double[fft_size / 2 + 1];
        for (int i = 0; i <= fft_size / 2; ++i) {
            tmp_real[i] = forward_real_fft.spectrum[i][0];
            tmp_imag[i] = forward_real_fft.spectrum[i][1];
        }

        for (int i = 0; i < fft_size; ++i) {
            forward_real_fft.waveform[i] *= i + 1.0;
        }
        Fft.fft_execute(forward_real_fft.forward_fft);
        for (int i = 0; i <= fft_size / 2; ++i) {
            centroid[i] = forward_real_fft.spectrum[i][0] * tmp_real[i]
                + tmp_imag[i] * forward_real_fft.spectrum[i][1];
        }
    }


    //-----------------------------------------------------------------------------
    // GetStaticCentroid() calculates the temporally static energy centroid.
    // Basic idea was proposed by H. Kawahara.
    //-----------------------------------------------------------------------------
    static void getStaticCentroid(final double[] x, int x_length, int fs,
        double current_f0, int fft_size, double current_position,
        final Common.ForwardRealFFT forward_real_fft, double[] static_centroid,
        MatlabFunctions.Random random
    ) {
        double[] centroid1 = new double[fft_size / 2 + 1];
        double[] centroid2 = new double[fft_size / 2 + 1];

        getCentroid(x, x_length, fs, current_f0, fft_size,
            current_position - 0.25 / current_f0, forward_real_fft, centroid1, random);
        getCentroid(x, x_length, fs, current_f0, fft_size,
            current_position + 0.25 / current_f0, forward_real_fft, centroid2, random);

        for (int i = 0; i <= fft_size / 2; ++i) {
            static_centroid[i] = centroid1[i] + centroid2[i];
        }

        Common.dcCorrection(static_centroid, current_f0, fs, fft_size, static_centroid);
    }


    //-----------------------------------------------------------------------------
    // GetSmoothedPowerSpectrum() calculates the smoothed power spectrum.
    // The parameters used for smoothing are optimized in davance.
    //-----------------------------------------------------------------------------
    static void getSmoothedPowerSpectrum(final double[] x, int x_length, int fs,
        double current_f0, int fft_size, double current_position,
        final Common.ForwardRealFFT forward_real_fft,
        double[] smoothed_power_spectrum, MatlabFunctions.Random random
    ) {
        for (int i = 0; i < fft_size; ++i) {
            forward_real_fft.waveform[i] = 0.0;
        }
        getWindowedWaveform(x, x_length, fs, current_f0,
            current_position, ConstantNumbers.kHanning, 4.0,
            forward_real_fft.waveform, random
        );

        Fft.fft_execute(forward_real_fft.forward_fft);
        for (int i = 0; i <= fft_size / 2; ++i) {
            smoothed_power_spectrum[i] =
                forward_real_fft.spectrum[i][0] * forward_real_fft.spectrum[i][0]
                + forward_real_fft.spectrum[i][1] * forward_real_fft.spectrum[i][1];
        }
        Common.dcCorrection(smoothed_power_spectrum, current_f0, fs, fft_size,
            smoothed_power_spectrum);
        Common.linearSmoothing(smoothed_power_spectrum, current_f0, fs, fft_size,
            smoothed_power_spectrum);
    }


    //-----------------------------------------------------------------------------
    // GetStaticGroupDelay() calculates the temporally static group delay.
    // This is the fundamental parameter in D4C.
    //-----------------------------------------------------------------------------
    static void getStaticGroupDelay(final double[] static_centroid,
        final double[] smoothed_power_spectrum, int fs, double f0,
        int fft_size, double[] static_group_delay
    ) {
        for (int i = 0; i <= fft_size / 2; ++i) {
            static_group_delay[i] = static_centroid[i] / smoothed_power_spectrum[i];
        }
        Common.linearSmoothing(
            static_group_delay, f0 / 2.0, fs, fft_size, static_group_delay
        );

        double[] smoothed_group_delay = new double[fft_size / 2 + 1];
        Common.linearSmoothing(
            static_group_delay, f0, fs, fft_size, smoothed_group_delay
        );

        for (int i = 0; i <= fft_size / 2; ++i) {
            static_group_delay[i] -= smoothed_group_delay[i];
        }
    }


    //-----------------------------------------------------------------------------
    // GetCoarseAperiodicity() calculates the aperiodicity in multiples of 3 kHz.
    // The upper limit is given based on the sampling frequency.
    //-----------------------------------------------------------------------------
    static void getCoarseAperiodicity(
            final double[] static_group_delay, int fs,
        int fft_size, int number_of_aperiodicities, final double[] window,
        int window_length, final Common.ForwardRealFFT forward_real_fft, 
        double[] coarse_aperiodicity
    ) {
        int boundary = (int) Math.round(fft_size * 8.0 / window_length);
        int half_window_length = window_length / 2;

        for (int i = 0; i < fft_size; ++i) {
            forward_real_fft.waveform[i] = 0.0;
        }

        double[] power_spectrum = new double[fft_size / 2 + 1];
        int center;
        for (int i = 0; i < number_of_aperiodicities; ++i) {
            center = (int) (ConstantNumbers.kFrequencyInterval * (i + 1) * fft_size / fs);
            for (int j = 0; j <= half_window_length * 2; ++j) {
                forward_real_fft.waveform[j] =
                    static_group_delay[center - half_window_length + j] * window[j];
            }
            Fft.fft_execute(forward_real_fft.forward_fft);
            for (int j = 0; j <= fft_size / 2; ++j) {
                power_spectrum[j] =
                    forward_real_fft.spectrum[j][0] * forward_real_fft.spectrum[j][0]
                    + forward_real_fft.spectrum[j][1] * forward_real_fft.spectrum[j][1];
            }
            // std::sort(power_spectrum, power_spectrum + fft_size / 2 + 1);
            Arrays.sort(power_spectrum, 0, fft_size / 2 + 1);
            for (int j = 1; j <= fft_size / 2; ++j) {
                power_spectrum[j] += power_spectrum[j - 1];
            }
            coarse_aperiodicity[i] =
                10 * Math.log10(power_spectrum[fft_size / 2 - boundary - 1]
                / power_spectrum[fft_size / 2]);
        }
    }


    static double d4cLoveTrainSub(final double[] x, int fs, int x_length,
        double current_f0, double current_position, int f0_length, int fft_size,
        int boundary0, int boundary1, int boundary2,
        Common.ForwardRealFFT forward_real_fft, MatlabFunctions.Random random
    ) {
        double[] power_spectrum = new double[fft_size];

        int window_length = (int) Math.round(1.5 * fs / current_f0) * 2 + 1;
        getWindowedWaveform(x, x_length, fs, current_f0, current_position,
            ConstantNumbers.kBlackman, 3.0, 
            forward_real_fft.waveform, random
        );

        for (int i = window_length; i < fft_size; ++i) {
            forward_real_fft.waveform[i] = 0.0;
        }
        Fft.fft_execute(forward_real_fft.forward_fft);

        for (int i = 0; i <= boundary0; ++i) {
            power_spectrum[i] = 0.0;
        }
        for (int i = boundary0 + 1; i < fft_size / 2 + 1; ++i) {
            power_spectrum[i] =
                forward_real_fft.spectrum[i][0] * forward_real_fft.spectrum[i][0]
                + forward_real_fft.spectrum[i][1] * forward_real_fft.spectrum[i][1];
        }
        for (int i = boundary0; i <= boundary2; ++i) {
            power_spectrum[i] += +power_spectrum[i - 1];
        }

        // double aperiodicity0 = power_spectrum[boundary1] / power_spectrum[boundary2];
        // return aperiodicity0;
        return power_spectrum[boundary1] / power_spectrum[boundary2];
    }


    //-----------------------------------------------------------------------------
    // D4CLoveTrain() determines the aperiodicity with VUV detection.
    // If a frame was determined as the unvoiced section, aperiodicity is set to
    // very high value as the safeguard.
    // If it was voiced section, the aperiodicity of 0 Hz is set to -60 dB.
    //-----------------------------------------------------------------------------
    static void d4cLoveTrain(final double[] x, int fs, int x_length,
        final double[] f0, int f0_length, final double[] temporal_positions,
        double[] aperiodicity0, MatlabFunctions.Random random
    ) {
        double lowest_f0 = 40.0;
        int fft_size = (int) (Math.pow(2.0,
            1.0 + (int) (Math.log(3.0 * fs / lowest_f0 + 1) / ConstantNumbers.kLog2)));
        Common.ForwardRealFFT forward_real_fft = new Common.ForwardRealFFT(fft_size);

        // Cumulative powers at 100, 4000, 7900 Hz are used for VUV identification.
        int boundary0 = (int) (Math.ceil(100.0 * fft_size / fs));
        int boundary1 = (int) (Math.ceil(4000.0 * fft_size / fs));
        int boundary2 = (int) (Math.ceil(7900.0 * fft_size / fs));
        for (int i = 0; i < f0_length; ++i) {
            if (f0[i] == 0.0) {
                aperiodicity0[i] = 0.0;
                continue;
            }
            aperiodicity0[i] = d4cLoveTrainSub(x, fs, x_length,
            Common.myMaxDouble(f0[i], lowest_f0), temporal_positions[i], f0_length,
            fft_size, boundary0, boundary1, boundary2, forward_real_fft, random);
        }
    }


    //-----------------------------------------------------------------------------
    // D4CGeneralBody() calculates a spectral envelope at a temporal
    // position. This function is only used in D4C().
    // Caution:
    //   forward_fft is allocated in advance to speed up the processing.
    //-----------------------------------------------------------------------------
    static void d4cGeneralBody(final double[] x, int x_length, int fs,
        double current_f0, int fft_size, double current_position,
        int number_of_aperiodicities, final double[] window, int window_length,
        final Common.ForwardRealFFT forward_real_fft, double[] coarse_aperiodicity,
        MatlabFunctions.Random random
    ) {
        double[] static_centroid = new double[fft_size / 2 + 1];
        double[] smoothed_power_spectrum = new double[fft_size / 2 + 1];
        double[] static_group_delay = new double[fft_size / 2 + 1];

        getStaticCentroid(x, x_length, fs, current_f0, fft_size, current_position,
            forward_real_fft, static_centroid, random);
        getSmoothedPowerSpectrum(x, x_length, fs, current_f0, fft_size,
            current_position, forward_real_fft, smoothed_power_spectrum, random);
        getStaticGroupDelay(static_centroid, smoothed_power_spectrum,
            fs, current_f0, fft_size, static_group_delay);

        getCoarseAperiodicity(static_group_delay, fs, fft_size,
            number_of_aperiodicities, window, window_length, forward_real_fft,
            coarse_aperiodicity);

        // Revision of the result based on the F0
        for (int i = 0; i < number_of_aperiodicities; ++i) {
            coarse_aperiodicity[i] = Common.myMinDouble(0.0,
                coarse_aperiodicity[i] + (current_f0 - 100) / 50.0);
        }
    }


    static void initializeAperiodicity(
        int f0_length, int fft_size, double[][] aperiodicity
    ) {
        for (int i = 0; i < f0_length; ++i) {
            for (int j = 0; j < fft_size / 2 + 1; ++j) {
                aperiodicity[i][j] = 1.0 - ConstantNumbers.kMySafeGuardMinimum;
            }
        }
    }


    static void getAperiodicity(final double[] coarse_frequency_axis,
        final double[] coarse_aperiodicity, int number_of_aperiodicities,
        final double[] frequency_axis, int fft_size, double[] aperiodicity
    ) {
        MatlabFunctions.interp1(coarse_frequency_axis, coarse_aperiodicity,
            number_of_aperiodicities + 2, frequency_axis, fft_size / 2 + 1,
            aperiodicity);
        for (int i = 0; i <= fft_size / 2; ++i) {
            aperiodicity[i] = Math.pow(10.0, aperiodicity[i] / 20.0);
        }
    }


    private static void d4cMain(final double[] x, int x_length, int fs,
        final double[] temporal_positions, final double[] f0, int f0_length,
        int fft_size, final D4C.Option option, double[][] aperiodicity
    ) { 
        MatlabFunctions.Random random = new MatlabFunctions.Random();
        initializeAperiodicity(f0_length, fft_size, aperiodicity);

        int fft_size_d4c = (int) (Math.pow(2.0,
            1.0 + (int) (Math.log(4.0 * fs / ConstantNumbers.kFloorF0D4C + 1)
            / ConstantNumbers.kLog2)));

        Common.ForwardRealFFT forward_real_fft = new Common.ForwardRealFFT(fft_size_d4c);

        int number_of_aperiodicities =
            (int) (Common.myMinDouble(ConstantNumbers.kUpperLimit,
            fs / 2.0 - ConstantNumbers.kFrequencyInterval)
            / ConstantNumbers.kFrequencyInterval);
        // Since the window function is common in D4CGeneralBody(),
        // it is designed here to speed up.
        int window_length =
            (int) (ConstantNumbers.kFrequencyInterval * fft_size_d4c / fs) * 2 + 1;
        double[] window =  new double[window_length];
        Common.nuttallWindow(window_length, window);

        // D4C Love Train (Aperiodicity of 0 Hz is given by the different algorithm)
        double[] aperiodicity0 = new double[f0_length];
        d4cLoveTrain(x, fs, x_length, f0, f0_length, temporal_positions,
            aperiodicity0, random);

        double[] coarse_aperiodicity = new double[number_of_aperiodicities + 2];
        coarse_aperiodicity[0] = -60.0;
        coarse_aperiodicity[number_of_aperiodicities + 1] =
            -ConstantNumbers.kMySafeGuardMinimum;

        double[] coarse_frequency_axis = new double[number_of_aperiodicities + 2];
        for (int i = 0; i <= number_of_aperiodicities; ++i) {
            coarse_frequency_axis[i] = i * ConstantNumbers.kFrequencyInterval;
            coarse_frequency_axis[number_of_aperiodicities + 1] = fs / 2.0;
        }

        double[] frequency_axis = new double[fft_size / 2 + 1];
        for (int i = 0; i <= fft_size / 2; ++i) {
            frequency_axis[i] = (double) (i) * fs / fft_size;
        }

        for (int i = 0; i < f0_length; ++i) {
            if (f0[i] == 0 || aperiodicity0[i] <= option.threshold) {
                continue;
            }

            double[] caTemp = Arrays.copyOfRange(
                coarse_aperiodicity, 1, coarse_aperiodicity.length);
            d4cGeneralBody(x, x_length, fs, Common.myMaxDouble(ConstantNumbers.kFloorF0D4C, f0[i]),
                fft_size_d4c, temporal_positions[i], number_of_aperiodicities, window,
                window_length, forward_real_fft, caTemp /* &coarse_aperiodicity[1] */, random);
            System.arraycopy(caTemp, 0, coarse_aperiodicity, 1, caTemp.length);

            // Linear interpolation to convert the coarse aperiodicity into its
            // spectral representation.
            getAperiodicity(coarse_frequency_axis, coarse_aperiodicity,
                number_of_aperiodicities, frequency_axis, fft_size, aperiodicity[i]);
        }
    }
}
