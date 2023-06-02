package jp.f_matano44.ja_world;

/**
 * F0 estimation based on instantaneous frequency.
 * This method is carried out by using the output of Dio().
 */
public class StoneMask {
    private StoneMask() {
        throw new IllegalStateException("StoneMask isn't allowed to create instance.");
    }


    /**
     * StoneMask() refines the estimated F0 by Dio.
     *
     * @param x Input signal
     * @param fs Sampling frequency
     * @param temporal_positions Temporal information
     * @param f0 f0 contour
     *
     * @return double[] refined_f0 : Refined F0
     */
    public static double[] refineF0(
        final double[] x, final double[] f0, final double[] temporal_positions, final int fs
    ) {
        double[] refined_f0 = new double[f0.length];
        stoneMaskMain(x.clone(), x.length, fs, temporal_positions.clone(),
            f0.clone(), f0.length, refined_f0);
        return refined_f0;
    }


    // -----------------------------------------------------------------------------
    // GetBaseIndex() calculates the temporal positions for windowing.
    // Since the result includes negative value and the value that exceeds the
    // length of the input signal, it must be modified appropriately.
    // -----------------------------------------------------------------------------
    private static void getBaseIndex(
        double current_position, final double[] base_time,
        int base_time_length, int fs, int[] index_raw
    ) {
        for (int i = 0; i < base_time_length; ++i) {
            index_raw[i] = (int) Math.round((current_position + base_time[i]) * fs);
        }
    }


    // -----------------------------------------------------------------------------
    // GetMainWindow() generates the window function.
    // -----------------------------------------------------------------------------
    private static void getMainWindow(
        double current_position, final int[] index_raw, int base_time_length,
        int fs, double window_length_in_time, double[] main_window
    ) {
        double tmp = 0.0;
        for (int i = 0; i < base_time_length; ++i) {
            tmp = (index_raw[i] - 1.0) / fs - current_position;
            main_window[i] = 0.42
                + 0.5 * Math.cos(2.0 * Math.PI * tmp / window_length_in_time)
                + 0.08 * Math.cos(4.0 * Math.PI * tmp / window_length_in_time);
        }
    }


    // -----------------------------------------------------------------------------
    // GetDiffWindow() generates the differentiated window.
    // Diff means differential.
    // -----------------------------------------------------------------------------
    private static void getDiffWindow(
        final double[] main_window, int base_time_length, double[] diff_window
    ) {
        diff_window[0] = -main_window[1] / 2.0;
        for (int i = 1; i < base_time_length - 1; ++i) {
            diff_window[i] = -(main_window[i + 1] - main_window[i - 1]) / 2.0;
        }
        diff_window[base_time_length - 1] = main_window[base_time_length - 2] / 2.0;
    }


    // -----------------------------------------------------------------------------
    // GetSpectra() calculates two spectra of the waveform windowed by windows
    // (main window and diff window).
    // -----------------------------------------------------------------------------
    private static void getSpectra(
        final double[] x, int x_length, int fft_size, final int[] index_raw,
        final double[] main_window, final double[] diff_window, int base_time_length,
        final Common.ForwardRealFFT forward_real_fft,
        /* fft_complex * */ double[][] main_spectrum,
        /* fft_complex * */ double[][] diff_spectrum
    ) {
        int[] index = new int[base_time_length];

        for (int i = 0; i < base_time_length; ++i) {
            index[i] = Common.myMaxInt(0, Common.myMinInt(x_length - 1, index_raw[i] - 1));
        }
        for (int i = 0; i < base_time_length; ++i) {
            forward_real_fft.waveform[i] = x[index[i]] * main_window[i];
        }
        for (int i = base_time_length; i < fft_size; ++i) {
            forward_real_fft.waveform[i] = 0.0;
        }

        Fft.fft_execute(forward_real_fft.forward_fft);
        for (int i = 0; i <= fft_size / 2; ++i) {
            main_spectrum[i][0] = forward_real_fft.spectrum[i][0];
            main_spectrum[i][1] = forward_real_fft.spectrum[i][1];
        }

        for (int i = 0; i < base_time_length; ++i) {
            forward_real_fft.waveform[i] = x[index[i]] * diff_window[i];
        }
        for (int i = base_time_length; i < fft_size; ++i) {
            forward_real_fft.waveform[i] = 0.0;
        }
        Fft.fft_execute(forward_real_fft.forward_fft);
        for (int i = 0; i <= fft_size / 2; ++i) {
            diff_spectrum[i][0] = forward_real_fft.spectrum[i][0];
            diff_spectrum[i][1] = forward_real_fft.spectrum[i][1];
        }
    }


    // -----------------------------------------------------------------------------
    // FixF0() fixed the F0 by instantaneous frequency.
    // -----------------------------------------------------------------------------
    private static double fixF0(
        final double[] power_spectrum, final double[] numerator_i,
        int fft_size, int fs, double initial_f0, int number_of_harmonics
    ) {
        double[] amplitude_list = new double[number_of_harmonics];
        double[] instantaneous_frequency_list = new double[number_of_harmonics];
        int index;
        for (int i = 0; i < number_of_harmonics; ++i) {
            index = Common.myMinInt(
                (int) Math.round(initial_f0 * fft_size / fs * (i + 1)),
                fft_size / 2);
            instantaneous_frequency_list[i] =
                power_spectrum[index] == 0.0 ? 0.0 : (double) (index) * fs / fft_size
                + numerator_i[index] / power_spectrum[index] * fs / 2.0 / Math.PI;
            amplitude_list[i] = Math.sqrt(power_spectrum[index]);
        }

        double denominator = 0.0;
        double numerator = 0.0;
        for (int i = 0; i < number_of_harmonics; ++i) {
            numerator += amplitude_list[i] * instantaneous_frequency_list[i];
            denominator += amplitude_list[i] * (i + 1);
        }
        return numerator / (denominator + ConstantNumbers.kMySafeGuardMinimum);
    }


    // -----------------------------------------------------------------------------
    // GetTentativeF0() calculates the F0 based on the instantaneous frequency.
    // -----------------------------------------------------------------------------
    private static double getTentativeF0(
        final double[] power_spectrum, final double[] numerator_i,
        int fft_size, int fs, double initial_f0
    ) {
        double tentative_f0 = fixF0(power_spectrum, numerator_i, fft_size, fs, initial_f0, 2);

        // If the fixed value is too large, the result will be rejected.
        if (tentative_f0 <= 0.0 || tentative_f0 > initial_f0 * 2) {
            return 0.0;
        }

        return fixF0(power_spectrum, numerator_i, fft_size, fs, tentative_f0, 6);
    }


    // -----------------------------------------------------------------------------
    // GetMeanF0() calculates the instantaneous frequency.
    // -----------------------------------------------------------------------------
    private static double getMeanF0(final double[] x, int x_length, int fs,
        double current_position, double initial_f0, int fft_size,
        double window_length_in_time, final double[] base_time, int base_time_length
    ) {
        int[] index_raw = new int[base_time_length];
        double[] main_window = new double[base_time_length];
        double[] diff_window = new double[base_time_length];
        getBaseIndex(current_position, base_time, base_time_length, fs, index_raw);
        getMainWindow(current_position, index_raw, base_time_length, fs,
                window_length_in_time, main_window);
        getDiffWindow(main_window, base_time_length, diff_window);

        Common.ForwardRealFFT forward_real_fft = new Common.ForwardRealFFT(fft_size);
        /* fft_complex * */ double[][] main_spectrum = new double[fft_size][2];
        /* fft_complex * */ double[][] diff_spectrum = new double[fft_size][2];
        getSpectra(x, x_length, fft_size, index_raw, main_window, diff_window,
                base_time_length, forward_real_fft, main_spectrum, diff_spectrum);

        double[] power_spectrum = new double[fft_size / 2 + 1];
        double[] numerator_i = new double[fft_size / 2 + 1];
        for (int j = 0; j <= fft_size / 2; ++j) {
            numerator_i[j] = main_spectrum[j][0] * diff_spectrum[j][1]
                - main_spectrum[j][1] * diff_spectrum[j][0];
            power_spectrum[j] = main_spectrum[j][0] * main_spectrum[j][0]
                + main_spectrum[j][1] * main_spectrum[j][1];
        }

        double tentative_f0 = getTentativeF0(power_spectrum, numerator_i,
                fft_size, fs, initial_f0);

        return tentative_f0;
    }


    // -----------------------------------------------------------------------------
    // GetRefinedF0() fixes the F0 estimated by Dio(). This function uses
    // instantaneous frequency.
    // -----------------------------------------------------------------------------
    private static double getRefinedF0(
        final double[] x, int x_length, int fs, double current_position, double initial_f0
    ) {
        if (initial_f0 <= ConstantNumbers.kFloorF0StoneMask || initial_f0 > fs / 12.0) {
            return 0.0;
        }
    
        int half_window_length = (int) (1.5 * fs / initial_f0 + 1.0);
        double window_length_in_time = (2.0 * half_window_length + 1.0) / fs;
        double[] base_time = new double[half_window_length * 2 + 1];
        for (int i = 0; i < half_window_length * 2 + 1; i++) {
            base_time[i] = (double) (-half_window_length + i) / fs;
        }
        int fft_size = (int) (Math.pow(2.0,
            2.0 + (int) (Math.log(half_window_length * 2.0 + 1.0) / ConstantNumbers.kLog2)));
    
        double mean_f0 = getMeanF0(x, x_length, fs, current_position,
            initial_f0, fft_size, window_length_in_time, base_time,
            half_window_length * 2 + 1);
    
        // If amount of correction is overlarge (20 %), initial F0 is employed.
        if (Math.abs(mean_f0 - initial_f0) > initial_f0 * 0.2) {
            mean_f0 = initial_f0;
        }

        return mean_f0;
    }


    //-----------------------------------------------------------------------------
    // StoneMask() refines the estimated F0 by Dio()
    //
    // Input:
    //   x                      : Input signal
    //   x_length               : Length of the input signal
    //   fs                     : Sampling frequency
    //   time_axis              : Temporal information
    //   f0                     : f0 contour
    //   f0_length              : Length of f0
    //
    // Output:
    //   refined_f0             : Refined F0
    //-----------------------------------------------------------------------------
    private static void stoneMaskMain(
        final double[] x, int x_length, int fs, final double[] temporal_positions,
        final double[] f0, int f0_length, double[] refined_f0
    ) {
        for (int i = 0; i < f0_length; i++) {
            refined_f0[i] = getRefinedF0(x, x_length, fs, temporal_positions[i], f0[i]);
        }
    }
}
