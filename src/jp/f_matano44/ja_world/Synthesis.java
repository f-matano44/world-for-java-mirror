package jp.f_matano44.ja_world;

/**
 * Voice synthesis based on f0, spectrogram and aperiodicity.
 *  forward_real_fft, inverse_real_fft and minimum_phase are used to speed up.
 */
public class Synthesis {
    /**
     * Synthesis: generate signal.
     *
     * @param f0 f0
     * @param spectrogram spectrogram
     * @param aperiodicity aperiodicity
     * @param fs sampling rate
     *
     * @return double[] signal
     * 
     * @throws IllegalArgumentException
     *      f0.length != sp.length != ap.length
     *      and spectrogram[n].length != aperiodicity[n].length
     */
    public static double[] getSignal(
        final double[] f0, final double[][] spectrogram, 
        final double[][] aperiodicity, int fs
    ) throws IllegalArgumentException {
        return getSignal(f0, spectrogram, aperiodicity, fs, ConstantNumbers.framePeriod);
    }


    /**
     * Synthesis: generate signal.
     *
     * @param f0 f0
     * @param spectrogram spectrogram
     * @param aperiodicity aperiodicity
     * @param fs sampling rate
     * @param frame_period frame_period
     *
     * @return double[] signal
     * 
     * @throws IllegalArgumentException
     *      f0.length == sp.length == ap.length
     *      and spectrogram[n].length != aperiodicity[n].length
     */
    public static double[] getSignal(
        final double[] f0, final double[][] spectrogram, 
        final double[][] aperiodicity, final int fs, final double frame_period
    ) throws IllegalArgumentException {
        if (
            f0.length != spectrogram.length || f0.length != aperiodicity.length
        ) {
            throw new IllegalArgumentException(
                "Mismatched number of frames between F0 (" + f0.length + "), "
                + "spectrogram (" + spectrogram.length
                + ") and aperiodicty (" + aperiodicity.length + ")."
            );

        } else if (
            spectrogram[0].length != aperiodicity[0].length
        ) {
            throw new IllegalArgumentException(
                "Mismatched dimensionality (spec size) between "
                + "spectrogram (" + spectrogram[0].length + ") and "
                + "aperiodicity (" + aperiodicity[0].length + ")"
            );

        } else {
            final int fft_size = (spectrogram[0].length - 1) * 2;
            final int y_length = (int) ((f0.length - 1) * frame_period / 1000.0 * fs) + 1;
            double[] y = new double[y_length];

            synthesis(f0.clone(), f0.length, spectrogram.clone(), aperiodicity.clone(), 
                fft_size, frame_period, fs, y_length, y);

            return y;
        }
    }


    private static void getNoiseSpectrum(
        int noise_size, int fft_size, final Common.ForwardRealFFT forward_real_fft,
        MatlabFunctions.Random random
    ) {
        double average = 0.0;
        for (int i = 0; i < noise_size; ++i) {
            forward_real_fft.waveform[i] = random.randn();
            average += forward_real_fft.waveform[i];
        }

        average /= noise_size;
        for (int i = 0; i < noise_size; ++i) {
            forward_real_fft.waveform[i] -= average;
        }
        for (int i = noise_size; i < fft_size; ++i) {
            forward_real_fft.waveform[i] = 0.0;
        }
        Fft.fft_execute(forward_real_fft.forward_fft);
    }


    //-----------------------------------------------------------------------------
    // GetAperiodicResponse() calculates an aperiodic response.
    //-----------------------------------------------------------------------------
    private static void getAperiodicResponse(int noise_size, int fft_size,
        final double[] spectrum, final double[] aperiodic_ratio, double current_vuv,
        final Common.ForwardRealFFT forward_real_fft,
        final Common.InverseRealFFT inverse_real_fft,
        final Common.MinimumPhaseAnalysis minimum_phase,
        double[] aperiodic_response, MatlabFunctions.Random random
    ) {
        getNoiseSpectrum(noise_size, fft_size, forward_real_fft, random);

        if (current_vuv != 0.0) {
            for (int i = 0; i <= minimum_phase.fft_size / 2; ++i) {
                minimum_phase.log_spectrum[i] =
                    Math.log(spectrum[i] * aperiodic_ratio[i]) / 2.0;
            }
        } else {
            for (int i = 0; i <= minimum_phase.fft_size / 2; ++i) {
                minimum_phase.log_spectrum[i] = Math.log(spectrum[i]) / 2.0;
            }
        }
        Common.GetMinimumPhaseSpectrum(minimum_phase);

        for (int i = 0; i <= fft_size / 2; ++i) {
            inverse_real_fft.spectrum[i][0] =
                minimum_phase.minimum_phase_spectrum[i][0]
                * forward_real_fft.spectrum[i][0]
                - minimum_phase.minimum_phase_spectrum[i][1]
                * forward_real_fft.spectrum[i][1];
            inverse_real_fft.spectrum[i][1] =
                minimum_phase.minimum_phase_spectrum[i][0]
                * forward_real_fft.spectrum[i][1]
                + minimum_phase.minimum_phase_spectrum[i][1]
                * forward_real_fft.spectrum[i][0];
        }
        Fft.fft_execute(inverse_real_fft.inverse_fft);
        MatlabFunctions.fftshift(inverse_real_fft.waveform, fft_size, aperiodic_response);
    }


    //-----------------------------------------------------------------------------
    // RemoveDCComponent()
    //-----------------------------------------------------------------------------
    private static void removeDCComponent(
        final double[] periodic_response, int fft_size,
        final double[] dc_remover, double[] new_periodic_response
    ) {
        double dc_component = 0.0;
        for (int i = fft_size / 2; i < fft_size; ++i) {
            dc_component += periodic_response[i];
        }
        for (int i = 0; i < fft_size / 2; ++i) {
            new_periodic_response[i] = -dc_component * dc_remover[i];
        }
        for (int i = fft_size / 2; i < fft_size; ++i) {
            new_periodic_response[i] -= dc_component * dc_remover[i];
        }
    }


    //-----------------------------------------------------------------------------
    // GetSpectrumWithFractionalTimeShift() calculates a periodic spectrum with
    // the fractional time shift under 1/fs.
    //-----------------------------------------------------------------------------
    static void getSpectrumWithFractionalTimeShift(
        final int fft_size, final double coefficient, 
        final Common.InverseRealFFT inverse_real_fft
    ) {
        for (int i = 0; i <= fft_size / 2; ++i) {
            final double re = inverse_real_fft.spectrum[i][0];
            final double im = inverse_real_fft.spectrum[i][1];
            final double re2 = Math.cos(coefficient * i);
            final double im2 = Math.sqrt(1.0 - re2 * re2);  // sin(pshift)
        
            inverse_real_fft.spectrum[i][0] = re * re2 + im * im2;
            inverse_real_fft.spectrum[i][1] = im * re2 - re * im2;
        }
    }


    //-----------------------------------------------------------------------------
    // GetPeriodicResponse() calculates a periodic response.
    //-----------------------------------------------------------------------------
    static void getPeriodicResponse(int fft_size, final double[] spectrum,
        final double[] aperiodic_ratio, double current_vuv,
        final Common.InverseRealFFT inverse_real_fft,
        final Common.MinimumPhaseAnalysis minimum_phase, final double[] dc_remover,
        double fractional_time_shift, int fs, double[] periodic_response,
        MatlabFunctions.Random random
    ) {
        if (current_vuv <= 0.5 || aperiodic_ratio[0] > 0.999) {
            for (int i = 0; i < fft_size; ++i) {
                periodic_response[i] = 0.0;
            }
            return;
        }
        
        for (int i = 0; i <= minimum_phase.fft_size / 2; ++i) {
            minimum_phase.log_spectrum[i] =
                Math.log(spectrum[i] * (1.0 - aperiodic_ratio[i])
                + ConstantNumbers.kMySafeGuardMinimum) / 2.0;
        }
        Common.GetMinimumPhaseSpectrum(minimum_phase);
        
        for (int i = 0; i <= fft_size / 2; ++i) {
            inverse_real_fft.spectrum[i][0] =
                minimum_phase.minimum_phase_spectrum[i][0];
            inverse_real_fft.spectrum[i][1] =
                minimum_phase.minimum_phase_spectrum[i][1];
        }
        
        // apply fractional time delay of fractional_time_shift seconds
        // using linear phase shift
        double coefficient =
            2.0 * Math.PI * fractional_time_shift * fs / fft_size;
        getSpectrumWithFractionalTimeShift(fft_size, coefficient, inverse_real_fft);
        
        Fft.fft_execute(inverse_real_fft.inverse_fft);
        MatlabFunctions.fftshift(inverse_real_fft.waveform, fft_size, periodic_response);
        removeDCComponent(
            periodic_response, fft_size, dc_remover, periodic_response);
    }


    private static void getSpectralEnvelope(
        final double current_time, final double frame_period, final int f0_length,
        final double[][] spectrogram, final int fft_size, final double[] spectral_envelope
    ) {
        final int current_frame_floor = Common.MyMinInt(f0_length - 1,
            (int) (Math.floor(current_time / frame_period)));
        final int current_frame_ceil = Common.MyMinInt(f0_length - 1,
            (int) (Math.ceil(current_time / frame_period)));
        final double interpolation = current_time / frame_period - current_frame_floor;
    
        if (current_frame_floor == current_frame_ceil) {
            for (int i = 0; i <= fft_size / 2; ++i) {
                spectral_envelope[i] = Math.abs(spectrogram[current_frame_floor][i]);
            }
        } else {
            for (int i = 0; i <= fft_size / 2; ++i) {
                spectral_envelope[i] =
                    (1.0 - interpolation) * Math.abs(spectrogram[current_frame_floor][i])
                    + interpolation * Math.abs(spectrogram[current_frame_ceil][i]);
            }
        }
    }


    private static void getAperiodicRatio(double current_time, double frame_period,
        int f0_length, final double[][] aperiodicity, int fft_size,
        double[] aperiodic_spectrum
    ) {
        int current_frame_floor = Common.MyMinInt(f0_length - 1,
            (int) (Math.floor(current_time / frame_period)));
        int current_frame_ceil = Common.MyMinInt(f0_length - 1,
            (int) (Math.ceil(current_time / frame_period)));
        double interpolation = current_time / frame_period - current_frame_floor;
    
        if (current_frame_floor == current_frame_ceil) {
            for (int i = 0; i <= fft_size / 2; ++i) {
                aperiodic_spectrum[i] = Math.pow(
                    Common.GetSafeAperiodicity(aperiodicity[current_frame_floor][i]), 
                    2.0
                    );
            }
        } else {
            for (int i = 0; i <= fft_size / 2; ++i) {
                aperiodic_spectrum[i] = Math.pow((1.0 - interpolation)
                    * Common.GetSafeAperiodicity(aperiodicity[current_frame_floor][i])
                    + interpolation
                    * Common.GetSafeAperiodicity(aperiodicity[current_frame_ceil][i]),
                    2.0);
            }
        }
    }


    //-----------------------------------------------------------------------------
    // GetOneFrameSegment() calculates a periodic and aperiodic response at a time.
    //-----------------------------------------------------------------------------
    private static void getOneFrameSegment(double current_vuv, int noise_size,
        final double[][] spectrogram, int fft_size,
        final double[][] aperiodicity, int f0_length, double frame_period,
        double current_time, double fractional_time_shift, int fs,
        Common.ForwardRealFFT forward_real_fft,
        final Common.InverseRealFFT inverse_real_fft,
        final Common.MinimumPhaseAnalysis minimum_phase, final double[] dc_remover,
        double[] response, MatlabFunctions.Random random
    ) {
        double[] spectral_envelope = new double[fft_size];
        double[] aperiodic_ratio = new double[fft_size];
        getSpectralEnvelope(current_time, frame_period, f0_length, spectrogram,
            fft_size, spectral_envelope);
        getAperiodicRatio(current_time, frame_period, f0_length, aperiodicity,
            fft_size, aperiodic_ratio);

        // Synthesis of the periodic response
        double[] periodic_response = new double[fft_size];
        getPeriodicResponse(fft_size, spectral_envelope, aperiodic_ratio,
            current_vuv, inverse_real_fft, minimum_phase, dc_remover,
            fractional_time_shift, fs, periodic_response, random);


        // Synthesis of the aperiodic response
        double[] aperiodic_response = new double[fft_size];
        getAperiodicResponse(noise_size, fft_size, spectral_envelope,
            aperiodic_ratio, current_vuv, forward_real_fft,
            inverse_real_fft, minimum_phase, aperiodic_response, random);

        double sqrt_noise_size = Math.sqrt((double) (noise_size));
        for (int i = 0; i < fft_size; ++i) {
            response[i] =
                (periodic_response[i] * sqrt_noise_size + aperiodic_response[i])
                / fft_size;
        }
    }


    private static void getTemporalParametersForTimeBase(
        final double[] f0, int f0_length,
        int fs, int y_length, double frame_period, double lowest_f0,
        double[] time_axis, double[] coarse_time_axis, double[] coarse_f0,
        double[] coarse_vuv
    ) {
        for (int i = 0; i < y_length; ++i) {
            time_axis[i] = i / (double) (fs);
        }
        // the array 'coarse_time_axis' is supposed to have 'f0_length + 1' positions
        for (int i = 0; i < f0_length; ++i) {
            coarse_time_axis[i] = i * frame_period;
            coarse_f0[i] = f0[i] < lowest_f0 ? 0.0 : f0[i];
            coarse_vuv[i] = coarse_f0[i] == 0.0 ? 0.0 : 1.0;
        }
        coarse_time_axis[f0_length] = f0_length * frame_period;
        coarse_f0[f0_length] = coarse_f0[f0_length - 1] * 2 
            - coarse_f0[f0_length - 2];
        coarse_vuv[f0_length] = coarse_vuv[f0_length - 1] * 2
            - coarse_vuv[f0_length - 2];
    }


    private static int getPulseLocationsForTimeBase(
        final double[] interpolated_f0,
        final double[] time_axis, int y_length, int fs, double[] pulse_locations,
        int[] pulse_locations_index, double[] pulse_locations_time_shift
    ) {
        double[] total_phase = new double[y_length];
        double[] wrap_phase = new double[y_length];
        double[] wrap_phase_abs = new double[y_length - 1];
        total_phase[0] = 2.0 * Math.PI * interpolated_f0[0] / fs;
        wrap_phase[0] = fmod(total_phase[0], 2.0 * Math.PI);
        for (int i = 1; i < y_length; ++i) {
            total_phase[i] = total_phase[i - 1]
                + 2.0 * Math.PI * interpolated_f0[i] / fs;
            wrap_phase[i] = fmod(total_phase[i], 2.0 * Math.PI);
            wrap_phase_abs[i - 1] = Math.abs(wrap_phase[i] - wrap_phase[i - 1]);
        }

        int number_of_pulses = 0;
        for (int i = 0; i < y_length - 1; ++i) {
            if (wrap_phase_abs[i] > Math.PI) {
                pulse_locations[number_of_pulses] = time_axis[i];
                pulse_locations_index[number_of_pulses] = i;

                // calculate the time shift in seconds between exact fractional pulse
                // position and the integer pulse position (sample i)
                // as we don't have access to the exact pulse position, we infer it
                // from the point between sample i and sample i + 1 where the
                // accummulated phase cross a multiple of 2pi
                // this point is found by solving y1 + x * (y2 - y1) = 0 for x, where y1
                // and y2 are the phases corresponding to sample i and i + 1, offset so
                // they cross zero; x >= 0
                double y1 = wrap_phase[i] - 2.0 * Math.PI;
                double y2 = wrap_phase[i + 1];
                double x = -y1 / (y2 - y1);
                pulse_locations_time_shift[number_of_pulses] = x / fs;
        
                ++number_of_pulses;
            }
        }
        return number_of_pulses;
    }


    private static int getTimeBase(
        final double[] f0, int f0_length, int fs,
        double frame_period, int y_length, double lowest_f0,
        double[] pulse_locations, int[] pulse_locations_index,
        double[] pulse_locations_time_shift, double[] interpolated_vuv
    ) {
        double[] time_axis = new double[y_length];
        double[] coarse_time_axis = new double[f0_length + 1];
        double[] coarse_f0 = new double[f0_length + 1];
        double[] coarse_vuv = new double[f0_length + 1];
        getTemporalParametersForTimeBase(f0, f0_length, fs, y_length, frame_period,
            lowest_f0, time_axis, coarse_time_axis, coarse_f0, coarse_vuv);
        double[] interpolated_f0 = new double[y_length];
        MatlabFunctions.interp1(coarse_time_axis, coarse_f0, f0_length + 1,
            time_axis, y_length, interpolated_f0);
        MatlabFunctions.interp1(coarse_time_axis, coarse_vuv, f0_length + 1,
            time_axis, y_length, interpolated_vuv);

        for (int i = 0; i < y_length; ++i) {
            interpolated_vuv[i] = interpolated_vuv[i] > 0.5 ? 1.0 : 0.0;
            interpolated_f0[i] =
            interpolated_vuv[i] == 0.0 ? ConstantNumbers.kDefaultF0 : interpolated_f0[i];
        }

        int number_of_pulses = getPulseLocationsForTimeBase(interpolated_f0,
            time_axis, y_length, fs, pulse_locations, pulse_locations_index,
            pulse_locations_time_shift);

        return number_of_pulses;
    }

    private static double[] getDCRemover(int fft_size) {
        final double[] dc_remover = new double[fft_size];
        double dc_component = 0.0;

        for (int i = 0; i < fft_size / 2; ++i) {
            dc_remover[i] = 0.5
                - 0.5 * Math.cos(2.0 * Math.PI * (i + 1.0) / (1.0 + fft_size));
            dc_remover[fft_size - i - 1] = dc_remover[i];
            dc_component += dc_remover[i] * 2.0;
        }
        for (int i = 0; i < fft_size / 2; ++i) {
            dc_remover[i] /= dc_component;
            dc_remover[fft_size - i - 1] = dc_remover[i];
        }

        return dc_remover;
    }


    private static void synthesis(final double[] f0, int f0_length,
        final double[][] spectrogram, final double[][] aperiodicity,
        final int fft_size, double frame_period, int fs, int y_length, double[] y
    ) {
        MatlabFunctions.Random random = new MatlabFunctions.Random();
        double[] impulse_response = new double[fft_size];

        Common.MinimumPhaseAnalysis minimum_phase = new Common.MinimumPhaseAnalysis(fft_size);
        Common.InverseRealFFT inverse_real_fft = new Common.InverseRealFFT(fft_size);
        Common.ForwardRealFFT forward_real_fft = new Common.ForwardRealFFT(fft_size);

        double[] pulse_locations = new double[y_length];
        int[] pulse_locations_index = new int[y_length];
        double[] pulse_locations_time_shift = new double[y_length];
        double[] interpolated_vuv = new double[y_length];
        final int number_of_pulses = getTimeBase(f0, f0_length, fs, frame_period / 1000.0,
            y_length, fs / fft_size + 1.0, pulse_locations, pulse_locations_index,
            pulse_locations_time_shift, interpolated_vuv);

        double[] dc_remover = getDCRemover(fft_size);

        frame_period /= 1000.0;
        for (int i = 0; i < number_of_pulses; ++i) {
            final int noise_size = 
                pulse_locations_index[Common.MyMinInt(number_of_pulses - 1, i + 1)]
                - pulse_locations_index[i];

            getOneFrameSegment(interpolated_vuv[pulse_locations_index[i]], noise_size,
                spectrogram, fft_size, aperiodicity, f0_length, frame_period,
                pulse_locations[i], pulse_locations_time_shift[i], fs,
                forward_real_fft, inverse_real_fft, minimum_phase, dc_remover,
                impulse_response, random);

            final int offset = pulse_locations_index[i] - fft_size / 2 + 1;
            final int lower_limit = Common.MyMaxInt(0, -offset);
            final int upper_limit = Common.MyMinInt(fft_size, y_length - offset);
            for (int j = lower_limit; j < upper_limit; ++j) {
                final int index = j + offset;
                y[index] += impulse_response[j];
            }
        }
    }


    // JA-WORLD original function
    private static double fmod(double a, double b) {
        return a - b * Math.floor(a / b);
    }
}
