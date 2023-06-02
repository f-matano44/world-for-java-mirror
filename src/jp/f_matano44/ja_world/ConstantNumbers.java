package jp.f_matano44.ja_world;

/**
 * This file only defines constant numbers used for several function.
 */
final class ConstantNumbers {
    private ConstantNumbers() {
        throw new IllegalStateException("ConstantNumbers isn't allowed to create instance.");
    }

    // for Dio()
    static final double kCutOff = 50.0;

    // for StoneMask()
    static final double kFloorF0StoneMask = 40.0;

    // static final double kPi = 3.1415926535897932384;
    static final double kMySafeGuardMinimum = 0.000000000001;
    static final double kEps = 0.00000000000000022204460492503131;
    static final double kFloorF0 = 71.0;
    static final double kCeilF0 = 800.0;
    static final double kDefaultF0 = 500.0;
    static final double kLog2 = 0.69314718055994529;
    // Maximum standard deviation not to be selected as a best f0.
    static final double kMaximumValue = 100000.0;

    // Note to me (fs: 48000)
    // 71 Hz is the limit to maintain the FFT size at 2048.
    // If we use 70 Hz as FLOOR_F0, the FFT size of 4096 is required.

    // for D4C()
    static final int kHanning = 1;
    static final int kBlackman = 2;
    static final double kFrequencyInterval = 3000.0;
    static final double kUpperLimit = 15000.0;
    static final double kThreshold = 0.85;
    static final double kFloorF0D4C = 47.0;

    // for Codec (Mel scale)
    // S. Stevens & J. Volkmann,
    // The Relation of Pitch to Frequency: A Revised Scale,
    // American Journal of Psychology, vol. 53, no. 3, pp. 329-353, 1940.
    static final double kM0 = 1127.01048;
    static final double kF0 = 700.0;
    static final double kFloorFrequency = 40.0;
    static final double kCeilFrequency = 20000.0;

    // other default value
    static final double framePeriod = 5.0;
}
