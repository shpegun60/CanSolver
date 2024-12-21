```cpp

  struct Result {
        u32 bitrate;             // result bit rate [bit/sec]
        u32 prescaler;           // The clock prescaler value used for generating the CAN timing.
        u32 tq_per_bit;          // Number of Time Quanta (TQ) per bit in the CAN bit timing.
        u32 tBS1;                // Time Segment 1 (in TQ), which includes Propagation and Phase Segment 1.
        u32 tBS2;                // Time Segment 2 (in TQ), used for Phase Segment 2.
        u32 sync_jump_width;     // Synchronization Jump Width (SJW) in TQ, defining the maximum adjustment for resynchronization.
        float sample_point;      // The sample point location as a percentage [%] of the bit period.
        float freq_diff;         // The frequency deviation [%] between the calculated and desired bit rate.
        bool solutionFinded;	 // Solution finded indicator (only the best structure)

        // Class method for comparison
        bool operator==(const Result& other) const {
            // Compare all fields to determine equality between two results.
            return std::tie(prescaler, tq_per_bit, tBS1, tBS2, sync_jump_width, sample_point, freq_diff) ==
                   std::tie(other.prescaler, other.tq_per_bit, other.tBS1, other.tBS2, other.sync_jump_width, other.sample_point, other.freq_diff);
        }
    };


    struct Input {
        u32 clock;               		// The FDCAN peripheral's input clock frequency in Hz.
        u32 clock_divider;      		// Divider applied to the input clock to determine the CAN clock frequency.
        u32 rate;               		// Desired nominal bit rate in bits per second (bps).
        float sampling_point_min;       // Minimum allowed sampling point percentage [%] (0–100).
        float sampling_point_max;  		// Maximum allowed sampling point percentage [%] (0–100).

        u32 prescaler_min;              // Minimum value of the prescaler for clock division.
        u32 prescaler_max;              // Maximum value of the prescaler for clock division.

        u32 time_seg1_min;              // Minimum duration of Time Segment 1 in Time Quanta (TQ).
        u32 time_seg1_max;              // Maximum duration of Time Segment 1 in Time Quanta (TQ).

        u32 time_seg2_min;              // Minimum duration of Time Segment 2 in Time Quanta (TQ).
        u32 time_seg2_max;              // Maximum duration of Time Segment 2 in Time Quanta (TQ).

        u32 sjw_min;                    // Minimum Synchronization Jump Width (SJW) in Time Quanta (TQ).
        u32 sjw_max;                    // Maximum Synchronization Jump Width (SJW) in Time Quanta (TQ).

        float baudrate_tolerance = 1.0f; // Allowed tolerance for the deviation of the actual bit rate from the desired rate, in percent.
    };

  /*
     * Estimated cost = (prescaler_max - prescaler_min) * (sampling_point_min - sampling_point_max) / 0.1  --> FAST
	 * from https://github.com/phryniszak/stm32g-fdcan
	 */
    const Result& calculate1(const Input& in);          	// Computes all possible CAN timing configurations based on the input parameters.

    /*
	 * Estimated cost = (time_seg1_max - time_seg1_min) * (time_seg2_max - time_seg2_min)
	 *
	 * from https://github.com/waszil/can_bit_timing_calculator/tree/master
	 */
    const Result& calculate2(const Input& in);			// Computes all possible CAN timing configurations based on the input parameters.

    /*
     * Estimated cost = ???
     * from https://github.com/rhyttr/SocketCAN/blob/master/can-utils/can-calc-bit-timing.c
	 */
    const Result& calculate3(const Input& in);			// Computes all possible CAN timing configurations based on the input parameters.

    /*
     * Estimated cost = ??? ULTRA - FAST
     * from https://github.com/shpegun60/MCP251XFD/blob/master/Core/MCP251XFD/MCP251XFD.c
     */
    const Result& calculate4(const Input &in);

    /*
     * Methods for microcontrollers (based on calculate4)
     */
    inline static Result calculate4_can(const Input &head);
    inline static std::array<Result, 2> calculate4_fdcan(const Input &head, const Input &data);
```
