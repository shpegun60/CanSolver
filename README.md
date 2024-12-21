```cpp
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
