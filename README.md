## Sites for calculation Can/FDCan Prescalars
* [Bit Timing Calculator for CAN FD :: kvaser](https://kvaser.com/support/calculators/can-fd-bit-timing-calculator/)
* [CAN Bus Bit Timing Calculator :: kvaser](https://kvaser.com/support/calculators/bit-timing-calculator/)
* [CAN Bit Time Calculation :: bittiming](http://www.bittiming.can-wiki.info/)
* [STM32G4 FDCAN :: Create by : phryniszak](https://phryniszak.github.io/stm32g-fdcan/)

All algorithms is tested on those site`s
## Algoritms
* **calulate0** - Brute Force Method - The most difficult method!!!! Use only when all methods have failed
* [calculate1](https://github.com/phryniszak/stm32g-fdcan)
* [calculate2](https://github.com/waszil/can_bit_timing_calculator/tree/master)
* [calculate3](https://github.com/rhyttr/SocketCAN/blob/master/can-utils/can-calc-bit-timing.c)
* [calculate4](https://github.com/shpegun60/MCP251XFD/blob/master/Core/MCP251XFD/MCP251XFD.c)

## Arguments for enable code parts:
 * **Qt in .pro file**: DEFINES += CANSOLVER_TEST=1 CANSOLVER_PRINT=1 CANSOLVER_MINIMUM=0
 * **Eclipse**: need add preprocessor variables to: compiler settings/preprocessor CANSOLVER_TEST=1 CANSOLVER_PRINT=1 CANSOLVER_MINIMUM=0
 * **CMake**: add_definitions(-DCANSOLVER_PRINT=1 -DCANSOLVER_TEST=1 -DCANSOLVER_MINIMUM=0)
 * **clang**: g++ -DCANSOLVER_PRINT=1 -DCANSOLVER_TEST=1 -DCANSOLVER_MINIMUM=0
 * **Microcontroller**: nothing to add, all hard function is disabled



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
     * Brute Force Method - The most difficult method!!!! Use only when all methods have failed
     * Estimated cost = (prescaler_max - prescaler_min) * (time_seg1_max - time_seg1_min) * (time_seg2_max - time_seg2_min)
     * From Shpegun60
     */
    const Result& calculate0(const Input& in);

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
## Test

```cpp
void CanSolver::test()
{
    CanSolver nom;
    CanSolver::Input nominal {
        100000000, // clock
        1,        // clock_divider
        1000000,  // rate
        75,       // sampling_point_min
        80,       // sampling_point_max
        1,        // prescaler_min
        512,      // prescaler_max
        1,        // time_seg1_min
        256,      // time_seg1_max
        1,        // time_seg2_min
        128,      // time_seg2_max
        1,        // sjw_min
        128,      // sjw_max
        0.1f      // baudrate_tolerance
    };


    {
        nom.calculate1(nominal);
        nom.remove_duplicates();
        nom.sort();
        std::cout << "Results for First Algititm:\n";
        nom.print_results();
        std::cout << "Benchmark: " << nom.getBenchmark() << "\n";
        std::cout << "\n";
    }

    {
        nom.calculate2(nominal);
        nom.remove_duplicates();
        nom.sort();
        std::cout << "Results for Second Algititm:\n";
        nom.print_results();
        std::cout << "Benchmark: " << nom.getBenchmark() << "\n";
        std::cout << "\n";
    }

    {
        nom.calculate3(nominal);
        nom.remove_duplicates();
        nom.sort();
        std::cout << "Results for Third Algititm:\n";
        nom.print_results();
        std::cout << "Benchmark: " << nom.getBenchmark() << "\n";
        std::cout << "\n";
    }

    {
        nom.calculate4(nominal);
        nom.remove_duplicates();
        nom.sort();
        std::cout << "Results for Four Algititm:\n";
        nom.print_results();
        std::cout << "Benchmark: " << nom.getBenchmark() << "\n";
        std::cout << "\n";
    }

    {
        std::cout << "Results for Simple CAN:\n";
        CanSolver can;
        auto res = CanSolver::calculate4_can(nominal);
        can.print_results(res);

        std::cout << "\nPassed: "<< (nom.getBest() == res)  << "\n\n";
    }

    {
        std::cout << "Results for FDCAN:\n";
        CanSolver can;
        auto res = CanSolver::calculate4_fdcan(nominal, nominal);
        can.print_results(res[0]);
        can.print_results(res[1]);

        std::cout << "\nHead Passed: "<< (nom.getBest() == res[0])  << "\n";
        std::cout << "Data Passed: "<< (nom.getBest() == res[0])  << "\n\n";
    }

}
```

## Output:
```
Results for First Algititm:
Rate[Mb]| Prescaler | TQ/bit | TimeSeg1 | TimeSeg2 | SJW  | SampleP[%] | Diff[%]
------------------------------------------------------------------------------- 
    1 |        25 |       4 |        2 |        1 |     1 |         75 |      0 
    1 |         5 |      20 |       14 |        5 |     5 |         75 |      0 
    1 |         1 |     100 |       74 |       25 |    25 |         75 |      0 
    1 |         4 |      25 |       18 |        6 |     6 |         76 |      0 
    1 |         2 |      50 |       37 |       12 |    12 |         76 |      0 
    1 |         1 |     100 |       75 |       24 |    24 |         76 |      0 
    1 |         1 |     100 |       76 |       23 |    23 |         77 |      0 
    1 |         2 |      50 |       38 |       11 |    11 |         78 |      0 
    1 |         1 |     100 |       77 |       22 |    22 |         78 |      0 
    1 |         1 |     100 |       78 |       21 |    21 |         79 |      0 
    1 |        20 |       5 |        3 |        1 |     1 |         80 |      0 
    1 |        10 |      10 |        7 |        2 |     2 |         80 |      0 
    1 |         5 |      20 |       15 |        4 |     4 |         80 |      0 
    1 |         4 |      25 |       19 |        5 |     5 |         80 |      0 
    1 |         2 |      50 |       39 |       10 |    10 |         80 |      0 
    1 |         1 |     100 |       79 |       20 |    20 |         80 |      0 

------
The best:
    1 |         1 |     100 |       74 |       25 |    25 |         75 |      0 

Benchmark: 1912

Results for Second Algititm:
Rate[Mb]| Prescaler | TQ/bit | TimeSeg1 | TimeSeg2 | SJW  | SampleP[%] | Diff[%]
------------------------------------------------------------------------------- 
    1 |        25 |       4 |        2 |        1 |     1 |         75 |      0
    1 |         5 |      20 |       14 |        5 |     5 |         75 |      0
    1 |         1 |     100 |       74 |       25 |    25 |         75 |      0
    1 |         4 |      25 |       18 |        6 |     6 |         76 |      0
    1 |         2 |      50 |       37 |       12 |    12 |         76 |      0
    1 |         1 |     100 |       75 |       24 |    24 |         76 |      0
    1 |         1 |     100 |       76 |       23 |    23 |         77 |      0
    1 |         2 |      50 |       38 |       11 |    11 |         78 |      0
    1 |         1 |     100 |       77 |       22 |    22 |         78 |      0
    1 |         1 |     100 |       78 |       21 |    21 |         79 |      0 
    1 |        20 |       5 |        3 |        1 |     1 |         80 |      0 
    1 |        10 |      10 |        7 |        2 |     2 |         80 |      0 
    1 |         5 |      20 |       15 |        4 |     4 |         80 |      0 
    1 |         4 |      25 |       19 |        5 |     5 |         80 |      0 
    1 |         2 |      50 |       39 |       10 |    10 |         80 |      0 
    1 |         1 |     100 |       79 |       20 |    20 |         80 |      0 

------
The best:
    1 |         1 |     100 |       74 |       25 |    25 |         75 |      0 

Benchmark: 32768

Results for Third Algititm:
Rate[Mb]| Prescaler | TQ/bit | TimeSeg1 | TimeSeg2 | SJW  | SampleP[%] | Diff[%]
------------------------------------------------------------------------------- 
    1 |        25 |       4 |        2 |        1 |     1 |         75 |      0
    1 |         5 |      20 |       14 |        5 |     5 |         75 |      0
    1 |         1 |     100 |       74 |       25 |    25 |         75 |      0
    1 |         4 |      25 |       18 |        6 |     6 |         76 |      0
    1 |         2 |      50 |       37 |       12 |    12 |         76 |      0
    1 |         1 |     100 |       75 |       24 |    24 |         76 |      0 
    1 |         1 |     100 |       76 |       23 |    23 |         77 |      0 
    1 |         2 |      50 |       38 |       11 |    11 |         78 |      0 
    1 |         1 |     100 |       77 |       22 |    22 |         78 |      0 
    1 |         1 |     100 |       78 |       21 |    21 |         79 |      0
    1 |        20 |       5 |        3 |        1 |     1 |         80 |      0 
    1 |        10 |      10 |        7 |        2 |     2 |         80 |      0 
    1 |         5 |      20 |       15 |        4 |     4 |         80 |      0 
    1 |         4 |      25 |       19 |        5 |     5 |         80 |      0 
    1 |         2 |      50 |       39 |       10 |    10 |         80 |      0 
    1 |         1 |     100 |       79 |       20 |    20 |         80 |      0 

------
The best:
    1 |         1 |     100 |       74 |       25 |    25 |         75 |      0 

Benchmark: 39066

Results for Four Algititm:
Rate[Mb]| Prescaler | TQ/bit | TimeSeg1 | TimeSeg2 | SJW  | SampleP[%] | Diff[%]
------------------------------------------------------------------------------- 
    1 |        25 |       4 |        2 |        1 |     1 |         75 |      0 
    1 |   20 |       5 |        3 |        1 |     1 |         80 |      0
    1 |        10 |      10 |        7 |        2 |     2 |         80 |      0 
    1 |         5 |      20 |       15 |        4 |     4 |         80 |      0 
    1 |         4 |      25 |       19 |        5 |     5 |         80 |      0 
    1 |         2 |      50 |       39 |       10 |    10 |         80 |      0 
    1 |         1 |     100 |       79 |       20 |    20 |         80 |      0 

------
The best:
    1 |         1 |     100 |       79 |       20 |    20 |         80 |      0 

Benchmark: 486

Results for Simple CAN:
Rate[Mb]| Prescaler | TQ/bit | TimeSeg1 | TimeSeg2 | SJW  | SampleP[%] | Diff[%]
------------------------------------------------------------------------------- 

------
The best:
    1 |         1 |     100 |       79 |       20 |    20 |         80 |      0 


Passed: 1

Results for FDCAN:
Rate[Mb]| Prescaler | TQ/bit | TimeSeg1 | TimeSeg2 | SJW  | SampleP[%] | Diff[%]
------------------------------------------------------------------------------- 

------
The best:
    1 |         1 |     100 |       79 |       20 |    20 |         80 |      0 

Rate[Mb]| Prescaler | TQ/bit | TimeSeg1 | TimeSeg2 | SJW  | SampleP[%] | Diff[%]
------------------------------------------------------------------------------- 

------
The best:
    1 |         1 |     100 |       79 |       20 |    20 |         80 |      0 


Head Passed: 1
Data Passed: 1


Process exited with code: 0
```
