/*
 * CanSolver.h
 *
 *  Created on: Dec 10, 2024
 *      Author: Shpegun60
 *
 *
 *	Insights:
 *		remove_duplicates: Uses std::unordered_set for efficient duplicate elimination, relying on the custom std::hash and operator==.
 *		Sorting Functions: Clearly implement domain-specific sorting logic for easy adjustments or extensions.
 *		Test Function: CanSolverTest provides a well-structured example of how to use the solver.
 *
 *
 *		CAN BitTiming calculation.
 *
 * F_in:
 *     given input clock frequency of the CAN protocol controller
 * Prescaler:
 *     scales down the input clock for the internal time base:
 *
 * TimeQuantum = 1/F_in * Prescaler
 *     e.g.:
 *         F_in        = 80 MHz
 *         Prescaler   = 8
 *         TimeQuantum = 1/80000000*8 = 0.0000001 s = 100 ns
 *
 * Bit Time:
 *     Each bit is divided into 4 segments:
 *
 *     <---------------- T_bit ------------------->
 *     --------------------------------------------
 *     |  SYNC   |  PROP   |  PHASE1   !  PHASE2  |
 *     --------------------------------------------
 *               <-------TSEG1---------><--TSEG2-->
 *                                      <-SJW->
 *
 *     SYNC:   Synchronization segment
 *                 always 1 TimeQuantum long by the CAN standard
 *     PROP:   Propagation segment
 *     PHASE1: Phase segment1
 *     PHASE2: Phase segment2
 *     !:      Sample point: always at the end of Phase segment 1. Can be expressed as a
 *                 position in the BitTime in percentage.
 *     TSEG1:  Time segment 1: usually this can be configured in a CAN controller
 *                 minimum 1 TimeQuantum long
 *     TSEG2:  Time segment 2: usually this can be configured in a CAN controller
 *                 minimum 1 TimeQuantum long
 *     SJW:    (Re)Synchronization Jump Width: The controller can lengthen TSEG1 or
 *             shorten TSEG2 by this value at resynchronization.
 *
 * TimeQuanta/BitTime:
 *     The number of TimeQuanta in a bit time:
 *         TQ/BT = SYNC + TSEG1 + TSEG2 = 1 + TSEG1 + TSEG2
 *
 * SamplePoint:
 *     Position of the sample point in the bit time in percentage:
 *         SP = (SYNC + TSEG1) / TQ/BT * 100 %
 *
 * BaudRate:
 *     Comes from the TimeQuantum, and the TQ/BT
 *                          1
 *     BaudRate = ---------------------
 *                 TimeQuantum * TQ/BT
 *
 * So with the wbove example and with TSEG1 = 15 and TSEG2 = 4:
 *     TQ/BT = 1 + 15 + 4 = 20
 *     BaudRate = 1 / (20 * 100ns) = 500000 = 500 kbps
 *     SamplePoint = (1 + 15) / 20 * 100% = 80%
 *
 */

#ifndef STM32_TOOLS_CAN_CANSOLVER_H_
#define STM32_TOOLS_CAN_CANSOLVER_H_

#include "basic_types.h"
#include <vector>
#include <array>
#include <tuple>

#define CANSOLVER_PRINT 1
#define CANSOLVER_TEST 1

class CanSolver
{
public:
    CanSolver() = default;
    ~CanSolver() = default;

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

private:
    // method need for calculate3
    static int calculate3_helper(const Input &in, u8 SyncSegment, int sampl_pt, int tseg, int *tseg1, int *tseg2);
    Result calculate3_helper_2(const Input &in, double freq, u8 SyncSegment, long best_error, long error, int sampl_pt, int best_tseg, int *tseg1, int *tseg2, int best_brp);
    // method need for calculate4
    Result calculate4_helper(const Input &in, double freq, u8 SyncSegment, u32 BRP, u32 TQbits);

    // method need for calculate4_can or calculate4_fdcan
    static Result calculate_simple(const Input &in);

public:
    inline Result getBest() const { return m_best; }
    void remove_duplicates();                 			// Removes duplicate results to ensure the results list contains only unique configurations.
    void sort_prescalar();                    			// Sorts the results by prescaler value and frequency deviation.
    void sort_percent();                      			// Sorts the results by the sampling point percentage.
    void sort_diff();                         			// Sorts the results by the frequency deviation.
    void sort();                         				// Sorts the results

    inline reg getBenchmark() const { return counts; }

    static constexpr float cia_sample_point(const u32 bitrate)
    {
        float sampl_pt = 50.0;

        if (bitrate > 800000) {
            sampl_pt = 75.0;
        } else if (bitrate > 500000) {
            sampl_pt = 80.0;
        } else {
            sampl_pt = 87.5;
        }

        return sampl_pt;
    }

#if CANSOLVER_PRINT
    void print_results(int max_rows = 50);              // Prints the calculated results to the console, limited by max_rows. Displays the best result as well.
    void print_results(const Result);                   // Prints the calculated results to the console, limited by max_rows. Displays the best result as well.
#endif /* CANSOLVER_PRINT */


public:
    std::vector<Result> m_results;  // Container holding all valid results of the CAN timing calculations.
    Result m_best = {};             // Stores the best result based on predefined criteria (e.g., lowest prescaler and smallest frequency deviation).
    reg counts = 0;
};


CanSolver::Result CanSolver::calculate4_can(const Input &head)
{
    const Result result = calculate_simple(head);
    return result;
}

std::array<CanSolver::Result, 2> CanSolver::calculate4_fdcan(const Input &head, const Input &data)
{
    std::array<Result, 2> result = {};
    result[0] = calculate_simple(head);
    result[1] = calculate_simple(data);
    return result;
}

#if CANSOLVER_TEST
//provides a well-structured example of how to use the solver.
void CanSolverTest();
#endif /* CANSOLVER_TEST */

#endif /* STM32_TOOLS_CAN_CANSOLVER_H_ */
