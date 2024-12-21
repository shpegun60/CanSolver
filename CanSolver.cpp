/*
 * CanSolver.cpp
 *
 *  Created on: Dec 10, 2024
 *      Author: Shpegun60
 */

#include "CanSolver.h"

#if CANSOLVER_PRINT || CANSOLVER_TEST
#include <iostream>
#include <iomanip>
#endif /* CANSOLVER_PRINT || CANSOLVER_TEST */

#include <algorithm>
#include <cmath>
#include <unordered_set>

const CanSolver::Result& CanSolver::calculate1(const Input &in)
{
    const double tfdcan_tq_clk = static_cast<double>(in.clock) / static_cast<double>(in.clock_divider);
    const u8 SyncSegment = 1;

    m_results.clear();
    m_best = {}; // Initialize with an empty value.
    counts = 0;

    if (std::isnan(tfdcan_tq_clk) || std::isnan(in.sampling_point_min) || std::isnan(in.sampling_point_max) || std::isnan(in.baudrate_tolerance)) {
        return m_best;
    }

    // main calc
    {
        const int sampleP_from = static_cast<int>(in.sampling_point_min * 10.0f);
        const int sampleP_to = static_cast<int>(in.sampling_point_max * 10.0f);

        for (u32 prescaler = in.prescaler_min; prescaler <= in.prescaler_max; ++prescaler) {
            const double ftq = tfdcan_tq_clk / prescaler; 		// Calculate frequency of Time Quanta (TQ) clock.
            const u32 tq_per_bit = static_cast<u32>(std::round(ftq / in.rate)); 	// Calculate the number of TQs per bit.

            if (std::isnan(ftq)) {
                ++counts;
                continue;
            }

            // can not find working solution with less then 4 tq per bit
            // NominalTimeSeg1Min(1) + NominalTimeSeg2Min(1) + tSyncSeg(1)
            // Minimum number of Time Quanta per bit (SyncSeg(1) + TSeg1(min1) + TSeg2(min1) = 4)
            if (tq_per_bit < 4) {
                ++counts;
                continue; // Skip if it does not meet the minimum requirement.
            }

            // Iterate over the sampling point range specified in the input
            for (int sampling_point = sampleP_from; sampling_point <= sampleP_to; ++sampling_point) {

                ++counts;
                if (std::isnan(sampling_point)) {
                    continue;
                }

                u32 tBS1 = static_cast<u32>(std::round((static_cast<double>(tq_per_bit) * static_cast<double>(sampling_point)) / 1000.0));
                const u32 tBS2 = tq_per_bit - tBS1; // Calculate Time Segment 2.

                // reduce by tSyncSeg
                tBS1 -= SyncSegment; // Subtract SyncSeg (1 TQ).

                // Ensure that the calculated TimeSeg1 and TimeSeg2 are within their respective limits.
                if (tBS1 < in.time_seg1_min || tBS1 > in.time_seg1_max || tBS2 < in.time_seg2_min || tBS2 > in.time_seg2_max) {
                    continue; // Skip invalid configurations.
                }

                // Calculate the Synchronization Jump Width (SJW).
                const u32 sync_jump_width = std::max(in.sjw_min, std::min(tBS2, in.sjw_max));

                // A maximum node-to-node oscillator variation of 1.7% is allowed.
                // https://www.calculatorsoup.com/calculators/algebra/percent-difference-calculator.php
                // Check for deviation in baud rate tolerance.
                const double tfq_real = tq_per_bit * in.rate; // Actual TQ frequency.
                const double freq_diff = ((ftq - tfq_real) * 2.0 * 100.0) / (ftq + tfq_real);

                // is baudrate in limits
                if (std::isnan(tfq_real) || std::isnan(freq_diff) || std::abs(freq_diff) > in.baudrate_tolerance) {
                    continue; // Skip if outside tolerance limits.
                }

                // Add the result to the list of valid results.
                const double sample_point_real = ((tBS1 + 1) * 100.0) / tq_per_bit;

                if (std::isnan(sample_point_real)) {
                    continue;
                }

                const u32 bitrate = static_cast<u32>(tfdcan_tq_clk) / (prescaler * (tBS1 + tBS2 + SyncSegment));

                Result result = {bitrate, prescaler, tq_per_bit, tBS1, tBS2, sync_jump_width, static_cast<float>(sample_point_real), static_cast<float>(freq_diff)};
                m_results.push_back(result);

                // Update the best result if:
                // - It's the first result or
                // - The new result has a smaller prescaler or
                // - The same prescaler but a smaller frequency deviation.
                if (m_best.prescaler == 0 ||
                    result.prescaler < m_best.prescaler ||
                    (result.prescaler == m_best.prescaler && std::abs(result.freq_diff) < std::abs(m_best.freq_diff))) {
                    m_best = result;
                    m_best.solutionFinded = true;
                }
            }
        }
    }

    return m_best;
}


const CanSolver::Result& CanSolver::calculate2(const Input &in)
{
    const double tfdcan_tq_clk = static_cast<double>(in.clock) / static_cast<double>(in.clock_divider);
    const u8 SyncSegment = 1;

    m_results.clear();
    m_best = {}; // Initialize with an empty value.
    counts = 0;

    if (std::isnan(tfdcan_tq_clk) || std::isnan(in.sampling_point_min) || std::isnan(in.sampling_point_max) || std::isnan(in.baudrate_tolerance)) {
        return m_best;
    }

    // main calc
    {
        for (u32 TS1 = in.time_seg1_min; TS1 <= in.time_seg1_max; ++TS1) {
            for (u32 TS2 = in.time_seg2_min; TS2 <= in.time_seg2_max; ++TS2) {
                const u32 time_quanta_per_bit_time = SyncSegment + TS1 + TS2;
                const double prescaler = std::round(tfdcan_tq_clk / (in.rate * time_quanta_per_bit_time));
                ++counts;

                // Minimum number of Time Quanta per bit (SyncSeg(1) + TSeg1(min1) + TSeg2(min1) = 4)
                if (time_quanta_per_bit_time < 4) {
                    continue; // Skip if it does not meet the minimum requirement.
                }

                if (std::isnan(prescaler) || prescaler < in.prescaler_min || prescaler > in.prescaler_max) {
                    continue;
                }

                const u32 prescaler_int 		= static_cast<u32>(prescaler);
                const double sample_point 		= (static_cast<double>(SyncSegment + TS1) / time_quanta_per_bit_time) * 100.0;

                if (std::isnan(sample_point) || sample_point < in.sampling_point_min || sample_point > in.sampling_point_max) {
                    continue;
                }


                // Check for deviation in baud rate tolerance.
                const double ftq = tfdcan_tq_clk / prescaler; 					// Calculate frequency of Time Quanta (TQ) clock.
                const double tfq_real = time_quanta_per_bit_time * in.rate; 	// Actual TQ frequency.
                const double freq_diff = ((ftq - tfq_real) * 2.0 * 100.0) / (ftq + tfq_real);

                if (std::isnan(ftq) || std::isnan(tfq_real) || std::isnan(freq_diff) || std::abs(freq_diff) > in.baudrate_tolerance) {
                    continue; // Skip if outside tolerance limits.
                }

                const u32 act_sjw = std::max(in.sjw_min, std::min(std::min(TS1, TS2), in.sjw_max));
                const u32 bitrate = static_cast<u32>(tfdcan_tq_clk) / (prescaler_int * (TS1 + TS2 + SyncSegment));

                const Result result = {bitrate, prescaler_int, time_quanta_per_bit_time, TS1, TS2, act_sjw,  static_cast<float>(sample_point), static_cast<float>(freq_diff)};
                m_results.push_back(result);

                // Update the best result if:
                // - It's the first result or
                // - The new result has a smaller prescaler or
                // - The same prescaler but a smaller frequency deviation.
                if (m_best.prescaler == 0 ||
                    result.prescaler < m_best.prescaler ||
                    (result.prescaler == m_best.prescaler && std::abs(result.freq_diff) < std::abs(m_best.freq_diff))) {
                    m_best = result;
                    m_best.solutionFinded = true;
                }
            }
        }
    }

    return m_best;
}

const CanSolver::Result& CanSolver::calculate3(const Input &in)
{
    const double tfdcan_tq_clk = static_cast<double>(in.clock) / static_cast<double>(in.clock_divider);
    const u8 SyncSegment = 1;
    m_results.clear();
    m_best = {}; // Initialize with an empty value.
    counts = 0;

    if (std::isnan(tfdcan_tq_clk) || std::isnan(in.sampling_point_min) || std::isnan(in.sampling_point_max) || std::isnan(in.baudrate_tolerance)) {
        return m_best;
    }

    // main calc
    {
        // constants
        const int sampleP_from = static_cast<int>(in.sampling_point_min * 10.0f);
        const int sampleP_to = static_cast<int>(in.sampling_point_max * 10.0f);

        /* tseg even = round down, odd = round up */
        const int tseg_from = (in.time_seg1_max + in.time_seg2_max) * 2 + SyncSegment;
        const int tseg_to =  (in.time_seg1_min + in.time_seg2_min) * 2;

        for (int sampl_pt = sampleP_from; sampl_pt <= sampleP_to; ++sampl_pt) {

            long rate, best_rate = 0;
            long best_error = 1000000000, error = 0;
            int best_tseg = 0, best_brp = 0, brp = 0;
            int tsegall, tseg = 0, tseg1 = 0, tseg2 = 0;
            int spt_error = 1000, spt = 0;

            for (tseg = tseg_from; tseg >= tseg_to; --tseg) {

                ++counts;
                tsegall = SyncSegment + tseg / 2;

                if (tsegall < 4) {
                    continue; // Skip if it does not meet the minimum requirement.
                }

                /* Compute all possible tseg choices (tseg=tseg1+tseg2) */
                brp = static_cast<int>(tfdcan_tq_clk) / (tsegall * in.rate) + tseg % 2;
                /* chose brp step which is possible in system */
                if ((brp < static_cast<int>(in.prescaler_min)) || (brp > static_cast<int>(in.prescaler_max))) {
                    continue;
                }

                rate = static_cast<int>(tfdcan_tq_clk) / (brp * tsegall);
                error = in.rate - rate;
                /* tseg brp biterror */
                if (error < 0) {
                    error = -error;
                }

                if (error > best_error) {
                    continue;
                }

                best_error = error;

                if (error == 0) {
                    spt = calculate3_helper(in, SyncSegment, sampl_pt, tseg / 2, &tseg1, &tseg2);

                    error = sampl_pt - spt;

                    if (error < 0) {
                        error = -error;
                    }

                    if (error > spt_error) {
                        continue;
                    }

                    spt_error = error;
                }

                best_tseg = tseg / 2;
                best_brp = brp;
                best_rate = rate;


                const Result result = calculate3_helper_2(in, tfdcan_tq_clk, SyncSegment, best_error, error, sampl_pt, best_tseg, &tseg1, &tseg2, best_brp);

                // Update the best result if:
                // - It's the first result or
                // - The new result has a smaller prescaler or
                // - The same prescaler but a smaller frequency deviation.
                if (m_best.prescaler == 0 ||
                    result.prescaler < m_best.prescaler ||
                    (result.prescaler == m_best.prescaler && std::abs(result.freq_diff) < std::abs(m_best.freq_diff))) {
                    m_best = result;
                    m_best.solutionFinded = true;
                }

                // if (error == 0) {
                //     break;
                // }
            }

            (void)best_rate;
        }
    }
    return m_best;
}

int CanSolver::calculate3_helper(const Input &in, u8 SyncSegment, int sampl_pt, int tseg, int *tseg1, int *tseg2)
{
    *tseg2 = tseg + SyncSegment - (sampl_pt * (tseg + SyncSegment)) / 1000;

    if (*tseg2 < static_cast<int>(in.time_seg2_min)) {
        *tseg2 = static_cast<int>(in.time_seg2_min);
    }

    if (*tseg2 > static_cast<int>(in.time_seg2_max)) {
        *tseg2 = static_cast<int>(in.time_seg2_max);
    }

    *tseg1 = tseg - *tseg2;

    if (*tseg1 > static_cast<int>(in.time_seg1_max)) {
        *tseg1 = static_cast<int>(in.time_seg1_max);
        *tseg2 = tseg - *tseg1;
    }

    return 1000 * (tseg + SyncSegment - *tseg2) / (tseg + 1);
}

CanSolver::Result CanSolver::calculate3_helper_2(const Input &in, double freq, u8 SyncSegment, long best_error, long error, int sampl_pt, int best_tseg, int *tseg1, int *tseg2, int best_brp)
{
    double freq_diff = 0.0;

    if (best_error) {
        /* Error in one-tenth of a percent */
        error = (best_error * 1000) / in.rate;
        freq_diff = error / 10.0f;
    }

    // is baudrate in limits
    if (std::isnan(freq_diff) || std::abs(freq_diff) > in.baudrate_tolerance) {
        return m_best; // Skip if outside tolerance limits.
    }

    /* real sample point */
    /*const u32 sample_point = */calculate3_helper(in, SyncSegment, sampl_pt, best_tseg, tseg1, tseg2);
    const u32 tq = *tseg1 + *tseg2 + SyncSegment;
    //const u32 prop_seg = *tseg1 / 2;
    // const u32 phase_seg1 = *tseg1 - prop_seg;
    // const u32 phase_seg2 = *tseg2;
    const u32 sjw = std::max(in.sjw_min, std::min(static_cast<u32>(std::min(*tseg1, *tseg2)), in.sjw_max));
    const u32 new_brp = best_brp;

    /* real bit-rate */
    const u32 bitrate = static_cast<u32>(freq) / (new_brp * (*tseg1 + *tseg2 + 1));
    const double sample_point_real = (static_cast<double>(SyncSegment + (*tseg1)) / tq) * 100.0;

    // Check for deviation in baud rate tolerance.
    const double ftq = freq / new_brp;                          // Calculate frequency of Time Quanta (TQ) clock.
    const double tfq_real = tq * in.rate;                       // Actual TQ frequency.
    const double freq_diff1 = ((ftq - tfq_real) * 2.0 * 100.0) / (ftq + tfq_real);

    // is baudrate in limits
    if (std::isnan(freq_diff1) || std::abs(freq_diff1) > in.baudrate_tolerance) {
        return m_best; // Skip if outside tolerance limits.
    }

    const Result result = {bitrate, new_brp, tq, static_cast<u32>(*tseg1), static_cast<u32>(*tseg2), sjw,  static_cast<float>(sample_point_real), static_cast<float>(std::max(freq_diff, freq_diff1))};
    m_results.push_back(result);

    return result;
}


const CanSolver::Result& CanSolver::calculate4(const Input &in)
{
    const double tfdcan_tq_clk = static_cast<double>(in.clock) / static_cast<double>(in.clock_divider);
    const u8 SyncSegment = 1;

    m_results.clear();
    m_best = {}; // Initialize with an empty value.
    counts = 0;

    if (std::isnan(tfdcan_tq_clk) || std::isnan(in.sampling_point_min) || std::isnan(in.sampling_point_max) || std::isnan(in.baudrate_tolerance)) {
        return m_best;
    }

    // main calc
    {
        // constants
        const u32 TQBIT_MIN = std::max(static_cast<u32>(SyncSegment + in.time_seg1_min + in.time_seg2_min), static_cast<u32>(4));
        const u32 TQBIT_MAX = SyncSegment + in.time_seg1_max + in.time_seg2_max;

        //--- Declaration -----------------------------------------
        u32 ErrorTQ, ErrorNTQ;
        u32 BestBRP = in.prescaler_max, BestNTQbits = TQBIT_MAX;

        //--- Calculate Nominal & Data Bit Time parameter ---------
        u32 MinErrorBR = UINT32_MAX;
        u32 BRP = in.prescaler_max;                                                 // Select the worst BRP value. Here all value from max to min will be tested to get the best tuple of NBRP and DBRP, identical TQ in both phases prevents quantization errors during bit rate switching

        while(--BRP >= in.prescaler_min) {
            const u32 NTQbits = static_cast<u32>(tfdcan_tq_clk) / in.rate / BRP;    // Calculate the NTQbits according to BRP and the desired Nominal Bitrate

            if ((NTQbits < TQBIT_MIN) || (NTQbits > TQBIT_MAX)) {
                ++counts;
                continue;   // This TQbits count is not possible with this BRP, then do the next BRP value
            }

            // NTQ & DTQ bits count
            ErrorNTQ = (static_cast<int>(tfdcan_tq_clk) - (in.rate * NTQbits * BRP));     // Calculate NTQ error
            ErrorTQ = ErrorNTQ;

            if (ErrorTQ <= MinErrorBR) {    // If better error then
                // Save best parameters
                MinErrorBR      = ErrorTQ;
                BestBRP         = BRP;
                BestNTQbits     = NTQbits;
                calculate4_helper(in, tfdcan_tq_clk, SyncSegment, BestBRP, BestNTQbits);
            }

            // NTQ+1 & DTQ bits count
            if (NTQbits < TQBIT_MAX) {
                ErrorNTQ = ((in.rate * (NTQbits+1) * BRP) - static_cast<u32>(tfdcan_tq_clk)); // Calculate NTQ error with NTQbits+1
                ErrorTQ = ErrorNTQ;

                if (ErrorTQ <= MinErrorBR) { // If better error then
                    // Save best parameters
                    MinErrorBR = ErrorTQ;
                    BestBRP = BRP;
                    BestNTQbits = NTQbits+1;
                    calculate4_helper(in, tfdcan_tq_clk, SyncSegment, BestBRP, BestNTQbits);
                }
            }
        }

        if (MinErrorBR == UINT32_MAX) {
            return m_best;          // Impossible to find a good BRP
        }

        //--- Calculate segments --------------------------
        m_best = calculate4_helper(in, tfdcan_tq_clk, SyncSegment, BestBRP, BestNTQbits);
        m_best.solutionFinded = true;
    }

    return m_best;
}


CanSolver::Result CanSolver::calculate_simple(const Input &in)
{
    const u32 fsysclk = in.clock / in.clock_divider;
    const u8 SyncSegment = 1;

    // main calc
    {
        // Time quanta boarders
        const u32 TQBIT_MIN = std::max(static_cast<u32>(SyncSegment + in.time_seg1_min + in.time_seg2_min), static_cast<u32>(4));
        const u32 TQBIT_MAX = SyncSegment + in.time_seg1_max + in.time_seg2_max;

        // Bit rate prescalar boarders
        const u32 BRP_MIN = in.prescaler_min;
        const u32 BRP_MAX = in.prescaler_max;

        // Time seq 1 boarders
        const u32 TSEG1_MIN = in.time_seg1_min;
        const u32 TSEG1_MAX = in.time_seg1_max;

        // Time seq 2 boarders
        const u32 TSEG2_MIN = in.time_seg2_min;
        const u32 TSEG2_MAX = in.time_seg2_max;

        // Time sjw boarders
        const u32 SJW_MIN = in.sjw_min;
        const u32 SJW_MAX = in.sjw_max;

        // BitRate
        const u32 desiredNominalBitrate = in.rate;


        //--- Declaration -----------------------------------------
        u32 ErrorTQ, ErrorNTQ;
        u32 BestBRP = BRP_MAX, BestNTQbits = TQBIT_MAX;

        //--- Calculate Nominal & Data Bit Time parameter ---------
        u32 MinErrorBR = UINT32_MAX;

        // Select the worst BRP value. Here all value from max to min will be tested to get the best tuple of NBRP and DBRP, identical TQ in both phases prevents quantization errors during bit rate switching
        u32 BRP = BRP_MAX;
        while (--BRP >= BRP_MIN) {

            // Calculate the NTQbits according to BRP and the desired Nominal Bitrate
            u32 NTQbits = fsysclk / desiredNominalBitrate / BRP;
            if((NTQbits < TQBIT_MIN) || (NTQbits > TQBIT_MAX)) {
                // This TQbits count is not possible with this BRP, then do the next BRP value
                continue;
            }

            // NTQ & DTQ bits count
            ErrorNTQ = (fsysclk - (desiredNominalBitrate * NTQbits * BRP)); // Calculate NTQ error
            ErrorTQ = ErrorNTQ;

            // If better error then
            if (ErrorTQ <= MinErrorBR) {
                // Save best parameters
                MinErrorBR = ErrorTQ;
                BestBRP = BRP;
                BestNTQbits = NTQbits;
            }

            // NTQ+1 bits count
            if (NTQbits < TQBIT_MAX) {
                // Calculate NTQ error with NTQbits+1
                ErrorNTQ = ((desiredNominalBitrate * (NTQbits+1) * BRP) - fsysclk);
                ErrorTQ = ErrorNTQ;

                // If better error then
                if (ErrorTQ <= MinErrorBR) {
                    // Save best parameters
                    MinErrorBR = ErrorTQ;
                    BestBRP = BRP;
                    BestNTQbits = NTQbits+1;
                }
            }
        }

        if (MinErrorBR == UINT32_MAX) {
            return Result();          // Impossible to find a good BRP
        }

        //--- Calculate Nominal segments --------------------------
        {
            const u32 brp = BestBRP;                                            // ** Save the best NBRP in the configuration **
            u32 NTSEG2 = BestNTQbits / 5;                                       // The Nominal Sample Point must be close to 80% (5x20%) of NTQ per bits so NTSEG2 should be 20% of NTQbits
            if ((BestNTQbits % 5) > 2) NTSEG2++;                                // To be as close as possible to 80%
            if (NTSEG2 < TSEG2_MIN) NTSEG2 = TSEG2_MIN;                         // Correct NTSEG2 if < 1
            if (NTSEG2 > TSEG2_MAX) NTSEG2 = TSEG2_MAX;                         // Correct NTSEG2 if > 128
            const u32 tseg2 = NTSEG2;                                           // ** Save the NTSEG2 in the configuration **
            u32 NTSEG1 = BestNTQbits - NTSEG2 - SyncSegment;                    // NTSEG1  = NTQbits - NTSEG2 - 1 (NSYNC)
            if (NTSEG1 < TSEG1_MIN) NTSEG1 = TSEG1_MIN;                         // Correct NTSEG1 if < 1
            if (NTSEG1 > TSEG1_MAX) NTSEG1 = TSEG1_MAX;                         // Correct NTSEG1 if > 256
            const u32 tseg1 = NTSEG1;                                           // ** Save the NTSEG1 in the configuration **
            u32 NSJW = NTSEG2;                                                  // Normally NSJW = NTSEG2, maximizing NSJW lessens the requirement for the oscillator tolerance
            if (NTSEG1 < NTSEG2) NSJW = NTSEG1;                                 // But NSJW = min(NPHSEG1, NPHSEG2)
            if (NSJW < SJW_MIN) NSJW = SJW_MIN;                                 // Correct NSJW if < 1
            if (NSJW > SJW_MAX) NSJW = SJW_MAX;                                 // Correct NSJW if > 128
            const u32 sjw = NSJW;                                               // ** Save the NSJW in the configuration **

            const u32 time_quanta_per_bit_time = SyncSegment + tseg1 + tseg2;
            const u32 bitrate = fsysclk / (brp * time_quanta_per_bit_time);
            const float sample_point = (static_cast<float>(SyncSegment + tseg1) / static_cast<float>(time_quanta_per_bit_time)) * 100.0f;

            // Check for deviation in baud rate tolerance.
            const float ftq = static_cast<float>(fsysclk / brp);                                               // Calculate frequency of Time Quanta (TQ) clock.
            const float tfq_real = static_cast<float>(time_quanta_per_bit_time * desiredNominalBitrate);       // Actual TQ frequency.
            const float freq_diff = ((ftq - tfq_real) * 2.0f * 100.0f) / (ftq + tfq_real);

            const Result result = {bitrate, brp, time_quanta_per_bit_time, tseg1, tseg2, sjw, sample_point, freq_diff, true};
            return result;
        }
    }
}


CanSolver::Result CanSolver::calculate4_helper(const Input &in, double freq, u8 SyncSegment, u32 BRP, u32 TQbits)
{
    //--- Calculate segments --------------------------
    const u32 brp = BRP /*- 1*/;                                        // ** Save the best NBRP in the configuration **
    u32 NTSEG2 = TQbits / 5;                                            // The Nominal Sample Point must be close to 80% (5x20%) of NTQ per bits so NTSEG2 should be 20% of NTQbits
    if ((TQbits % 5) > 2) NTSEG2++;                                     // To be as close as possible to 80%
    if (NTSEG2 < in.time_seg2_min) NTSEG2 = in.time_seg2_min;           // Correct NTSEG2 if < 1
    if (NTSEG2 > in.time_seg2_max) NTSEG2 = in.time_seg2_max;           // Correct NTSEG2 if > 128
    const u32 tseg2 = NTSEG2/* - 1*/;                                   // ** Save the NTSEG2 in the configuration **
    u32 NTSEG1 = TQbits - NTSEG2 - SyncSegment;                         // NTSEG1  = NTQbits - NTSEG2 - 1 (NSYNC)
    if (NTSEG1 < in.time_seg1_min) NTSEG1 = in.time_seg1_min;           // Correct NTSEG1 if < 1
    if (NTSEG1 > in.time_seg1_max) NTSEG1 = in.time_seg1_max;           // Correct NTSEG1 if > 256
    const u32 tseg1 = NTSEG1/* - 1*/;                                   // ** Save the NTSEG1 in the configuration **
    u32 NSJW = NTSEG2;                                                  // Normally NSJW = NTSEG2, maximizing NSJW lessens the requirement for the oscillator tolerance
    if (NTSEG1 < NTSEG2) NSJW = NTSEG1;                                 // But NSJW = min(NPHSEG1, NPHSEG2)
    if (NSJW < in.sjw_min) NSJW = in.sjw_min;                           // Correct NSJW if < 1
    if (NSJW > in.sjw_max) NSJW = in.sjw_max;                           // Correct NSJW if > 128
    const u32 sjw = NSJW/* - 1*/;                                       // ** Save the NSJW in the configuration **

    const u32 time_quanta_per_bit_time = SyncSegment + tseg1 + tseg2;
    const u32 bitrate = static_cast<u32>(freq) / (brp * time_quanta_per_bit_time);
    const double sample_point = (static_cast<double>(SyncSegment + tseg1) / time_quanta_per_bit_time) * 100.0;

    // Check for deviation in baud rate tolerance.
    const double ftq = freq / brp;                                                      // Calculate frequency of Time Quanta (TQ) clock.
    const double tfq_real = time_quanta_per_bit_time * in.rate;                         // Actual TQ frequency.
    const double freq_diff = ((ftq - tfq_real) * 2.0 * 100.0) / (ftq + tfq_real);

    const Result result = {bitrate, brp, time_quanta_per_bit_time, tseg1, tseg2, sjw, static_cast<float>(sample_point), static_cast<float>(freq_diff)};
    m_results.push_back(result);
    return result;
}




void CanSolver::sort_prescalar()
{
    if (m_results.empty()) {
        return; // No results to sort.
    }

    // Sort the results based on prescaler and sample point.
    std::sort(m_results.begin(), m_results.end(), [](const Result& a, const Result& b) {
        if (a.prescaler != b.prescaler) {
            return a.prescaler < b.prescaler; // Prefer smaller prescaler values.
        }
        return std::abs(a.sample_point) < std::abs(b.sample_point); // Prefer smaller sample point deviation.
    });
}


void CanSolver::sort_percent()
{
    if (m_results.empty()) {
        return; // No results to sort.
    }

    // Sort by sample point deviation.
    std::sort(m_results.begin(), m_results.end(), [](const Result& a, const Result& b) {
        return std::abs(a.sample_point) < std::abs(b.sample_point);
    });
}


void CanSolver::sort_diff()
{
    if (m_results.empty()) {
        return; // No results to sort.
    }

    // Sort by frequency deviation.
    std::sort(m_results.begin(), m_results.end(), [](const Result& a, const Result& b) {
        return std::abs(a.freq_diff) < std::abs(b.freq_diff);
    });
}

void CanSolver::sort()
{
    if (m_results.empty()) {
        return; // No results to sort.
    }

    // Sort
    std::sort(m_results.begin(), m_results.end(), [](const Result& a, const Result& b) {
        if (std::abs(a.sample_point - b.sample_point) < 0.001) {
            return a.prescaler > b.prescaler;
        }
        return a.sample_point < b.sample_point;
    });
}




namespace std {
template <>
struct hash<CanSolver::Result> {
    std::size_t operator()(const CanSolver::Result& res) const noexcept {
        std::size_t h1 = std::hash<u32>{}(res.prescaler);
        std::size_t h2 = std::hash<u32>{}(res.tq_per_bit);
        std::size_t h3 = std::hash<u32>{}(res.tBS1);
        std::size_t h4 = std::hash<u32>{}(res.tBS2);
        std::size_t h5 = std::hash<u32>{}(res.sync_jump_width);
        std::size_t h6 = std::hash<float>{}(res.sample_point);
        std::size_t h7 = std::hash<float>{}(res.freq_diff);

        // Комбінуємо хеші:
        return (((((h1 ^ (h2 << 1)) ^ (h3 << 2)) ^ (h4 << 3)) ^ (h5 << 4)) ^ (h6 << 5) ^ (h7 << 6));
    }
};
}

void CanSolver::remove_duplicates()
{
    std::unordered_set<Result> seen; 		// Track unique results.
    std::vector<Result> unique_results; 	// Store unique results.

    for (const auto& result : m_results) {
        if (seen.insert(result).second) { 		// Insert and check if it's new.
            unique_results.push_back(result); 	// Add to unique results if new.
        }
    }

    m_results = std::move(unique_results); // Replace original results with unique results.
}

#if CANSOLVER_PRINT

void CanSolver::print_results(int max_rows)
{
    std::cout << "Rate[Mb]| Prescaler | TQ/bit | TimeSeg1 | TimeSeg2 | SJW  | SampleP[%] | Diff[%]\n";
    std::cout << "-------------------------------------------------------------------------------\n";

    int rows = 0;
    for (const auto& res : m_results) {
        if (rows++ >= max_rows) break; // Stop if the row limit is reached.

        std::cout << std::setw(5) << static_cast<float>(res.bitrate)/1000000.0
                  << " | " << std::setw(9) << res.prescaler
                  << " | " << std::setw(7) << res.tq_per_bit
                  << " | " << std::setw(8) << res.tBS1
                  << " | " << std::setw(8) << res.tBS2
                  << " | " << std::setw(5) << res.sync_jump_width
                  << " | " << std::setw(10) << res.sample_point
                  << " | " << std::setw(6) << res.freq_diff << "\n";
    }

    std::cout << std::endl << "------" << std::endl;

    // Print the best result.
    std::cout << "The best:\n";
    std::cout << std::setw(5) <<static_cast<float>(m_best.bitrate)/1000000
              << " | " << std::setw(9) << m_best.prescaler
              << " | " << std::setw(7) << m_best.tq_per_bit
              << " | " << std::setw(8) << m_best.tBS1
              << " | " << std::setw(8) << m_best.tBS2
              << " | " << std::setw(5) << m_best.sync_jump_width
              << " | " << std::setw(10) << m_best.sample_point
              << " | " << std::setw(6) << m_best.freq_diff << "\n";
    std::cout << "\n";
}

void CanSolver::print_results(Result result)
{
    m_best = result;
    print_results();
}

#endif /* CANSOLVER_PRINT */


#if CANSOLVER_TEST

void CanSolverTest()
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

#endif /* CANSOLVER_TEST */




