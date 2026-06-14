/*
 * GEMStripMap usage
 *
 * Initialize one strip map for each GEM detector:
 *
 *     GEMStripMap stripMap(0); // detectorId = 0
 *
 * Find a strip by its zero-based strip number:
 *
 *     GEMStrip* xStrip = stripMap.GetXStrip(10);
 *     GEMStrip* yStrip = stripMap.GetYStrip(20);
 *     GEMStrip* strip =
 *         stripMap.GetStrip(GEMStripDirection::X, 10);
 *
 * Find a strip by its X or Y position in mm:
 *
 *     GEMStrip* xStrip = stripMap.FindXStrip(12.4);
 *     GEMStrip* yStrip = stripMap.FindYStrip(-30.0);
 *     GEMStrip* strip =
 *         stripMap.FindStrip(GEMStripDirection::Y, -30.0);
 *
 * Position lookup returns the strip whose interval contains the given
 * coordinate. The valid range is [map_low, map_high); a position outside
 * the map returns nullptr:
 *
 *     if (xStrip != nullptr) {
 *         xStrip->pulse[0] = 100.;
 *         Double_t center = xStrip->position;
 *     }
 *
 * Find the X/Y strip region around a position. Num_half is the number of
 * strips included on each side of the center strip. For example,
 * Num_half = 6 normally gives 13 strips in each direction. The region is
 * clipped at the map boundaries:
 *
 *     GEMStripMap::GEMStripRegion region =
 *         stripMap.findRegion(stripMap, x, y, 6);
 *
 *     if (region.IsValid()) {
 *         Int_t firstXStrip = region.lowX_strip;
 *         Int_t lastXStrip = region.upX_strip;
 *         Double_t leftBoundary = region.xLow;
 *         Double_t rightBoundary = region.xHigh;
 *     }
 *
 * GEMStripRegion parameters:
 *
 *     lowX_strip, upX_strip first/last X strip number, inclusive
 *     lowY_strip, upY_strip first/last Y strip number, inclusive
 *     xLow, xHigh             physical X boundaries in mm
 *     yLow, yHigh             physical Y boundaries in mm
 *
 * Clear all event signal information while keeping the strip geometry:
 *
 *     stripMap.Clear();
 *
 * GEMStrip parameters:
 *
 *     detectorId            GEM detector ID
 *     direction             X or Y strip direction
 *     stripNumber           zero-based strip number
 *     position              strip center position in mm
 *     activeRegionHalfWidth active region is position +/- this value
 *     charge                total integrated charge on the strip
 *     pulse                 ADC values of the six time samples
 *     peakADC               peak ADC value
 *     meanADC               mean ADC value of the six samples
 *     peakSample            zero-based sample index of peakADC
 *
 * Clear() resets charge, pulse, peakADC, meanADC, and peakSample for all
 * strips. This header organizes strip geometry, event data, and lookup
 * regions; it does not calculate those signal quantities.
 */

#ifndef GEM_STRIP_MAP_H
#define GEM_STRIP_MAP_H

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include "Rtypes.h"

#include "GEMdigVariables.h"

enum class GEMStripDirection {
    X,
    Y
};

struct GEMStrip {
    Int_t detectorId = -1;
    GEMStripDirection direction = GEMStripDirection::X;
    Int_t stripNumber = -1;
    Double_t position = 0.; // Center position in the X or Y direction, mm
    Double_t activeRegionHalfWidth = stripWidth / 2.; // Active region: position +/- this value, mm

    Double_t charge = 0.; // Total charge on the strip
    std::array<Double_t, samplingNum> pulse = {};
    Double_t peakADC = 0.;
    Double_t meanADC = 0.;
    Int_t peakSample = -1; // Zero-based index in pulse
};

class GEMStripMap {
public:
    explicit GEMStripMap(Int_t detector = -1)
        : detectorId(detector)
    {
        BuildStrips();
    }

    GEMStrip* GetStrip(GEMStripDirection direction, Int_t stripNumber)
    {
        std::vector<GEMStrip>& strips = SelectStrips(direction);
        if (stripNumber < 0 ||
            static_cast<std::size_t>(stripNumber) >= strips.size()) {
            return nullptr;
        }
        return &strips[stripNumber];
    }

    const GEMStrip* GetStrip(GEMStripDirection direction,
                             Int_t stripNumber) const
    {
        const std::vector<GEMStrip>& strips = SelectStrips(direction);
        if (stripNumber < 0 ||
            static_cast<std::size_t>(stripNumber) >= strips.size()) {
            return nullptr;
        }
        return &strips[stripNumber];
    }

    GEMStrip* FindStrip(GEMStripDirection direction, Double_t position)
    {
        return GetStrip(direction, PositionToStripNumber(direction, position));
    }

    GEMStrip* GetXStrip(Int_t stripNumber)
    {
        return GetStrip(GEMStripDirection::X, stripNumber);
    }

    GEMStrip* GetYStrip(Int_t stripNumber)
    {
        return GetStrip(GEMStripDirection::Y, stripNumber);
    }

    GEMStrip* FindXStrip(Double_t xPosition)
    {
        return FindStrip(GEMStripDirection::X, xPosition);
    }

    GEMStrip* FindYStrip(Double_t yPosition)
    {
        return FindStrip(GEMStripDirection::Y, yPosition);
    }

    const std::vector<GEMStrip>& GetXStrips() const
    {
        return xStrips;
    }

    const std::vector<GEMStrip>& GetYStrips() const
    {
        return yStrips;
    }

    void Clear()
    {
        ClearStrips(xStrips);
        ClearStrips(yStrips);
    }

    struct GEMStripRegion {
        Int_t lowX_strip = -1;
        Int_t lowY_strip = -1;
        Int_t upX_strip = -1;
        Int_t upY_strip = -1;

        Double_t xLow = 0.;
        Double_t xHigh = 0.;
        Double_t yLow = 0.;
        Double_t yHigh = 0.;

        bool IsValid() const
        {
            return lowX_strip >= 0 && lowY_strip >= 0 && upX_strip >= lowX_strip && upY_strip >= lowY_strip &&
                xHigh > xLow && yHigh > yLow;
        }
    };

    inline GEMStripRegion findRegion(
        GEMStripMap& gemStripMap, Double_t x, Double_t y, Int_t Num_half)
    {
        if(Num_half < 0) return {};

        const GEMStrip* centerXStrip = gemStripMap.FindXStrip(x);
        const GEMStrip* centerYStrip = gemStripMap.FindYStrip(y);
        if(centerXStrip == nullptr || centerYStrip == nullptr) return {};

        GEMStripRegion region;
        region.lowX_strip = centerXStrip->stripNumber - Num_half;
        region.lowY_strip = centerYStrip->stripNumber - Num_half;
        region.upX_strip = centerXStrip->stripNumber + Num_half;
        region.upY_strip = centerYStrip->stripNumber + Num_half;

        if(region.lowX_strip < 0) region.lowX_strip = 0;
        if(region.lowY_strip < 0) region.lowY_strip = 0;
        if(region.upX_strip >= map_x_bins) region.upX_strip = map_x_bins - 1;
        if(region.upY_strip >= map_y_bins) region.upY_strip = map_y_bins - 1;

        const GEMStrip* leftStrip = gemStripMap.GetXStrip(region.lowX_strip);
        const GEMStrip* rightStrip = gemStripMap.GetXStrip(region.upX_strip);
        const GEMStrip* bottomStrip = gemStripMap.GetYStrip(region.lowY_strip);
        const GEMStrip* topStrip = gemStripMap.GetYStrip(region.upY_strip);
        if(leftStrip == nullptr || rightStrip == nullptr ||
        bottomStrip == nullptr || topStrip == nullptr) {
            return {};
        }

        region.xLow =
            leftStrip->position - stripWidth / 2.;
        region.xHigh =
            rightStrip->position + stripWidth / 2.;
        region.yLow =
            bottomStrip->position - stripWidth / 2.;
        region.yHigh =
            topStrip->position + stripWidth / 2.;
        return region;
    }

private:
    Int_t detectorId;
    std::vector<GEMStrip> xStrips;
    std::vector<GEMStrip> yStrips;

    void BuildStrips()
    {
        xStrips.resize(map_x_bins);
        for (Int_t stripNumber = 0; stripNumber < map_x_bins; ++stripNumber) {
            GEMStrip& strip = xStrips[stripNumber];
            strip.detectorId = detectorId;
            strip.direction = GEMStripDirection::X;
            strip.stripNumber = stripNumber;
            strip.position =
                map_x_low + (static_cast<Double_t>(stripNumber) + 0.5) * stripWidth;
        }

        yStrips.resize(map_y_bins);
        for (Int_t stripNumber = 0; stripNumber < map_y_bins; ++stripNumber) {
            GEMStrip& strip = yStrips[stripNumber];
            strip.detectorId = detectorId;
            strip.direction = GEMStripDirection::Y;
            strip.stripNumber = stripNumber;
            strip.position =
                map_y_low + (static_cast<Double_t>(stripNumber) + 0.5) * stripWidth;
        }
    }

    Int_t PositionToStripNumber(GEMStripDirection direction,
                                Double_t position) const
    {
        const Double_t low =
            direction == GEMStripDirection::X ? map_x_low : map_y_low;
        const Double_t high =
            direction == GEMStripDirection::X ? map_x_high : map_y_high;

        if (position < low || position >= high) {
            return -1;
        }

        return static_cast<Int_t>(
            std::floor((position - low) / stripWidth));
    }

    std::vector<GEMStrip>& SelectStrips(GEMStripDirection direction)
    {
        return direction == GEMStripDirection::X ? xStrips : yStrips;
    }

    const std::vector<GEMStrip>& SelectStrips(
        GEMStripDirection direction) const
    {
        return direction == GEMStripDirection::X ? xStrips : yStrips;
    }

    void ClearStrips(std::vector<GEMStrip>& strips)
    {
        for (GEMStrip& strip : strips) {
            strip.pulse.fill(0.);
            strip.charge = 0.;
            strip.peakADC = 0.;
            strip.meanADC = 0.;
            strip.peakSample = -1;
        }
    }
};

#endif
