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
 * A lookup outside the map returns nullptr. Check the pointer before
 * reading or writing the strip:
 *
 *     if (xStrip != nullptr) {
 *         xStrip->pulse[0] = 100.;
 *         Double_t center = xStrip->position;
 *     }
 *
 * GEMStrip parameters:
 *
 *     detectorId            GEM detector ID
 *     direction             X or Y strip direction
 *     stripNumber           zero-based strip number
 *     position              strip center position in mm
 *     activeRegionHalfWidth active region is position +/- this value
 *     pulse                 ADC values of the six time samples
 *     peakADC               peak ADC value
 *     meanADC               mean ADC value of the six samples
 *     peakSample            zero-based sample index of peakADC
 *
 * This header organizes the strip data and lookup map only. It does not
 * calculate peakADC, meanADC, or peakSample.
 */

#ifndef GEM_STRIP_MAP_H
#define GEM_STRIP_MAP_H

#include <array>
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
        return static_cast<Int_t>((position - low) / stripWidth);
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
};

#endif
