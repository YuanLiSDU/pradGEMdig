#ifndef GEM_PLANE_H
#define GEM_PLANE_H

#include <cstdint>
#include "ConfigParser.h"
#include "GEMStruct.h"

class GEMDetector;
class GEMCluster;

class GEMPlane
{
public:
    enum Type
    {
        Undefined_Type = -1,
        Plane_X = 0,
        Plane_Y,
        Max_Types,
    };
    // macro in ConfigParser.h
    ENUM_MAP(Type, 0, "X|Y");

public:
    // constructors
    GEMPlane(GEMDetector *det = nullptr);
    GEMPlane(const std::string &n, const int &t, const float &s,
            const int &d, const int &o = 0, GEMDetector *det = nullptr);

    // copy/move constructors
    GEMPlane(const GEMPlane &that);
    GEMPlane(GEMPlane &&that);

    // destructor
    virtual ~GEMPlane();

    // copy/move assignment operators
    GEMPlane &operator= (const GEMPlane &rhs);
    GEMPlane &operator= (GEMPlane &&rhs);

    // public member functions
    void AddStripHit(int strip, float charge, short timebin, const std::vector<float> &ts_adc, bool xtalk = false, int crate = 0, int mpd = 0, int adc = 0);
    void ClearStripHits();
    float GetStripPosition(const int &plane_strip) const;
    void FormClusters(GEMCluster *method);

    // set parameter
    void SetDetector(GEMDetector *det, bool force_set = false);
    void UnsetDetector(bool force_unset = false);
    void SetName(const std::string &n) {name = n;}
    void SetType(const Type &t) {type = t;}
    void SetSize(const float &s) {size = s;}
    void SetOrientation(const int &o) {orient = o;}
    void SetCapacity(int c);
    void SetOffset(double o){offset = o;}

    // get parameter
    GEMDetector *GetDetector() const {return detector;}
    const std::string &GetName() const {return name;}
    Type GetType() const {return type;}
    float GetSize() const {return size;}
    int GetOrientation() const {return orient;}
    double GetOfffset() const {return offset;}
    std::vector<StripHit> &GetStripHits() {return strip_hits;}
    const std::vector<StripHit> &GetStripHits() const {return strip_hits;}
    std::vector<StripCluster> &GetStripClusters() {return strip_clusters;}
    const std::vector<StripCluster> &GetStripClusters() const {return strip_clusters;};
    void PrintStatus();

private:
    GEMDetector *detector;
    std::string name;
    Type type;
    float size;
    int orient;
    int direction;

    // plane raw hits and clusters
    std::vector<StripHit> strip_hits;
    std::vector<StripCluster> strip_clusters;

    double offset = 0.;
};

#endif
