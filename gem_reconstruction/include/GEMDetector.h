#ifndef GEM_DETECTOR_H
#define GEM_DETECTOR_H

#include <string>
#include <vector>
#include "GEMStruct.h"
#include "GEMPlane.h"

// reserve space for faster filling of clusters
#define GEM_CLUSTERS_BUFFER 500

class GEMDetectorLayer;
class GEMCluster;

class GEMDetector
{
public:
    friend class GEMCluster;

public:
    // constructor
    GEMDetector(GEMDetectorLayer *layer,
                const std::string &readoutBoard,
                const std::string &detectorType,
                const std::string &detectorName,
                const int &detectorID,
                const int &layerID, // the id of the layer
                const int &pos_index); // the position index of this detector inside the layer

    // default constructor
    GEMDetector();
    // a convenient constructor
    GEMDetector(GEMDetectorLayer *layer,
                const std::string &detectorName,
                const int &detectorID,
                const int &layerID,
                const int &pos_index);

    // copy/move constructors
    GEMDetector(const GEMDetector &that);
    GEMDetector(GEMDetector &&that);

    // desctructor
    virtual ~GEMDetector();

    // copy/move assignment operators
    GEMDetector &operator =(const GEMDetector &rhs);
    GEMDetector &operator =(GEMDetector &&rhs);

    // public member functions
    void SetGEMLayer(GEMDetectorLayer *l);
    void SetResolution(double r) {res = r;}
    void SetID(int i){det_id = i;}
    void setLayerID(int i){layer_id = i;}
    void SetLayerPositionIndex(int i){layer_position_index = i;}
    void SetOffset(double x, double y);
    void UnsetSystem(bool false_unset = false);
    bool AddPlane(GEMPlane *plane);
    bool AddPlane(const int &plane_axis, const std::string &name, const double &size,
                  const int &dir);
    void RemovePlane(const int &plane_axis);
    void DisconnectPlane(const int &plane_axis, bool false_disconn = false);
    void ConnectPlanes();
    void Reconstruct(GEMCluster *c);
    void CollectHits();
    void ClearHits();
    void Reset();

    // get parameters
    GEMDetectorLayer *GetLayer() const {return gem_layer;}
    double GetResolution() const {return res;}
    const std::string &GetType() const {return type;}
    const std::string &GetReadoutBoard() const {return readout_board;}
    GEMPlane *GetPlane(const int &plane_axis) const;
    GEMPlane *GetPlane(const std::string &plane_axis) const;
    std::vector<GEMPlane*> GetPlaneList() const;
    std::vector<GEMHit> &GetHits() {return gem_hits;}
    const std::vector<GEMHit> &GetHits() const {return gem_hits;}
    int GetDetID() const {return det_id;}
    int GetLayerID() const {return layer_id;}
    int GetDetLayerPositionIndex() const {return layer_position_index;}
    double GetXOffset() const {return x_offset;}
    double GetYOffset() const {return y_offset;}
    const std::string &GetName() const {return det_name;}
    void PrintStatus();

private:
    GEMDetectorLayer *gem_layer;
    std::string det_name;
    int det_id = -1;
    int layer_id;
    int layer_position_index;
    std::string type;
    std::string readout_board;
    float res = 0.07; // units in mm
    std::vector<GEMPlane*> planes;
    std::vector<GEMHit> gem_hits;
    double x_offset = 0.;
    double y_offset = 0.;
};

#endif
