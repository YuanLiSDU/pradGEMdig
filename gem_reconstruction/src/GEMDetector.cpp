//============================================================================//
// GEM detector class                                                         //
// A detector has several planes (currently they are X, Y)                    //
// The planes are managed by detector, thus copy a detector will also copy    //
// the planes that consist of it                                              //
//                                                                            //
// Chao Peng     10/07/2016                                                  //
// Xinzhan Bai   12/02/2020                                                  //
//============================================================================//


#include "GEMDetector.h"
#include "GEMDetectorLayer.h"
#include "GEMCluster.h"

////////////////////////////////////////////////////////////////////////////////
// constructor

GEMDetector::GEMDetector(GEMDetectorLayer *layer,
                         const std::string &readoutBoard,
                         const std::string &detectorType,
                         const std::string &detector,
                         const int &detectorID,
                         const int &layerID,
                         const int &position_index)
: gem_layer(layer), det_name(detector), det_id(detectorID), layer_id(layerID), 
    layer_position_index(position_index),
    type(detectorType), readout_board(readoutBoard), res(0.07)
{
    planes.resize(GEMPlane::Max_Types, nullptr);

    gem_hits.reserve(GEM_CLUSTERS_BUFFER);
}

////////////////////////////////////////////////////////////////////////////////
// default constructor

GEMDetector::GEMDetector()
: gem_layer(nullptr), det_name("Undefined_GEM"), det_id(-1), layer_id(-1), 
    layer_position_index(-1),
    type("PRADSIMGEM"), readout_board("CARTESIAN"), res(0.07)
{
    planes.resize(GEMPlane::Max_Types, nullptr);

    gem_hits.reserve(GEM_CLUSTERS_BUFFER);
}

////////////////////////////////////////////////////////////////////////////////
// a convenient constructor

GEMDetector::GEMDetector(GEMDetectorLayer *layer,
                         const std::string &detectorName,
                         const int &detectorID,
                         const int &layerID,
                         const int &position_index)
: gem_layer(layer), det_name(detectorName), det_id(detectorID), layer_id(layerID), 
    layer_position_index(position_index),
    type("PRADSIMGEM"), readout_board("CARTESIAN"), res(0.07)
{
    planes.resize(GEMPlane::Max_Types, nullptr);

    gem_hits.reserve(GEM_CLUSTERS_BUFFER);
}

////////////////////////////////////////////////////////////////////////////////
// copy and move assignment will copy or move the planes because planes are 
// managed by detectors, but it won't copy the connection to gem system and so
// won't the id assigned by gem system
// copy constructor

GEMDetector::GEMDetector(const GEMDetector &that)
: det_name(that.det_name), det_id(that.det_id), layer_id(that.layer_id),
    layer_position_index(that.layer_position_index), type(that.type),
  readout_board(that.readout_board), gem_hits(that.gem_hits), res(that.res)
{
    for(auto &plane : that.planes)
    {
        if(plane != nullptr)
            planes.push_back(new GEMPlane(*plane));
        else
            planes.push_back(nullptr);
    }

    ConnectPlanes();
}

////////////////////////////////////////////////////////////////////////////////
// move constructor

GEMDetector::GEMDetector(GEMDetector &&that)
: det_name(std::move(that.det_name)), det_id(std::move(that.det_id)),
    layer_id(std::move(that.layer_id)),
    layer_position_index(std::move(that.layer_position_index)), type(std::move(that.type)),
  readout_board(std::move(that.readout_board)), planes(std::move(that.planes)),
  gem_hits(std::move(that.gem_hits)), res(that.res)
{
    // reset the planes' detector
    ConnectPlanes();
}

////////////////////////////////////////////////////////////////////////////////
// destructor

GEMDetector::~GEMDetector()
{
    for(auto &plane : planes)
    {
        if(plane)
            plane->UnsetDetector(true);

        delete plane;
    }
}


////////////////////////////////////////////////////////////////////////////////
// copy assignment operator

GEMDetector &GEMDetector::operator =(const GEMDetector &rhs)
{
    if(this == &rhs)
        return *this;

    GEMDetector that(rhs); // use copy constructor
    *this = std::move(that); // use move assignment
    return *this;
}

////////////////////////////////////////////////////////////////////////////////
// move assignment operator

GEMDetector &GEMDetector::operator =(GEMDetector &&rhs)
{
    if(this == &rhs)
        return *this;

    det_name = std::move(rhs.det_name);
    det_id = std::move(rhs.det_id);
    layer_id = std::move(rhs.layer_id);
    layer_position_index = std::move(rhs.layer_position_index);
    type = std::move(rhs.type);
    readout_board = std::move(rhs.readout_board);
    planes = std::move(rhs.planes);
    gem_hits = std::move(rhs.gem_hits);
    res = rhs.res;

    ConnectPlanes();

    return *this;
}

//============================================================================//
// Public Member Functions                                                    //
//============================================================================//

////////////////////////////////////////////////////////////////////////////////
// connect detector to gem layer

void GEMDetector::SetGEMLayer(GEMDetectorLayer *l)
{
    if(l == gem_layer)
        return;

    gem_layer = l;
}


////////////////////////////////////////////////////////////////////////////////
// add plane to the detector

bool GEMDetector::AddPlane(const int &plane_type,
                           const std::string &name,
                           const double &size,
                           const int &dir)
{
    GEMPlane *new_plane = new GEMPlane(name, plane_type, size, dir);

    // successfully added plane in
    if(AddPlane(new_plane))
        return true;

    delete new_plane;
    return false;
}

////////////////////////////////////////////////////////////////////////////////
// add plane to the detector

bool GEMDetector::AddPlane(GEMPlane *plane)
{
    if(plane == nullptr) {
        std::cerr<<"Error: adding nullptr plane pointers to detectors."<<std::endl;
        return false;
    }

    int idx = (int)plane->GetType();

    if(plane->GetDetector() != nullptr) {
        std::cerr << " GEM Detector Error: "
                  << "Trying to add plane " << plane->GetName()
                  << " to detector " << det_name
                  << ", but the plane is belong to " << plane->GetDetector()->GetName()
                  << ". Action aborted."
                  << std::endl;
        return false;
    }

    if(planes[idx] != nullptr) {
        // already added, do nothing
        // Do not use this error message, because this AddPlane() operation will be called
        //      for each APV, we are now using a different mapping file mechanism
        //std::cerr << " GEM Detector Error: "
        //          << "Trying to add multiple planes with the same plane type " << idx
        //          << "to detector " << det_name << ". Action aborted."
        //          << std::endl;
        return false;
    }

    plane->SetDetector(this);
    planes[idx] = plane;
    return true;
}

////////////////////////////////////////////////////////////////////////////////
// remove plane

void GEMDetector::RemovePlane(const int &plane_type)
{
    if((uint32_t)plane_type >= planes.size())
        return;

    auto &plane = planes[plane_type];

    if(plane) {
        plane->UnsetDetector(true);
        delete plane, plane = nullptr;
    }
}

////////////////////////////////////////////////////////////////////////////////
// remove plane

void GEMDetector::DisconnectPlane(const int &plane_type, bool force_disconn)
{
    if((uint32_t)plane_type >= planes.size())
        return;

    auto &plane = planes[plane_type];

    if(!plane)
        return;

    if(!force_disconn)
        plane->UnsetDetector(true);

    plane = nullptr;
}

////////////////////////////////////////////////////////////////////////////////
// Asure the connection to all the planes

void GEMDetector::ConnectPlanes()
{
    for(auto &plane : planes)
    {
        if(plane != nullptr)
            plane->SetDetector(this);
    }
}

////////////////////////////////////////////////////////////////////////////////
// set detector offset
void GEMDetector::SetOffset(double x_o, double y_o)
{ 
    if(planes.size() != 2) {
        std::cout<<"WARNING: Invalid Operation. You must add 2 planes to detector before setting offsets."<<std::endl
                 <<"         Default offset values: (x_offset=0, y_offset=0) will be used."<<std::endl;
        return;
    }

    x_offset = x_o;
    y_offset = y_o;

    planes[0] -> SetOffset(x_offset);
    planes[1] -> SetOffset(y_offset);
}

////////////////////////////////////////////////////////////////////////////////
// reconstruct hits on planes, need GEMCluster as an input

void GEMDetector::Reconstruct(GEMCluster *gem_recon)
{
    // group strip hits into clusters
    for(auto &plane : planes)
    {
        if(plane == nullptr)
            continue;
        plane->FormClusters(gem_recon);
    }

    // Cartesian reconstruction method
    // reconstruct event hits from clusters
    GEMPlane *plane_x = GetPlane(GEMPlane::Plane_X);
    GEMPlane *plane_y = GetPlane(GEMPlane::Plane_Y);
    // do not have these two planes, cannot reconstruct
    if(!plane_x || !plane_y)
        return;
    gem_recon->CartesianReconstruct(plane_x->GetStripClusters(),
                                    plane_y->GetStripClusters(),
                                    gem_hits,
                                    det_id,
                                    res);
}

////////////////////////////////////////////////////////////////////////////////
// clear all the hits on plane

void GEMDetector::ClearHits()
{
    for(auto &plane : planes)
    {
        if(plane != nullptr)
            plane->ClearStripHits();
    }
}

////////////////////////////////////////////////////////////////////////////////
// clear all the hits and clusters

void GEMDetector::Reset()
{
    ClearHits();
    gem_hits.clear();
}

////////////////////////////////////////////////////////////////////////////////
// get plane

GEMPlane *GEMDetector::GetPlane(const std::string &plane_type)
const
{
    uint32_t idx = (uint32_t) GEMPlane::str2Type(plane_type.c_str());

    if(idx >= planes.size())
        return nullptr;

    return planes.at(idx);
}

////////////////////////////////////////////////////////////////////////////////
// get the plane list

std::vector<GEMPlane*> GEMDetector::GetPlaneList()
const
{
    // since it allows nullptr in planes
    // for safety issue, only pack existing planes and return
    std::vector<GEMPlane*> result;

    for(auto &plane : planes)
    {
        if(plane != nullptr)
            result.push_back(plane);
    }

    return result;
}


////////////////////////////////////////////////////////////////////////////////
// get plane by type

GEMPlane *GEMDetector::GetPlane(const int &plane_type)
const
{
    return planes[plane_type];
}


////////////////////////////////////////////////////////////////////////////////
// get apv lists from a plane

void GEMDetector::PrintStatus()
{
    std::cout<<"detector id = "<<GetDetID()<<" contains "
        <<planes.size()<<" planes."<<std::endl;

    int nPlane = 0;
    for(auto &i: planes) {
        if(i != nullptr) {
            i -> PrintStatus();
            nPlane++;
        }
        else
            std::cout<<" with undefined plane "<<nPlane++<<std::endl;
    }
}
