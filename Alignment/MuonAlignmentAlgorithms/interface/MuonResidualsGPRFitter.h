// #ifndef Alignment_MuonAlignmentAlgorithms_MuonResidualsGPRFitter_H
// #define Alignment_MuonAlignmentAlgorithms_MuonResidualsGPRFitter_H
#ifndef MUON_RESIDUALS_GPR_FITTER_H
#define MUON_RESIDUALS_GPR_FITTER_H

#ifndef STANDALONE_FITTER
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "Alignment/CommonAlignment/interface/Alignable.h"
#endif

#include "Alignment/MuonAlignmentAlgorithms/interface/MuonResidualsTwoBin.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"

#include "TMinuit.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TMatrixDSym.h"

#include <cstdio>
#include <iostream>
#include <string>
#include <sstream>
#include <map>

// #ifdef STANDALONE_FITTER
// #include "MuonResidualsFitter.h"
// #else
// #include "Alignment/MuonAlignmentAlgorithms/interface/MuonResidualsFitter.h"
// #endif

class MuonResidualsGPRFitter {
public:
// enum for station 1,2,3 data
    enum class Data_6DOF {
        kResidX = 0,
        kResidY = 1,
        kResSlopeX = 2,
        kResSlopeY = 3,
        kPositionX = 4,
        kPositionY = 5,
        kAngleX = 6,
        kAngleY = 7,
        kRedChi2 = 8,
        kPz = 9,
        kPt = 10,
        kCharge = 11,
        kStation = 12,
        kWheel = 13,
        kSector = 14,
        kChambW = 15,
        kChambl = 16,
        kNData = 17
    };

    // enum for station 4 data
    enum class Data_5DOF {
        kResid = 0,
        kResSlope = 1,
        kPositionX = 2,
        kPositionY = 3,
        kAngleX = 4,
        kAngleY = 5,
        kRedChi2 = 6,
        kPz = 7,
        kPt = 8,
        kCharge = 9,
        kStation = 10,
        kWheel = 11,
        kSector = 12,
        kChambW = 13,
        kChambl = 14,
        kNData = 15
    };

    enum class PARAMS {
        kAlignX,
        kAlignY,
        kAlignZ,
        kAlignPhiX,
        kAlignPhiY,
        kAlignPhiZ,
        kResidXSigma,
        kResidYSigma,
        kResSlopeXSigma,
        kResSlopeYSigma,
        kCount // needed to count number of minuit parameters
    };

    explicit MuonResidualsGPRFitter(const DTGeometry* dt_Geometry);

    ~MuonResidualsGPRFitter() {}

    void setPrintLevel(int printLevel) { m_printLevel = printLevel; }
    void setStrategy(int strategy) { m_strategy = strategy; }

    // TMatrixDSym covarianceMatrix() const { return m_cov; }
    double loglikelihood() const { return m_loglikelihood; }

    //returns number of parameters to be fitted
    int npar() const { return static_cast<int>(PARAMS::kCount); }

    // method selecting a subset of a map

    // methods returning all residuals
    std::map<Alignable*, MuonResidualsTwoBin*>::const_iterator datamap_begin() const { return m_datamap.begin(); }
    std::map<Alignable*, MuonResidualsTwoBin*>::const_iterator datamap_end() const { return m_datamap.end(); }

    // method returning residuals of given alignable
    // better to change at to find
    std::vector<double *>::const_iterator selected_chamber_residualsPos_begin(Alignable* ali) const { return m_datamap.at(ali)->residualsPos_begin(); }
    std::vector<double *>::const_iterator selected_chamber_residualsPos_end(Alignable* ali) const { return m_datamap.at(ali)->residualsPos_end(); }
    std::vector<double *>::const_iterator selected_chamber_residualsNeg_begin(Alignable* ali) const { return m_datamap.at(ali)->residualsNeg_begin(); }
    std::vector<double *>::const_iterator selected_chamber_residualsNeg_end(Alignable* ali) const { return m_datamap.at(ali)->residualsNeg_begin(); }

    // method filling pairs (or const_iterator) to m_datamap
    void fill(std::map<Alignable*, MuonResidualsTwoBin*>::const_iterator ali_and_data);

    //returns number of all residuals
    int getSize() const { return m_datamap.size(); }

    // dt geometry getter
    const DTGeometry* getDTGeometry() const { return m_gpr_dtGeometry; }

    // function which is called to do a fit on a set of alignables
    /// implement version of this function to be able to fit a subset of DT system
    // wrapper-function only preparing stuff for dofit
    bool fit();

private:
    // map store all pairs alignable chamber - TwoBin with residuals for this chamber
    std::map<Alignable*, MuonResidualsTwoBin*> m_datamap;

    // pointer to DT geometry to access methods for coordinate conversion in FCN
    const DTGeometry* m_gpr_dtGeometry;

    int m_printLevel;
    int m_strategy;
    // bool m_weightAlignment;

    // ? for MINUIT's output ?
    std::vector<double> m_value;
    std::vector<double> m_error;
    // TMatrixDSym m_cov;
    double m_loglikelihood;

    void inform(TMinuit *tMinuit); // add this later

    // function actually calculating the shift; to be called from fit function above
    // I don not understand usage of parNum and parName yet
    // takes a function pointer to function FCN
    bool dofit(void (*fcn)(int &, double *, double &, double *, int),
               std::vector<int> &parNum,
               std::vector<std::string> &parName,
               std::vector<double> &start,
               std::vector<double> &step,
               std::vector<double> &low,
               std::vector<double> &high);

};


// Auxilliary class to get information into the fit function; Idk what its doing, copied from MuonResidualsfitter.h
class MuonResidualsGPRFitterFitInfo : public TObject {
public:
    MuonResidualsGPRFitterFitInfo(MuonResidualsGPRFitter *gpr_fitter) : m_gpr_fitter(gpr_fitter) {}
    MuonResidualsGPRFitter *gpr_fitter() { return m_gpr_fitter; } // is needed to get access to gpr fitter object in likelihod calc

private:
    MuonResidualsGPRFitter *m_gpr_fitter;
#ifdef STANDALONE_FITTER
    ClassDef(MuonResidualsGPRFitterFitInfo, 1);
#endif
};

#ifdef STANDALONE_FITTER
    ClassImp(MuonResidualsGPRFitterFitInfo);
#endif

double MuonResidualsGPRFitter_logPureGaussian(double residual, double center, double sigma);

// #endif  // Alignment_MuonAlignmentAlgorithms_MuonResidualsGPRFitter_H
#endif //MUON_RESIDUALS_GPR_FITTER_H
