// Rotate 1D spherical or 2D axisymmetric coords into a higher dimension

#include <memory>
#include <iostream>

#include "AMReX_VisMF.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_ParmParse.H"
#include "AMReX_StateDescriptor.H"

using namespace amrex;

const std::string CP_VERSION_STR("CheckPointVersion_1.0");

struct StateDataLike
{
    struct TimeInterval
    {
        Real start, stop;
    };
    const StateDescriptor *desc;
    Box domain;
    BoxArray grids;
    TimeInterval new_time;
    TimeInterval old_time;
    std::unique_ptr<MultiFab> new_data;
    std::unique_ptr<MultiFab> old_data;
    Vector< Vector<BCRec> > bc;
};

struct AmrLevelLike
{
    int level;
    Geometry geom;
    BoxArray grids;
    DistributionMapping dmap;
    IntVect crse_ratio;
    IntVect fine_ratio;
    Vector<StateDataLike> state;
    Vector<StateDataLike> new_state;
};

struct AmrLike
{
    int spdim;
    int finest_level;
    Real cumtime;
    Vector<Real> dt_level;
    Vector<int> level_steps;
    Vector<int> level_count;
    Vector<int> n_cycle;
    Vector<Real> dt_min;
    Vector<IntVect> ref_ratio;
    Vector<Geometry> geom;
    Vector< Vector<Real> > dx;
    Vector<AmrLevelLike> amrLevels;
};

class Rotator
{
public:
    
    Rotator();
    
private:
    
    void readCheckpoint(const std::string& filePath);
    
    AmrLike amrlike;
    
    std::string chkfile;
    int nghost = 0;
};

Rotator::Rotator()
{
    ParmParse pp;
    
    if(pp.contains("checkin"))
    {
        pp.get("checkin", chkfile);
    }
    else
    {
        Abort("No checkpoint file supplied!");
    }
    
    readCheckpoint(chkfile);
}

void Rotator::readCheckpoint(const std::string& filePath)
{
    // Preread FabArray headers
    std::map<std::string, Vector<char> > faHeaderMap;
    std::string faHeaderFile(filePath + "/FabArrayHeaders.txt");
    Vector<char> faHeaderFileChars;
    bool exitOnError = false;
    ParallelDescriptor::ReadAndBcastFile(faHeaderFile, faHeaderFileChars, exitOnError);
        
    if(faHeaderFileChars.size() > 0)
    {
        std::string charPtrString(faHeaderFileChars.dataPtr());
        std::istringstream fais(charPtrString, std::istringstream::in);
            
        while(!fais.eof())
        {
            std::string faHeaderName;
            fais >> faHeaderName;
                
            if(fais.eof()) break;
                
            std::string faHeaderFullName(filePath + '/' + faHeaderName + "_H");
            Vector<char> &tempCharArray = faHeaderMap[faHeaderFullName];
            ParallelDescriptor::ReadAndBcastFile(faHeaderFullName, tempCharArray); 
        }
    }
    
    // Create istringstream from Header file
    std::string headerFile(filePath + "/Header");
    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(headerFile, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);
    
    // Start reading
    std::string first_line;
    bool new_checkpoint_format = false;
    std::getline(is, first_line);
    
    // Old vs. new checkpoint format
    if(first_line == CP_VERSION_STR) 
    {
        new_checkpoint_format = true;
        is >> amrlike.spdim;
    } 
    else 
    {
        amrlike.spdim = atoi(first_line.c_str());
    }
    
    // Cumulative time and max level
    is >> amrlike.cumtime;
    int mx_lev;
    is >> mx_lev;
    is >> amrlike.finest_level;
    
    if(ParallelDescriptor::IOProcessor())
    {
        std::cout << "finest level in lower dim: " << amrlike.finest_level <<  std::endl;
    }
     
    // Resize vectors
    amrlike.geom.resize(mx_lev + 1);
    amrlike.ref_ratio.resize(mx_lev);
    amrlike.dt_level.resize(mx_lev + 1);
    amrlike.dt_min.resize(mx_lev + 1);
    amrlike.n_cycle.resize(mx_lev + 1);
    amrlike.level_steps.resize(mx_lev + 1);
    amrlike.level_count.resize(mx_lev + 1);
    amrlike.amrLevels.resize(mx_lev + 1);
    
    // Read info into Amr-like object
    for(int i = 0; i <= mx_lev; ++i) is >> amrlike.geom[i];
    for(int i = 0; i <  mx_lev; ++i) is >> amrlike.ref_ratio[i];
    for(int i = 0; i <= mx_lev; ++i) is >> amrlike.dt_level[i];

    if(new_checkpoint_format)
    {
        for(int i = 0; i <= mx_lev; ++i) is >> amrlike.dt_min[i];
    }
    else
    {
        for(int i = 0; i <= mx_lev; ++i) amrlike.dt_min[i] = amrlike.dt_level[i];
    }

    amrlike.n_cycle.resize(mx_lev + 1);
    for(int i = 0; i <= mx_lev; ++i) is >> amrlike.n_cycle[i];

    for(int i = 0; i <= mx_lev; ++i) is >> amrlike.level_steps[i];
    for(int i = 0; i <= mx_lev; ++i) is >> amrlike.level_count[i];

    // Read levels as in AmrLevel::restart
    for(int lev = 0; lev <= amrlike.finest_level; ++lev)
    {
        AmrLevelLike& amrlev = amrlike.amrLevels[lev];
        
        is >> amrlev.level;
        is >> amrlev.geom;
        
        amrlev.fine_ratio = IntVect::TheUnitVector();
        amrlev.fine_ratio.scale(-1);
        amrlev.crse_ratio = IntVect::TheUnitVector();
        amrlev.crse_ratio.scale(-1);
        
        if(amrlev.level > 0)
        {
            amrlev.crse_ratio = amrlike.ref_ratio[amrlev.level-1];
        }
        if(amrlev.level < mx_lev)
        {
            amrlev.fine_ratio = amrlike.ref_ratio[amrlev.level];
        }
        
        amrlev.grids.readFrom(is);
        
        int nstate;
        is >> nstate;
        int ndesc = nstate;
    
        // Need a distribution mapping for parallel rotations
        amrlev.dmap.define(amrlev.grids);
        
        // Managed pointer to factory
        std::unique_ptr< FabFactory<FArrayBox> > mfab_factory;
        mfab_factory.reset(new FArrayBoxFactory());
        
        amrlev.state.resize(ndesc);
        amrlev.new_state.resize(ndesc);
        
        // Read StateData as in StateData::restart
        for(int i = 0; i < ndesc; i++)
        {
            StateDataLike& sdlike = amrlev.state[i];
            
            is >> sdlike.domain;
            sdlike.grids.readFrom(is);
            
            is >> sdlike.old_time.start;
            is >> sdlike.old_time.stop;
            is >> sdlike.new_time.start;
            is >> sdlike.new_time.stop;
            
            // 1 for new data, 2 for both old and new
            int nsets;
            is >> nsets;
            
            BL_ASSERT(nsets < 3);
            
            sdlike.new_data.reset(new MultiFab(amrlev.grids, amrlev.dmap,
                nstate, nghost, MFInfo().SetTag("StateData"), *mfab_factory));
                
            sdlike.old_data.reset();
            if(nsets >= 2)
            {
                sdlike.old_data.reset(new MultiFab(amrlev.grids, amrlev.dmap,
                    nstate, nghost, MFInfo().SetTag("StateData"), *mfab_factory));
            }
            
            std::string mfab_name;
            std::string full_path;
            MultiFab* mfab = nullptr;
            
            // Load each MultiFab
            for(int i = 0; i < nsets; ++i)
            {
                if(i == 0) mfab = sdlike.new_data.get();
                else if(i == 1) mfab = sdlike.old_data.get();
                
                is >> mfab_name;
                full_path = filePath + "/" + mfab_name;
                
                // Check for preread header
                std::map<std::string, Vector<char> >::iterator mapIter;
                const char *faHeader = nullptr;
                mapIter = faHeaderMap.find(full_path + "_H");
                if(mapIter != faHeaderMap.end())
                {
                    faHeader = mapIter->second.dataPtr();
                }
                
                VisMF::Read(*mfab, full_path, faHeader);
            }
        }
    }
}

int main(int argc, char* argv[])
{
    std::cout << argv[0] << std::endl;
    std::cout << argv[1] << std::endl;
    std::cout << "AMReX init..." << std::endl;
    
    Initialize(argc, argv);
    
    std::cout << "Creating rotator..." << std::endl;
    
    // Create Rotator object from checkpoint file
    // Rotator rotator;
    
    return 0;
}
