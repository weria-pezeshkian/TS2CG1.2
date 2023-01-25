#if !defined(AFX_ARGUMENT_H_7F4A21B8_C13C_11D3_BF23_124095086234__INCLUDED_)
#define AFX_ARGUMENT_H_7F4A21B8_C13C_11D3_BF23_124095086234__INCLUDED_

#include "SimDef.h"


class Argument
{
public:
    
	Argument(std::vector <std::string> arg);
	 ~Argument();
    
public:
        inline const std::vector <std::string> GetArgumentString()      const {return m_Argument;}
    
        inline const int         GetArgCon()                            const {return m_ArgCon;}
        inline const std::string GetTrjFile()       			const {return m_TrjFile;}
        inline const int GetBStep()           			const {return m_IniStep;}
        inline const int GetEStep()            			const {return m_FinStep;}
        inline const int GetSeed()        	                const {return m_Seed;}
	inline const std::string GetGeneralOutputFilename()             	const {return m_GeneralOutputFilename;}
	inline const int GetJump()             	const {return m_Jump;}
    inline const std::string GetIndexFileName()                 const {return m_IndexFileName;}
    inline const std::string GetIndexFileName2()                 const {return m_IndexFileName2;}
    inline const std::string GetGroFileName()                 const {return m_GroFileName;}
    inline const std::string GetJobType()                 const {return m_AnalysisType;}
    inline const bool GetHealth()            const {return m_Health;}
    inline  double GetCutoff()            const {return m_Cutoff;}
    inline  double GetGrayAngle()            const {return m_GrayAngle;}
    inline  std::string GetMoreData()            const {return m_MoreData;}
    inline  std::string GetpdbFileName()            const {return m_pdbFileName;}
    inline  int GetBindsNumber()            const {return m_NoBins;}
    inline  int GetNoAverageFrame()            const {return m_TimeAverageFrames;}
    
    inline  std::string GetTargetGroupName()            const {return m_TargetGroupName;}
    inline  double GetBinLength()            const {return m_binL;}
    inline  std::string GetReferenceGroupName()            const {return m_ReferenceGroupName;}



    inline const double GetSoftwareVersion()                 const {return m_SoftWareVersion;}

    // // =================

private:
    std::vector <std::string> m_Argument;
    std::string m_excutablename;
    std::string m_GroFileName;
    std::string m_TrjFile;    //
    std::string m_MoreData;    //
    int m_IniStep;          // Initial step number: the is zero mostly except for rerun
    int m_FinStep;          // Final step number.
    std::string m_IndexFileName;
    std::string m_IndexFileName2;
    int m_ArgCon;

    bool m_Health;
    int m_Seed;
    std::string  m_GeneralOutputFilename;
    int m_Jump;
    double m_SoftWareVersion;
    std::string m_AnalysisType;
    double m_Cutoff;
    double m_GrayAngle;
    int m_NoBins;
    double m_binL;
    int m_TimeAverageFrames;
    std::string m_pdbFileName;
    
    std::string m_TargetGroupName;
    std::string m_ReferenceGroupName;


};

#endif
