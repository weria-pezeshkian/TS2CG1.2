#if !defined(AFX_ARGUMENT_H_7F4A21B8_C13C_11D3_BF23_124095086234__INCLUDED_)
#define AFX_ARGUMENT_H_7F4A21B8_C13C_11D3_BF23_124095086234__INCLUDED_

#include "Def.h"
#include <vector>
#include <string>

class Argument {
public:
    Argument(std::vector<std::string> arg);
    ~Argument();

    inline const std::vector<std::string> GetArgumentString() const { return m_Argument; }
    inline const int GetArgCon() const { return m_ArgCon; }
    inline const std::string GetIn_GroFileName() const { return m_In_GroFileName; }
    inline const std::string GetOut_GroFileName() const { return m_Out_GroFileName; }
    inline const std::string GetTemplateFileName() const { return m_Tem_GroFileName; }
    inline const double GetRCutOff() const { return m_RCutOff; }
    inline const int GetSeed() const { return m_Seed; }
    inline const double GetDB() const { return m_DB; }
    inline const double GetUCELLSize() const { return m_UCELLSize; }
    inline const std::vector<int> GetIon() const { return m_Ion; }
    inline std::string GetNegativeIonName() const { return m_NegName; }
    inline std::string GetPositiveIonName() const { return m_PosName; }

private:
    std::vector<std::string> m_Argument;
    std::string m_In_GroFileName;
    std::string m_Out_GroFileName;
    std::string m_Tem_GroFileName;
    int m_Seed;
    double m_RCutOff;
    int m_ArgCon;
    double m_DB;
    double m_UCELLSize;
    std::vector<int> m_Ion;

public:
    std::string m_NegName;
    std::string m_PosName;

};

#endif
