////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputVtkBase.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: VTK file format output base class.
//
////////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <string>

#include <boost/core/ignore_unused.hpp>
#include <boost/format.hpp>

#include <LibUtilities/BasicUtils/FileSystem.h>

#include "OutputVtkBase.h"

namespace Nektar
{
namespace FieldUtils
{

// Disable the base VTK factory if using the VTK library, this is so we can
// register the same extension for both.
#if !NEKTAR_USING_VTK
ModuleKey OutputVtkBase::m_className =
    GetModuleFactory().RegisterCreatorFunction(ModuleKey(eOutputModule, "vtu"),
                                               OutputVtkBase::create,
                                               "Writes a VTU file.");
#endif

OutputVtkBase::OutputVtkBase(FieldSharedPtr f) : OutputFileBase(f)
{
    m_requireEquiSpaced = true;
}

OutputVtkBase::~OutputVtkBase()
{
}

void OutputVtkBase::OutputFromPts(po::variables_map &vm)
{
    int i, j;
    LibUtilities::PtsFieldSharedPtr fPts = m_f->m_fieldPts;

    // Extract the output filename and extension
    std::string filename = PrepareOutput(vm);

    // Write solution.
    std::ofstream outfile(filename.c_str());
    WriteVtkHeader(outfile);
    int nfields = 0;
    int dim     = fPts->GetDim();

    int nvert   = 1;
    int vtktype = 1;
    switch (fPts->GetPtsType())
    {
        case LibUtilities::ePtsFile:
        case LibUtilities::ePtsLine:
        {
            NEKERROR(ErrorUtil::efatal,
                     "VTK output needs setting up for ePtsFile and ePtsLine");
            break;
        }
        case LibUtilities::ePtsPlane:
        {
            NEKERROR(ErrorUtil::efatal,
                     "VTK output needs setting up for PtsPlane");
            break;
        }
        case LibUtilities::ePtsBox:
        {
            NEKERROR(ErrorUtil::efatal,
                     "VTK output needs setting up for PtsBox");
            break;
        }
        case LibUtilities::ePtsSegBlock:
        {
            nvert   = 2;
            vtktype = 3;
            break;
        }
        case LibUtilities::ePtsTriBlock:
        {
            nvert   = 3;
            vtktype = 5;
            break;
        }
        case LibUtilities::ePtsTetBlock:
        {
            nvert   = 4;
            vtktype = 10;
            break;
        }
        default:
            NEKERROR(ErrorUtil::efatal, "ptsType not supported yet.");
    }

    std::vector<Array<OneD, int>> ptsConn;
    fPts->GetConnectivity(ptsConn);

    nfields = fPts->GetNFields();

    int nPts      = fPts->GetNpoints();
    int numBlocks = 0;
    for (i = 0; i < ptsConn.size(); ++i)
    {
        numBlocks += ptsConn[i].size() / nvert;
    }

    // write out pieces of data.
    outfile << "    <Piece NumberOfPoints=\"" << nPts << "\" NumberOfCells=\""
            << numBlocks << "\">" << std::endl;
    outfile << "      <Points>" << std::endl;
    outfile << "        <DataArray type=\"Float64\" "
            << "NumberOfComponents=\"" << 3 << "\" format=\"ascii\">" << std::endl;
    for (i = 0; i < nPts; ++i)
    {
        for (j = 0; j < dim; ++j)
        {
            outfile << "          " << std::setprecision(8) << std::scientific
                    << fPts->GetPointVal(j, i) << " ";
        }
        for (j = dim; j < 3; ++j)
        {
            // pack to 3D since paraview does not seem to handle 2D
            outfile << "          0.000000";
        }
        outfile << std::endl;
    }
    outfile << "        </DataArray>" << std::endl;
    outfile << "      </Points>" << std::endl;
    outfile << "      <Cells>" << std::endl;
    outfile << "        <DataArray type=\"Int32\" "
            << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

    // dump connectivity data if it exists
    outfile << "          ";
    int cnt = 1;
    for (i = 0; i < ptsConn.size(); ++i)
    {
        for (j = 0; j < ptsConn[i].size(); ++j)
        {
            outfile << ptsConn[i][j] << " ";
            if ((!(cnt % nvert)) && cnt)
            {
                outfile << std::endl;
                outfile << "          ";
            }
            cnt++;
        }
    }
    outfile << "        </DataArray>" << std::endl;
    outfile << "        <DataArray type=\"Int32\" "
            << "Name=\"offsets\" format=\"ascii\">" << std::endl;

    outfile << "          ";
    for (i = 0; i < numBlocks; ++i)
    {
        outfile << i * nvert + nvert << " ";
    }
    outfile << std::endl;
    outfile << "        </DataArray>" << std::endl;
    outfile << "        <DataArray type=\"UInt8\" "
            << "Name=\"types\" format=\"ascii\">" << std::endl;
    outfile << "          ";
    for (i = 0; i < numBlocks; ++i)
    {
        outfile << vtktype << " ";
    }
    outfile << std::endl;
    outfile << "        </DataArray>" << std::endl;
    outfile << "      </Cells>" << std::endl;
    outfile << "      <PointData>" << std::endl;

    // printing the fields
    for (j = 0; j < nfields; ++j)
    {
        outfile << "        <DataArray type=\"Float64\" Name=\""
                << m_f->m_variables[j] << "\">" << std::endl;
        outfile << "          ";
        for (i = 0; i < fPts->GetNpoints(); ++i)
        {
            outfile << fPts->GetPointVal(dim + j, i) << " ";
        }
        outfile << std::endl;
        outfile << "        </DataArray>" << std::endl;
    }

    outfile << "      </PointData>" << std::endl;
    outfile << "    </Piece>" << std::endl;

    WriteVtkFooter(outfile);
    std::cout << "Written file: " << filename << std::endl;

    // output parallel outline info if necessary
    if ((m_f->m_comm->GetRank() == 0) && (m_f->m_comm->GetSize() != 1))
    {
        WritePVtu(vm);
        std::cout << "Written file: " << filename << std::endl;
    }
}

void OutputVtkBase::OutputFromExp(po::variables_map &vm)
{
    int i, j;
    // Extract the output filename and extension
    std::string filename = PrepareOutput(vm);

    // Write solution.
    std::ofstream outfile(filename.c_str());
    WriteVtkHeader(outfile);
    int nfields = m_f->m_variables.size();

    int nstrips;
    m_f->m_session->LoadParameter("Strip_Z", nstrips, 1);

    // Homogeneous strip variant
    for (int s = 0; s < nstrips; ++s)
    {
        // For each field write out field data for each expansion.
        for (i = 0; i < m_f->m_exp[0]->GetNumElmts(); ++i)
        {
            m_f->m_exp[0]->WriteVtkPieceHeader(outfile, i, s);

            // For this expansion write out each field.
            for (j = 0; j < nfields; ++j)
            {
                m_f->m_exp[s * nfields + j]->WriteVtkPieceData(
                    outfile, i, m_f->m_variables[j]);
            }
            m_f->m_exp[0]->WriteVtkPieceFooter(outfile, i);
        }
    }

    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        WriteEmptyVtkPiece(outfile);
    }

    WriteVtkFooter(outfile);
    std::cout << "Written file: " << filename << std::endl;

    // output parallel outline info if necessary
    if ((m_f->m_comm->GetRank() == 0) && (m_f->m_comm->GetSize() != 1))
    {
        WritePVtu(vm);
    }
}

void OutputVtkBase::OutputFromData(po::variables_map &vm)
{
    boost::ignore_unused(vm);
    NEKERROR(ErrorUtil::efatal, "OutputVtk can't write using only FieldData.");
}

fs::path OutputVtkBase::GetPath(std::string &filename, po::variables_map &vm)
{
    boost::ignore_unused(vm);

    int nprocs = m_f->m_comm->GetSize();
    fs::path specPath;
    if (nprocs == 1)
    {
        specPath = fs::path(filename);
    }
    else
    {
        // replace .vtu by _vtu
        int dot     = filename.find_last_of('.');
        std::string path = filename.substr(0, dot) + "_vtu";
        specPath    = fs::path(path);
    }
    return fs::path(specPath);
}

fs::path OutputVtkBase::GetFullOutName(std::string &filename,
                                       po::variables_map &vm)
{
    int nprocs = m_f->m_comm->GetSize();

    fs::path fulloutname;
    if (nprocs == 1)
    {
        fulloutname = filename;
    }
    else
    {
        // Guess at filename that might belong to this process.
        boost::format pad("P%1$07d.%2$s");
        pad % m_f->m_comm->GetRank() % "vtu";

        // Generate full path name
        fs::path specPath = GetPath(filename, vm);
        fs::path poutfile(pad.str());
        fulloutname = specPath / poutfile;
    }
    return fulloutname;
}

void OutputVtkBase::WriteVtkHeader(std::ostream &outfile)
{
    outfile << "<?xml version=\"1.0\"?>" << std::endl;
    outfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
            << "byte_order=\"LittleEndian\">" << std::endl;
    outfile << "  <UnstructuredGrid>" << std::endl;
}

void OutputVtkBase::WriteVtkFooter(std::ostream &outfile)
{
    outfile << "  </UnstructuredGrid>" << std::endl;
    outfile << "</VTKFile>" << std::endl;
}

void OutputVtkBase::WriteEmptyVtkPiece(std::ofstream &outfile)
{
    // write out empty piece of data.
    outfile << "    <Piece NumberOfPoints=\"" << 0 << "\" NumberOfCells=\"" << 0
            << "\">" << std::endl;
    outfile << "      <Points>" << std::endl;
    outfile << "        <DataArray type=\"Float64\" "
            << "NumberOfComponents=\"" << 3 << "\" format=\"ascii\">" << std::endl;
    outfile << "        </DataArray>" << std::endl;
    outfile << "      </Points>" << std::endl;
    outfile << "      <Cells>" << std::endl;
    outfile << "        <DataArray type=\"Int32\" "
            << "Name=\"connectivity\" format=\"ascii\">" << std::endl;
    outfile << "        </DataArray>" << std::endl;
    outfile << "        <DataArray type=\"Int32\" "
            << "Name=\"offsets\" format=\"ascii\">" << std::endl;

    outfile << "          ";
    outfile << std::endl;
    outfile << "        </DataArray>" << std::endl;
    outfile << "        <DataArray type=\"UInt8\" "
            << "Name=\"types\" format=\"ascii\">" << std::endl;
    outfile << "          ";
    outfile << std::endl;
    outfile << "        </DataArray>" << std::endl;
    outfile << "      </Cells>" << std::endl;
    outfile << "      <PointData>" << std::endl;

    outfile << "      </PointData>" << std::endl;
    outfile << "    </Piece>" << std::endl;
}

void OutputVtkBase::WritePVtu(po::variables_map &vm)
{
    std::string filename = m_config["outfile"].as<std::string>();
    int dot         = filename.find_last_of('.');
    std::string body     = filename.substr(0, dot);
    filename        = body + ".pvtu";

    std::ofstream outfile(filename.c_str());

    int nprocs  = m_f->m_comm->GetSize();
    std::string path = LibUtilities::PortablePath(GetPath(filename, vm));

    outfile << "<?xml version=\"1.0\"?>" << std::endl;
    outfile << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" "
            << "byte_order=\"LittleEndian\">" << std::endl;
    outfile << "<PUnstructuredGrid GhostLevel=\"0\">" << std::endl;
    outfile << "<PPoints> " << std::endl;
    outfile << "<PDataArray type=\"Float64\" NumberOfComponents=\"" << 3
            << "\"/> " << std::endl;
    outfile << "</PPoints>" << std::endl;
    outfile << "<PCells>" << std::endl;
    outfile << "<PDataArray type=\"Int32\" Name=\"connectivity\" "
               "NumberOfComponents=\"1\"/>"
            << std::endl;
    outfile << "<PDataArray type=\"Int32\" Name=\"offsets\"      "
               "NumberOfComponents=\"1\"/>"
            << std::endl;
    outfile << "<PDataArray type=\"UInt8\" Name=\"types\"        "
               "NumberOfComponents=\"1\"/>"
            << std::endl;
    outfile << "</PCells>" << std::endl;
    outfile << "<PPointData Scalars=\"Material\">" << std::endl;
    for (int i = 0; i < m_f->m_variables.size(); ++i)
    {
        outfile << "<PDataArray type=\"Float64\" Name=\"" << m_f->m_variables[i]
                << "\"/>" << std::endl;
    }
    outfile << "</PPointData>" << std::endl;

    for (int i = 0; i < nprocs; ++i)
    {
        boost::format pad("P%1$07d.vtu");
        pad % i;
        outfile << "<Piece Source=\"" << path << "/" << pad.str() << "\"/>"
                << std::endl;
    }
    outfile << "</PUnstructuredGrid>" << std::endl;
    outfile << "</VTKFile>" << std::endl;

    std::cout << "Written file: " << filename << std::endl;
}

std::string OutputVtkBase::PrepareOutput(po::variables_map &vm)
{
    // Extract the output filename and extension
    std::string filename = m_config["outfile"].as<std::string>();

    fs::path specPath    = GetPath(filename, vm);
    fs::path fulloutname = GetFullOutName(filename, vm);
    filename             = LibUtilities::PortablePath(fulloutname);

    if (m_f->m_comm->GetSize() != 1)
    {
        if (m_f->m_comm->TreatAsRankZero())
        {
            try
            {
                fs::create_directory(specPath);
            }
            catch (fs::filesystem_error &e)
            {
                ASSERTL0(false, "Filesystem error: " + std::string(e.what()));
            }
            std::cout << "Writing files to directory: " << specPath << std::endl;
        }
        m_f->m_comm->Block();
    }
    else
    {
        std::cout << "Writing: " << specPath << std::endl;
    }
    return filename;
}

} // namespace FieldUtils
} // namespace Nektar
