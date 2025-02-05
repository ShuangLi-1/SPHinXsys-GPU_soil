#include "io_vtk_mesh.h"

namespace SPH
{
//=================================================================================================//
BodyStatesRecordingToMeshVtp::BodyStatesRecordingToMeshVtp(SPHBody &body, ANSYSMesh &ansys_mesh)
    : BodyStatesRecordingToVtp(body), node_coordinates_(ansys_mesh.node_coordinates_),
      elements_nodes_connection_(ansys_mesh.elements_nodes_connection_) {}
//=================================================================================================//
BodyStatesRecordingToMeshVtu::BodyStatesRecordingToMeshVtu(SPHBody &body, ANSYSMesh &ansys_mesh)
    : BodyStatesRecordingToVtp(body), node_coordinates_(ansys_mesh.node_coordinates_),
      elements_nodes_connection_(ansys_mesh.elements_nodes_connection_), bounds_(body) {};
//=================================================================================================//
void BodyStatesRecordingToMeshVtu::FileHeader(std::ofstream &out_file)
{
    out_file << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
    out_file << "<UnstructuredGrid>\n";
    out_file << "<FieldData>\n";
    out_file << "<DataArray type=\"Int32\" Name=\"ispatch\" NumberOfTuples=\"1\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"0\">\n";
    out_file << "0\n";
    out_file << "</DataArray>\n";
    out_file << "</FieldData>\n";
}
//=================================================================================================//
void BodyStatesRecordingToMeshVtu::FileNodeCoordinates(
    std::ofstream &out_file, StdLargeVec<Vecd> &node_coordinates_,
    StdLargeVec<StdVec<size_t>> &elements_nodes_connection_, SPHBody &bounds_, Real &rangemax)
{
    // Write point data
    out_file << "<Piece NumberOfPoints=\"" << node_coordinates_.size() << "\" NumberOfCells=\"" << elements_nodes_connection_.size() << "\">\n";
    out_file << "<PointData>\n";
    out_file << "</PointData>\n";
    out_file << "<CellData>\n";
    out_file << "</CellData>\n";
    out_file << "<Points>\n";
    BoundingBox bounds = bounds_.getSPHSystemBounds();
    Real first_max = bounds.first_.cwiseAbs().maxCoeff();
    Real second_max = bounds.second_.cwiseAbs().maxCoeff();
    rangemax = 1.03075 * (std::max(first_max, second_max));
    out_file << "<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"" << rangemax << "\">\n";

    size_t total_nodes = node_coordinates_.size();
    for (size_t node = 0; node != total_nodes; ++node)
    {
        out_file << node_coordinates_[node][0] << " " << node_coordinates_[node][1] << " " << node_coordinates_[node][2] << " \n";
    }
}
//=================================================================================================//
void BodyStatesRecordingToMeshVtu::FileInformationKey(std::ofstream &out_file, Real &rangemax)
{
    out_file << "<InformationKey name=\"L2_NORM_RANGE\" location=\"vtkDataArray\" length=\"2\">\n";
    out_file << "<Value index=\"0\">\n";
    out_file << "0\n";
    out_file << "</Value>\n";
    out_file << "<Value index=\"1\">\n";
    out_file << rangemax << " \n";
    out_file << "</Value>\n";
    out_file << "</InformationKey>\n";
    out_file << "<InformationKey name=\"L2_NORM_FINITE_RANGE\" location=\"vtkDataArray\" length=\"2\">\n";
    out_file << "<Value index=\"0\">\n";
    out_file << "0\n";
    out_file << "</Value>\n";
    out_file << "<Value index=\"1\">\n";
    out_file << rangemax << " \n";
    out_file << "</Value>\n";
    out_file << "</InformationKey>\n";
    out_file << "</DataArray>\n";
    out_file << "</Points>\n";
}
//=================================================================================================//
void BodyStatesRecordingToMeshVtu::FileCellConnectivity(
    std::ofstream &out_file,
    StdLargeVec<StdVec<size_t>> &elements_nodes_connection_, StdLargeVec<Vecd> &node_coordinates_)
{
    out_file << "<Cells>\n";
    // Write Cell data
    out_file << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"" << node_coordinates_.size() - 1 << "\">\n";

    for (const auto &cell : elements_nodes_connection_)
    {
        for (const auto &vertex : cell)
        {
            out_file << vertex << " ";
        }
        out_file << "\n";
    }

    out_file << "</DataArray>\n";
}
//=================================================================================================//
void BodyStatesRecordingToMeshVtu::FileOffsets(std::ofstream &out_file,
                                               StdLargeVec<StdVec<size_t>> &elements_nodes_connection_)
{
    out_file << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\" RangeMin=\"4\" RangeMax=\"" << 4 * elements_nodes_connection_.size() << "\">\n";

    size_t offset = 0;
    for (const auto &face : elements_nodes_connection_)
    {
        offset += face.size();
        out_file << offset << " ";
    }
    out_file << "\n</DataArray>\n";
}
//=================================================================================================//
void BodyStatesRecordingToMeshVtu::FileTypeOfCell(std::ofstream &out_file,
                                                  StdLargeVec<StdVec<size_t>> &elements_nodes_connection_)
{
    size_t type = 10;
    // Specifies type of cell 10 = tetrahedral
    out_file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\" RangeMin=\"10\" RangeMax=\"10\">\n";
    for (const auto &types : elements_nodes_connection_)
    {
        for (size_t i = 0; i < types.size(); ++i)
        {
            out_file << type << " ";
        }
        out_file << "\n";
    }
    // Write face attribute data
    out_file << "</DataArray>\n";
    out_file << "</Cells>\n";
}
//=================================================================================================//
} // namespace SPH
