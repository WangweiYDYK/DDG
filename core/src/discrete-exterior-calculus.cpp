// PLEASE READ:
//
// This file additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because we are
// "inside" the class, we no longer have to call
//
//          geometry->inputVertexPositions[v], etc.
//
// We can just call
//
//          this->inputVertexPositions[v], etc.
//
// or simply
//
//          inputVertexPositions[v], etc.
//
// In addition, we no longer access the corresponding surface mesh via
//
//          mesh->vertices(), etc.
//
// but instead <mesh> is not a pointer anymore, so we use
//
//          mesh.vertices(), etc.
//
// Functions in this file can be called from other projects simply by using geometry->buildHodgeStar0Form(), etc. where
// "geometry" is a pointer to a VertexPositionGeometry. This avoids having to declare a GeometryRoutines object in every
// project, and also mimics the way that geometry routines are normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {


/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {

    // TODO
    int n = mesh.nVertices();
    std::vector<Eigen::Triplet<double>> tripList;
    SparseMatrix<double> HodgeStar(n, n);
    for (const auto& v: mesh.vertices())
    {
        const auto vid = v.getIndex();
        double A = barycentricDualArea(v);
        tripList.emplace_back(Eigen::Triplet<double>(vid, vid, A));
    }
    HodgeStar.setFromTriplets(tripList.begin(), tripList.end());
    return HodgeStar; // placeholder
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {

    // TODO
    int n = mesh.nEdges();
    std::vector<Eigen::Triplet<double>> tripList;
    SparseMatrix<double> HodgeStar(n, n);
    for (const auto& e: mesh.edges())
    {
        const auto eid = e.getIndex();
        double cotValue = 0.0;
        cotValue += cotan(e.halfedge());
        cotValue += cotan(e.halfedge().twin());
        tripList.emplace_back(Eigen::Triplet<double>(eid, eid, 0.5 * cotValue));
    }
    HodgeStar.setFromTriplets(tripList.begin(), tripList.end());
    return HodgeStar; // placeholder
}

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {

    // TODO
    int n = mesh.nFaces();
    SparseMatrix<double> HodgeStar(n, n);
    for (const auto& f : mesh.faces())
    {
        HodgeStar.insert(f.getIndex(), f.getIndex()) = 1 / faceArea(f);
    }
    return HodgeStar; // placeholder
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {

    // TODO
    std::vector<Eigen::Triplet<double>> tripList;
    SparseMatrix<double> exteriorDerivative(mesh.nEdges(), mesh.nVertices());
    for (const auto& e: mesh.edges())
    {
        tripList.emplace_back(Eigen::Triplet<double>(e.getIndex(), e.firstVertex().getIndex(), -1));
        tripList.emplace_back(Eigen::Triplet<double>(e.getIndex(), e.secondVertex().getIndex(), 1));
    }
    exteriorDerivative.setFromTriplets(tripList.begin(), tripList.end());
    return exteriorDerivative; // placeholder
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {

    // TODO
    std::vector<Eigen::Triplet<double>> tripList;
    SparseMatrix<double> exteriorDerivative(mesh.nFaces(), mesh.nEdges());
    for (const auto& f: mesh.faces())
    {
        for (const auto& he : f.adjacentHalfedges())
        {
            double coeff = he.orientation() ? 1.0 : -1.0;
            tripList.emplace_back(Eigen::Triplet<double>(f.getIndex(), he.edge().getIndex(), coeff));
		}
    }
    exteriorDerivative.setFromTriplets(tripList.begin(), tripList.end());
    return exteriorDerivative; // placeholder
}

} // namespace surface
} // namespace geometrycentral