// PLEASE READ:
//
// This file implements additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because
// we are "inside" the class, we no longer have to call
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
// Functions in this file can be called from other projects simply by using geometry->cotan(he),
// geometry->barycentricDualArea(v), etc. where "geometry" is a pointer to a VertexPositionGeometry. This avoids having
// to declare a GeometryRoutines object in every project, and also mimics the way that geometry routines are normally
// called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.


#include "geometrycentral/surface/vertex_position_geometry.h"
#include <complex>

namespace geometrycentral {
namespace surface {

/*
 * Compute the Euler characteristic of the mesh.
 */
int VertexPositionGeometry::eulerCharacteristic() const {
    return (int)mesh.nVertices() - (int)mesh.nEdges() + (int)mesh.nFaces();
}

/*
 * Compute the mean length of all the edges in the mesh.
 *
 * Input:
 * Returns: The mean edge length.
 */
double VertexPositionGeometry::meanEdgeLength() const {

    double total = 0.0;
    for (Edge e : mesh.edges()) {
        total += edgeLength(e);
    }
    return total / mesh.nEdges();
}

/*
 * Compute the total surface area of the mesh.
 *
 * Input:
 * Returns: The surface area of the mesh.
 */
double VertexPositionGeometry::totalArea() const {

    double total = 0.0;
    for (Face f : mesh.faces()) {
        total += faceArea(f);
    }
    return total;
}

/*
 * Computes the cotangent of the angle opposite to a halfedge. (Do NOT use built-in function for this)
 *
 * Input: The halfedge whose cotan weight is to be computed.
 * Returns: The cotan of the angle opposite the given halfedge.
 */
double VertexPositionGeometry::cotan(Halfedge he) const {

    // TODO
    Halfedge he1 = he.next();
    Halfedge he2 = he1.next().twin();
    auto l1 = this->edgeLength(he1.edge());
    auto l2 = this->edgeLength(he2.edge());
    auto l3 = this->edgeLength(he.edge());

    auto v1 = this->halfedgeVector(he1);
    auto v2 = this->halfedgeVector(he2);

    double mode = abs(l1 * l2);

    double cos = dot(v1, v2) / mode;
    double sin = cross(v1, v2).norm() / mode;

    return (cos / sin);
}

/*
 * Computes the barycentric dual area of a vertex.
 *
 * Input: The vertex whose barycentric dual area is to be computed.
 * Returns: The barycentric dual area of the given vertex.
 */
double VertexPositionGeometry::barycentricDualArea(Vertex v) const {

    // TODO
    double barycentricDualArea = 0.0;
    for (auto f : v.adjacentFaces())
    {
        barycentricDualArea += faceArea(f) / 3.0;
    }
    return barycentricDualArea; // placeholder
}

/*
 * Computes the angle (in radians) at a given corner. (Do NOT use built-in function for this)
 *
 *
 * Input: The corner at which the angle needs to be computed.
 * Returns: The angle clamped between 0 and π.
 */
double VertexPositionGeometry::angle(Corner c) const {

    // TODO
    auto he1 = c.halfedge();
    if (he1.edge().isBoundary())
    {
        return 0;
    }
    auto l1 = halfedgeVector(he1);
    auto l2 = halfedgeVector(he1.next().next().twin());
    return acos(dot(l1, l2) / (norm(l1) * norm(l2))); // placeholder
}

/*
 * Computes the signed angle (in radians) between two adjacent faces. (Do NOT use built-in function for this)
 *
 * Input: The halfedge (shared by the two adjacent faces) on which the dihedral angle is computed.
 * Returns: The dihedral angle.
 */
double VertexPositionGeometry::dihedralAngle(Halfedge he) const {

    // TODO
    auto n1 = faceNormal(he.face());
    auto n2 = faceNormal(he.twin().face());
    auto e = halfedgeVector(he).normalize();

    return atan2(dot(e, cross(n1, n2)), dot(n1, n2));
}

/*
 * Computes the normal at a vertex using the "equally weighted" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "equally weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalEquallyWeighted(Vertex v) const {

    // TODO
    Vector3 vertexNormal = {0, 0, 0};
    for (const auto& f : v.adjacentFaces())
    {
        vertexNormal += faceNormal(f);
    }
    vertexNormal.normalize();
    return vertexNormal; // placeholder
}

/*
 * Computes the normal at a vertex using the "tip angle weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "tip angle weights" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAngleWeighted(Vertex v) const {

    // TODO
    Vector3 vertexNormal = {0, 0, 0};
    for (const auto& c : v.adjacentCorners())
    {
        vertexNormal += faceNormal(c.face()) * angle(c);
    }
    return vertexNormal.normalize(); // placeholder
}

/*
 * Computes the normal at a vertex using the "inscribed sphere" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "inscribed sphere" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalSphereInscribed(Vertex v) const {

    // TODO
    Vector3 vertexNormal = {0, 0, 0};
    for (const auto& c: v.adjacentCorners())
    {
        const auto he1 = c.halfedge();
        const auto he2 = c.halfedge().next().next().twin();
        const auto vec1 = halfedgeVector(he1);
        const auto vec2 = halfedgeVector(he2);
        vertexNormal += cross(vec1, vec2) / (norm2(vec1) * norm2(vec2));
    }
    return vertexNormal.normalize(); // placeholder
}

/*
 * Computes the normal at a vertex using the "face area weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "face area weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAreaWeighted(Vertex v) const {

    // TODO
    Vector3 vertexNormal = {0, 0, 0};
    for (const auto& f : v.adjacentFaces()) {
        vertexNormal += faceArea(f) * faceNormal(f);
    }
    return vertexNormal.normalize(); // placeholder
}

/*
 * Computes the normal at a vertex using the "Gauss curvature" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "Gauss curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalGaussianCurvature(Vertex v) const {

    // TODO
    Vector3 gaussianCurvatureNorm = {0, 0, 0};
    for (const auto& he: v.outgoingHalfedges()) 
    {
        gaussianCurvatureNorm += dihedralAngle(he) * halfedgeVector(he).normalize();
    }
    return gaussianCurvatureNorm.normalize(); // placeholder
}

/*
 * Computes the normal at a vertex using the "mean curvature" method (equivalent to the "area gradient" method).
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "mean curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalMeanCurvature(Vertex v) const {

    // TODO
    Vector3 meanCurvatureNorm = {0, 0, 0};
    for (const auto& he: v.outgoingHalfedges())
    {
        meanCurvatureNorm += halfedgeVector(he) * (cotan(he) + cotan(he.twin()));
    }
    return meanCurvatureNorm.normalize(); // placeholder
}

/*
 * Computes the angle defect at a vertex.
 *
 * Input: The vertex whose angle defect is to be computed.
 * Returns: The angle defect of the given vertex.
 */
double VertexPositionGeometry::angleDefect(Vertex v) const {

    // TODO
    double totalAngle = 0.0;
    for (const auto& c : v.adjacentCorners())
    {
        totalAngle += angle(c);
    }
    return 2 * PI - totalAngle; // placeholder
}

/*
 * Computes the total angle defect of the mesh.
 *
 * Input:
 * Returns: The total angle defect
 */
double VertexPositionGeometry::totalAngleDefect() const {

    // TODO
    double totalAngleDefect = 0.0;
    for (const auto& v : mesh.vertices())
    {
        totalAngleDefect += angleDefect(v);
    }
    return totalAngleDefect; // placeholder
}

/*
 * Computes the (integrated) scalar mean curvature at a vertex.
 *
 * Input: The vertex whose mean curvature is to be computed.
 * Returns: The mean curvature at the given vertex.
 */
double VertexPositionGeometry::scalarMeanCurvature(Vertex v) const {

    // TODO
    double scalarMeanCurvature = 0.0;
    for (const auto& he : v.outgoingHalfedges())
    {
        scalarMeanCurvature += edgeLength(he.edge()) * dihedralAngle(he); 
    }
    return scalarMeanCurvature /= 2;// placeholder
}

/*
 * Computes the circumcentric dual area of a vertex.
 *
 * Input: The vertex whose circumcentric dual area is to be computed.
 * Returns: The circumcentric dual area of the given vertex.
 */
double VertexPositionGeometry::circumcentricDualArea(Vertex v) const {

    // TODO
    double area = 0.0f;
    for (const auto& he : v.outgoingHalfedges())
    {
        const auto eLength = edgeLength(he.edge());
        const auto cot = cotan(he) + cotan(he.twin());
        area += eLength * eLength * cot;
    }
    return area / 8; // placeholder
}

/*
 * Computes the (pointwise) minimum and maximum principal curvature values at a vertex.
 *
 * Input: The vertex on which the principal curvatures need to be computed.
 * Returns: A std::pair containing the minimum and maximum principal curvature values at a vertex.
 */
std::pair<double, double> VertexPositionGeometry::principalCurvatures(Vertex v) const {

    // TODO
    double circumcentricDualArea = this->circumcentricDualArea(v);
    double H = this->scalarMeanCurvature(v) / circumcentricDualArea;
    double K = this->angleDefect(v) / circumcentricDualArea;
    double k1 = H - sqrt(H * H - K);
    double k2 = H + sqrt(H * H - K);
    return std::make_pair(k1, k2); // placeholder
}


/*
 * Builds the sparse POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace matrix,
 * multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse positive definite Laplace matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::laplaceMatrix() const {

    // TODO
    std::vector<Eigen::Triplet<double>> tripList;
    for (const auto& v: mesh.vertices())
    {
        double sumCot = 0.0f;
        for (const auto& he : v.outgoingHalfedges())
        {
            double cot = edgeCotanWeight(he.edge());
            sumCot += cot;
            tripList.emplace_back(Eigen::Triplet<double>(he.tailVertex().getIndex(), he.tipVertex().getIndex(), -cot));
        }
        tripList.emplace_back(Eigen::Triplet<double>(v.getIndex(), v.getIndex(), sumCot + 1e-8));
    }
    SparseMatrix<double> laplaceMatrix(mesh.nVertices(), mesh.nVertices());
    laplaceMatrix.setFromTriplets(tripList.begin(), tripList.end());
    return laplaceMatrix; // placeholder
}

/*
 * Builds the sparse diagonal mass matrix containing the barycentric dual area of each vertex.
 *
 * Input:
 * Returns: Sparse mass matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::massMatrix() const {

    // TODO
    std::vector<Eigen::Triplet<double>> tripList;
    for (const auto& v: mesh.vertices())
    {
        tripList.emplace_back(Eigen::Triplet<double>(v.getIndex(), v.getIndex(), barycentricDualArea(v)));
    }
    SparseMatrix<double> massMatrix(mesh.nVertices(), mesh.nVertices());
    massMatrix.setFromTriplets(tripList.begin(), tripList.end());
    return massMatrix; // placeholder
}

/*
 * Builds the sparse complex POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace
 * matrix, multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse complex positive definite Laplace matrix for the mesh.
 */
SparseMatrix<std::complex<double>> VertexPositionGeometry::complexLaplaceMatrix() const {

    // TODO
    std::vector<Eigen::Triplet<double>> tripList;
    for (const auto& v: mesh.vertices())
    {
        std::complex<double> value(0.0, 0.0);
        for (const auto& he : v.outgoingHalfedges())
        {
            double cot = edgeCotanWeight(he.edge());
            value += cot;
            tripList.emplace_back(Eigen::Triplet<std::complex<double>>(v.getIndex(), v.getIndex(), std::complex<double>(-cot, 0.0));
        }
        tripList.emplace_back(Eigen::Triplet<std::complex<double>>(v.getIndex(), v.getIndex(), value + 1e-8));
    }
    SparseMatrix<std::complex<double>> result(mesh.nVertices(), mesh.nVertices());
    result.setFromTriplets(tripList.begin(), tripList.end());
    return result; // placeholder
}

/*
 * Compute the center of mass of a mesh.
 */
Vector3 VertexPositionGeometry::centerOfMass() const {

    // Compute center of mass.
    Vector3 center = {0.0, 0.0, 0.0};
    for (Vertex v : mesh.vertices()) {
        center += inputVertexPositions[v];
    }
    center /= mesh.nVertices();

    return center;
}

/*
 * Centers a mesh about the origin.
 * Also rescales the mesh to unit radius if <rescale> == true.
 */
void VertexPositionGeometry::normalize(const Vector3& origin, bool rescale) {

    // Compute center of mass.
    Vector3 center = centerOfMass();

    // Translate to origin [of original mesh].
    double radius = 0;
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] -= center;
        radius = std::max(radius, inputVertexPositions[v].norm());
    }

    // Rescale.
    if (rescale) {
        for (Vertex v : mesh.vertices()) {
            inputVertexPositions[v] /= radius;
        }
    }

    // Translate to origin [of original mesh].
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] += origin;
    }
}

} // namespace surface
} // namespace geometrycentral