#pragma once
#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/heat_method_distance.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"

#include <complex>
#include <memory>
#include <tuple>
#include <vector>

using namespace geometrycentral;
using namespace geometrycentral::surface;

class VectorHeatMethodSolver {
public:
    VectorHeatMethodSolver(IntrinsicGeometryInterface& geom);

    std::unique_ptr<PositiveDefiniteSolver<double>> scalarHeatSolver;
    std::unique_ptr<LinearSolver<std::complex<double>>> vectorHeatSolver;
    std::unique_ptr<PositiveDefiniteSolver<double>> poissonSolver;

    // method
    // extend the scalar
    VertexData<double> extendScalar(const std::vector<std::tuple<SurfacePoint, double>>& sources);

    // transport the tangent vector
    // if the sources is not single, then need to use extendScalar to get the scalar
    VertexData<Vector2> transportTangentVectors(const std::vector<std::tuple<SurfacePoint, Vector2>>& sources);

    // compute the log map
    // log map is a map for M -> TpM
    VertexData<Vector2> computeLogMap(const Vertex& sourceVert, double vertexDistanceShift = 0.);
    void buildRadialRHS(Vertex vert, Vector<std::complex<double>>& distGradRHS);

private:
    SurfaceMesh& mesh;
    IntrinsicGeometryInterface& geom;

    // matrix for solve the equation
    SparseMatrix<double> massMatrix;
    SparseMatrix<double> LC;
    SparseMatrix<double> scalarOp;
    SparseMatrix<double> L;
    SparseMatrix<std::complex<double>> Lconn;
    SparseMatrix<std::complex<double>> vectorOp;
    double timestep;
};