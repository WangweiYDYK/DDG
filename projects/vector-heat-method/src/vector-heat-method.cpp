#include "vector-heat-method.h"

VectorHeatMethodSolver::VectorHeatMethodSolver(IntrinsicGeometryInterface& geom) : mesh(geom.mesh), geom(geom) {
    // init the timestep by caculate the mean length of edges
    geom.requireEdgeLengths();
    double meanEdgeLength = 0;
    for (const auto& e : mesh.edges()) {
        meanEdgeLength += geom.edgeLengths[e];
    }
    meanEdgeLength /= mesh.nEdges();
    this->timestep = meanEdgeLength * meanEdgeLength;
    geom.unrequireEdgeLengths();

    // init the mass matrix
    geom.requireVertexLumpedMassMatrix();
    this->massMatrix = geom.vertexLumpedMassMatrix;
    geom.unrequireVertexLumpedMassMatrix();

    // init the cot laplace matrix
    geom.requireCotanLaplacian();
    this->LC = geom.cotanLaplacian;
    geom.unrequireCotanLaplacian();

    // init the laplace matrix
    this->scalarOp = this->massMatrix + timestep * this->LC;
    scalarHeatSolver.reset(new PositiveDefiniteSolver<double>(this->scalarOp));
    poissonSolver.reset(new PositiveDefiniteSolver<double>(this->LC));

    geom.requireVertexConnectionLaplacian();
    this->Lconn = geom.vertexConnectionLaplacian;
    this->vectorOp = this->massMatrix.cast<std::complex<double>>() + timestep * this->Lconn;
    vectorHeatSolver.reset(new SquareSolver<std::complex<double>>(this->vectorOp));
    geom.unrequireVertexConnectionLaplacian();
}

VertexData<double> VectorHeatMethodSolver::extendScalar(const std::vector<std::tuple<SurfacePoint, double>>& sources) 
{
    if (sources.empty()) {
        return VertexData<double>(mesh, std::numeric_limits<double>::quiet_NaN());
    }

    // define the RHS
    // deltaRHS which indicate whether the point is a source
    Vector<double> deltaRHS = Vector<double>::Zero(mesh.nVertices());
    // scalarValueRHS is the scalarValue of the point to solve heat diffusion equation
    Vector<double> scalarValueRHS = Vector<double>::Zero(mesh.nVertices());

    // first step
    // get the RHS of the two equation
    geom.requireVertexIndices();
    for (const auto& tup : sources) {
        auto pt = std::get<0>(tup);
        auto scalarValue = std::get<1>(tup);

        auto facePoint = pt.inSomeFace();
        
        {
            Halfedge he = facePoint.edge.halfedge();
            size_t index = geom.vertexIndices[he.vertex()];
            double w = facePoint.faceCoords.x;
            scalarValueRHS[index] = w * scalarValue;
            deltaRHS[index] = w;
            he = he.next();
            // caculate the coef of second adjacent vertex
            index = geom.vertexIndices[he.vertex()];
            w = facePoint.faceCoords.y;
            scalarValueRHS[index] = w * scalarValue;
            deltaRHS[index] = w;
            he = he.next();
            // caculate the coef of third adjacent vertex
            index = geom.vertexIndices[he.vertex()];
            w = facePoint.faceCoords.z;
            scalarValueRHS[index] = w * scalarValue;
            deltaRHS[index] = w;
        }
    }
    geom.unrequireVertexIndices();

    // second step
    // solve the two equation
    Vector<double> scalarValueSolution = scalarHeatSolver->solve(scalarValueRHS);
    Vector<double> deltaSolution = scalarHeatSolver->solve(deltaRHS);

    // third step
    // get the result
    Vector<double> interpResult = scalarValueSolution.array() / deltaSolution.array();
    VertexData<double> result(mesh, interpResult);

    return result;
}

VertexData<Vector2>
VectorHeatMethodSolver::transportTangentVectors(const std::vector<std::tuple<SurfacePoint, Vector2>>& sources) {
    assert(sources.empty());
    const bool singleVec = sources.size() == 1;

    // first step
    // Integrate the vector heat flow
    Vector<std::complex<double>> dirRHS = Vector<std::complex<double>>::Zero(mesh.nVertices());

    std::vector<std::tuple<SurfacePoint, double>> magnitudeSources;

    for (const auto& tup : sources) {
        SurfacePoint pt = std::get<0>(tup);
        Vector2 vec = std::get<1>(tup);

        magnitudeSources.emplace_back(pt, vec.norm());
        std::complex<double> vecComplex = Vector2::fromComplex(vec).normalize();
        SurfacePoint facePoint = pt.inSomeFace();
        Halfedge he = facePoint.edge.halfedge();

        // get the dirRHS of the facePoint
        {
            // caculate the coef of first adjacent vertex
            size_t index = geom.vertexIndices[he.vertex()];
            double w = facePoint.faceCoords.x;
            dirRHS[index] = w * vecComplex;
            he = he.next();
            // caculate the coef of second adjacent vertex
            index = geom.vertexIndices[he.vertex()];
            w = facePoint.faceCoords.y;
            dirRHS[index] = w * vecComplex;
            he = he.next();
            // caculate the coef of third adjacent vertex
            index = geom.vertexIndices[he.vertex()];
            w = facePoint.faceCoords.z;
            dirRHS[index] = w * vecComplex;
        }
    }

    // solve the vector diffusion equation
    Vector<std::complex<double>> vecSolution = vectorHeatSolver->solve(dirRHS);

    // second step
    // get the magnitude right

    VertexData<Vector2> result(mesh);

    geom.requireVertexIndices();
    if (singleVec) {
        double sourceNorm = std::get<1>(sources[0]).norm();
        vecSolution = vecSolution.array() / vecSolution.array().abs() * sourceNorm;
        for (const auto& v : mesh.vertices())
        {
            result[v] = Vector2::fromComplex(vecSolution[geom.vertexIndices[v]]);
        }
    } else {
        // For multiple sources, need to interpolate magnitudes

        // === Perform scalar interpolation
        VertexData<double> interpMags = extendScalar(magnitudeSources);

        // Scale and copy to result
        for (Vertex v : mesh.vertices()) {
            Vector2 dir = Vector2::fromComplex(vecSolution[geom.vertexIndices[v]]).normalize();
            result[v] = dir * interpMags[v];
        }
    }

    geom.unrequireVertexIndices();
    return result;
}

VertexData<Vector2> VectorHeatMethodSolver::computeLogMap(const Vertex& sourceVert, double vertexDistanceShift = 0.) 
{

    // to compute the log map of descrete surface, we need to solve the poisson equation
    // first need to get the descrete radialRHS

    Vector<std::complex<double>> radialRHS = Vector<std::complex<double>>::Zero(mesh.nVertices());
    buildRadialRHS();

    Vector<std::complex<double>> radialSol = vectorHeatSolver->solve(radialRHS);

}