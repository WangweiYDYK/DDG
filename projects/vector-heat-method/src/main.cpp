#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "vector-heat-method.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr;
std::unique_ptr<VertexPositionGeometry> geometry_uptr;

// so we can more easily pass these to different classes
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;
std::string MESHNAME;

// Some global variables
std::unique_ptr<VectorHeatMethodSolver> solver;
Vector<double> DELTA;                      // sources
polyscope::SurfaceGraphQuantity* currVert; // currently active vertex
Vector<double> SOLUTION;
polyscope::SurfaceVertexColorQuantity* solnColors;
polyscope::SurfaceGraphQuantity* isolines;
double maxPhi = 0.0;
double vertexRadius;
double isolinesRadius;

struct SourceVertex {
	Vertex vertex;
    float scalarValue = 1.;

    // the position of the vertex can be defined by the complex number
    // by transform it into the local coordinate system of the vertex and it's one-ring-neighbor
    float vectorMagnitude = 1.;
    float vectorAngleRad = 0.;
};
std::vector<SourceVertex> sourcePoints;

struct SiteVert {
    Vertex vertex;
    float weight = 1.;
};
std::vector<SiteVert> siteVerts;

bool vizFirstRun = true;
void updateSourceSetViz() {

    // Scalar balls around sources
    std::vector<std::pair<size_t, double>> sourcePairs;
    for (SourceVertex& s : sourcePoints) {
        size_t ind = geometry->vertexIndices[s.vertex];
        sourcePairs.emplace_back(ind, s.scalarValue);
    }
    auto scalarQ = polyscope::getSurfaceMesh()->addVertexIsolatedScalarQuantity("source scalars", sourcePairs);
    scalarQ->setPointRadius(scalarQ->getPointRadius() * 2.f);
    scalarQ
        ->setColorMap("reds");

    if (vizFirstRun) {
        scalarQ->setEnabled(true);
    }

    // Vectors at sources
    VertexData<Vector2> sourceVectors(*mesh, Vector2::zero());
    for (SourceVertex& s : sourcePoints) {
        sourceVectors[s.vertex] = Vector2::fromAngle(s.vectorAngleRad) * s.vectorMagnitude;
    }
    auto vectorQ = polyscope::getSurfaceMesh()->addVertexIntrinsicVectorQuantity("source vectors", sourceVectors);
    vectorQ->setVectorLengthScale(.05);
    vectorQ->setVectorRadius(.005);
    vectorQ->setVectorColor(glm::vec3{227 / 255., 52 / 255., 28 / 255.});
    if (vizFirstRun) {
        vectorQ->setEnabled(true);
    }

    vizFirstRun = false;
}

bool vizFirstRunSite = true;
void updateSiteSetViz() {

    // Scalar balls around sources
    std::vector<std::pair<size_t, double>> sourcePairs;
    for (SiteVert& s : siteVerts) {
        size_t ind = geometry->vertexIndices[s.vertex];
        sourcePairs.emplace_back(ind, s.weight);
    }
    auto scalarQ = polyscope::getSurfaceMesh()->addVertexIsolatedScalarQuantity("averaging sites", sourcePairs);
    scalarQ->setPointRadius(scalarQ->getPointRadius() * 2.);
    scalarQ->setColorMap("blues");
    if (vizFirstRunSite) {
        scalarQ->setEnabled(true);
    }

    vizFirstRunSite = false;
}

void addSourceVertex(size_t index) 
{
    Vertex v = mesh->vertex(index);
    for (const auto& s : sourcePoints) {
        if (s.vertex == v) {
			polyscope::warning("vertex already a source");
			return;
		}
    }
    SourceVertex sv;
    sv.vertex = v;
    sourcePoints.emplace_back(sv);
    updateSourceSetViz();
}

void addVertexSite(size_t ind) {
    Vertex v = mesh->vertex(ind);

    // Make sure not already used
    for (SiteVert& s : siteVerts) {
                if (s.vertex == v) {
                        std::stringstream ss;
                        ss << "Vertex " << v;
                        std::string vStr = ss.str();
                        polyscope::warning("Vertex " + vStr + " is already a site");
                        return;
                }
    }

    SiteVert newV;
    newV.vertex = v;
    newV.weight = 1.0;
    siteVerts.push_back(newV);
    updateSiteSetViz();
}

void extendScalar() 
{
}
 
void vectorTransport() 
{
    if (solver == nullptr) {
                solver.reset(new VectorHeatMethodSolver(*geometry));
    }
    if (sourcePoints.size() == 0) {
                polyscope::warning("no source points set");
                return;
    }
    std::vector<std::tuple<SurfacePoint, Vector2>> points;
    for (const auto& p:sourcePoints)
    {
                points.emplace_back(std::make_tuple(p.vertex, Vector2::fromAngle(p.vectorAngleRad) * p.vectorMagnitude));
    }
    VertexData<Vector2> result = solver->transportTangentVectors(points);
    auto psVec = psMesh->addVertexIntrinsicVectorQuantity("vector extension", result);
    psVec->setEnabled(true);
}

void computeLogMap() {}

void computeCenter() {}

void flipZ() {
    // Rotate mesh 180 deg about up-axis on startup
    glm::mat4x4 rot = glm::rotate(glm::mat4x4(1.0f), static_cast<float>(PI), glm::vec3(0, 1, 0));
    for (Vertex v : mesh->vertices()) {
        Vector3 vec = geometry->inputVertexPositions[v];
        glm::vec4 rvec = {vec[0], vec[1], vec[2], 1.0};
        rvec = rot * rvec;
        geometry->inputVertexPositions[v] = {rvec[0], rvec[1], rvec[2]};
    }
    psMesh->updateVertexPositions(geometry->inputVertexPositions);
}

/*
 * Map solution to mesh colors.
 */
//std::vector<std::array<double, 3>> computeColors(const Vector<double>& sol) {
//
//     Determine maximum-magnitude element for scaling purposes
//    maxPhi = 0;
//    for (size_t i = 0; i < mesh->nVertices(); i++) {
//        maxPhi = std::max(maxPhi, sol[i]); // distances should be >= 0
//    }
//
//    std::vector<std::array<double, 3>> colors;
//    for (Vertex v : mesh->vertices()) {
//        colors.push_back(mapToColor(maxPhi - sol[v.getIndex()], 0, maxPhi, "hot")); // invert colormap
//    }
//    return colors;
//}

// Set mesh color.
void setColors(const std::vector<std::array<double, 3>>& colors) {
    solnColors = psMesh->addVertexColorQuantity("Solution", colors);
    solnColors->setEnabled(true);
}


/*
 * Display isolines.
 */
void showIsolines() {

    std::vector<Vector3> positions;
    std::vector<std::array<size_t, 2>> edgeInds;
    double distBetweenLines = maxPhi / 20.0; // enforce spacing
    for (Face f : mesh->faces()) {
        std::vector<Vector3> pos;
        for (Halfedge he : f.adjacentHalfedges()) {
            double vs = SOLUTION[he.tailVertex().getIndex()];
            double vd = SOLUTION[he.tipVertex().getIndex()];
            int region1 = floor(vs / distBetweenLines);
            int region2 = floor(vd / distBetweenLines);
            if (region1 != region2) {
                double val = region2 * distBetweenLines;
                if (region1 > region2) {
                    val = region1 * distBetweenLines;
                }
                double t = (val - vs) / (vd - vs);
                Vector3 ps = geometry->inputVertexPositions[he.tailVertex()];
                Vector3 pd = geometry->inputVertexPositions[he.tipVertex()];
                Vector3 p = ps + t * (pd - ps);
                pos.push_back(p);
            }
        }
        if (pos.size() == 2) {
            positions.push_back(pos[0]);
            positions.push_back(pos[1]);
            edgeInds.push_back({positions.size() - 2, positions.size() - 1});
        }
    }
    isolines = psMesh->addSurfaceGraphQuantity("Isolines", positions, edgeInds);
    isolines->setEnabled(true);
    isolines->setRadius(isolinesRadius);
    isolines->setColor({0.0, 0.0, 0.0});
}

/*
 * Show selected vertices.
 * This function gets called every time an element is selected on-screen.
 */
void showSelected() {

    // Show selected vertices in yellow
    std::vector<Vector3> vertPos;
    std::vector<std::array<size_t, 2>> vertInd;
    DELTA = Vector<double>::Zero(mesh->nVertices());
    for (std::set<size_t>::iterator it = polyscope::state::subset.vertices.begin();
         it != polyscope::state::subset.vertices.end(); ++it) {
        vertPos.push_back(geometry->inputVertexPositions[*it]);
        DELTA[*it] = 1;
    }
    polyscope::SurfaceGraphQuantity* showVerts = psMesh->addSurfaceGraphQuantity("selected vertices", vertPos, vertInd);
    showVerts->setEnabled(true);
    showVerts->setRadius(vertexRadius);
    showVerts->setColor({1.0, 0.65, 0.0});

    // Show the currently selected vertex in red.
    int currIdx = polyscope::state::currVertexIndex;
    if (currIdx != -1) {
        std::vector<Vector3> pos = {geometry->inputVertexPositions[currIdx]};
        currVert = psMesh->addSurfaceGraphQuantity("current vertex", pos, std::vector<std::array<size_t, 2>>());
        currVert->setEnabled(true);
        currVert->setRadius(vertexRadius);
        currVert->setColor({1.0, 0.0, 0.0});
    } else {
        currVert->setEnabled(false);
    }
}

void redraw() {
    showSelected();
    polyscope::requestRedraw();
}

void buildPointsMenu() 
{
    bool anyChanged = false;

    ImGui::PushItemWidth(200);

    int id = 0;
    int eraseInd = -1;
    for (SourceVertex& s : sourcePoints) {
        std::stringstream ss;
        ss << "Vertex " << s.vertex;
        std::string vStr = ss.str();
        ImGui::PushID(vStr.c_str());

        ImGui::TextUnformatted(vStr.c_str());

        ImGui::SameLine();
        if (ImGui::Button("delete")) {
            eraseInd = id;
            anyChanged = true;
        }
        ImGui::Indent();

        if (ImGui::InputFloat("scalar value", &s.scalarValue)) anyChanged = true;
        if (ImGui::InputFloat("vector mag", &s.vectorMagnitude)) anyChanged = true;
        if (ImGui::SliderAngle("vector angle", &s.vectorAngleRad)) anyChanged = true;

        ImGui::Unindent();
        ImGui::PopID();
    }
    ImGui::PopItemWidth();

    // actually do erase, if requested
    if (eraseInd != -1) {
        sourcePoints.erase(sourcePoints.begin() + eraseInd);
    }

    if (ImGui::Button("add point")) {
        long long int pickVert = polyscope::getSurfaceMesh()->selectVertex();
        if (pickVert >= 0) {
            addSourceVertex(pickVert);
            anyChanged = true;
        }
    }

    if (anyChanged) {
        //updateSourceSetViz();
    }
}

void buildSitesMenu() {

    bool anyChanged = false;

    ImGui::PushItemWidth(200);

    int id = 0;
    int eraseInd = -1;
    for (SiteVert& s : siteVerts) {
        std::stringstream ss;
        ss << "Vertex " << s.vertex;
        std::string vStr = ss.str();
        ImGui::PushID(vStr.c_str());

        ImGui::TextUnformatted(vStr.c_str());

        ImGui::SameLine();
        if (ImGui::Button("delete")) {
            eraseInd = id;
            anyChanged = true;
        }
        ImGui::Indent();

        if (ImGui::InputFloat("weight", &s.weight)) anyChanged = true;

        ImGui::Unindent();
        ImGui::PopID();
    }
    ImGui::PopItemWidth();

    // actually do erase, if requested
    if (eraseInd != -1) {
        siteVerts.erase(siteVerts.begin() + eraseInd);
    }

    if (ImGui::Button("add site")) {
        long long int pickVert = polyscope::getSurfaceMesh()->selectVertex();
        if (pickVert >= 0) {
            addVertexSite(pickVert);
            anyChanged = true;
        }
    }

    if (anyChanged) {
        updateSiteSetViz();
    }
}

void functionCallback() {

    if (ImGui::Button("Solve")) {
        //if (polyscope::state::subset.vertices.size() > 0) {
        //    SOLUTION = HM.compute(DELTA);
        //    if (SOLUTION.norm() > 0) {
        //        setColors(computeColors(SOLUTION));
        //        showIsolines();
        //        redraw();
        //    }
        //}
    }
    if (ImGui::Button("Reset")) {
        polyscope::state::subset.vertices.clear();
        polyscope::state::currVertexIndex = -1;
        psMesh->setSurfaceColor({1.0, 0.45, 0.0});
        solnColors->setEnabled(false);
        isolines->setEnabled(false);
        redraw();
    }

    if (ImGui::TreeNode("select source points"))
    {
        //buildPointsMenu();
        ImGui::TreePop();
    }

    if (ImGui::Button("Vector Transport")) {
        vectorTransport();
		redraw();
    }
}

int main(int argc, char** argv) {

    // Configure the argument parser
    args::ArgumentParser parser("15-458 HW5");
    args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

    // Parse args
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help&) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    // If a mesh name was not given, use default mesh.
    std::string filepath = "../../../input/bunny.obj";
    if (inputFilename) {
        filepath = args::get(inputFilename);
    }

    MESHNAME = polyscope::guessNiceNameFromPath(filepath);

    // Load mesh
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();

    // Get indices for element picking
    polyscope::state::facePickIndStart = mesh->nVertices();
    polyscope::state::edgePickIndStart = polyscope::state::facePickIndStart + mesh->nFaces();
    polyscope::state::halfedgePickIndStart = polyscope::state::edgePickIndStart + mesh->nEdges();

    // Initialize polyscope
    polyscope::init();

    // Set the callback function
    polyscope::state::userCallback = functionCallback;

    // Add mesh to GUI
    psMesh = polyscope::registerSurfaceMesh(MESHNAME, geometry->inputVertexPositions, mesh->getFaceVertexList(),
                                            polyscopePermutations(*mesh));

    // Initalize
    flipZ();
    double lengthScale = geometry->meanEdgeLength();
    vertexRadius = lengthScale * 0.2;
    isolinesRadius = lengthScale * 0.05;
    psMesh->setSmoothShade(true);
    psMesh->setSurfaceColor({1.0, 0.45, 0.0}); // orange
    currVert =
        psMesh->addSurfaceGraphQuantity("current vertex", std::vector<Vector3>(), std::vector<std::array<size_t, 2>>());
    // initialize to something in case "Reset" is pressed before anything happens
    solnColors = psMesh->addVertexColorQuantity("Solution", std::vector<std::array<double, 3>>(mesh->nVertices()));
    isolines =
        psMesh->addSurfaceGraphQuantity("Isolines", std::vector<Vector3>(), std::vector<std::array<size_t, 2>>());
    //HM = VectorHeatMethodSolver(geometry, 1.0f);
    DELTA = Vector<double>::Zero(mesh->nVertices());

    // Give control to the polyscope gui
    polyscope::show();

    delete mesh;
    delete geometry;

    return EXIT_SUCCESS;
}