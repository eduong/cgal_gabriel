#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <boost/chrono.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/pending/disjoint_sets.hpp>

#include "GraphDefs.h"
#include "Timer.h"
#include "GraphUtil.h"
#include "GraphParse.h"
#include "GraphGen.h"

#include <fstream>
#include <string>

#define SHOW_DEBUG false

void printGraph(const char* title, VertexVector* vertices, EdgeVector* edges, bool printVertices = false) {
	std::cout << std::endl << "=== " << title << std::endl;

	if (printVertices) {
		// Iterate through the vertices and print them out
		for (int i = 0; i < vertices->size(); i++) {
			CGALPoint* v = (*vertices)[i];
			std::cout << "index: " << v << " (" << (*v) << ")" << std::endl;
		}
	}

	// Iterate through the edges and print them out
	std::cout << std::endl << "Edges: " << std::endl;
	for (int i = 0; i < edges->size(); i++) {
		SimpleEdge* e = (*edges)[i];
		CGALPoint* src = (*vertices)[e->u];
		CGALPoint* tar = (*vertices)[e->v];
		//std::cout << edgeWeightMap->at(e).weight << " (" << (*g)[src].pt << ") (" << (*g)[tar].pt << ")" << std::endl;
		std::cout << e->weight << " (" << (*src) << ") (" << (*tar) << ")" << std::endl;
	}

	std::cout << std::endl;
}

CDT* computeCdt(VertexVector* vertices, EdgeVector* edges, boost::unordered_map<TriVertexHandle, VertexIndex>** handlesToIndex) {
	CDT* cdt = new CDT();

	(*handlesToIndex) = new boost::unordered_map<TriVertexHandle, VertexIndex>();
	boost::unordered_map<VertexIndex, TriVertexHandle> vertexHandles;
	for (int i = 0; i < vertices->size(); i++) {
		CGALPoint* pt = (*vertices)[i];
		TriVertexHandle vHandle = cdt->insert(*pt);
		vertexHandles.emplace(i, vHandle);
		(*handlesToIndex)->emplace(vHandle, i);
	}

	// Insert constraint edges
	for (int i = 0; i < edges->size(); i++) {
		SimpleEdge* edge = (*edges)[i];
		VertexIndex u = edge->u;
		VertexIndex v = edge->v;
		TriVertexHandle uH = vertexHandles[u];
		TriVertexHandle vH = vertexHandles[v];
		cdt->insert_constraint(uH, vH);
	}

	assert(cdt->is_valid());

	return cdt;
}

// Given 3 colinear points, if the 2 furthest points are constrained, how is this delaunay triangulated?
// It seems the outcome is implementation dependent. CGAL will split the constraint into 2 edges.
// Other implementations might allow for overlapping, collinear edges. Either way, the current plan
// is to assume that an input forest, F, with collinear edges will be replaced by smaller constraint edges.
EdgeVector* newConstraintSetFromCdt(CDT* cdt, VertexVector* originalVertices) {
	EdgeVector* newEdgeVector = new EdgeVector();

	boost::unordered_map<CGALPoint, VertexIndex> vertexIndex;
	for (int i = 0; i < originalVertices->size(); i++) {
		vertexIndex[*(*originalVertices)[i]] = i;
	}

	// Add edges to graph
	for (CDT::Edge_iterator eit = cdt->edges_begin(); eit != cdt->edges_end(); ++eit) {
		CDT::Edge cgal_e = *eit;

		if (cdt->is_constrained(cgal_e)) {
			// Assumes point coord are unique
			CGALSegment segement = cdt->segment(cgal_e);
			CGALPoint cgal_u = segement.point(0);
			CGALPoint cgal_v = segement.point(1);
			VertexIndex u = vertexIndex[cgal_u];
			VertexIndex v = vertexIndex[cgal_v];

			SimpleEdge* edge = new SimpleEdge(u, v, 0);
			newEdgeVector->push_back(edge);
		}
	}
	return newEdgeVector;
}

TriVertexHandle OppositeOfEdge(TriVertexHandle ev0, TriVertexHandle ev1, TriFaceHandle f) {
	TriVertexHandle v0 = f->vertex(0);
	if (v0 != ev0 && v0 != ev1) {
		return v0;
	}

	TriVertexHandle v1 = f->vertex(1);
	if (v1 != ev0 && v1 != ev1) {
		return v1;
	}

	TriVertexHandle v2 = f->vertex(2);
	if (v2 != ev0 && v2 != ev1) {
		return v2;
	}

	assert(false); // Unreachable
	return NULL;
}

bool IsInsideCircle(CGALPoint* p, CGALPoint* q, CGALPoint* t) {
	if (SHOW_DEBUG) {
		std::cout << "(" << (*p) << ") (" << (*q) << ") (" << (*t) << ")" << std::endl;
	}
	CGALCircle c(*p, *q);
	CGAL::Bounded_side side = c.bounded_side(*t);

	// v is outside or on the circumcircle of f
	return side == CGAL::Bounded_side::ON_BOUNDED_SIDE
		|| side == CGAL::Bounded_side::ON_BOUNDARY;
}

void computeCgg(
	VertexVector* vertices,
	EdgeVector* edges,
	boost::unordered_map<TriVertexHandle, VertexIndex>** handlesToIndex,
	EdgeVector** NewEdges,
	EdgeVector** S_Edges) {

	CDT* cdt = computeCdt(vertices, edges, handlesToIndex);
	(*NewEdges) = newConstraintSetFromCdt(cdt, vertices);

	boost::unordered_set<TriEdge>* S = new boost::unordered_set<TriEdge>();

	TriVertexHandle infiniteVertex = cdt->infinite_vertex();
	int edgeCount = 0;

	boost::chrono::high_resolution_clock::time_point startTotal = boost::chrono::high_resolution_clock::now();

	for (FiniteEdgeIter iter = cdt->finite_edges_begin(); iter != cdt->finite_edges_end(); ++iter) {
		edgeCount++;

		// typedef std::pair<Face_handle, int> Edge;
		TriEdge e = *iter;
		int eIndex = e.second;

		// Edge shared by faces f0 and f1
		TriFaceHandle f0 = e.first;
		TriFaceHandle f1 = e.first->neighbor(eIndex);

		// Vertex opposite of edge e in f0, f1
		TriVertexHandle opp0 = f0->vertex(eIndex);

		// Vertex of edge e
		TriVertexHandle e0 = f0->vertex(f0->cw(eIndex));
		TriVertexHandle e1 = f0->vertex(f0->ccw(eIndex));

		TriVertexHandle opp1 = OppositeOfEdge(e0, e1, f1);

		if (SHOW_DEBUG) {
			// Vertex endpoint v0, v1 of edge e for f0 (doesn't seem to be possible to identify the same edge for f1 in a similar way)
			std::cout << eIndex << std::endl;
			std::cout << "(" << *e0 << ") (" << *e1 << ")" << std::endl;
			std::cout << "(" << *opp0 << ")" << std::endl;
		}

		if (e0 != infiniteVertex && e1 != infiniteVertex) {
			CGALPoint p = e0->point();
			CGALPoint q = e1->point();

			bool addToConstraint = false;

			if (opp0 != infiniteVertex) {
				CGALPoint t = opp0->point();
				if (IsInsideCircle(&p, &q, &t)) {
					addToConstraint = true;
				}
			}

			if (opp1 != infiniteVertex) {
				CGALPoint t = opp1->point();
				if (IsInsideCircle(&p, &q, &t)) {
					addToConstraint = true;
				}
			}

			if (addToConstraint) {
				S->emplace(e);
			}
		}
	}

	(*S_Edges) = new EdgeVector();
	for (boost::unordered_set<TriEdge>::iterator iter = S->begin(); iter != S->end(); ++iter) {
		TriEdge e = *iter;
		int eIndex = e.second;
		TriFaceHandle f0 = e.first;
		TriVertexHandle e0 = f0->vertex(f0->cw(eIndex));
		TriVertexHandle e1 = f0->vertex(f0->ccw(eIndex));

		VertexIndex u = (**handlesToIndex)[e0];
		VertexIndex v = (**handlesToIndex)[e1];

		(*S_Edges)->push_back(new SimpleEdge(u, v, 0));
	}

	boost::chrono::high_resolution_clock::time_point endTotal = boost::chrono::high_resolution_clock::now();
	boost::chrono::milliseconds total = (boost::chrono::duration_cast<boost::chrono::milliseconds>(endTotal - startTotal));
	printDuration("Total duration", total);

	delete S;
	delete cdt;
}

EdgeVector* intersectInputSetWithConstraintSet(EdgeVector* inputSet, EdgeVector* constraintSet) {
	boost::unordered_set<SimpleEdge>* NewEdges_Hashset = new boost::unordered_set<SimpleEdge>();
	for (int i = 0; i < inputSet->size(); i++) {
		SimpleEdge* e = (*inputSet)[i];
		NewEdges_Hashset->emplace((*e));
	}

	EdgeVector* intersect = new EdgeVector();
	for (int i = 0; i < constraintSet->size(); i++) {
		SimpleEdge* se = (*constraintSet)[i];
		if (NewEdges_Hashset->count(*se) > 0) {
			intersect->push_back(new SimpleEdge(se->u, se->v, se->weight));
		}
	}

	return intersect;
}

bool containsEdge(boost::unordered_set<SimpleEdge>* edgeSet, SimpleEdge* edge) {
	SimpleEdge se(edge->u, edge->v, 0);
	return edgeSet->count(se) > 0;
}

// True if A a subgraph of B
bool isSubgraph(VertexVector* vertices, EdgeVector* a, EdgeVector* b) {
	boost::unordered_set<SimpleEdge>* bEdgeSet = createSimpleEdgeSet(b);

	// Iterate through the edges
	for (int i = 0; i < a->size(); i++) {
		SimpleEdge* edge = (*a)[i];
		if (!containsEdge(bEdgeSet, edge)) {
			CGALPoint* u = (*vertices)[edge->u];
			CGALPoint* v = (*vertices)[edge->v];
			EdgeWeight weight = CGAL::squared_distance(*u, *v);
			std::cout << "b contains edge not in a: " << weight << " (" << *u << ") (" << *v << ")" << std::endl;
			return false;
		}
	}

	return true;
}

EdgeVector* convertCdtToGraph(VertexVector* vertices, CDT* cdt) {
	EdgeVector* edgeVec = new EdgeVector();
	// edgeVec->reserve(cdt->number_of_faces() * 3); // Upper bound on number of edges (typically too much)

	// Map CGALPoint -> VertexIndex
	boost::unordered_map<CGALPoint, VertexIndex> vertexIndex;
	for (int i = 0; i < vertices->size(); i++) {
		vertexIndex[*(*vertices)[i]] = i;
	}

	// Add edges to graph
	for (CDT::Edge_iterator eit = cdt->edges_begin(); eit != cdt->edges_end(); ++eit) {
		CDT::Edge cgal_e = *eit;
		CGALSegment segement = cdt->segment(cgal_e);
		CGALPoint cgal_u = segement.point(0);
		CGALPoint cgal_v = segement.point(1);
		VertexIndex u = vertexIndex[cgal_u];
		VertexIndex v = vertexIndex[cgal_v];

		SimpleEdge* edge = new SimpleEdge(u, v, 0);
		edgeVec->push_back(edge);
	}

	return edgeVec;
}

bool isCdtSubgraph(VertexVector* vertices, EdgeVector* edgesF, EdgeVector* edgesS) {
	boost::unordered_map<TriVertexHandle, VertexIndex>* handlesToIndex;
	CDT* cdtS = computeCdt(vertices, edgesS, &handlesToIndex); // CDT of mimimum edge constraint
	EdgeVector* ev_cdtS = convertCdtToGraph(vertices, cdtS);
	bool res = isSubgraph(vertices, edgesF, ev_cdtS);

	delete ev_cdtS;
	delete cdtS;

	return res;
}

bool isCggSubgraph(VertexVector* vertices, EdgeVector* edgesF, EdgeVector* edgesS) {
	boost::unordered_map<TriVertexHandle, VertexIndex>* handlesToIndex;
	EdgeVector* newEdgesS;
	EdgeVector* cggS;
	computeCgg(vertices, edgesS, &handlesToIndex, &newEdgesS, &cggS);
	bool res = isSubgraph(vertices, newEdgesS, cggS);

	delete cggS;
	delete newEdgesS;
	delete handlesToIndex;

	return res;
}

int main(int argc, char* argv[]) {
	const char* vertFile = (argc > 2) ? argv[1] : NULL;
	const char* edgeFile = (argc > 2) ? argv[2] : NULL;

	VertexVector* vertices = NULL;
	EdgeVector* edges = NULL;

	if (vertFile == NULL || edgeFile == NULL) {
		// Random graph
		createRandomPlaneForest(1000, 1000, 1000, &vertices, &edges);
	}
	else {
		// Load graph from file
		parseGraph(vertFile, edgeFile, &vertices, &edges);
	}

	// Compute CGG(V, E)
	boost::unordered_map<TriVertexHandle, VertexIndex>* handlesToIndex;
	EdgeVector* NewEdges;
	EdgeVector* cggS;
	computeCgg(vertices, edges, &handlesToIndex, &NewEdges, &cggS);
	EdgeVector* S = intersectInputSetWithConstraintSet(NewEdges, cggS);
	std::cout << "Edges in E: " << NewEdges->size() << " Edges in S: " << S->size() << " Ratio: " << (double)((double)S->size() / (double)NewEdges->size()) << std::endl;

	// Validatation:
	// S ⊆ E
	if (!isSubgraph(vertices, S, NewEdges)) {
		std::cout << "Error: isSubgraph is false" << std::endl;
	}

	// F ⊆ CDT(V, S) (Note: CGG is contained in CDT)
	if (!isCdtSubgraph(vertices, NewEdges, S)) {
		std::cout << "Error: isCdtSubgraph is false" << std::endl;
	}

	// F ⊆ CGG(V, S)
	if (!isCggSubgraph(vertices, NewEdges, S)) {
		std::cout << "Error: isCggSubgraph is false" << std::endl;
	}

	deleteEdgeVector(S);
	deleteEdgeVector(cggS);
	deleteEdgeVector(NewEdges);
	delete handlesToIndex;
	deleteEdgeVector(edges);
	deleteVerticesVector(vertices);

	return 0;
}