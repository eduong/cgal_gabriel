#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <boost/chrono.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/pending/disjoint_sets.hpp>

#include "GraphDefs.h"
#include "timer.h"
#include "verify.h"

#include <fstream>
#include <string>

#define SHOW_DEBUG false

void printGraph(const char* title, BoostGraph* g) {
	std::cout << std::endl << "=== " << title << std::endl;

	// Iterate through the vertices and print them out
	std::cout << "Vertices: " << std::endl;
	std::pair<VertexIter, VertexIter> vp;
	for (vp = boost::vertices(*g); vp.first != vp.second; ++vp.first) {
		Vertex v = *vp.first;
		std::cout << "index: " << v << " (" << (*g)[v].pt << ")" << std::endl;
	}

	// Iterate through the edges and print them out
	std::cout << std::endl << "Edges: " << std::endl;
	std::pair<EdgeIter, EdgeIter> ep;
	EdgeIter ei, ei_end;
	for (tie(ei, ei_end) = boost::edges(*g); ei != ei_end; ++ei) {
		Edge e = *ei;
		Vertex src = source(e, *g);
		Vertex tar = target(e, *g);
		//std::cout << edgeWeightMap->at(e).weight << " (" << (*g)[src].pt << ") (" << (*g)[tar].pt << ")" << std::endl;
		std::cout << (*g)[e].weight << " (" << (*g)[src].pt << ") (" << (*g)[tar].pt << ")" << std::endl;
	}
	std::cout << std::endl;
}

/**
* Naive linear time intersection
* Returns true if edge (u, v) intersects an edge in g, otherwise false
**/
bool DoesIntersect(BoostGraph* g, Vertex u, Vertex v) {
	CGALPoint uPt = (*g)[u].pt;
	CGALPoint vPt = (*g)[v].pt;
	CGALSegment segUV(uPt, vPt);

	std::pair<EdgeIter, EdgeIter> ep;
	EdgeIter ei, ei_end;
	for (tie(ei, ei_end) = boost::edges(*g); ei != ei_end; ++ei) {
		Edge e = *ei;
		Vertex src = source(e, *g);
		Vertex tar = target(e, *g);
		CGALSegment seg((*g)[src].pt, (*g)[tar].pt);

		CGAL::cpp11::result_of<CGALIntersect(CGALSegment, CGALSegment)>::type result = intersection(seg, segUV);
		if (result) {
			if (const CGALSegment* s = boost::get<CGALSegment>(&*result)) {
				//std::cout << *s << std::endl;
				return true;
			}
			else if (const CGALPoint* p = boost::get<CGALPoint >(&*result)) {
				//std::cout << " i " << *p;
				// Ignore intersection at segment endpoints
				if (*p != uPt && *p != vPt) {
					return true;
				}
			}
		}
	}
	return false;
}

BoostGraph* CreateRandomPlaneForest(int numVertices, int radius, int edgeRolls, double edgeProbability) {
	BoostGraph* g = new BoostGraph();

	CGAL::Random_points_on_circle_2<TriPoint, Creator> randPts(radius);

	// Generate vertices with random coordinated within bounds
	for (int i = 0; i < numVertices; i++) {
		Vertex v = add_vertex(*g);
		(*g)[v].pt = (*randPts++);
	}

	printGraph("Random vert graph", g);

	// Define edge random gen
	boost::uniform_real<> edgeRange(0, 1);
	boost::variate_generator<boost::minstd_rand, boost::uniform_real<>> edgeDice(gen, edgeRange);
	edgeDice.engine().seed(static_cast<unsigned int>(std::time(0)));

	boost::uniform_int<> vertexRange(0, numVertices - 1);
	boost::variate_generator<boost::minstd_rand, boost::uniform_int<>> vertexDice(gen, vertexRange);
	vertexDice.engine().seed(static_cast<unsigned int>(std::time(0)));

	std::vector<int> rank(numVertices);
	std::vector<int> parent(numVertices);
	boost::disjoint_sets<int*, int*, boost::find_with_full_path_compression> ds(&rank[0], &parent[0]);
	for (int i = 0; i < rank.size(); i++) {
		rank[i] = i;
		parent[i] = i;
	}

	// Select random vertices u, v for edgeRolls number of times
	// An edge connects u, v:
	//		1. u != v
	//		2. roll <= edgeProbability
	//		3. adding edge(u, v) does not create a cycle
	//		4. edge(u, v) does not intersect any other edge
	for (int i = 0; i < edgeRolls; i++) {
		Vertex u = vertexDice();
		Vertex v = vertexDice();
		double d = edgeDice();
		if (u != v
			&& d <= edgeProbability
			&& ds.find_set(u) != ds.find_set(v)
			&& !DoesIntersect(g, u, v)) {

			// Add edge(u, v)
			std::pair<Edge, bool> result = add_edge(u, v, *g);
			assert(result.second);
			ds.link(u, v);
			//std::cout << " - added";
		}
	}

	return g;
}

CDT* computeCdt(BoostGraph* g) {
	CDT* cdt = new CDT();

	// Insert vertices
	boost::unordered_map<CGALPoint, TriVertexHandle> vertexHandles;
	std::pair<VertexIter, VertexIter> vp;
	for (vp = boost::vertices(*g); vp.first != vp.second; ++vp.first) {
		Vertex v = *vp.first;
		CGALPoint pt = (*g)[v].pt;
		TriVertexHandle vHandle = cdt->insert(pt);
		vertexHandles.emplace(pt, vHandle);
	}

	// Insert constraint edges
	std::pair<EdgeIter, EdgeIter> ep;
	EdgeIter ei, ei_end;
	for (tie(ei, ei_end) = boost::edges(*g); ei != ei_end; ++ei) {
		Edge e = *ei;
		Vertex u = source(e, *g);
		Vertex v = target(e, *g);
		TriVertexHandle uH = vertexHandles[(*g)[u].pt];
		TriVertexHandle vH = vertexHandles[(*g)[v].pt];
		cdt->insert_constraint(uH, vH);
	}

	assert(cdt->is_valid());

	return cdt;
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

int main(int argc, char* argv[]) {
	BoostGraph* f = CreateRandomPlaneForest(10, 10000, 5000, 0.50);

	CDT* cdt = computeCdt(f);

	boost::unordered_set<TriEdge>* s = new boost::unordered_set<TriEdge>();

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
				s->emplace(e);
			}
		}
	}

	boost::chrono::high_resolution_clock::time_point endTotal = boost::chrono::high_resolution_clock::now();
	boost::chrono::milliseconds total = (boost::chrono::duration_cast<boost::chrono::milliseconds>(endTotal - startTotal));
	printDuration("Total duration", total);

	std::cout << "Edges in T: " << edgeCount << " Edges in S: " << s->size() << " Ratio: " << (double)((double)s->size() / (double)edgeCount) << std::endl;

	delete s;

	return 0;
}