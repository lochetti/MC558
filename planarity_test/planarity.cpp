/*
    Planarity test
    Copyright (C) 2012  Renato Tadeu Lochetti

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/adaptors.h>
#include <lemon/connectivity.h>
#include <lemon/bfs.h>
#include <queue>
#include <lemon/path.h>
#include <fstream>
#include <string>
#define INFINITY 999999;
using namespace std;
typedef enum color { WHITE, GRAY, BLACK } color;

typedef struct attributes {
	color c; // color
	int d; // distance
	lemon::ListGraph::Node p; // predecessor
} attributes;

void readGraph(string filename, lemon::ListGraph& g){

	// Read file.
	ifstream file(filename.c_str());
	stringstream input;
	if(file){
		input << file.rdbuf();
		file.close();
	} else {
		cout << "File not found: " << filename << endl;
		exit(1);
	}
	
	// Parser.
	lemon::ListGraph::NodeMap<string> rotulo(g);
	lemon::GraphReader<lemon::ListGraph>(g, input).run();
}

bool eulersFormulaVerification(lemon::ListGraph& g){
	int nodes = countNodes(g);
	int edges = countEdges(g);
	return edges <= 3*nodes - 6;
}

void dfsRecursive(lemon::ListGraph& g, lemon::SubGraph<lemon::ListGraph>& gi, lemon::ListGraph::NodeMap<attributes> &dfs_attrs, lemon::ListGraph::Node s){
	dfs_attrs[s].c = GRAY;
	for(lemon::ListGraph::IncEdgeIt e(g, s); e != lemon::INVALID; ++e){
		lemon::ListGraph::Node otherNode = g.oppositeNode(s, e);
		if(dfs_attrs[s].p == lemon::INVALID ||g.id(otherNode) != g.id(dfs_attrs[s].p)){
			if(dfs_attrs[otherNode].c == GRAY){ // found one cicle. Go back to your father does not count!
				gi.enable(otherNode);
				gi.enable(e);
				lemon::ListGraph::Node aux = s;
				lemon::ListGraph::Edge auxE;
				do{
					gi.enable(aux);
					auxE = findEdge(g, aux, dfs_attrs[aux].p);
					gi.enable(auxE);
					aux = dfs_attrs[aux].p;
				}while(g.id(aux) != g.id(otherNode));
				return;
			}else{
				dfs_attrs[otherNode].p = s;
				dfsRecursive(g, gi, dfs_attrs, otherNode);
				if(countNodes(gi) > 0){ //we already have a cicle!
					return;
				}
			}
		}
	}
	dfs_attrs[s].c = BLACK;
}

void dfs(lemon::ListGraph& g, lemon::SubGraph<lemon::ListGraph>& gi, lemon::ListGraph::NodeMap<attributes> &dfs_attrs){
	for (lemon::ListGraph::NodeIt n(g); n!=lemon::INVALID ; ++n){
		dfs_attrs[n].c = WHITE;
		dfs_attrs[n].p = lemon::INVALID;
	}
	for(lemon::ListGraph::NodeIt n(g); n!=lemon::INVALID; ++n){
		if(dfs_attrs[n].c == WHITE){
			dfsRecursive(g, gi, dfs_attrs, n);
			if(countNodes(gi) > 0){ //we already have a cicle!
				return;
			}
		}
	}
}

void produceFragments(lemon::ListGraph& g, lemon::SubGraph<lemon::ListGraph>& gi, lemon::ListGraph::NodeMap<attributes> &dfs_attrs, 
						vector<lemon::ListGraph::Edge>& edgeFragments, vector< vector<lemon::ListGraph::Node> >& otherFragments,
						vector< vector<lemon::ListGraph::Node> >& edgeFragmentsAttach, vector< vector<lemon::ListGraph::Node> >& otherFragmentsAttach){
	
	//edges fragments first.
	for(lemon::ListGraph::NodeIt n(g); n!=lemon::INVALID ; ++n){
		if(gi.status(n)){
			for(lemon::ListGraph::IncEdgeIt e(g, n); e != lemon::INVALID; ++e){
				lemon::ListGraph::Node otherNode = g.oppositeNode(n, e);
				if(gi.status(otherNode) && !(dfs_attrs[n].p == otherNode) && !(dfs_attrs[otherNode].p == n) && !gi.status(e)){
					edgeFragments.push_back(e); // edge 'e' is a fragment!
					vector<lemon::ListGraph::Node> attchments;
					attchments.push_back(otherNode);
					attchments.push_back(n);
					edgeFragmentsAttach.push_back(attchments); 
				}
			}
		}
	}

	//other fragments.
	lemon::ListGraph::NodeMap<int> components(g);	
	lemon::ListGraph::NodeMap<bool> node_filterL(g, true);
	lemon::ListGraph::EdgeMap<bool> edge_filterL(g, true);
	lemon::SubGraph<lemon::ListGraph> gMinusGi(g, node_filterL, edge_filterL);
	//compute G - Gi
	for(lemon::ListGraph::NodeIt n(g); n!=lemon::INVALID ; ++n){
		if(gi.status(n)){
			gMinusGi.disable(n);
		}
	}
	for(lemon::ListGraph::EdgeIt e(g); e!=lemon::INVALID; ++e){
		if(gi.status(e)){
			gMinusGi.disable(e);
		}
	}
	//calculate the connectedComponenets of GMinusGi. each connected component will be a fragment.
	int connectedParts = connectedComponents(gMinusGi, components);
	for(int i = 0; i < connectedParts; i++){
		vector<lemon::ListGraph::Node> empty;
		vector<lemon::ListGraph::Node> empty2;		
		otherFragments.push_back(empty);
		otherFragmentsAttach.push_back(empty2);
	}
	for(lemon::ListGraph::NodeIt n(g); n!=lemon::INVALID ; ++n){
		if(!gi.status(n)){
			int ccNumber = components[n];
			otherFragments.at(ccNumber).push_back(n);
		}
	}
	for(int i = 0; i < otherFragments.size(); i++){
		vector<lemon::ListGraph::Node> fragmentAttachment;
		for(int j = 0; j < otherFragments[i].size(); j++){
			lemon::ListGraph::Node n = otherFragments.at(i).at(j);
			for(lemon::ListGraph::IncEdgeIt e(g, n); e != lemon::INVALID; ++e){
				lemon::ListGraph::Node possibleAttach = g.oppositeNode(n, e);
				if(!gi.status(n) && gi.status(possibleAttach) && !gi.status(e) && find(fragmentAttachment.begin(), fragmentAttachment.end(), possibleAttach) == fragmentAttachment.end()){
					fragmentAttachment.push_back(possibleAttach);
				}
			}
		}
		otherFragmentsAttach.at(i) = fragmentAttachment;
	}
}
bool fragmentsFacesBuildAndVerification(lemon::ListGraph& g, vector< vector<lemon::ListGraph::Node> > faces,
										vector< vector<lemon::ListGraph::Node> > edgeFragmentsAttach, vector< vector<lemon::ListGraph::Node> > otherFragmentsAttach,
										vector< vector<int> >& edgeFragmentsFaces, vector< vector<int> >& otherFragmentsFaces){
										
	//edges first
	for(int i = 0; i < edgeFragmentsAttach.size(); i++){
		vector<int> empty;
		edgeFragmentsFaces.push_back(empty);
		for(int j = 0; j < faces.size(); j++){
			bool thisFaceHasAllNodes = true;
			for(int k = 0; k < edgeFragmentsAttach.at(i).size(); k++){
				lemon::ListGraph::Node attachNode = edgeFragmentsAttach.at(i).at(k);
				if(find(faces.at(j).begin(), faces.at(j).end(), attachNode) == faces.at(j).end()){
					thisFaceHasAllNodes = false;
					break;
				}
			}
			if(thisFaceHasAllNodes){
				edgeFragmentsFaces.at(i).push_back(j);
			}
		}
		if(edgeFragmentsFaces.at(i).size() == 0 || edgeFragmentsFaces.at(i).empty()){
			cout << "nao" << endl;
			exit(0);
		}
	}
	
	//other fragments.
		for(int i = 0; i < otherFragmentsAttach.size(); i++){
		vector<int> empty;
		otherFragmentsFaces.push_back(empty);
		for(int j = 0; j < faces.size(); j++){
			bool thisFaceHasAllNodes = true;
			for(int k = 0; k < otherFragmentsAttach.at(i).size(); k++){
				lemon::ListGraph::Node attachNode = otherFragmentsAttach.at(i).at(k);
				if(find(faces.at(j).begin(), faces.at(j).end(), attachNode) == faces.at(j).end()){
					thisFaceHasAllNodes = false;
					break;
				}
			}
			if(thisFaceHasAllNodes){
				otherFragmentsFaces.at(i).push_back(j);
			}
		}
		if(otherFragmentsFaces.at(i).size() == 0 || otherFragmentsFaces.at(i).empty()){
			cout << "nao" << endl;
			exit(0);
		}
	}
	return true;
}

void bfs(lemon::ListGraph &g, lemon::SubGraph<lemon::ListGraph>& gi, lemon::ListGraph::NodeMap<attributes>& a, lemon::ListGraph::Node s, lemon::ListGraph::Node d) {
	for (lemon::SubGraph<lemon::ListGraph>::NodeIt n(gi); n != lemon::INVALID; ++n) {
		a[n].c = WHITE;
		a[n].p = lemon::INVALID;
	}
	a[s].c = GRAY;

	queue<lemon::ListGraph::Node> q;
	q.push(s);
	lemon::ListGraph::Node u, v;

	while (!q.empty()) {
		u = q.front();
		q.pop();
		for (lemon::SubGraph<lemon::ListGraph>::IncEdgeIt e(gi, u); e != lemon::INVALID; ++e) {
			v = gi.oppositeNode(u, e);
			if (a[v].c == WHITE) {
				a[v].c = GRAY;       
				a[v].p = u;           
				q.push(v);            
			}
			if(g.id(d) == g.id(v)){
				return;
			}
		}
		a[u].c = BLACK;
	}
}

void findPathThroughFragmentOther(lemon::ListGraph& g, lemon::SubGraph<lemon::ListGraph>& gi, vector<lemon::ListGraph::Node> attachments, vector<lemon::ListGraph::Node>& pathVector, vector<lemon::ListGraph::Node> frag){
	lemon::ListGraph::NodeMap<bool> node_filterR(g, false);
	lemon::ListGraph::EdgeMap<bool> edge_filterR(g, false);
	lemon::SubGraph<lemon::ListGraph> giPlusFrag(g, node_filterR, edge_filterR);

	lemon::ListGraph::Node source = attachments.at(0);
	lemon::ListGraph::Node destination = attachments.at(1);
	giPlusFrag.enable(source);
	giPlusFrag.enable(destination);

	for(int i = 0; i < frag.size(); i++){
		lemon::ListGraph::Node n = frag.at(i);
		giPlusFrag.enable(n);
	}
	for(int i = 0; i < frag.size(); i++){
		lemon::ListGraph::Node n = frag.at(i);
		for(lemon::ListGraph::IncEdgeIt e(g, n); e != lemon::INVALID; ++e){
			lemon::ListGraph::Node v = g.oppositeNode(n, e);
			if(giPlusFrag.status(v)){
				giPlusFrag.enable(e);
			}
		}
	}
	
	lemon::ListGraph::NodeMap<attributes> bfs_attrs(g);
	
	//runs a bfs in order to find a path from Source to Destination
	bfs(g, giPlusFrag, bfs_attrs, source, destination);
	pathVector.push_back(destination);
	lemon::ListGraph::Node aux = destination;
	while(bfs_attrs[aux].p != lemon::INVALID){
		aux = bfs_attrs[aux].p;
		pathVector.insert(pathVector.begin(), aux);
	}
}


int main(int argc, char *argv[]){
	
	if(argc != 2){
		cout << "uso: " << argv[0] << " <filename>" << endl;
		return 1;
	}

	// Read graph.
	string filename = argv[1];	
	lemon::ListGraph g;
	readGraph(filename, g);
	
	//Euler's formula validation
	if(!eulersFormulaVerification(g)){
		cout << "nao" << endl;
		return 0;
	}
	lemon::ListGraph::NodeMap<attributes> dfs_attrs(g);
	lemon::ListGraph::NodeMap<bool> node_filter(g, false);
	lemon::ListGraph::EdgeMap<bool> edge_filter(g, false);
	lemon::SubGraph<lemon::ListGraph> gi(g, node_filter, edge_filter);
	
	//Runs a modified dfs in order to find a cicle in the inicial graph 'g' and put him in 'gi'.
	dfs(g, gi, dfs_attrs);
	//do we have a cicle?
	if(!countNodes(gi) > 0){
		cout << "sim" << endl;
		return 0;
	}
	
	//Initially we have 2 faces. One internal and other external. Both faces have the same nodes.
	vector< vector<lemon::ListGraph::Node> > faces;
	vector<lemon::ListGraph::Node> initialFace;
	vector<lemon::ListGraph::Edge> edgesAtTheInicialFace;
	lemon::ListGraph::Node firstNodeOfTheFace;
	lemon::ListGraph::Node actualNodeOfTheFace;
	for(lemon::ListGraph::NodeIt n(g); n!=lemon::INVALID ; ++n){
		if(gi.status(n)){
			firstNodeOfTheFace = n;
			actualNodeOfTheFace = n;
			break;
		}
	}
	
	do{ //calculate the cicle of the face.
		for(lemon::ListGraph::IncEdgeIt e(g, actualNodeOfTheFace); e != lemon::INVALID; ++e){
			if(gi.status(e) && find(edgesAtTheInicialFace.begin(), edgesAtTheInicialFace.end(), e) == edgesAtTheInicialFace.end()){
				lemon::ListGraph::Node nodeAux = g.oppositeNode(actualNodeOfTheFace, e);
				if(gi.status(nodeAux)){
					edgesAtTheInicialFace.push_back(e);
					initialFace.push_back(nodeAux);
					actualNodeOfTheFace = nodeAux;
					break;
				}
			}
		}
	
	}while(g.id(firstNodeOfTheFace) != g.id(actualNodeOfTheFace));;
	faces.push_back(initialFace); //external
	faces.push_back(initialFace); //internal
	
	//main looping
	while(countNodes(gi) < countNodes(g)){
		vector<lemon::ListGraph::Edge> edgeFragments; //will represent the edge type fragments
		vector< vector<lemon::ListGraph::Node> > edgeFragmentsAttach; //will represent the attachment nodes of each edge type fragment
		vector< vector<lemon::ListGraph::Node> > otherFragments; // will represent the other fragments.
		vector< vector<lemon::ListGraph::Node> > otherFragmentsAttach; //will represent the attachment nodes of each other type fragments.
		vector< vector<int> > edgeFragmentsFaces; //will represent which faces each fragment can be embbeded
		vector< vector<int> > otherFragmentsFaces; //will represent which faces each fragment can be embbeded
		//produces all gi-fragments of the graph g.
		produceFragments(g, gi, dfs_attrs, edgeFragments, otherFragments, edgeFragmentsAttach, otherFragmentsAttach);
		// Compute the faces that each fragment can be embbeded and verify if this number of faces is not null
		if(!fragmentsFacesBuildAndVerification(g, faces, edgeFragmentsAttach, otherFragmentsAttach, edgeFragmentsFaces, otherFragmentsFaces)){
			cout << "nao" << endl;
			exit(0);
		}
		
		//Choose one fragment! 
		int minNumberOfFaces = INFINITY;
		int indexOfTheFragment = INFINITY;
		int typeOfTheFragment = INFINITY; //0 = edge. 1 = other.

		for(int i = 0; i < edgeFragmentsFaces.size(); i++){
			int numberOfFaces = edgeFragmentsFaces.at(i).size();
			if(numberOfFaces < minNumberOfFaces){
				minNumberOfFaces = numberOfFaces;
				indexOfTheFragment = i;
				typeOfTheFragment = 0;
			}
		}
		for(int i = 0; i < otherFragmentsFaces.size(); i++){
			int numberOfFaces = otherFragmentsFaces.at(i).size();
			if(numberOfFaces < minNumberOfFaces){
				minNumberOfFaces = numberOfFaces;
				indexOfTheFragment = i;
				typeOfTheFragment = 1;
			}
		}

		vector<lemon::ListGraph::Node> attachments = typeOfTheFragment == 0 ? edgeFragmentsAttach.at(indexOfTheFragment) : otherFragmentsAttach.at(indexOfTheFragment);
		vector<lemon::ListGraph::Node> path;
		//Produces de path of the fragment that will be embbeded in the Gi graph
		if(typeOfTheFragment == 0){
			path.push_back(attachments.at(0));
			path.push_back(attachments.at(1));
		}else{
			findPathThroughFragmentOther(g, gi, attachments, path, otherFragments.at(indexOfTheFragment));
		}
		
		//re-define the faces.
		vector<int> fragmentFaces = typeOfTheFragment == 0 ? edgeFragmentsFaces.at(indexOfTheFragment) : otherFragmentsFaces.at(indexOfTheFragment);
		int faceSwitch = 1;
		int choosenFace = fragmentFaces.at(0);
		vector<lemon::ListGraph::Node> face = faces.at(choosenFace);
		faces.erase(faces.begin()+choosenFace);
		vector<lemon::ListGraph::Node> newFace1;
		vector<lemon::ListGraph::Node> newFace2;
		lemon::ListGraph::Node startNode = face.at(0);
		lemon::ListGraph::Node attach1 = attachments.at(0);
		lemon::ListGraph::Node attach2 = attachments.at(1);
		lemon::ListGraph::Node actualNode = startNode;
		int faceIterator = 0;
		vector<lemon::ListGraph::Edge> usedEdges;
		while(faceIterator < face.size()){
			actualNode = face.at(faceIterator);
			if(faceSwitch == 1){
				newFace1.push_back(actualNode);
				if(g.id(actualNode) == g.id(attach1) || g.id(actualNode) == g.id(attach2)){
					faceSwitch*=-1;
					newFace2.push_back(actualNode);
					if(g.id(actualNode) != g.id(path.at(0))){
						reverse(path.begin(),path.end());
					}
					for(int i = 0; i < path.size(); i++){
						lemon::ListGraph::Node n = path.at(i);
						if(find(face.begin(), face.end(), n) == face.end()){
							newFace1.push_back(n);
						}
					}
				}
			}else{
				newFace2.push_back(actualNode);
				if(g.id(actualNode) == g.id(attach1) || g.id(actualNode) == g.id(attach2)){
					faceSwitch*=-1;
					newFace1.push_back(actualNode);
					if(g.id(actualNode) != g.id(path.at(0))){
						reverse(path.begin(),path.end());
					}
					for(int i = 0; i < path.size(); i++){
						lemon::ListGraph::Node n = path.at(i);
						if(find(face.begin(), face.end(), n) == face.end()){
							newFace2.push_back(n);
						}
					}
				}
			}
			faceIterator++;
		}
		//new faces !
		faces.push_back(newFace1);
		faces.push_back(newFace2);
		
		//put the path together if the gi graph
		for(int i = 0; i < (path.size() - 1); i++){
			lemon::ListGraph::Node node1 = path.at(i);
			int j = i+1;
			lemon::ListGraph::Node node2 = path.at(j);
			lemon::ListGraph::Edge edge = findEdge(g, node1, node2);
			gi.enable(node1);
			gi.enable(edge);
			gi.enable(node2);
		}
	}
	
	//if we're here so far, the graph is planar.
	cout << "sim" << endl;
	return 0;
}