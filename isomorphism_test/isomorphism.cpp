/*
    Isomorphism test with no regular graphs
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
#include <vector>
#include <fstream>
#include <string>

/* 	returns true if the first argument goes before 
	the second argument in the specific strict 
	weak ordering it defines, and false otherwise. */
bool compint(int i, int j){return i<j;}
bool compVector(std::vector<int> vectorA, std::vector<int> vectorB){
		if(vectorA.size() < vectorB.size())
			return true;
		else if(vectorB.size() < vectorA.size())
			return false;
		else
			return vectorA < vectorB;
}

int nodeCount(lemon::ListGraph& g){
	int count=0;
    for (lemon::ListGraph::NodeIt n(g); n!=lemon::INVALID; ++n) 
		++count;
	return count;
}	
int edgeCount(lemon::ListGraph& g){
	int count=0;
	for(lemon::ListGraph::EdgeIt e(g); e!=lemon::INVALID; ++e)
		++count;
	return count;
}
int nodeDegree(lemon::ListGraph& g, lemon::ListGraph::Node n){
	int count=0;
	for(lemon::ListGraph::IncEdgeIt e(g, n); e != lemon::INVALID; ++e)
		++count;
	return count;
}
std::vector< std::vector<int> > neighborhoodDegree(lemon::ListGraph& g){
	std::vector< std::vector<int> > retorno;
	for (lemon::ListGraph::NodeIt v1(g); v1!=lemon::INVALID; ++v1){
		std::vector<int> listaVizinhos;
		int grauV1 = nodeDegree(g, v1);
		listaVizinhos.push_back(grauV1);
		for(lemon::ListGraph::IncEdgeIt e(g, v1); e != lemon::INVALID; ++e){
			lemon::ListGraph::Node vizinho = g.oppositeNode(v1, e);
			int grauVizinho = nodeDegree(g, vizinho);
			listaVizinhos.push_back(grauVizinho);
		}
		sort(listaVizinhos.begin(), listaVizinhos.end(), compint);
		retorno.push_back(listaVizinhos);
	}
	//sort(retorno.begin(), retorno.end(), compVector);
	return retorno;
}
	

std::string degreeSeq(lemon::ListGraph& g){

	int i=0;
		
	/* iniciar seq */
	std::vector<int> seq;
	
	for(lemon::ListGraph::NodeIt v(g); v != lemon::INVALID; ++v){
		int cont = 0;
		for(lemon::ListGraph::IncEdgeIt e(g, v); e != lemon::INVALID; ++e){
			cont++;
		}
		
		/* incluir cont em seq */
		seq.push_back(cont);
	}
	
	/* ordenar seq */
	sort(seq.begin(), seq.end(), compint);
	
	/* retornar seq */
	std::stringstream result;
	std::copy(seq.begin(), seq.end(), std::ostream_iterator<int>(result, " "));
	return result.str();
}

void readGraph(std::string filename, lemon::ListGraph& g){

	// Le arquivo.
	std::ifstream file(filename.c_str());
	std::stringstream input;
	if(file){
		input << file.rdbuf();
		file.close();
	} else {
		std::cout << "File not found: " << filename << std::endl;
		exit(1);
	}
	
	// Parser.
	lemon::ListGraph::NodeMap<std::string> rotulo(g);
	lemon::GraphReader<lemon::ListGraph>(g, input).
		run();
}

int main(int argc, char *argv[]){
	if(argc != 3){
		std::cout << "uso: " << argv[0] << " <filename> <filename>" << std::endl;
		return 1;
	}

	// Read first graph
	std::string filenameG = argv[1];	
	lemon::ListGraph G;
	readGraph(filenameG, G);

	// Read second graph
	std::string filenameH = argv[2];
	lemon::ListGraph H;
	readGraph(filenameH, H);
	
	bool isomorfos = false;
	
	int countNodeG = nodeCount(G);
	int countNodeH = nodeCount(H);
	if(countNodeG != countNodeH){
		std::cout << "NAO" << std::endl;
		return 0;
	}
	int countEdgeG = edgeCount(G);
	int countEdgeH = edgeCount(H);
	if(countEdgeG != countEdgeH){
		std::cout << "NAO" << std::endl;
		return 0;
	}		
	std::string listaG = degreeSeq(G);
	std::string listaH = degreeSeq(H);
	if(listaG != listaH){
		std::cout << "NAO" << std::endl;
		return 0;
	}
	std::vector< std::vector<int> > listaVizinhosG = neighborhoodDegree(G);
	std::vector< std::vector<int> > listaVizinhosH = neighborhoodDegree(H);
	bool existeElementosDuplicados = false;
	for(std::vector< std::vector<int> >::size_type i = 0; i < listaVizinhosG.size() - 1; i++) {
		if(std::equal(listaVizinhosG[i].begin(), listaVizinhosG[i].end(), listaVizinhosG[i+1].begin()) || 
			std::equal(listaVizinhosH[i].begin(), listaVizinhosH[i].end(), listaVizinhosH[i+1].begin())){
			existeElementosDuplicados = true;
			break;
		}
	}
	bool listasIguais = true;
	if(existeElementosDuplicados){
		for(std::vector< std::vector<int> >::size_type i = 0; i < listaVizinhosG.size(); i++) {
			if(!std::equal(listaVizinhosG[i].begin(), listaVizinhosG[i].end(), listaVizinhosH[i].begin()))
				listasIguais = false;
				break;
		}
		if(listasIguais){
			std::cout << "SIM" << std::endl; // In fact I'm not sure that they're isomorphics. I'm guessing that they're!! 
			return 0;
		}else{
			std::cout << "NAO" << std::endl;
			return 0;
		}
	}else{
		for(std::vector< std::vector<int> >::size_type i = 0; i < listaVizinhosG.size(); i++) {
			if(!std::equal(listaVizinhosG[i].begin(), listaVizinhosG[i].end(), listaVizinhosH[i].begin()))
				listasIguais = false;
				break;
		}
		if(listasIguais){
			std::cout << "SIM" << std::endl;
			return 0;
		}else{
			std::cout << "NAO" << std::endl;
			return 0;
		}
	}
	

}
